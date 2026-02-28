//! Gaia DR3 life-map ingestion lane (streaming).
//!
//! Reads a Gaia-style CSV export, applies derived habitability + entropy stage gates
//! per star, and emits a map + ranked signal targets.

use gutoe_physics::{
    classify_stage, derive_thresholds_and_multipliers, evaluate_abiogenesis_gate,
    evaluate_great_filter,
    habitability_score, infer_component_from_position, is_habitable, main_sequence_lifetime_gyr,
    stage_entropy_multiplier, AbiogenesisWindows, CivilizationStage, GalacticLifeSeed,
    GreatFilterCivilizationInput, GreatFilterWindows, KAUFFMAN_CLOSURE_THRESHOLD,
    MILKY_WAY_STELLAR_COUNT_ESTIMATE, SOLAR_GALACTIC_RADIUS_LY,
};
use image::{Rgb, RgbImage};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::json;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;

const PC_TO_LY: f64 = 3.26156;
const Z_SOLAR: f64 = 0.0142;

#[derive(Debug, Clone)]
struct GaiaColumnMap {
    source_id: Option<usize>,
    x_ly: Option<usize>,
    y_ly: Option<usize>,
    z_ly: Option<usize>,
    x_pc: Option<usize>,
    y_pc: Option<usize>,
    z_pc: Option<usize>,
    ra_deg: Option<usize>,
    dec_deg: Option<usize>,
    parallax_mas: Option<usize>,
    distance_pc: Option<usize>,
    age_gyr: Option<usize>,
    age_years: Option<usize>,
    log_age_years: Option<usize>,
    metallicity_z: Option<usize>,
    metallicity_dex: Option<usize>,
    mass_solar: Option<usize>,
    teff_k: Option<usize>,
}

#[derive(Debug, Clone, Copy)]
struct ProcessedStar {
    source_id: u64,
    seed: GalacticLifeSeed,
    ra_deg: f64,
    dec_deg: f64,
    distance_ly: f64,
    stage: CivilizationStage,
    entropy_multiplier: f64,
    habitable: bool,
    habitability_score: f64,
    local_n_times_p: f64,
    signal: bool,
}

#[derive(Debug, Clone, Copy)]
struct RenderPoint {
    x_ly: f64,
    y_ly: f64,
    z_ly: f64,
    stage: CivilizationStage,
    habitable: bool,
    signal: bool,
}

#[derive(Debug, Clone, Copy)]
struct EpochAgg {
    delta_gyr: f64,
    habitable: usize,
    signal: usize,
    type_i: usize,
    type_ii: usize,
    type_iii: usize,
    entropy_sum: f64,
}

#[derive(Debug, Clone, Copy)]
struct SignalTarget {
    source_id: u64,
    distance_from_sun_ly: f64,
    ra_deg: f64,
    dec_deg: f64,
    x_ly: f64,
    y_ly: f64,
    z_ly: f64,
    stage: CivilizationStage,
    habitability_score: f64,
    entropy_multiplier: f64,
    local_n_times_p: f64,
    survival_fraction: f64,
    strict_pass_fraction: f64,
    energy_pass_fraction: f64,
    conflict_pass_fraction: f64,
    self_destruction_pass_fraction: f64,
    environment_pass_fraction: f64,
    stellar_stability_likelihood: f64,
    orbital_architecture_likelihood: f64,
    metallicity_band_likelihood: f64,
    galactic_environment_likelihood: f64,
    transition_likelihood: f64,
    self_destruction_likelihood: f64,
}

fn normalize_header(s: &str) -> String {
    s.trim()
        .trim_matches('"')
        .trim_matches('\'')
        .to_ascii_lowercase()
}

fn parse_f64_at(fields: &[&str], idx: Option<usize>) -> Option<f64> {
    let i = idx?;
    let raw = fields.get(i)?.trim().trim_matches('"');
    if raw.is_empty() {
        return None;
    }
    raw.parse::<f64>().ok().filter(|v| v.is_finite())
}

fn parse_u64_at(fields: &[&str], idx: Option<usize>) -> Option<u64> {
    let i = idx?;
    let raw = fields.get(i)?.trim().trim_matches('"');
    if raw.is_empty() {
        return None;
    }
    raw.parse::<u64>().ok()
}

fn header_idx(headers: &[String], keys: &[&str]) -> Option<usize> {
    for k in keys {
        let kk = normalize_header(k);
        if let Some(i) = headers.iter().position(|h| h == &kk) {
            return Some(i);
        }
    }
    None
}

fn build_column_map(headers: &[String]) -> GaiaColumnMap {
    GaiaColumnMap {
        source_id: header_idx(headers, &["source_id", "dr3_source_id", "id"]),
        x_ly: header_idx(headers, &["x_ly", "gal_x_ly"]),
        y_ly: header_idx(headers, &["y_ly", "gal_y_ly"]),
        z_ly: header_idx(headers, &["z_ly", "gal_z_ly"]),
        x_pc: header_idx(headers, &["x_pc", "gal_x_pc"]),
        y_pc: header_idx(headers, &["y_pc", "gal_y_pc"]),
        z_pc: header_idx(headers, &["z_pc", "gal_z_pc"]),
        ra_deg: header_idx(headers, &["ra", "ra_deg", "ra_icrs"]),
        dec_deg: header_idx(headers, &["dec", "dec_deg", "de_icrs"]),
        parallax_mas: header_idx(headers, &["parallax", "parallax_mas", "plx"]),
        distance_pc: header_idx(
            headers,
            &["distance_pc", "dist_pc", "r_est", "distance_gspphot"],
        ),
        age_gyr: header_idx(headers, &["age_gyr", "age_flame_gyr", "age_flame", "age"]),
        age_years: header_idx(headers, &["age_years", "age_yr", "age_flame_years"]),
        log_age_years: header_idx(headers, &["log_age", "log_age_yr", "log10_age_years"]),
        metallicity_z: header_idx(headers, &["metallicity", "z", "metallicity_z"]),
        metallicity_dex: header_idx(headers, &["mh_gspphot", "feh", "mh", "metallicity_dex"]),
        mass_solar: header_idx(headers, &["mass_solar", "mass_flame", "mass"]),
        teff_k: header_idx(headers, &["teff_gspphot", "teff", "effective_temperature"]),
    }
}

fn estimate_mass_solar(fields: &[&str], cols: &GaiaColumnMap) -> Option<f64> {
    if let Some(m) = parse_f64_at(fields, cols.mass_solar).filter(|v| *v > 0.0) {
        return Some(m.clamp(0.08, 60.0));
    }
    if let Some(teff) = parse_f64_at(fields, cols.teff_k).filter(|v| *v > 1500.0) {
        let m = (teff / 5772.0).powf(1.7);
        return Some(m.clamp(0.08, 3.5));
    }
    None
}

fn estimate_age_gyr(fields: &[&str], cols: &GaiaColumnMap) -> Option<f64> {
    if let Some(v) = parse_f64_at(fields, cols.age_gyr).filter(|v| *v > 0.0) {
        return Some(if v > 200.0 { v / 1.0e3 } else { v });
    }
    if let Some(v) = parse_f64_at(fields, cols.age_years).filter(|v| *v > 0.0) {
        return Some(v / 1.0e9);
    }
    if let Some(v) = parse_f64_at(fields, cols.log_age_years).filter(|v| *v > 0.0) {
        return Some(10.0_f64.powf(v) / 1.0e9);
    }
    None
}

fn estimate_metallicity_z(fields: &[&str], cols: &GaiaColumnMap) -> Option<f64> {
    if let Some(z) = parse_f64_at(fields, cols.metallicity_z).filter(|v| *v > 0.0) {
        return Some(z.clamp(1.0e-4, 0.2));
    }
    if let Some(dex) = parse_f64_at(fields, cols.metallicity_dex) {
        return Some((Z_SOLAR * 10.0_f64.powf(dex)).clamp(1.0e-4, 0.2));
    }
    None
}

fn ra_dec_to_galactic_xyz_ly(ra_deg: f64, dec_deg: f64, distance_ly: f64) -> (f64, f64, f64) {
    let ra = ra_deg.to_radians();
    let dec = dec_deg.to_radians();
    let x_eq = dec.cos() * ra.cos();
    let y_eq = dec.cos() * ra.sin();
    let z_eq = dec.sin();

    // ICRS -> Galactic rotation (J2000).
    let x_gal = -0.0548755604 * x_eq - 0.8734370902 * y_eq - 0.4838350155 * z_eq;
    let y_gal = 0.4941094279 * x_eq - 0.4448296300 * y_eq + 0.7469822445 * z_eq;
    let z_gal = -0.8676661490 * x_eq - 0.1980763734 * y_eq + 0.4559837762 * z_eq;

    (distance_ly * x_gal, distance_ly * y_gal, distance_ly * z_gal)
}

fn estimate_position(fields: &[&str], cols: &GaiaColumnMap) -> Option<(f64, f64, f64, f64, f64, f64)> {
    // Returns (x_ly, y_ly, z_ly, ra_deg, dec_deg, distance_ly)
    if let (Some(x), Some(y), Some(z)) = (
        parse_f64_at(fields, cols.x_ly),
        parse_f64_at(fields, cols.y_ly),
        parse_f64_at(fields, cols.z_ly),
    ) {
        let d = (x * x + y * y + z * z).sqrt();
        return Some((x, y, z, f64::NAN, f64::NAN, d));
    }
    if let (Some(x), Some(y), Some(z)) = (
        parse_f64_at(fields, cols.x_pc),
        parse_f64_at(fields, cols.y_pc),
        parse_f64_at(fields, cols.z_pc),
    ) {
        let x_ly = x * PC_TO_LY;
        let y_ly = y * PC_TO_LY;
        let z_ly = z * PC_TO_LY;
        let d = (x_ly * x_ly + y_ly * y_ly + z_ly * z_ly).sqrt();
        return Some((x_ly, y_ly, z_ly, f64::NAN, f64::NAN, d));
    }

    let ra = parse_f64_at(fields, cols.ra_deg)?;
    let dec = parse_f64_at(fields, cols.dec_deg)?;
    let distance_pc = if let Some(p) = parse_f64_at(fields, cols.parallax_mas).filter(|v| *v > 0.0)
    {
        1000.0 / p
    } else {
        parse_f64_at(fields, cols.distance_pc).filter(|v| *v > 0.0)?
    };
    let d_ly = distance_pc * PC_TO_LY;
    let (x, y, z) = ra_dec_to_galactic_xyz_ly(ra, dec, d_ly);
    Some((x, y, z, ra, dec, d_ly))
}

fn stage_color(stage: CivilizationStage) -> Rgb<u8> {
    match stage {
        CivilizationStage::BareRock => Rgb([80, 80, 90]),
        CivilizationStage::PrebioticChemistry => Rgb([70, 120, 255]),
        CivilizationStage::AutocatalyticLife => Rgb([80, 220, 120]),
        CivilizationStage::PhotosyntheticBiosphere => Rgb([175, 235, 80]),
        CivilizationStage::MulticellularEcosystem => Rgb([255, 175, 70]),
        CivilizationStage::TechnologicalIntelligence => Rgb([70, 240, 255]),
        CivilizationStage::KardashevTypeI => Rgb([175, 120, 255]),
        CivilizationStage::KardashevTypeII => Rgb([255, 95, 150]),
        CivilizationStage::KardashevTypeIII => Rgb([255, 255, 255]),
    }
}

fn scale_rgb(c: Rgb<u8>, k: f64) -> Rgb<u8> {
    let kk = k.clamp(0.0, 3.0);
    Rgb([
        (c[0] as f64 * kk).clamp(0.0, 255.0) as u8,
        (c[1] as f64 * kk).clamp(0.0, 255.0) as u8,
        (c[2] as f64 * kk).clamp(0.0, 255.0) as u8,
    ])
}

fn put_px(img: &mut RgbImage, x: i32, y: i32, c: Rgb<u8>) {
    if x >= 0 && y >= 0 && (x as u32) < img.width() && (y as u32) < img.height() {
        img.put_pixel(x as u32, y as u32, c);
    }
}

fn draw_dot(img: &mut RgbImage, x: i32, y: i32, r: i32, c: Rgb<u8>) {
    for dy in -r..=r {
        for dx in -r..=r {
            if dx * dx + dy * dy <= r * r {
                put_px(img, x + dx, y + dy, c);
            }
        }
    }
}

fn draw_line(img: &mut RgbImage, x0: i32, y0: i32, x1: i32, y1: i32, c: Rgb<u8>) {
    let mut x = x0;
    let mut y = y0;
    let dx = (x1 - x0).abs();
    let sx = if x0 < x1 { 1 } else { -1 };
    let dy = -(y1 - y0).abs();
    let sy = if y0 < y1 { 1 } else { -1 };
    let mut err = dx + dy;
    loop {
        put_px(img, x, y, c);
        if x == x1 && y == y1 {
            break;
        }
        let e2 = 2 * err;
        if e2 >= dy {
            err += dy;
            x += sx;
        }
        if e2 <= dx {
            err += dx;
            y += sy;
        }
    }
}

fn project(x: f64, y: f64, z: f64, w: u32, h: u32) -> (i32, i32, f64) {
    let yaw = -35.0_f64.to_radians();
    let pitch = 22.0_f64.to_radians();
    let cy = yaw.cos();
    let sy = yaw.sin();
    let cp = pitch.cos();
    let sp = pitch.sin();

    let x1 = x * cy - y * sy;
    let y1 = x * sy + y * cy;
    let z1 = z;
    let y2 = y1 * cp - z1 * sp;
    let z2 = y1 * sp + z1 * cp;

    let d = 120_000.0;
    let persp = d / (d + y2 + 1.0e-6);
    let s = 0.0108;
    let sx = w as f64 * 0.5 + x1 * s * persp;
    let sy = h as f64 * 0.58 - z2 * s * persp;
    (sx.round() as i32, sy.round() as i32, y2)
}

fn render_map(points: &[RenderPoint], out_path: &PathBuf) -> Result<(), String> {
    let w = 1600u32;
    let h = 1000u32;
    let mut img = RgbImage::new(w, h);
    for px in img.pixels_mut() {
        *px = Rgb([0, 0, 0]);
    }

    let mut proj = Vec::with_capacity(points.len());
    for p in points {
        let (sx, sy, depth) = project(p.x_ly, p.y_ly, p.z_ly, w, h);
        proj.push((sx, sy, depth, *p));
    }
    proj.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal));

    for (sx, sy, depth, p) in proj {
        let base = stage_color(p.stage);
        let depth_mod = (0.72 + 0.35 * (1.0 - ((depth + 45_000.0) / 90_000.0).clamp(0.0, 1.0)))
            .clamp(0.45, 1.2);
        let hab = if p.habitable { 1.0 } else { 0.35 };
        let c = scale_rgb(base, depth_mod * hab);
        draw_dot(&mut img, sx, sy, if p.signal { 2 } else { 1 }, c);
    }

    let (ox, oy, _) = project(0.0, 0.0, 0.0, w, h);
    let (xx, xy, _) = project(40_000.0, 0.0, 0.0, w, h);
    let (yx, yy, _) = project(0.0, 40_000.0, 0.0, w, h);
    let (zx, zy, _) = project(0.0, 0.0, 1_600.0, w, h);
    draw_line(&mut img, ox, oy, xx, xy, Rgb([110, 110, 110]));
    draw_line(&mut img, ox, oy, yx, yy, Rgb([110, 110, 110]));
    draw_line(&mut img, ox, oy, zx, zy, Rgb([110, 110, 110]));
    draw_dot(&mut img, ox, oy, 3, Rgb([220, 220, 220]));

    let (sx, sy, _) = project(SOLAR_GALACTIC_RADIUS_LY, 0.0, 30.0, w, h);
    draw_dot(&mut img, sx, sy, 4, Rgb([0, 170, 255]));
    draw_dot(&mut img, sx, sy, 2, Rgb([255, 255, 255]));

    img.save(out_path)
        .map_err(|e| format!("save {}: {e}", out_path.display()))
}

fn reservoir_push(rng: &mut StdRng, seen: usize, cap: usize, reservoir: &mut Vec<RenderPoint>, p: RenderPoint) {
    if reservoir.len() < cap {
        reservoir.push(p);
        return;
    }
    let j = rng.gen_range(0..seen);
    if j < cap {
        reservoir[j] = p;
    }
}

fn push_nearest(nearest: &mut Vec<SignalTarget>, target: SignalTarget, k: usize) {
    if nearest.len() < k {
        nearest.push(target);
        return;
    }
    let mut worst_i = 0usize;
    let mut worst_d = nearest[0].distance_from_sun_ly;
    for (i, t) in nearest.iter().enumerate().skip(1) {
        if t.distance_from_sun_ly > worst_d {
            worst_d = t.distance_from_sun_ly;
            worst_i = i;
        }
    }
    if target.distance_from_sun_ly < worst_d {
        nearest[worst_i] = target;
    }
}

fn main() {
    let input = std::env::var("GUTOE_GAIA_DR3_CSV")
        .expect("set GUTOE_GAIA_DR3_CSV to Gaia DR3 CSV path");
    let out_dir = std::env::var("GUTOE_GAIA_LIFE_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/gaia_life_map".to_string());
    let max_rows = std::env::var("GUTOE_GAIA_MAX_ROWS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok());
    let reservoir_cap = std::env::var("GUTOE_GAIA_RENDER_SAMPLE")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(200_000);
    let nearest_k = std::env::var("GUTOE_GAIA_NEAREST_K")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(256);
    let gf_trials = std::env::var("GUTOE_GF_TRIALS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(512);
    let gf_seed_salt = std::env::var("GUTOE_GF_SEED_SALT")
        .ok()
        .and_then(|s| s.parse::<u64>().ok())
        .unwrap_or(9_813_311);
    let write_all_signals = std::env::var("GUTOE_GAIA_WRITE_SIGNAL_CATALOG")
        .ok()
        .map(|s| s == "1" || s.eq_ignore_ascii_case("true"))
        .unwrap_or(false);

    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);
    let txt_path = out.join("gaia_life_report.txt");
    let json_path = out.join("gaia_life_report.json");
    let png_path = out.join("gaia_life_map.png");
    let targets_csv_path = out.join("gaia_signal_targets.csv");
    let all_signals_path = out.join("gaia_signal_catalog.csv");

    let file = File::open(&input).expect("open Gaia csv");
    let mut reader = BufReader::new(file);
    let mut header_line = String::new();
    reader.read_line(&mut header_line).expect("read header");
    let headers = header_line
        .trim_end()
        .split(',')
        .map(normalize_header)
        .collect::<Vec<_>>();
    let cols = build_column_map(&headers);

    let (universe_age_gyr, thresholds, multipliers) = derive_thresholds_and_multipliers();
    let abiogenesis = evaluate_abiogenesis_gate(AbiogenesisWindows::default(), 298.15);
    let baseline_n_times_p = abiogenesis.closure.n_times_p;
    let gf_windows = GreatFilterWindows {
        trials_per_civilization: gf_trials,
        ..GreatFilterWindows::default()
    };

    let mut rng = StdRng::seed_from_u64(2_026_033_1);
    let mut render_reservoir = Vec::<RenderPoint>::new();
    let mut nearest_targets = Vec::<SignalTarget>::new();
    let mut stage_counts = [0usize; 9];
    let mut epochs = [
        EpochAgg {
            delta_gyr: 0.0,
            habitable: 0,
            signal: 0,
            type_i: 0,
            type_ii: 0,
            type_iii: 0,
            entropy_sum: 0.0,
        },
        EpochAgg {
            delta_gyr: 0.5,
            habitable: 0,
            signal: 0,
            type_i: 0,
            type_ii: 0,
            type_iii: 0,
            entropy_sum: 0.0,
        },
        EpochAgg {
            delta_gyr: 1.0,
            habitable: 0,
            signal: 0,
            type_i: 0,
            type_ii: 0,
            type_iii: 0,
            entropy_sum: 0.0,
        },
        EpochAgg {
            delta_gyr: 2.0,
            habitable: 0,
            signal: 0,
            type_i: 0,
            type_ii: 0,
            type_iii: 0,
            entropy_sum: 0.0,
        },
    ];

    let mut all_signals_writer = if write_all_signals {
        let mut w = BufWriter::new(File::create(&all_signals_path).expect("create all signal csv"));
        writeln!(
            w,
            "source_id,ra_deg,dec_deg,distance_ly,x_ly,y_ly,z_ly,stage,habitability_score,entropy_multiplier,local_n_times_p,survival_fraction,strict_pass_fraction,energy_pass_fraction,conflict_pass_fraction,self_destruction_pass_fraction,environment_pass_fraction,stellar_stability_likelihood,orbital_architecture_likelihood,metallicity_band_likelihood,galactic_environment_likelihood,transition_likelihood,self_destruction_likelihood"
        )
        .expect("write all signal header");
        Some(w)
    } else {
        None
    };

    let mut seen = 0usize;
    let mut used = 0usize;
    let mut skipped_missing = 0usize;
    let mut skipped_bad_coords = 0usize;
    let mut skipped_bad_physics = 0usize;
    let mut gf_signal_count = 0usize;
    let mut gf_survival_sum = 0.0;
    let mut gf_strict_sum = 0.0;
    let mut gf_energy_sum = 0.0;
    let mut gf_conflict_sum = 0.0;
    let mut gf_self_sum = 0.0;
    let mut gf_env_sum = 0.0;
    let mut gf_stellar_sum = 0.0;
    let mut gf_orbital_sum = 0.0;
    let mut gf_metal_band_sum = 0.0;
    let mut gf_galactic_sum = 0.0;
    let mut gf_transition_sum = 0.0;
    let mut gf_self_likelihood_sum = 0.0;

    let mut line = String::new();
    loop {
        line.clear();
        let nread = reader.read_line(&mut line).expect("read line");
        if nread == 0 {
            break;
        }
        seen += 1;
        if let Some(maxr) = max_rows {
            if seen > maxr {
                break;
            }
        }
        if line.trim().is_empty() {
            skipped_missing += 1;
            continue;
        }
        let fields = line.trim_end().split(',').collect::<Vec<_>>();

        let (x_ly, y_ly, z_ly, mut ra_deg, mut dec_deg, distance_ly) =
            if let Some(v) = estimate_position(&fields, &cols) {
                v
            } else {
                skipped_bad_coords += 1;
                continue;
            };
        if !(distance_ly.is_finite() && distance_ly > 0.0) {
            skipped_bad_coords += 1;
            continue;
        }

        if !ra_deg.is_finite() || !dec_deg.is_finite() {
            let r_xy = (x_ly * x_ly + y_ly * y_ly).sqrt().max(1.0e-12);
            ra_deg = y_ly.atan2(x_ly).to_degrees();
            dec_deg = (z_ly / distance_ly).asin().to_degrees();
            if !ra_deg.is_finite() {
                ra_deg = 0.0;
            }
            if !dec_deg.is_finite() {
                dec_deg = (z_ly / (r_xy.hypot(z_ly))).asin().to_degrees();
            }
        }

        let mass = if let Some(v) = estimate_mass_solar(&fields, &cols) {
            v
        } else {
            skipped_missing += 1;
            continue;
        };
        let lifetime = main_sequence_lifetime_gyr(mass);
        let age = if let Some(v) = estimate_age_gyr(&fields, &cols) {
            v.min(0.98 * lifetime).max(0.01)
        } else {
            skipped_missing += 1;
            continue;
        };
        let metallicity = if let Some(v) = estimate_metallicity_z(&fields, &cols) {
            v
        } else {
            skipped_missing += 1;
            continue;
        };

        let source_id = parse_u64_at(&fields, cols.source_id).unwrap_or(seen as u64);
        let component = infer_component_from_position(x_ly, y_ly, z_ly);
        let seed = GalacticLifeSeed {
            id: source_id,
            component,
            x_ly,
            y_ly,
            z_ly,
            galactic_radius_ly: (x_ly * x_ly + y_ly * y_ly).sqrt(),
            mass_solar: mass,
            age_gyr: age,
            metallicity,
            main_sequence_lifetime_gyr: lifetime,
        };

        let local_n_times_p = baseline_n_times_p * (metallicity / Z_SOLAR).sqrt();
        let kauffman_ok = local_n_times_p >= KAUFFMAN_CLOSURE_THRESHOLD;
        let h_score = habitability_score(seed);
        let hab = is_habitable(seed, h_score);
        if !age.is_finite() || !metallicity.is_finite() {
            skipped_bad_physics += 1;
            continue;
        }

        let mut stage = classify_stage(age, thresholds);
        if stage.rank() >= CivilizationStage::AutocatalyticLife.rank() && !kauffman_ok {
            stage = CivilizationStage::PrebioticChemistry;
        }
        let entropy = stage_entropy_multiplier(stage, multipliers);
        let signal = hab && kauffman_ok && stage.is_signal();

        let star = ProcessedStar {
            source_id,
            seed,
            ra_deg,
            dec_deg,
            distance_ly,
            stage,
            entropy_multiplier: entropy,
            habitable: hab,
            habitability_score: h_score,
            local_n_times_p,
            signal,
        };
        used += 1;
        stage_counts[stage.rank() as usize] += 1;

        for e in &mut epochs {
            let aged = star.seed.age_gyr + e.delta_gyr;
            if aged > 0.98 * star.seed.main_sequence_lifetime_gyr {
                continue;
            }
            let aged_seed = GalacticLifeSeed {
                age_gyr: aged,
                ..star.seed
            };
            let aged_score = habitability_score(aged_seed);
            let aged_hab = is_habitable(aged_seed, aged_score);
            if !aged_hab {
                continue;
            }
            let local_np = baseline_n_times_p * (star.seed.metallicity / Z_SOLAR).sqrt();
            let mut st = classify_stage(aged, thresholds);
            if st.rank() >= CivilizationStage::AutocatalyticLife.rank()
                && local_np < KAUFFMAN_CLOSURE_THRESHOLD
            {
                st = CivilizationStage::PrebioticChemistry;
            }
            let ent = stage_entropy_multiplier(st, multipliers);
            e.habitable += 1;
            e.entropy_sum += ent;
            if st.is_signal() {
                e.signal += 1;
            }
            if st.rank() >= CivilizationStage::KardashevTypeI.rank() {
                e.type_i += 1;
            }
            if st.rank() >= CivilizationStage::KardashevTypeII.rank() {
                e.type_ii += 1;
            }
            if st.rank() >= CivilizationStage::KardashevTypeIII.rank() {
                e.type_iii += 1;
            }
        }

        reservoir_push(
            &mut rng,
            used.max(1),
            reservoir_cap,
            &mut render_reservoir,
            RenderPoint {
                x_ly: star.seed.x_ly,
                y_ly: star.seed.y_ly,
                z_ly: star.seed.z_ly,
                stage: star.stage,
                habitable: star.habitable,
                signal: star.signal,
            },
        );

        if star.signal {
            let gf = evaluate_great_filter(
                GreatFilterCivilizationInput {
                    stage_rank: star.stage.rank(),
                    mass_solar: star.seed.mass_solar,
                    metallicity_z: star.seed.metallicity,
                    age_gyr: star.seed.age_gyr,
                    habitability_score: star.habitability_score,
                    entropy_multiplier: star.entropy_multiplier,
                    local_n_times_p: star.local_n_times_p,
                    galactic_radius_ly: star.seed.galactic_radius_ly,
                    galactic_z_ly: star.seed.z_ly,
                },
                gf_windows,
                star.source_id ^ gf_seed_salt,
            );
            gf_signal_count += 1;
            gf_survival_sum += gf.survival_fraction;
            gf_strict_sum += gf.strict_pass_fraction;
            gf_energy_sum += gf.energy_pass_fraction;
            gf_conflict_sum += gf.conflict_pass_fraction;
            gf_self_sum += gf.self_destruction_pass_fraction;
            gf_env_sum += gf.environment_pass_fraction;
            gf_stellar_sum += gf.stellar_stability_likelihood;
            gf_orbital_sum += gf.orbital_architecture_likelihood;
            gf_metal_band_sum += gf.metallicity_band_likelihood;
            gf_galactic_sum += gf.galactic_environment_likelihood;
            gf_transition_sum += gf.transition_likelihood;
            gf_self_likelihood_sum += gf.self_destruction_likelihood;

            if let Some(w) = &mut all_signals_writer {
                writeln!(
                    w,
                    "{},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9}",
                    star.source_id,
                    star.ra_deg,
                    star.dec_deg,
                    star.distance_ly,
                    star.seed.x_ly,
                    star.seed.y_ly,
                    star.seed.z_ly,
                    star.stage.as_str(),
                    star.habitability_score,
                    star.entropy_multiplier,
                    star.local_n_times_p,
                    gf.survival_fraction,
                    gf.strict_pass_fraction,
                    gf.energy_pass_fraction,
                    gf.conflict_pass_fraction,
                    gf.self_destruction_pass_fraction,
                    gf.environment_pass_fraction,
                    gf.stellar_stability_likelihood,
                    gf.orbital_architecture_likelihood,
                    gf.metallicity_band_likelihood,
                    gf.galactic_environment_likelihood,
                    gf.transition_likelihood,
                    gf.self_destruction_likelihood
                )
                .expect("write all signals row");
            }
            let dx = star.seed.x_ly - SOLAR_GALACTIC_RADIUS_LY;
            let dy = star.seed.y_ly;
            let dz = star.seed.z_ly;
            let d = (dx * dx + dy * dy + dz * dz).sqrt();
            push_nearest(
                &mut nearest_targets,
                SignalTarget {
                    source_id: star.source_id,
                    distance_from_sun_ly: d,
                    ra_deg: star.ra_deg,
                    dec_deg: star.dec_deg,
                    x_ly: star.seed.x_ly,
                    y_ly: star.seed.y_ly,
                    z_ly: star.seed.z_ly,
                    stage: star.stage,
                    habitability_score: star.habitability_score,
                    entropy_multiplier: star.entropy_multiplier,
                    local_n_times_p: star.local_n_times_p,
                    survival_fraction: gf.survival_fraction,
                    strict_pass_fraction: gf.strict_pass_fraction,
                    energy_pass_fraction: gf.energy_pass_fraction,
                    conflict_pass_fraction: gf.conflict_pass_fraction,
                    self_destruction_pass_fraction: gf.self_destruction_pass_fraction,
                    environment_pass_fraction: gf.environment_pass_fraction,
                    stellar_stability_likelihood: gf.stellar_stability_likelihood,
                    orbital_architecture_likelihood: gf.orbital_architecture_likelihood,
                    metallicity_band_likelihood: gf.metallicity_band_likelihood,
                    galactic_environment_likelihood: gf.galactic_environment_likelihood,
                    transition_likelihood: gf.transition_likelihood,
                    self_destruction_likelihood: gf.self_destruction_likelihood,
                },
                nearest_k,
            );
        }
    }

    if let Some(w) = &mut all_signals_writer {
        w.flush().expect("flush all signal csv");
    }

    render_map(&render_reservoir, &png_path).expect("render map");

    nearest_targets.sort_by(|a, b| {
        a.distance_from_sun_ly
            .partial_cmp(&b.distance_from_sun_ly)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let mut target_csv = BufWriter::new(File::create(&targets_csv_path).expect("create targets csv"));
    writeln!(
        target_csv,
        "source_id,distance_from_sun_ly,ra_deg,dec_deg,x_ly,y_ly,z_ly,stage,habitability_score,entropy_multiplier,local_n_times_p,survival_fraction,strict_pass_fraction,energy_pass_fraction,conflict_pass_fraction,self_destruction_pass_fraction,environment_pass_fraction,stellar_stability_likelihood,orbital_architecture_likelihood,metallicity_band_likelihood,galactic_environment_likelihood,transition_likelihood,self_destruction_likelihood"
    )
    .expect("targets header");
    for t in &nearest_targets {
        writeln!(
            target_csv,
            "{},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9}",
            t.source_id,
            t.distance_from_sun_ly,
            t.ra_deg,
            t.dec_deg,
            t.x_ly,
            t.y_ly,
            t.z_ly,
            t.stage.as_str(),
            t.habitability_score,
            t.entropy_multiplier,
            t.local_n_times_p,
            t.survival_fraction,
            t.strict_pass_fraction,
            t.energy_pass_fraction,
            t.conflict_pass_fraction,
            t.self_destruction_pass_fraction,
            t.environment_pass_fraction,
            t.stellar_stability_likelihood,
            t.orbital_architecture_likelihood,
            t.metallicity_band_likelihood,
            t.galactic_environment_likelihood,
            t.transition_likelihood,
            t.self_destruction_likelihood
        )
        .expect("targets row");
    }
    target_csv.flush().expect("flush targets");

    let mut txt = BufWriter::new(File::create(&txt_path).expect("create txt"));
    let signal_count_present = epochs[0].signal;
    let habitable_count_present = epochs[0].habitable;
    let signal_fraction = if used > 0 {
        signal_count_present as f64 / used as f64
    } else {
        0.0
    };
    let predicted_signal_count_milky_way_present = signal_fraction * MILKY_WAY_STELLAR_COUNT_ESTIMATE;
    let gf_survival_fraction_conditional = if gf_signal_count > 0 {
        gf_survival_sum / gf_signal_count as f64
    } else {
        0.0
    };
    let gf_strict_pass_mean = if gf_signal_count > 0 {
        gf_strict_sum / gf_signal_count as f64
    } else {
        0.0
    };
    let gf_energy_pass_mean = if gf_signal_count > 0 {
        gf_energy_sum / gf_signal_count as f64
    } else {
        0.0
    };
    let gf_conflict_pass_mean = if gf_signal_count > 0 {
        gf_conflict_sum / gf_signal_count as f64
    } else {
        0.0
    };
    let gf_self_pass_mean = if gf_signal_count > 0 {
        gf_self_sum / gf_signal_count as f64
    } else {
        0.0
    };
    let gf_environment_pass_mean = if gf_signal_count > 0 {
        gf_env_sum / gf_signal_count as f64
    } else {
        0.0
    };
    let gf_stellar_likelihood_mean = if gf_signal_count > 0 {
        gf_stellar_sum / gf_signal_count as f64
    } else {
        0.0
    };
    let gf_orbital_likelihood_mean = if gf_signal_count > 0 {
        gf_orbital_sum / gf_signal_count as f64
    } else {
        0.0
    };
    let gf_metal_band_likelihood_mean = if gf_signal_count > 0 {
        gf_metal_band_sum / gf_signal_count as f64
    } else {
        0.0
    };
    let gf_galactic_likelihood_mean = if gf_signal_count > 0 {
        gf_galactic_sum / gf_signal_count as f64
    } else {
        0.0
    };
    let gf_transition_likelihood_mean = if gf_signal_count > 0 {
        gf_transition_sum / gf_signal_count as f64
    } else {
        0.0
    };
    let gf_self_likelihood_mean = if gf_signal_count > 0 {
        gf_self_likelihood_sum / gf_signal_count as f64
    } else {
        0.0
    };
    let surviving_signal_count_present_expected =
        signal_count_present as f64 * gf_survival_fraction_conditional;
    let surviving_signal_fraction_present = if used > 0 {
        surviving_signal_count_present_expected / used as f64
    } else {
        0.0
    };
    let surviving_signal_count_milky_way_expected =
        predicted_signal_count_milky_way_present * gf_survival_fraction_conditional;
    writeln!(txt, "[gaia_dr3_life_map]").expect("write");
    writeln!(txt, "input_csv = {}", input).expect("write");
    writeln!(txt, "rows_seen = {}", seen).expect("write");
    writeln!(txt, "rows_used = {}", used).expect("write");
    writeln!(txt, "skipped_missing = {}", skipped_missing).expect("write");
    writeln!(txt, "skipped_bad_coords = {}", skipped_bad_coords).expect("write");
    writeln!(txt, "skipped_bad_physics = {}", skipped_bad_physics).expect("write");
    writeln!(txt, "universe_age_gyr = {:.9}", universe_age_gyr).expect("write");
    writeln!(
        txt,
        "kauffman_baseline_n_times_p = {:.9}",
        baseline_n_times_p
    )
    .expect("write");
    writeln!(txt, "habitable_count_present = {}", habitable_count_present).expect("write");
    writeln!(txt, "signal_count_present = {}", signal_count_present).expect("write");
    writeln!(txt, "signal_fraction_present = {:.12}", signal_fraction).expect("write");
    writeln!(
        txt,
        "predicted_signal_count_milky_way_present = {:.6e}",
        predicted_signal_count_milky_way_present
    )
    .expect("write");
    writeln!(txt, "great_filter_trials_per_civilization = {}", gf_trials).expect("write");
    writeln!(
        txt,
        "great_filter_survival_fraction_conditional = {:.9}",
        gf_survival_fraction_conditional
    )
    .expect("write");
    writeln!(
        txt,
        "surviving_signal_count_present_expected = {:.6}",
        surviving_signal_count_present_expected
    )
    .expect("write");
    writeln!(
        txt,
        "surviving_signal_fraction_present = {:.12}",
        surviving_signal_fraction_present
    )
    .expect("write");
    writeln!(
        txt,
        "surviving_signal_count_milky_way_expected = {:.6e}",
        surviving_signal_count_milky_way_expected
    )
    .expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[stage_counts_present]").expect("write");
    for (i, n) in stage_counts.iter().enumerate() {
        let st = match i {
            0 => CivilizationStage::BareRock,
            1 => CivilizationStage::PrebioticChemistry,
            2 => CivilizationStage::AutocatalyticLife,
            3 => CivilizationStage::PhotosyntheticBiosphere,
            4 => CivilizationStage::MulticellularEcosystem,
            5 => CivilizationStage::TechnologicalIntelligence,
            6 => CivilizationStage::KardashevTypeI,
            7 => CivilizationStage::KardashevTypeII,
            _ => CivilizationStage::KardashevTypeIII,
        };
        writeln!(txt, "{} = {}", st.as_str(), n).expect("write");
    }
    writeln!(txt).expect("write");
    writeln!(txt, "[forward_epochs]").expect("write");
    for e in &epochs {
        writeln!(
            txt,
            "delta={:.3}Gyr habitable={} signal={} typeI={} typeII={} typeIII={} mean_entropy={:.9}",
            e.delta_gyr,
            e.habitable,
            e.signal,
            e.type_i,
            e.type_ii,
            e.type_iii,
            if e.habitable > 0 {
                e.entropy_sum / e.habitable as f64
            } else {
                0.0
            }
        )
        .expect("write");
    }
    writeln!(txt).expect("write");
    writeln!(txt, "[great_filter_gate_means]").expect("write");
    writeln!(txt, "strict_pass_mean = {:.9}", gf_strict_pass_mean).expect("write");
    writeln!(txt, "energy_pass_mean = {:.9}", gf_energy_pass_mean).expect("write");
    writeln!(txt, "conflict_pass_mean = {:.9}", gf_conflict_pass_mean).expect("write");
    writeln!(txt, "self_destruction_pass_mean = {:.9}", gf_self_pass_mean).expect("write");
    writeln!(txt, "environment_pass_mean = {:.9}", gf_environment_pass_mean).expect("write");
    writeln!(
        txt,
        "stellar_stability_likelihood_mean = {:.9}",
        gf_stellar_likelihood_mean
    )
    .expect("write");
    writeln!(
        txt,
        "orbital_architecture_likelihood_mean = {:.9}",
        gf_orbital_likelihood_mean
    )
    .expect("write");
    writeln!(
        txt,
        "metallicity_band_likelihood_mean = {:.9}",
        gf_metal_band_likelihood_mean
    )
    .expect("write");
    writeln!(
        txt,
        "galactic_environment_likelihood_mean = {:.9}",
        gf_galactic_likelihood_mean
    )
    .expect("write");
    writeln!(
        txt,
        "transition_likelihood_mean = {:.9}",
        gf_transition_likelihood_mean
    )
    .expect("write");
    writeln!(
        txt,
        "self_destruction_likelihood_mean = {:.9}",
        gf_self_likelihood_mean
    )
    .expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[nearest_signal_targets]").expect("write");
    for t in &nearest_targets {
        writeln!(
            txt,
            "source_id={} d_ly={:.6} ra={:.9} dec={:.9} stage={} np={:.9} survive={:.6}",
            t.source_id,
            t.distance_from_sun_ly,
            t.ra_deg,
            t.dec_deg,
            t.stage.as_str(),
            t.local_n_times_p,
            t.survival_fraction
        )
        .expect("write");
    }
    txt.flush().expect("flush txt");

    let json_payload = json!({
        "summary": {
            "input_csv": input,
            "rows_seen": seen,
            "rows_used": used,
            "skipped_missing": skipped_missing,
            "skipped_bad_coords": skipped_bad_coords,
            "skipped_bad_physics": skipped_bad_physics,
            "universe_age_gyr": universe_age_gyr,
            "kauffman_baseline_n_times_p": baseline_n_times_p,
            "habitable_count_present": habitable_count_present,
            "signal_count_present": signal_count_present,
            "signal_fraction_present": signal_fraction,
            "predicted_signal_count_milky_way_present": predicted_signal_count_milky_way_present,
            "surviving_signal_count_present_expected": surviving_signal_count_present_expected,
            "surviving_signal_fraction_present": surviving_signal_fraction_present,
            "surviving_signal_count_milky_way_expected": surviving_signal_count_milky_way_expected
        },
        "thresholds_gyr": {
            "prebiotic": thresholds.prebiotic_age_gyr,
            "autocatalytic": thresholds.autocatalytic_age_gyr,
            "photosynthetic": thresholds.photosynthetic_age_gyr,
            "multicellular": thresholds.multicellular_age_gyr,
            "intelligence": thresholds.intelligence_age_gyr,
            "kardashev_i": thresholds.kardashev_i_age_gyr,
            "kardashev_ii": thresholds.kardashev_ii_age_gyr,
            "kardashev_iii": thresholds.kardashev_iii_age_gyr
        },
        "stage_counts_present": [
            {"stage":"bare_rock","count":stage_counts[0]},
            {"stage":"prebiotic_chemistry","count":stage_counts[1]},
            {"stage":"autocatalytic_life","count":stage_counts[2]},
            {"stage":"photosynthetic_biosphere","count":stage_counts[3]},
            {"stage":"multicellular_ecosystem","count":stage_counts[4]},
            {"stage":"technological_intelligence","count":stage_counts[5]},
            {"stage":"kardashev_type_i","count":stage_counts[6]},
            {"stage":"kardashev_type_ii","count":stage_counts[7]},
            {"stage":"kardashev_type_iii","count":stage_counts[8]}
        ],
        "forward_epochs": epochs.iter().map(|e| json!({
            "delta_gyr": e.delta_gyr,
            "habitable": e.habitable,
            "signal": e.signal,
            "type_i": e.type_i,
            "type_ii": e.type_ii,
            "type_iii": e.type_iii,
            "mean_entropy_multiplier_habitable": if e.habitable > 0 { e.entropy_sum / e.habitable as f64 } else { 0.0 }
        })).collect::<Vec<_>>(),
        "great_filter": {
            "trials_per_civilization": gf_trials,
            "signal_count_modeled": gf_signal_count,
            "survival_fraction_conditional": gf_survival_fraction_conditional,
            "strict_pass_mean": gf_strict_pass_mean,
            "energy_pass_mean": gf_energy_pass_mean,
            "conflict_pass_mean": gf_conflict_pass_mean,
            "self_destruction_pass_mean": gf_self_pass_mean,
            "environment_pass_mean": gf_environment_pass_mean,
            "stellar_stability_likelihood_mean": gf_stellar_likelihood_mean,
            "orbital_architecture_likelihood_mean": gf_orbital_likelihood_mean,
            "metallicity_band_likelihood_mean": gf_metal_band_likelihood_mean,
            "galactic_environment_likelihood_mean": gf_galactic_likelihood_mean,
            "transition_likelihood_mean": gf_transition_likelihood_mean,
            "self_destruction_likelihood_mean": gf_self_likelihood_mean
        },
        "nearest_signal_targets": nearest_targets.iter().map(|t| json!({
            "source_id": t.source_id,
            "distance_from_sun_ly": t.distance_from_sun_ly,
            "ra_deg": t.ra_deg,
            "dec_deg": t.dec_deg,
            "x_ly": t.x_ly,
            "y_ly": t.y_ly,
            "z_ly": t.z_ly,
            "stage": t.stage.as_str(),
            "habitability_score": t.habitability_score,
            "entropy_multiplier": t.entropy_multiplier,
            "local_n_times_p": t.local_n_times_p,
            "survival_fraction": t.survival_fraction,
            "strict_pass_fraction": t.strict_pass_fraction,
            "energy_pass_fraction": t.energy_pass_fraction,
            "conflict_pass_fraction": t.conflict_pass_fraction,
            "self_destruction_pass_fraction": t.self_destruction_pass_fraction,
            "environment_pass_fraction": t.environment_pass_fraction,
            "stellar_stability_likelihood": t.stellar_stability_likelihood,
            "orbital_architecture_likelihood": t.orbital_architecture_likelihood,
            "metallicity_band_likelihood": t.metallicity_band_likelihood,
            "galactic_environment_likelihood": t.galactic_environment_likelihood,
            "transition_likelihood": t.transition_likelihood,
            "self_destruction_likelihood": t.self_destruction_likelihood
        })).collect::<Vec<_>>(),
        "artifacts": {
            "map_png": png_path.display().to_string(),
            "targets_csv": targets_csv_path.display().to_string(),
            "report_txt": txt_path.display().to_string(),
            "all_signals_csv": if write_all_signals { Some(all_signals_path.display().to_string()) } else { None }
        }
    });
    let mut jf = BufWriter::new(File::create(&json_path).expect("create json"));
    writeln!(
        jf,
        "{}",
        serde_json::to_string_pretty(&json_payload).expect("serialize")
    )
    .expect("write json");
    jf.flush().expect("flush json");

    println!("wrote {}", png_path.display());
    println!("wrote {}", targets_csv_path.display());
    if write_all_signals {
        println!("wrote {}", all_signals_path.display());
    }
    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "gaia_dr3_life_map: seen={} used={} habitable_now={} signal_now={} survive_cond={:.4} survive_MW={:.3e} nearest_targets={}",
        seen,
        used,
        habitable_count_present,
        signal_count_present,
        gf_survival_fraction_conditional,
        surviving_signal_count_milky_way_expected,
        nearest_targets.len()
    );
}
