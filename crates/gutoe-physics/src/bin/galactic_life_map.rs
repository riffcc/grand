//! Galactic life map report from entropy-progression stage classification.

use gutoe_physics::{
    evaluate_galactic_life_map, CivilizationStage, GalacticComponent, GalacticLifeMapConfig,
    SOLAR_GALACTIC_RADIUS_LY,
};
use image::{Rgb, RgbImage};
use serde_json::json;
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

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

fn component_luma(component: GalacticComponent) -> f64 {
    match component {
        GalacticComponent::Disk => 1.0,
        GalacticComponent::Bulge => 0.88,
        GalacticComponent::Halo => 0.75,
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
    // Camera yaw/pitch for a 3D-looking disk projection.
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

fn render_png(
    score: &gutoe_physics::GalacticLifeMapScorecard,
    out_path: &PathBuf,
) -> Result<(), String> {
    let w = 1600u32;
    let h = 1000u32;
    let mut img = RgbImage::new(w, h);
    for px in img.pixels_mut() {
        *px = Rgb([0, 0, 0]);
    }

    // Project + depth-sort points (far -> near).
    let mut proj = Vec::with_capacity(score.points.len());
    for p in &score.points {
        let (sx, sy, depth) = project(p.seed.x_ly, p.seed.y_ly, p.seed.z_ly, w, h);
        proj.push((sx, sy, depth, *p));
    }
    proj.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal));

    for (sx, sy, depth, p) in proj {
        let base = stage_color(p.stage);
        let depth_mod = (0.72 + 0.35 * (1.0 - ((depth + 45_000.0) / 90_000.0).clamp(0.0, 1.0)))
            .clamp(0.45, 1.2);
        let comp = component_luma(p.seed.component);
        let hab = if p.habitable { 1.0 } else { 0.35 };
        let c = scale_rgb(base, depth_mod * comp * hab);
        let r = if p.signal { 2 } else { 1 };
        draw_dot(&mut img, sx, sy, r, c);
    }

    // Axes from galactic center.
    let (ox, oy, _) = project(0.0, 0.0, 0.0, w, h);
    let (xx, xy, _) = project(40_000.0, 0.0, 0.0, w, h);
    let (yx, yy, _) = project(0.0, 40_000.0, 0.0, w, h);
    let (zx, zy, _) = project(0.0, 0.0, 1_600.0, w, h);
    draw_line(&mut img, ox, oy, xx, xy, Rgb([110, 110, 110]));
    draw_line(&mut img, ox, oy, yx, yy, Rgb([110, 110, 110]));
    draw_line(&mut img, ox, oy, zx, zy, Rgb([110, 110, 110]));
    draw_dot(&mut img, ox, oy, 3, Rgb([220, 220, 220]));

    // Solar reference marker near local arm.
    let (sx, sy, _) = project(SOLAR_GALACTIC_RADIUS_LY, 0.0, 30.0, w, h);
    draw_dot(&mut img, sx, sy, 4, Rgb([0, 170, 255]));
    draw_dot(&mut img, sx, sy, 2, Rgb([255, 255, 255]));

    img.save(out_path)
        .map_err(|e| format!("save {}: {e}", out_path.display()))
}

fn main() {
    let n = std::env::var("GUTOE_GALACTIC_MAP_N")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(120_000);
    let seed = std::env::var("GUTOE_GALACTIC_MAP_SEED")
        .ok()
        .and_then(|s| s.parse::<u64>().ok())
        .unwrap_or(7_031_337);
    let cfg = GalacticLifeMapConfig {
        sample_count: n,
        rng_seed: seed,
        ..GalacticLifeMapConfig::default()
    };
    let score = evaluate_galactic_life_map(cfg);

    let out_dir = std::env::var("GUTOE_GALACTIC_LIFE_MAP_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/galactic_life_map".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);
    let csv_path = out.join("galactic_life_catalog.csv");
    let json_path = out.join("galactic_life_report.json");
    let txt_path = out.join("galactic_life_report.txt");
    let png_path = out.join("galactic_life_map.png");

    let _ = render_png(&score, &png_path);

    let mut csv = File::create(&csv_path).expect("create csv");
    writeln!(
        csv,
        "id,component,x_ly,y_ly,z_ly,r_ly,mass_solar,age_gyr,metallicity,lifetime_gyr,habitable,habitability_score,stage,entropy_multiplier,signal"
    )
    .expect("write csv header");
    for p in &score.points {
        writeln!(
            csv,
            "{},{},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{},{:.6},{},{:.6},{}",
            p.seed.id,
            match p.seed.component {
                GalacticComponent::Disk => "disk",
                GalacticComponent::Bulge => "bulge",
                GalacticComponent::Halo => "halo",
            },
            p.seed.x_ly,
            p.seed.y_ly,
            p.seed.z_ly,
            p.seed.galactic_radius_ly,
            p.seed.mass_solar,
            p.seed.age_gyr,
            p.seed.metallicity,
            p.seed.main_sequence_lifetime_gyr,
            p.habitable,
            p.habitability_score,
            p.stage.as_str(),
            p.entropy_multiplier,
            p.signal
        )
        .expect("write csv row");
    }

    let mut signal_points: Vec<_> = score
        .points
        .iter()
        .filter(|p| p.signal)
        .map(|p| {
            let dx = p.seed.x_ly - SOLAR_GALACTIC_RADIUS_LY;
            let dy = p.seed.y_ly;
            let dz = p.seed.z_ly;
            let d = (dx * dx + dy * dy + dz * dz).sqrt();
            (d, p)
        })
        .collect();
    signal_points.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    let top_signals = signal_points.into_iter().take(24).collect::<Vec<_>>();

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[galactic_life_map]").expect("write");
    writeln!(txt, "sample_count = {}", score.config.sample_count).expect("write");
    writeln!(txt, "rng_seed = {}", score.config.rng_seed).expect("write");
    writeln!(txt, "universe_age_gyr = {:.9}", score.universe_age_gyr).expect("write");
    writeln!(
        txt,
        "habitable_count_present = {}",
        score.habitable_count_present
    )
    .expect("write");
    writeln!(txt, "signal_count_present = {}", score.signal_count_present).expect("write");
    writeln!(
        txt,
        "signal_fraction_present = {:.12}",
        score.present_signal_fraction()
    )
    .expect("write");
    writeln!(
        txt,
        "predicted_signal_count_milky_way_present = {:.6e}",
        score.predicted_signal_count_milky_way_present
    )
    .expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[stage_counts_present]").expect("write");
    for (stage, n) in &score.stage_counts_present {
        writeln!(txt, "{} = {}", stage.as_str(), n).expect("write");
    }
    writeln!(txt).expect("write");
    writeln!(txt, "[forecast]").expect("write");
    for f in &score.forecasts {
        writeln!(
            txt,
            "age={:.6} Gyr (Δ={:.3}): habitable={} signal={} typeI={} typeII={} typeIII={} mean_entropy={:.6}",
            f.epoch_age_gyr,
            f.delta_gyr,
            f.habitable_count,
            f.signal_count,
            f.type_i_count,
            f.type_ii_count,
            f.type_iii_count,
            f.mean_entropy_multiplier_habitable
        )
        .expect("write");
    }
    writeln!(txt).expect("write");
    writeln!(txt, "[closest_signals]").expect("write");
    for (d, p) in &top_signals {
        writeln!(
            txt,
            "id={} distance_ly={:.3} pos=({:.3},{:.3},{:.3}) stage={} entropy={:.6}",
            p.seed.id,
            d,
            p.seed.x_ly,
            p.seed.y_ly,
            p.seed.z_ly,
            p.stage.as_str(),
            p.entropy_multiplier
        )
        .expect("write");
    }

    let payload = json!({
        "summary": {
            "sample_count": score.config.sample_count,
            "rng_seed": score.config.rng_seed,
            "universe_age_gyr": score.universe_age_gyr,
            "habitable_count_present": score.habitable_count_present,
            "signal_count_present": score.signal_count_present,
            "signal_fraction_present": score.present_signal_fraction(),
            "predicted_signal_count_milky_way_present": score.predicted_signal_count_milky_way_present
        },
        "thresholds_gyr": {
            "prebiotic": score.thresholds.prebiotic_age_gyr,
            "autocatalytic": score.thresholds.autocatalytic_age_gyr,
            "photosynthetic": score.thresholds.photosynthetic_age_gyr,
            "multicellular": score.thresholds.multicellular_age_gyr,
            "intelligence": score.thresholds.intelligence_age_gyr,
            "kardashev_i": score.thresholds.kardashev_i_age_gyr,
            "kardashev_ii": score.thresholds.kardashev_ii_age_gyr,
            "kardashev_iii": score.thresholds.kardashev_iii_age_gyr
        },
        "entropy_multipliers": {
            "bare_rock": score.multipliers.bare_rock,
            "prebiotic": score.multipliers.prebiotic,
            "autocatalytic": score.multipliers.autocatalytic,
            "photosynthetic": score.multipliers.photosynthetic,
            "multicellular": score.multipliers.multicellular,
            "intelligence": score.multipliers.intelligence,
            "kardashev_i": score.multipliers.kardashev_i,
            "kardashev_ii": score.multipliers.kardashev_ii,
            "kardashev_iii": score.multipliers.kardashev_iii
        },
        "stage_counts_present": score.stage_counts_present.iter().map(|(s,n)| {
            json!({"stage": s.as_str(), "count": n})
        }).collect::<Vec<_>>(),
        "forecast": score.forecasts.iter().map(|f| json!({
            "epoch_age_gyr": f.epoch_age_gyr,
            "delta_gyr": f.delta_gyr,
            "habitable_count": f.habitable_count,
            "signal_count": f.signal_count,
            "type_i_count": f.type_i_count,
            "type_ii_count": f.type_ii_count,
            "type_iii_count": f.type_iii_count,
            "mean_entropy_multiplier_habitable": f.mean_entropy_multiplier_habitable
        })).collect::<Vec<_>>(),
        "closest_signals": top_signals.iter().map(|(d, p)| json!({
            "id": p.seed.id,
            "distance_ly": d,
            "x_ly": p.seed.x_ly,
            "y_ly": p.seed.y_ly,
            "z_ly": p.seed.z_ly,
            "stage": p.stage.as_str(),
            "entropy_multiplier": p.entropy_multiplier
        })).collect::<Vec<_>>(),
        "artifacts": {
            "map_png": png_path.display().to_string(),
            "catalog_csv": csv_path.display().to_string(),
            "report_txt": txt_path.display().to_string()
        }
    });
    let mut json_file = File::create(&json_path).expect("create json");
    writeln!(
        json_file,
        "{}",
        serde_json::to_string_pretty(&payload).expect("serialize")
    )
    .expect("write json");

    println!("wrote {}", png_path.display());
    println!("wrote {}", csv_path.display());
    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "galactic_life_map: sample={} habitable={} signal={} predicted_MW={:.3e}",
        score.config.sample_count,
        score.habitable_count_present,
        score.signal_count_present,
        score.predicted_signal_count_milky_way_present
    );
}
