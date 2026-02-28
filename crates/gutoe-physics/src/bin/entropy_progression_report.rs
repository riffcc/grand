//! Entropy-progression report and heatmap for the simulated known universe.

use gutoe_physics::{
    evaluate_entropy_progression_gate, DissipativeStage, EntropyProgressionWindows,
    UniverseAssumptions, UniverseSimulationDepth, UniverseWindows,
};
use image::{ImageBuffer, Rgb, RgbImage};
use serde_json::json;
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

fn inferno_like(v: f64) -> Rgb<u8> {
    let x = v.clamp(0.0, 1.0);
    // Compact inferno-like piecewise polynomial blend.
    let r = (255.0 * (0.25 + 0.75 * x.powf(0.85))).clamp(0.0, 255.0) as u8;
    let g = (255.0 * (x.powf(1.35) * (1.0 - 0.25 * x))).clamp(0.0, 255.0) as u8;
    let b = (255.0 * ((1.0 - x).powf(1.15) * (0.55 + 0.45 * x))).clamp(0.0, 255.0) as u8;
    Rgb([r, g, b])
}

fn channel_rows(
    score: &gutoe_physics::EntropyProgressionScorecard,
) -> Vec<(&'static str, Vec<f64>)> {
    let mut baseline = Vec::with_capacity(score.samples.len());
    let mut prebiotic = Vec::with_capacity(score.samples.len());
    let mut autocatalytic = Vec::with_capacity(score.samples.len());
    let mut photosynthetic = Vec::with_capacity(score.samples.len());
    let mut multicellular = Vec::with_capacity(score.samples.len());
    let mut intelligence = Vec::with_capacity(score.samples.len());
    let mut total = Vec::with_capacity(score.samples.len());

    for s in &score.samples {
        baseline.push(s.baseline_per_area_w_m2_k);
        prebiotic.push(s.prebiotic_per_area_w_m2_k);
        autocatalytic.push(s.autocatalytic_per_area_w_m2_k);
        photosynthetic.push(s.photosynthetic_per_area_w_m2_k);
        multicellular.push(s.multicellular_per_area_w_m2_k);
        intelligence.push(s.intelligence_per_area_w_m2_k);
        total.push(s.total_per_area_w_m2_k);
    }

    vec![
        ("baseline", baseline),
        ("prebiotic", prebiotic),
        ("autocatalytic", autocatalytic),
        ("photosynthetic", photosynthetic),
        ("multicellular", multicellular),
        ("intelligence", intelligence),
        ("total", total),
    ]
}

fn render_heatmap_png(
    score: &gutoe_physics::EntropyProgressionScorecard,
    out_path: &PathBuf,
) -> Result<(), String> {
    let rows = channel_rows(score);
    let width = score.samples.len() as u32;
    let row_h = 24u32;
    let gap = 2u32;
    let height = rows.len() as u32 * row_h + (rows.len() as u32 - 1) * gap;
    let mut img: RgbImage = ImageBuffer::new(width, height);

    let mut all_vals = Vec::new();
    for (_, r) in &rows {
        all_vals.extend(r.iter().copied().filter(|v| v.is_finite() && *v > 0.0));
    }
    let eps = 1.0e-24;
    let min_log = all_vals
        .iter()
        .map(|v| (v + eps).log10())
        .fold(f64::INFINITY, f64::min);
    let max_log = all_vals
        .iter()
        .map(|v| (v + eps).log10())
        .fold(f64::NEG_INFINITY, f64::max);
    let span = (max_log - min_log).max(1.0e-12);

    for (ri, (_, row)) in rows.iter().enumerate() {
        let y0 = ri as u32 * (row_h + gap);
        for (x, v) in row.iter().enumerate() {
            let t = ((v.max(0.0) + eps).log10() - min_log) / span;
            let c = inferno_like(t);
            for yy in y0..(y0 + row_h) {
                img.put_pixel(x as u32, yy, c);
            }
        }
    }

    // Stage activation overlays.
    for act in &score.stage_activations {
        let mut ix = 0usize;
        let mut best = f64::INFINITY;
        for (i, s) in score.samples.iter().enumerate() {
            let d = (s.age_gyr - act.activation_age_gyr).abs();
            if d < best {
                best = d;
                ix = i;
            }
        }
        let x = ix as u32;
        let line = match act.stage {
            DissipativeStage::PrebioticChemistry => Rgb([220, 220, 220]),
            DissipativeStage::AutocatalyticLife => Rgb([120, 250, 120]),
            DissipativeStage::PhotosyntheticBiosphere => Rgb([120, 220, 255]),
            DissipativeStage::MulticellularEcosystem => Rgb([255, 200, 80]),
            DissipativeStage::TechnologicalIntelligence => Rgb([255, 90, 90]),
            DissipativeStage::BareRock => Rgb([255, 255, 255]),
        };
        for y in 0..height {
            img.put_pixel(x, y, line);
        }
    }

    img.save(out_path)
        .map_err(|e| format!("save {}: {e}", out_path.display()))
}

fn main() {
    let assumptions = UniverseAssumptions::default();
    let universe_windows = UniverseWindows::default();
    let depth = UniverseSimulationDepth {
        history_points: 768,
        history_z_max: 1.0e9,
        integral_z_max: gutoe_physics::Z_INTEGRAL_MAX,
    };
    let windows = EntropyProgressionWindows::default();

    let score = evaluate_entropy_progression_gate(assumptions, universe_windows, depth, windows);

    let out_dir = std::env::var("GUTOE_ENTROPY_PROGRESSION_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/entropy_progression".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);
    let txt_path = out.join("entropy_progression_report.txt");
    let json_path = out.join("entropy_progression_report.json");
    let png_path = out.join("entropy_progression_heatmap.png");

    let _ = render_heatmap_png(&score, &png_path);

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[entropy_progression]").expect("write");
    writeln!(txt, "universe_age_gyr = {:.9}", score.universe_age_gyr).expect("write");
    writeln!(txt, "h0_km_s_mpc = {:.9}", score.h0_km_s_mpc).expect("write");
    writeln!(txt, "hubble_radius_m = {:.12e}", score.hubble_radius_m).expect("write");
    writeln!(
        txt,
        "hubble_surface_area_m2 = {:.12e}",
        score.hubble_surface_area_m2
    )
    .expect("write");
    writeln!(txt, "local_maxima_count = {}", score.local_maxima_count).expect("write");
    writeln!(txt, "local_minima_count = {}", score.local_minima_count).expect("write");
    writeln!(
        txt,
        "max_positive_step_age_gyr = {:.9}",
        score.max_positive_step_age_gyr
    )
    .expect("write");
    writeln!(
        txt,
        "max_positive_step_w_m2_k = {:.12e}",
        score.max_positive_step_w_m2_k
    )
    .expect("write");
    writeln!(
        txt,
        "monotone_stage_plateaus = {}",
        score.monotone_stage_plateaus
    )
    .expect("write");
    writeln!(
        txt,
        "intelligence_step_dominant = {}",
        score.intelligence_step_dominant
    )
    .expect("write");
    writeln!(txt, "extrema_present = {}", score.extrema_present).expect("write");
    writeln!(txt, "passes_all = {}", score.passes_all()).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[stage_activations]").expect("write");
    for a in &score.stage_activations {
        writeln!(
            txt,
            "{}: activation_age_gyr={:.9}, incremental_gain={:.9}",
            a.stage.as_str(),
            a.activation_age_gyr,
            a.incremental_gain
        )
        .expect("write");
    }
    writeln!(txt).expect("write");
    writeln!(txt, "[stage_plateaus]").expect("write");
    for p in &score.stage_plateaus {
        writeln!(
            txt,
            "{}: [{:.9}, {:.9}) mean_total_per_area_w_m2_k={:.12e}, mean_effective_multiplier={:.9}",
            p.stage.as_str(),
            p.age_start_gyr,
            p.age_end_gyr,
            p.mean_total_per_area_w_m2_k,
            p.mean_effective_multiplier
        )
        .expect("write");
    }

    let last = score.samples.last().copied().expect("sample");
    let payload = json!({
        "gate": {
            "monotone_stage_plateaus": score.monotone_stage_plateaus,
            "intelligence_step_dominant": score.intelligence_step_dominant,
            "extrema_present": score.extrema_present,
            "passes_all": score.passes_all()
        },
        "summary": {
            "universe_age_gyr": score.universe_age_gyr,
            "h0_km_s_mpc": score.h0_km_s_mpc,
            "hubble_radius_m": score.hubble_radius_m,
            "hubble_surface_area_m2": score.hubble_surface_area_m2,
            "local_maxima_count": score.local_maxima_count,
            "local_minima_count": score.local_minima_count,
            "max_positive_step_age_gyr": score.max_positive_step_age_gyr,
            "max_positive_step_w_m2_k": score.max_positive_step_w_m2_k,
            "final_total_per_area_w_m2_k": last.total_per_area_w_m2_k,
            "final_total_universe_w_k": last.total_universe_w_k
        },
        "stage_activations": score.stage_activations.iter().map(|a| json!({
            "stage": a.stage.as_str(),
            "activation_age_gyr": a.activation_age_gyr,
            "incremental_gain": a.incremental_gain
        })).collect::<Vec<_>>(),
        "stage_plateaus": score.stage_plateaus.iter().map(|p| json!({
            "stage": p.stage.as_str(),
            "age_start_gyr": p.age_start_gyr,
            "age_end_gyr": p.age_end_gyr,
            "mean_total_per_area_w_m2_k": p.mean_total_per_area_w_m2_k,
            "mean_effective_multiplier": p.mean_effective_multiplier
        })).collect::<Vec<_>>(),
        "artifacts": {
            "heatmap_png": png_path.display().to_string()
        }
    });

    let mut json_file = File::create(&json_path).expect("create json");
    writeln!(
        json_file,
        "{}",
        serde_json::to_string_pretty(&payload).expect("serialize")
    )
    .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!("wrote {}", png_path.display());
    println!(
        "Entropy progression: pass={} maxima={} minima={} final_total={:.3e} W/m^2/K",
        score.passes_all(),
        score.local_maxima_count,
        score.local_minima_count,
        last.total_per_area_w_m2_k
    );
}
