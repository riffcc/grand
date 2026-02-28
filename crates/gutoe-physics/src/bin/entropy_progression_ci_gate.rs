//! CI gate for entropy-progression thermodynamic lane.

use gutoe_physics::{
    evaluate_entropy_progression_gate, EntropyProgressionWindows, UniverseAssumptions,
    UniverseSimulationDepth, UniverseWindows, Z_INTEGRAL_MAX,
};
use serde_json::json;
use std::fs::{self, File};
use std::io::Write;
use std::process;

fn main() {
    let assumptions = UniverseAssumptions::default();
    let universe_windows = UniverseWindows::default();
    let depth = UniverseSimulationDepth {
        history_points: 768,
        history_z_max: 1.0e9,
        integral_z_max: Z_INTEGRAL_MAX,
    };
    let windows = EntropyProgressionWindows::default();
    let score = evaluate_entropy_progression_gate(assumptions, universe_windows, depth, windows);
    let overall_pass = score.passes_all();

    let out_dir = std::env::var("GUTOE_ENTROPY_PROGRESSION_GATE_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let json_path = format!("{out_dir}/entropy_progression_ci_gate.json");

    let last = score.samples.last().copied().expect("sample");
    let payload = json!({
        "overall_pass": overall_pass,
        "gate": {
            "monotone_stage_plateaus": score.monotone_stage_plateaus,
            "intelligence_step_dominant": score.intelligence_step_dominant,
            "extrema_present": score.extrema_present
        },
        "summary": {
            "universe_age_gyr": score.universe_age_gyr,
            "h0_km_s_mpc": score.h0_km_s_mpc,
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
        })).collect::<Vec<_>>()
    });

    let mut json_file = File::create(&json_path).expect("create gate json");
    writeln!(
        json_file,
        "{}",
        serde_json::to_string_pretty(&payload).expect("serialize")
    )
    .expect("write gate json");

    println!(
        "Entropy progression gate: pass={} (monotone={} intelligence_dominant={} extrema={} maxima={} minima={} final={:.3e} W/m^2/K)",
        overall_pass,
        score.monotone_stage_plateaus,
        score.intelligence_step_dominant,
        score.extrema_present,
        score.local_maxima_count,
        score.local_minima_count,
        last.total_per_area_w_m2_k
    );
    println!("wrote {json_path}");

    if !overall_pass {
        eprintln!(
            "FAIL: monotone={} intelligence_dominant={} extrema={} maxima={} minima={}",
            score.monotone_stage_plateaus,
            score.intelligence_step_dominant,
            score.extrema_present,
            score.local_maxima_count,
            score.local_minima_count,
        );
        process::exit(2);
    }
}
