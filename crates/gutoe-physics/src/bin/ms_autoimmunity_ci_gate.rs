//! MS autoimmunity lane CI gate.

use gutoe_physics::{
    default_ms_mimicry_input, default_natalizumab_proxy, default_ocrelizumab_proxy,
    default_targeted_blocker_input, evaluate_molecular_mimicry, evaluate_targeted_blocker,
    evaluate_therapy_effect,
};
use serde_json::json;
use std::fs::{self, File};
use std::io::Write;
use std::process;

fn env_f64(key: &str, default: f64) -> f64 {
    std::env::var(key)
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(default)
}

fn main() {
    let mimicry = evaluate_molecular_mimicry(default_ms_mimicry_input());
    let therapy = evaluate_therapy_effect(
        mimicry.misrecognition_risk_index,
        default_ocrelizumab_proxy(),
        default_natalizumab_proxy(),
    );
    let blocker =
        evaluate_targeted_blocker(mimicry.activation_excess_kj_mol, default_targeted_blocker_input());

    let gap_max = env_f64("GUTOE_MS_MIMICRY_GAP_MAX_KJ_MOL", 2.0);
    let activation_excess_min = env_f64("GUTOE_MS_ACTIVATION_EXCESS_MIN_KJ_MOL", 0.2);
    let activation_excess_max = env_f64("GUTOE_MS_ACTIVATION_EXCESS_MAX_KJ_MOL", 3.0);
    let therapy_reduction_min = env_f64("GUTOE_MS_THERAPY_REDUCTION_MIN", 0.5);
    let blocker_required_occ_max = env_f64("GUTOE_MS_BLOCKER_REQUIRED_OCC_MAX", 1.0);

    let gap_ok = mimicry.mimicry_gap_kj_mol <= gap_max;
    let activation_ok = mimicry.activation_excess_kj_mol >= activation_excess_min
        && mimicry.activation_excess_kj_mol <= activation_excess_max;
    let therapy_ok = therapy.relative_drive_reduction_fraction >= therapy_reduction_min;
    let blocker_ok = blocker.required_occupancy_fraction <= blocker_required_occ_max
        && blocker.feasible_at_given_concentration;
    let overall_pass = gap_ok && activation_ok && therapy_ok && blocker_ok;

    let out_dir = std::env::var("GUTOE_MS_AUTOIMMUNE_GATE_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let json_path = format!("{out_dir}/ms_autoimmunity_ci_gate.json");

    let payload = json!({
        "overall_pass": overall_pass,
        "windows": {
            "mimicry_gap_max_kj_mol": gap_max,
            "activation_excess_min_kj_mol": activation_excess_min,
            "activation_excess_max_kj_mol": activation_excess_max,
            "therapy_reduction_min": therapy_reduction_min,
            "blocker_required_occupancy_max": blocker_required_occ_max
        },
        "score": {
            "mimicry_gap_kj_mol": mimicry.mimicry_gap_kj_mol,
            "activation_excess_kj_mol": mimicry.activation_excess_kj_mol,
            "misrecognition_risk_index": mimicry.misrecognition_risk_index,
            "therapy_relative_drive_reduction_fraction": therapy.relative_drive_reduction_fraction,
            "blocker_required_occupancy_fraction": blocker.required_occupancy_fraction,
            "blocker_feasible": blocker.feasible_at_given_concentration
        },
        "gate": {
            "mimicry_gap_ok": gap_ok,
            "activation_window_ok": activation_ok,
            "therapy_reduction_ok": therapy_ok,
            "blocker_feasible_ok": blocker_ok
        }
    });

    let mut json_file = File::create(&json_path).expect("create gate json");
    writeln!(
        json_file,
        "{}",
        serde_json::to_string_pretty(&payload).expect("serialize")
    )
    .expect("write gate json");

    println!(
        "ms_autoimmunity_gate: pass={} (gap={:.3}, activation_excess={:.3}, therapy_reduction={:.3}, blocker_required_occ={:.3})",
        overall_pass,
        mimicry.mimicry_gap_kj_mol,
        mimicry.activation_excess_kj_mol,
        therapy.relative_drive_reduction_fraction,
        blocker.required_occupancy_fraction
    );
    println!("wrote {json_path}");

    if !overall_pass {
        eprintln!(
            "FAIL: gap_ok={} activation_ok={} therapy_ok={} blocker_ok={}",
            gap_ok, activation_ok, therapy_ok, blocker_ok
        );
        process::exit(2);
    }
}
