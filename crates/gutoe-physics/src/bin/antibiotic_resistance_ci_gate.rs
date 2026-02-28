//! CI gate for beta-lactamase inhibitor resistance ranking lane.

use gutoe_physics::default_antibiotic_resistance_panel;
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

fn env_usize(key: &str, default: usize) -> usize {
    std::env::var(key)
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(default)
}

fn main() {
    let temperature_k = env_f64("GUTOE_BETA_LACTAMASE_TEMP_K", 310.15);
    let panel = default_antibiotic_resistance_panel(temperature_k);

    let min_pair_count = env_usize("GUTOE_BETA_LACTAMASE_MIN_PAIR_COUNT", 15);
    let max_mean_abs_log10_error = env_f64("GUTOE_BETA_LACTAMASE_MAX_MEAN_ABS_LOG10_ERROR", 1.25);
    let max_ndm_occ_1u_m = env_f64("GUTOE_BETA_LACTAMASE_MAX_NDM_OCC_1UM", 0.10);

    let count_ok = panel.rows.len() >= min_pair_count;
    let error_ok = panel.mean_abs_log10_error <= max_mean_abs_log10_error;
    let ndm_ok = panel.ndm_max_predicted_occupancy_at_1u_m <= max_ndm_occ_1u_m;

    let tem_pred = panel.best_by_enzyme.iter().find(|b| b.enzyme_name == "TEM-1");
    let kpc_pred = panel.best_by_enzyme.iter().find(|b| b.enzyme_name == "KPC");

    let tem_winner_ok = tem_pred
        .map(|b| b.by_predicted_inhibitor == "avibactam")
        .unwrap_or(false);
    let kpc_winner_ok = kpc_pred
        .map(|b| b.by_predicted_inhibitor == "avibactam" || b.by_predicted_inhibitor == "vaborbactam")
        .unwrap_or(false);

    let overall_pass = count_ok && error_ok && ndm_ok && tem_winner_ok && kpc_winner_ok;

    let out_dir = std::env::var("GUTOE_BETA_LACTAMASE_GATE_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let json_path = format!("{out_dir}/antibiotic_resistance_ci_gate.json");

    let payload = json!({
        "overall_pass": overall_pass,
        "windows": {
            "temperature_k": temperature_k,
            "min_pair_count": min_pair_count,
            "max_mean_abs_log10_error": max_mean_abs_log10_error,
            "max_ndm_occupancy_at_1uM": max_ndm_occ_1u_m
        },
        "summary": {
            "pair_count": panel.rows.len(),
            "mean_abs_log10_error_pred_vs_anchor": panel.mean_abs_log10_error,
            "ndm_max_predicted_occupancy_at_1uM": panel.ndm_max_predicted_occupancy_at_1u_m,
            "tem_predicted_winner": tem_pred.map(|b| b.by_predicted_inhibitor).unwrap_or("missing"),
            "kpc_predicted_winner": kpc_pred.map(|b| b.by_predicted_inhibitor).unwrap_or("missing")
        },
        "gate": {
            "count_ok": count_ok,
            "error_ok": error_ok,
            "ndm_occupancy_ok": ndm_ok,
            "tem_winner_ok": tem_winner_ok,
            "kpc_winner_ok": kpc_winner_ok
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
        "antibiotic_resistance_gate: pass={} pairs={} mean_abs_log10_err={:.3} ndm_max_occ_1uM={:.3}",
        overall_pass,
        panel.rows.len(),
        panel.mean_abs_log10_error,
        panel.ndm_max_predicted_occupancy_at_1u_m
    );
    println!("wrote {json_path}");

    if !overall_pass {
        eprintln!(
            "FAIL: count_ok={} error_ok={} ndm_ok={} tem_winner_ok={} kpc_winner_ok={}",
            count_ok, error_ok, ndm_ok, tem_winner_ok, kpc_winner_ok
        );
        process::exit(2);
    }
}
