//! CI gate for phage-host matching lane.

use gutoe_physics::default_phage_host_matching_panel;
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
    let temperature_k = env_f64("GUTOE_PHAGE_MATCH_TEMP_K", 310.15);
    let panel = default_phage_host_matching_panel(temperature_k);

    let min_pair_count = env_usize("GUTOE_PHAGE_MATCH_MIN_PAIR_COUNT", 16);
    let min_mean_best_lysis = env_f64("GUTOE_PHAGE_MATCH_MIN_MEAN_BEST_LYSIS", 0.20);
    let max_probe_abs_delta = env_f64("GUTOE_PHAGE_MATCH_MAX_PROBE_ABS_DELTA", 1.0e-12);
    let min_ndm_best_lysis = env_f64("GUTOE_PHAGE_MATCH_MIN_NDM_BEST_LYSIS", 0.20);

    let count_ok = panel.rows.len() >= min_pair_count;
    let mean_lysis_ok = panel.mean_best_lysis_score >= min_mean_best_lysis;
    let probe_ok = panel.resistance_independence_probe_abs_delta <= max_probe_abs_delta;

    let ndm_best = panel
        .best_by_strain
        .iter()
        .find(|b| b.strain_name == "kp_ndm1_clinical");
    let ndm_match_ok = ndm_best
        .map(|b| b.best_phage_name == "phi_kp_omp")
        .unwrap_or(false);
    let ndm_lysis_ok = ndm_best
        .map(|b| b.best_lysis_score >= min_ndm_best_lysis)
        .unwrap_or(false);

    let overall_pass = count_ok && mean_lysis_ok && probe_ok && ndm_match_ok && ndm_lysis_ok;

    let out_dir = std::env::var("GUTOE_PHAGE_MATCH_GATE_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let json_path = format!("{out_dir}/phage_host_matching_ci_gate.json");

    let payload = json!({
        "overall_pass": overall_pass,
        "windows": {
            "temperature_k": temperature_k,
            "min_pair_count": min_pair_count,
            "min_mean_best_lysis_score": min_mean_best_lysis,
            "max_resistance_independence_probe_abs_delta": max_probe_abs_delta,
            "min_ndm_best_lysis_score": min_ndm_best_lysis
        },
        "summary": {
            "pair_count": panel.rows.len(),
            "mean_best_lysis_score": panel.mean_best_lysis_score,
            "resistance_independence_probe_abs_delta": panel.resistance_independence_probe_abs_delta,
            "ndm_best_phage": ndm_best
                .map(|b| b.best_phage_name.clone())
                .unwrap_or_else(|| "missing".to_string()),
            "ndm_best_lysis_score": ndm_best.map(|b| b.best_lysis_score).unwrap_or(0.0)
        },
        "gate": {
            "count_ok": count_ok,
            "mean_best_lysis_ok": mean_lysis_ok,
            "probe_ok": probe_ok,
            "ndm_match_ok": ndm_match_ok,
            "ndm_lysis_ok": ndm_lysis_ok
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
        "phage_host_matching_gate: pass={} pairs={} mean_best_lysis={:.3} probe_abs_delta={:.3e}",
        overall_pass,
        panel.rows.len(),
        panel.mean_best_lysis_score,
        panel.resistance_independence_probe_abs_delta
    );
    println!("wrote {json_path}");

    if !overall_pass {
        eprintln!(
            "FAIL: count_ok={} mean_lysis_ok={} probe_ok={} ndm_match_ok={} ndm_lysis_ok={}",
            count_ok, mean_lysis_ok, probe_ok, ndm_match_ok, ndm_lysis_ok
        );
        process::exit(2);
    }
}
