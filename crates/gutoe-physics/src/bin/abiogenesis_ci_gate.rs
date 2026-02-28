//! Abiogenesis CI gate: Kauffman closure inevitability check.

use gutoe_physics::{evaluate_abiogenesis_gate, AbiogenesisWindows};
use serde_json::json;
use std::fs::{self, File};
use std::io::Write;
use std::process;

const DEFAULT_TEMP_K: f64 = 298.15;

fn main() {
    let temperature_k = std::env::var("GUTOE_ABIOGENESIS_TEMP_K")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(DEFAULT_TEMP_K);
    let windows = AbiogenesisWindows::default();
    let score = evaluate_abiogenesis_gate(windows, temperature_k);
    let overall_pass = score.passes_all();

    let out_dir = std::env::var("GUTOE_ABIOGENESIS_GATE_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let json_path = format!("{out_dir}/abiogenesis_ci_gate.json");

    let payload = json!({
        "overall_pass": overall_pass,
        "temperature_k": temperature_k,
        "windows": {
            "closure_threshold": windows.closure_threshold,
            "catalytic_probability_min": windows.catalytic_probability_min,
            "robust_margin_min": windows.robust_margin_min
        },
        "prebiotic": {
            "feedstock_species": score.prebiotic.feedstock_species,
            "amino_acid_pool_left": score.prebiotic.amino_acid_pool_left,
            "nucleotide_pool": score.prebiotic.nucleotide_pool,
            "catalytic_probability_lower_bound": score.prebiotic.catalytic_probability_lower_bound
        },
        "closure": {
            "n_times_p": score.closure.n_times_p,
            "threshold": score.closure.threshold,
            "closure_excess": score.closure.closure_excess
        },
        "inevitability": {
            "n_times_p_lower_3sigma": score.inevitability.n_times_p_lower_3sigma,
            "robust_margin": score.inevitability.robust_margin,
            "pved_delta_e_ev": score.inevitability.pved_delta_e_ev
        },
        "gate": {
            "prebiotic_ok": score.prebiotic_ok,
            "closure_ok": score.closure_ok,
            "inevitability_ok": score.inevitability_ok
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
        "Abiogenesis gate: pass={} (N*p={:.6}, lower_3σ={:.6}, margin={:.6}, p_min={:.6}, threshold={:.6})",
        overall_pass,
        score.closure.n_times_p,
        score.inevitability.n_times_p_lower_3sigma,
        score.inevitability.robust_margin,
        score.prebiotic.catalytic_probability_lower_bound,
        score.closure.threshold,
    );
    println!("wrote {json_path}");

    if !overall_pass {
        eprintln!(
            "FAIL: prebiotic_ok={} closure_ok={} inevitability_ok={}",
            score.prebiotic_ok, score.closure_ok, score.inevitability_ok
        );
        process::exit(2);
    }
}
