//! Multi-cannabinoid panel CI gate.

use gutoe_physics::{default_cannabinoid_specs, evaluate_cannabinoid_panel, NeuronCouplingInput};
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
    let specs = default_cannabinoid_specs();
    let rows = evaluate_cannabinoid_panel(&specs, 310.15, NeuronCouplingInput::default());

    let min_count = env_f64("GUTOE_CANNABINOID_MIN_COUNT", 10.0) as usize;
    let explained_min = env_f64("GUTOE_CANNABINOID_MEAN_EXPLAINED_MIN", 0.10);
    let abs_resid_err_max = env_f64("GUTOE_CANNABINOID_MEAN_ABS_RESIDUAL_ERROR_MAX", 8.0);

    let mean_explained = rows
        .iter()
        .map(|r| r.explained_fraction_of_abs_delta_g)
        .sum::<f64>()
        / rows.len().max(1) as f64;
    let mean_abs_resid_error = rows
        .iter()
        .map(|r| r.residual_closure_error_kj_mol.abs())
        .sum::<f64>()
        / rows.len().max(1) as f64;

    let count_ok = rows.len() >= min_count;
    let explained_ok = mean_explained >= explained_min;
    let residual_ok = mean_abs_resid_error <= abs_resid_err_max;

    // Sanity: high-affinity THC should exceed CBD occupancy at 100 nM.
    let thc = rows.iter().find(|r| r.name == "delta9_thc");
    let cbd = rows.iter().find(|r| r.name == "cbd");
    let potency_order_ok = if let (Some(thc), Some(cbd)) = (thc, cbd) {
        thc.occupancy_100nm > cbd.occupancy_100nm
    } else {
        false
    };

    let overall_pass = count_ok && explained_ok && residual_ok && potency_order_ok;

    let out_dir = std::env::var("GUTOE_CANNABINOID_GATE_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let json_path = format!("{out_dir}/cannabinoid_panel_ci_gate.json");

    let payload = json!({
        "overall_pass": overall_pass,
        "windows": {
            "min_count": min_count,
            "mean_explained_min": explained_min,
            "mean_abs_residual_error_max": abs_resid_err_max
        },
        "summary": {
            "count": rows.len(),
            "mean_explained": mean_explained,
            "mean_abs_residual_error": mean_abs_resid_error
        },
        "gate": {
            "count_ok": count_ok,
            "mean_explained_ok": explained_ok,
            "mean_abs_residual_ok": residual_ok,
            "potency_order_ok": potency_order_ok
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
        "cannabinoid_panel_gate: pass={} (count={}, mean_explained={:.3}, mean_abs_resid_err={:.3}, potency_order_ok={})",
        overall_pass,
        rows.len(),
        mean_explained,
        mean_abs_resid_error,
        potency_order_ok
    );
    println!("wrote {json_path}");

    if !overall_pass {
        eprintln!(
            "FAIL: count_ok={} explained_ok={} residual_ok={} potency_order_ok={}",
            count_ok, explained_ok, residual_ok, potency_order_ok
        );
        process::exit(2);
    }
}
