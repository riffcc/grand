//! Cardiovascular binding CI gate.
//!
//! Ensures atorvastatin/HMGCR thermodynamic transduction remains in a
//! physically plausible band and the electrostatic floor remains non-trivial.

use gutoe_physics::{
    decompose_non_electrostatic_residual, evaluate_atorvastatin_hmgcr_binding,
    BindingBenchmarkInput, ElectrostaticProxyInput, ResidualProxyInput,
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
    let benchmark = BindingBenchmarkInput {
        ki_nanomolar: env_f64("GUTOE_CARDIO_KI_NM", BindingBenchmarkInput::default().ki_nanomolar),
        temperature_k: env_f64(
            "GUTOE_CARDIO_TEMP_K",
            BindingBenchmarkInput::default().temperature_k,
        ),
    };
    let proxy = ElectrostaticProxyInput::default();
    let score = evaluate_atorvastatin_hmgcr_binding(benchmark, proxy);
    let residual = decompose_non_electrostatic_residual(
        score.residual_required_kj_mol,
        ResidualProxyInput::default(),
    );

    let abs_exp = score.experimental_delta_g_kj_mol.abs();
    let abs_floor = score.qed_floor_total_kj_mol.abs();
    let abs_residual = score.residual_required_kj_mol.abs();
    let abs_residual_error = residual.closure_error_kj_mol.abs();

    let exp_abs_min = env_f64("GUTOE_CARDIO_EXP_ABS_MIN_KJ_MOL", 40.0);
    let exp_abs_max = env_f64("GUTOE_CARDIO_EXP_ABS_MAX_KJ_MOL", 55.0);
    let floor_abs_min = env_f64("GUTOE_CARDIO_QED_FLOOR_ABS_MIN_KJ_MOL", 20.0);
    let explained_min = env_f64("GUTOE_CARDIO_EXPLAINED_MIN", 0.50);
    let residual_abs_max = env_f64("GUTOE_CARDIO_RESIDUAL_ABS_MAX_KJ_MOL", 35.0);
    let residual_closure_abs_max = env_f64("GUTOE_CARDIO_RESIDUAL_CLOSURE_ABS_MAX_KJ_MOL", 3.0);

    let exp_band_ok = abs_exp >= exp_abs_min && abs_exp <= exp_abs_max;
    let qed_floor_ok = abs_floor >= floor_abs_min && score.qed_floor_total_kj_mol < 0.0;
    let explained_ok = score.explained_fraction_of_abs_delta_g >= explained_min;
    let residual_ok = score.residual_required_kj_mol < 0.0 && abs_residual <= residual_abs_max;
    let residual_closure_ok = abs_residual_error <= residual_closure_abs_max;
    let overall_pass =
        exp_band_ok && qed_floor_ok && explained_ok && residual_ok && residual_closure_ok;

    let out_dir = std::env::var("GUTOE_CARDIO_BINDING_GATE_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let json_path = format!("{out_dir}/cardiovascular_binding_ci_gate.json");

    let payload = json!({
        "overall_pass": overall_pass,
        "benchmark": {
            "ki_nanomolar": benchmark.ki_nanomolar,
            "temperature_k": benchmark.temperature_k
        },
        "windows": {
            "experimental_abs_min_kj_mol": exp_abs_min,
            "experimental_abs_max_kj_mol": exp_abs_max,
            "qed_floor_abs_min_kj_mol": floor_abs_min,
            "explained_min": explained_min,
            "residual_abs_max_kj_mol": residual_abs_max,
            "residual_closure_abs_max_kj_mol": residual_closure_abs_max
        },
        "score": {
            "experimental_delta_g_kj_mol": score.experimental_delta_g_kj_mol,
            "qed_floor_total_kj_mol": score.qed_floor_total_kj_mol,
            "residual_required_kj_mol": score.residual_required_kj_mol,
            "explained_fraction_of_abs_delta_g": score.explained_fraction_of_abs_delta_g,
            "residual_modeled_total_kj_mol": residual.modeled_residual_total_kj_mol,
            "residual_closure_error_kj_mol": residual.closure_error_kj_mol
        },
        "gate": {
            "experimental_band_ok": exp_band_ok,
            "qed_floor_ok": qed_floor_ok,
            "explained_ok": explained_ok,
            "residual_ok": residual_ok,
            "residual_closure_ok": residual_closure_ok
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
        "Cardiovascular binding gate: pass={} (|ΔG_exp|={:.3}, |QED_floor|={:.3}, |residual|={:.3}, |residual_error|={:.3}, explained={:.3})",
        overall_pass,
        abs_exp,
        abs_floor,
        abs_residual,
        abs_residual_error,
        score.explained_fraction_of_abs_delta_g
    );
    println!("wrote {json_path}");

    if !overall_pass {
        eprintln!(
            "FAIL: exp_band_ok={} qed_floor_ok={} explained_ok={} residual_ok={} residual_closure_ok={}",
            exp_band_ok, qed_floor_ok, explained_ok, residual_ok, residual_closure_ok
        );
        process::exit(2);
    }
}
