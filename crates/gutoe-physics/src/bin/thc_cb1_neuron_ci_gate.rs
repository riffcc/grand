//! THC-CB1 neuron CI gate.

use gutoe_physics::{
    decompose_thc_cb1_non_electrostatic_residual, evaluate_thc_cb1_binding,
    simulate_thc_cb1_neuron_response, NeuronCouplingInput, ThcCb1BindingInput,
    ThcElectrostaticProxyInput, ThcResidualProxyInput,
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
    let binding = ThcCb1BindingInput::default();
    let electro = ThcElectrostaticProxyInput::default();
    let residual_proxy = ThcResidualProxyInput::default();
    let coupling = NeuronCouplingInput::default();

    let score = evaluate_thc_cb1_binding(binding, electro);
    let residual =
        decompose_thc_cb1_non_electrostatic_residual(score.residual_required_kj_mol, residual_proxy);
    let sweep = [0.0, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0];
    let points = simulate_thc_cb1_neuron_response(binding, coupling, &sweep);

    let exp_abs_min = env_f64("GUTOE_THC_EXP_ABS_MIN_KJ_MOL", 38.0);
    let exp_abs_max = env_f64("GUTOE_THC_EXP_ABS_MAX_KJ_MOL", 50.0);
    let floor_abs_min = env_f64("GUTOE_THC_QED_FLOOR_ABS_MIN_KJ_MOL", 3.0);
    let floor_abs_max = env_f64("GUTOE_THC_QED_FLOOR_ABS_MAX_KJ_MOL", 15.0);
    let residual_closure_abs_max = env_f64("GUTOE_THC_RESIDUAL_CLOSURE_ABS_MAX_KJ_MOL", 5.0);

    let abs_exp = score.experimental_delta_g_kj_mol.abs();
    let abs_floor = score.qed_floor_total_kj_mol.abs();
    let abs_residual_closure = residual.closure_error_kj_mol.abs();

    let exp_band_ok = abs_exp >= exp_abs_min && abs_exp <= exp_abs_max;
    let floor_band_ok =
        abs_floor >= floor_abs_min && abs_floor <= floor_abs_max && score.qed_floor_total_kj_mol < 0.0;
    let residual_closure_ok = abs_residual_closure <= residual_closure_abs_max;

    let mut monotone_ok = true;
    for i in 1..points.len() {
        if points[i].occupancy_fraction + 1.0e-12 < points[i - 1].occupancy_fraction {
            monotone_ok = false;
            break;
        }
        if points[i].release_probability > points[i - 1].release_probability + 1.0e-12 {
            monotone_ok = false;
            break;
        }
        if points[i].firing_rate_hz > points[i - 1].firing_rate_hz + 1.0e-12 {
            monotone_ok = false;
            break;
        }
    }

    let overall_pass = exp_band_ok && floor_band_ok && residual_closure_ok && monotone_ok;

    let out_dir = std::env::var("GUTOE_THC_CB1_GATE_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let json_path = format!("{out_dir}/thc_cb1_neuron_ci_gate.json");

    let payload = json!({
        "overall_pass": overall_pass,
        "windows": {
            "experimental_abs_min_kj_mol": exp_abs_min,
            "experimental_abs_max_kj_mol": exp_abs_max,
            "qed_floor_abs_min_kj_mol": floor_abs_min,
            "qed_floor_abs_max_kj_mol": floor_abs_max,
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
            "qed_floor_band_ok": floor_band_ok,
            "residual_closure_ok": residual_closure_ok,
            "monotone_neuron_curve_ok": monotone_ok
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
        "THC-CB1 neuron gate: pass={} (|ΔG_exp|={:.3}, |QED_floor|={:.3}, |residual_error|={:.3}, monotone={})",
        overall_pass,
        abs_exp,
        abs_floor,
        abs_residual_closure,
        monotone_ok
    );
    println!("wrote {json_path}");

    if !overall_pass {
        eprintln!(
            "FAIL: exp_band_ok={} floor_band_ok={} residual_closure_ok={} monotone_ok={}",
            exp_band_ok, floor_band_ok, residual_closure_ok, monotone_ok
        );
        process::exit(2);
    }
}
