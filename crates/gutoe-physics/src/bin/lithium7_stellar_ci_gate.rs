//! Pop-II stellar lithium depletion CI gate.

use gutoe_physics::{
    evaluate_lithium7_stellar_depletion_default, lithium7_stellar_closure_pass,
    LI7_STELLAR_CLOSURE_DELTA_ABS_MAX,
};
use std::fs::{self, File};
use std::io::Write;
use std::process;

fn main() {
    let report = evaluate_lithium7_stellar_depletion_default();
    let best = report.best_match;
    let closure_delta_abs = best.closure_delta.abs();
    let overall_pass = lithium7_stellar_closure_pass(&report);

    let out_dir =
        std::env::var("GUTOE_LI7_STELLAR_GATE_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let json_path = format!("{out_dir}/lithium7_stellar_ci_gate.json");
    let mut json = File::create(&json_path).expect("create gate json");

    writeln!(
        json,
        "{{\n  \"overall_pass\": {},\n  \"closure_delta_abs_max\": {:.12},\n  \"target_from_bbn\": {{\"eta10\": {:.12}, \"required_survival_factor\": {:.12}, \"required_depletion_percent\": {:.12}}},\n  \"best_match\": {{\"label\": \"{}\", \"mass_solar\": {:.6}, \"metallicity_z\": {:.8}, \"survival_factor\": {:.12}, \"depletion_percent\": {:.12}, \"closure_delta\": {:.12}, \"closure_delta_abs\": {:.12}}},\n  \"agreement_with_required\": {:.12}\n}}",
        overall_pass,
        LI7_STELLAR_CLOSURE_DELTA_ABS_MAX,
        report.eta10,
        report.required_survival_factor,
        report.required_depletion_percent,
        best.input.label,
        best.input.mass_solar,
        best.input.metallicity_z,
        best.survival_factor,
        best.depletion_percent,
        best.closure_delta,
        closure_delta_abs,
        report.agreement_with_required
    )
    .expect("write gate json");

    println!(
        "Li-7 stellar gate: pass={} (required={:.12}, best={} -> {:.12}, |Δ|={:.12}, max={:.12})",
        overall_pass,
        report.required_survival_factor,
        best.input.label,
        best.survival_factor,
        closure_delta_abs,
        LI7_STELLAR_CLOSURE_DELTA_ABS_MAX
    );
    println!("wrote {json_path}");

    if !overall_pass {
        eprintln!(
            "FAIL: best_case={} closure_delta_abs={:.12} exceeds max={:.12}",
            best.input.label, closure_delta_abs, LI7_STELLAR_CLOSURE_DELTA_ABS_MAX
        );
        process::exit(2);
    }
}
