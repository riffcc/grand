//! GRAND-353 CI gate for uncertainty propagation lane.

use gutoe_physics::{
    evaluate_uncertainty_gate, UncertaintyAssumptions, UncertaintyWindows, UniverseAssumptions,
    UniverseWindows,
};
use std::fs::{self, File};
use std::io::Write;
use std::process;

fn main() {
    let ua = UncertaintyAssumptions::default();
    let uw = UncertaintyWindows::default();
    let gate = evaluate_uncertainty_gate(
        UniverseAssumptions::default(),
        UniverseWindows::default(),
        ua,
        uw,
    );

    let s = &gate.summary;
    let out_dir = std::env::var("GUTOE_UNCERTAINTY_GATE_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let json_path = format!("{out_dir}/uncertainty_ci_gate.json");
    let mut json = File::create(&json_path).expect("create gate json");

    writeln!(
        json,
        "{{\n  \"overall_pass\": {},\n  \"windows\": {{\"pass_fraction_min\": {:.9}, \"h0_p95_rel_span_max\": {:.9}, \"theta_star_p95_rel_span_max\": {:.9}, \"yp_network_p95_span_max\": {:.9}}},\n  \"score\": {{\"valid_samples\": {}, \"pass_fraction\": {:.9}, \"h0_p95_rel_span\": {:.9}, \"theta_star_p95_rel_span\": {:.9}, \"yp_network_p95_span\": {:.9}, \"pass_fraction_ok\": {}, \"h0_span_ok\": {}, \"theta_star_span_ok\": {}, \"yp_span_ok\": {}, \"passes_all\": {}}}\n}}",
        gate.passes_all(),
        uw.pass_fraction_min,
        uw.h0_p95_rel_span_max,
        uw.theta_star_p95_rel_span_max,
        uw.yp_network_p95_span_max,
        s.valid_samples,
        s.pass_fraction,
        s.h0_km_s_mpc.rel_span95(),
        s.theta_star_rad.rel_span95(),
        s.yp_network.abs_span95(),
        gate.pass_fraction_ok,
        gate.h0_span_ok,
        gate.theta_star_span_ok,
        gate.yp_span_ok,
        gate.passes_all(),
    )
    .expect("write gate json");

    println!(
        "Uncertainty gate: pass={} (valid={}, pass_frac={:.3}, H0_span95={:.3}, theta*_span95={:.3}, Yp_span95={:.4})",
        gate.passes_all(),
        s.valid_samples,
        s.pass_fraction,
        s.h0_km_s_mpc.rel_span95(),
        s.theta_star_rad.rel_span95(),
        s.yp_network.abs_span95(),
    );
    println!("wrote {json_path}");

    if !gate.passes_all() {
        eprintln!(
            "FAIL: pass_fraction_ok={} h0_span_ok={} theta_star_span_ok={} yp_span_ok={}",
            gate.pass_fraction_ok, gate.h0_span_ok, gate.theta_star_span_ok, gate.yp_span_ok,
        );
        process::exit(2);
    }
}
