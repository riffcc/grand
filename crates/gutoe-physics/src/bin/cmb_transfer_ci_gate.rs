//! GRAND-351 CI gate for CMB/BAO transfer lane.

use gutoe_physics::constants::{
    lambda_cosmological_full_candidate, C, DARK_TO_VISIBLE_GEOMETRIC_RATIO,
};
use gutoe_physics::dark_matter_falsification::OMEGA_BARYON_OBS;
use gutoe_physics::{
    evaluate_inflation_gate, evaluate_transfer_gate, InflationWindows, TransferAssumptions,
    TransferWindows,
};
use std::fs::{self, File};
use std::io::Write;
use std::process;

fn h0_from_lambda_and_omega_lambda(lambda: f64, omega_lambda: f64) -> f64 {
    let meter_per_mpc = 3.085_677_581_491_367e22;
    let h0_s_inv = C * (lambda / (3.0 * omega_lambda)).sqrt();
    h0_s_inv * meter_per_mpc / 1_000.0
}

fn main() {
    let inflation = evaluate_inflation_gate(InflationWindows::default());
    let omega_b0 = OMEGA_BARYON_OBS;
    let omega_dm0 = OMEGA_BARYON_OBS * DARK_TO_VISIBLE_GEOMETRIC_RATIO;
    let omega_m0 = omega_b0 + omega_dm0;
    let omega_r0 = 9.0e-5;
    let omega_k0 = 0.0;
    let omega_lambda0 = 1.0 - omega_m0 - omega_r0 - omega_k0;
    let h0 = h0_from_lambda_and_omega_lambda(lambda_cosmological_full_candidate(), omega_lambda0);

    let a = TransferAssumptions {
        h0_km_s_mpc: h0,
        omega_b0,
        omega_m0,
        omega_r0,
        omega_k0,
        omega_lambda0,
        n_s: inflation.n_s,
        a_s: inflation.a_s,
    };
    let w = TransferWindows::default();
    let s = evaluate_transfer_gate(a, w);

    let out_dir =
        std::env::var("GUTOE_CMB_GATE_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let json_path = format!("{out_dir}/cmb_transfer_ci_gate.json");
    let mut json = File::create(&json_path).expect("create gate json");

    writeln!(
        json,
        "{{\n  \"overall_pass\": {},\n  \"windows\": {{\"rs_rel_max\": {:.12}, \"theta_star_rel_max\": {:.12}, \"l1_rel_max\": {:.12}, \"l2_rel_max\": {:.12}}},\n  \"score\": {{\"rs_drag_mpc\": {:.9}, \"theta_star_rad\": {:.9e}, \"l_peak1\": {:.9}, \"l_peak2\": {:.9}, \"rs_rel_error\": {:.12}, \"theta_star_rel_error\": {:.12}, \"l1_rel_error\": {:.12}, \"l2_rel_error\": {:.12}, \"rs_ok\": {}, \"theta_star_ok\": {}, \"l1_ok\": {}, \"l2_ok\": {}, \"transfer_positive_ok\": {}, \"passes_all\": {}}}\n}}",
        s.passes_all(),
        w.rs_rel_max,
        w.theta_star_rel_max,
        w.l1_rel_max,
        w.l2_rel_max,
        s.rs_drag_mpc,
        s.theta_star_rad,
        s.l_peak1,
        s.l_peak2,
        s.rs_rel_error,
        s.theta_star_rel_error,
        s.l1_rel_error,
        s.l2_rel_error,
        s.rs_ok,
        s.theta_star_ok,
        s.l1_ok,
        s.l2_ok,
        s.transfer_positive_ok,
        s.passes_all(),
    ).expect("write gate json");

    println!(
        "CMB transfer gate: pass={} (r_s={:.2} Mpc, theta*={:.6e}, l1={:.2}, l2={:.2})",
        s.passes_all(),
        s.rs_drag_mpc,
        s.theta_star_rad,
        s.l_peak1,
        s.l_peak2,
    );
    println!("wrote {json_path}");

    if !s.passes_all() {
        eprintln!(
            "FAIL: rs_ok={} theta_star_ok={} l1_ok={} l2_ok={} transfer_positive_ok={}",
            s.rs_ok, s.theta_star_ok, s.l1_ok, s.l2_ok, s.transfer_positive_ok,
        );
        process::exit(2);
    }
}
