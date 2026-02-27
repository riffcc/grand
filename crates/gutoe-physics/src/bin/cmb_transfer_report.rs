//! GRAND-351: CMB/BAO transfer report from derived GUTOE inputs.

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

    let out_dir = std::env::var("GUTOE_CMB_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let txt_path = format!("{out_dir}/cmb_transfer_report.txt");
    let json_path = format!("{out_dir}/cmb_transfer_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[transfer_inputs]").expect("write");
    writeln!(txt, "H0_km_s_mpc = {:.9}", a.h0_km_s_mpc).expect("write");
    writeln!(txt, "omega_b0 = {:.12}", a.omega_b0).expect("write");
    writeln!(txt, "omega_m0 = {:.12}", a.omega_m0).expect("write");
    writeln!(txt, "omega_r0 = {:.12}", a.omega_r0).expect("write");
    writeln!(txt, "omega_lambda0 = {:.12}", a.omega_lambda0).expect("write");
    writeln!(txt, "n_s = {:.12}", a.n_s).expect("write");
    writeln!(txt, "A_s = {:.12e}", a.a_s).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[transfer_outputs]").expect("write");
    writeln!(txt, "z_drag = {:.9}", s.z_drag).expect("write");
    writeln!(txt, "r_s_drag_mpc = {:.9}", s.rs_drag_mpc).expect("write");
    writeln!(txt, "theta_star_rad = {:.9e}", s.theta_star_rad).expect("write");
    writeln!(txt, "l_peak1 = {:.9}", s.l_peak1).expect("write");
    writeln!(txt, "l_peak2 = {:.9}", s.l_peak2).expect("write");
    writeln!(txt, "growth_z0 = {:.9}", s.growth_z0).expect("write");
    writeln!(txt, "growth_z1 = {:.9}", s.growth_z1).expect("write");
    writeln!(txt, "pk_pivot_z0 = {:.12e}", s.pk_pivot_z0).expect("write");
    writeln!(txt, "pk_pivot_z1 = {:.12e}", s.pk_pivot_z1).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[residuals]").expect("write");
    writeln!(txt, "rs_rel_error = {:.12}", s.rs_rel_error).expect("write");
    writeln!(txt, "theta_star_rel_error = {:.12}", s.theta_star_rel_error).expect("write");
    writeln!(txt, "l1_rel_error = {:.12}", s.l1_rel_error).expect("write");
    writeln!(txt, "l2_rel_error = {:.12}", s.l2_rel_error).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[gate]").expect("write");
    writeln!(txt, "rs_ok = {}", s.rs_ok).expect("write");
    writeln!(txt, "theta_star_ok = {}", s.theta_star_ok).expect("write");
    writeln!(txt, "l1_ok = {}", s.l1_ok).expect("write");
    writeln!(txt, "l2_ok = {}", s.l2_ok).expect("write");
    writeln!(txt, "transfer_positive_ok = {}", s.transfer_positive_ok).expect("write");
    writeln!(txt, "passes_all = {}", s.passes_all()).expect("write");

    let mut json = File::create(&json_path).expect("create json");
    writeln!(
        json,
        "{{\n  \"inputs\": {{\"h0_km_s_mpc\": {:.9}, \"omega_b0\": {:.12}, \"omega_m0\": {:.12}, \"omega_r0\": {:.12}, \"omega_k0\": {:.12}, \"omega_lambda0\": {:.12}, \"n_s\": {:.12}, \"a_s\": {:.12e}}},\n  \"windows\": {{\"rs_rel_max\": {:.12}, \"theta_star_rel_max\": {:.12}, \"l1_rel_max\": {:.12}, \"l2_rel_max\": {:.12}}},\n  \"score\": {{\"h\": {:.9}, \"omega_b_h2\": {:.9}, \"omega_m_h2\": {:.9}, \"z_drag\": {:.9}, \"z_recomb\": {:.9}, \"rs_drag_mpc\": {:.9}, \"dm_recomb_mpc\": {:.9}, \"theta_star_rad\": {:.9e}, \"acoustic_scale_la\": {:.9}, \"l_peak1\": {:.9}, \"l_peak2\": {:.9}, \"growth_z0\": {:.9}, \"growth_z1\": {:.9}, \"pk_pivot_z0\": {:.12e}, \"pk_pivot_z1\": {:.12e}, \"rs_rel_error\": {:.12}, \"theta_star_rel_error\": {:.12}, \"l1_rel_error\": {:.12}, \"l2_rel_error\": {:.12}, \"rs_ok\": {}, \"theta_star_ok\": {}, \"l1_ok\": {}, \"l2_ok\": {}, \"transfer_positive_ok\": {}, \"passes_all\": {}}}\n}}",
        a.h0_km_s_mpc,
        a.omega_b0,
        a.omega_m0,
        a.omega_r0,
        a.omega_k0,
        a.omega_lambda0,
        a.n_s,
        a.a_s,
        w.rs_rel_max,
        w.theta_star_rel_max,
        w.l1_rel_max,
        w.l2_rel_max,
        s.h,
        s.omega_b_h2,
        s.omega_m_h2,
        s.z_drag,
        s.z_recomb,
        s.rs_drag_mpc,
        s.dm_recomb_mpc,
        s.theta_star_rad,
        s.acoustic_scale_la,
        s.l_peak1,
        s.l_peak2,
        s.growth_z0,
        s.growth_z1,
        s.pk_pivot_z0,
        s.pk_pivot_z1,
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
    ).expect("write json");

    println!("wrote {txt_path}");
    println!("wrote {json_path}");
    println!(
        "CMB transfer: r_s={:.3} Mpc, theta*={:.6e}, l1={:.2}, l2={:.2}, pass={}",
        s.rs_drag_mpc,
        s.theta_star_rad,
        s.l_peak1,
        s.l_peak2,
        s.passes_all(),
    );
}
