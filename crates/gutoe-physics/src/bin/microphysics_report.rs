//! GRAND-352: explicit BBN/recombination microphysics report.

use gutoe_physics::{evaluate_microphysics_gate, MicrophysicsAssumptions, MicrophysicsWindows};
use std::fs::{self, File};
use std::io::Write;

fn baseline() -> MicrophysicsAssumptions {
    MicrophysicsAssumptions {
        h0_km_s_mpc: 68.0163,
        omega_b0: 0.0493,
        omega_m0: 0.3182,
        omega_r0: 9.0e-5,
        omega_k0: 0.0,
        omega_lambda0: 1.0 - 0.3182 - 9.0e-5,
        eta10: 5.938,
    }
}

fn main() {
    let a = baseline();
    let w = MicrophysicsWindows::default();
    let s = evaluate_microphysics_gate(a, w);

    let out_dir =
        std::env::var("GUTOE_MICRO_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let txt_path = format!("{out_dir}/microphysics_report.txt");
    let json_path = format!("{out_dir}/microphysics_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[microphysics_inputs]").expect("write");
    writeln!(txt, "H0_km_s_mpc = {:.9}", a.h0_km_s_mpc).expect("write");
    writeln!(txt, "omega_b0 = {:.12}", a.omega_b0).expect("write");
    writeln!(txt, "omega_m0 = {:.12}", a.omega_m0).expect("write");
    writeln!(txt, "omega_r0 = {:.12}", a.omega_r0).expect("write");
    writeln!(txt, "omega_lambda0 = {:.12}", a.omega_lambda0).expect("write");
    writeln!(txt, "eta10 = {:.9}", a.eta10).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[microphysics_outputs]").expect("write");
    writeln!(txt, "Y_p_network = {:.9}", s.yp_network).expect("write");
    writeln!(txt, "D_H_network = {:.12e}", s.dh_network).expect("write");
    writeln!(txt, "He3_H_network = {:.12e}", s.he3h_network).expect("write");
    writeln!(
        txt,
        "BBN_freezeout_seconds = {:.6}",
        s.bbn_freezeout_seconds
    )
    .expect("write");
    writeln!(txt, "z_visibility_peak = {:.6}", s.z_visibility_peak).expect("write");
    writeln!(txt, "tau_recomb = {:.9e}", s.tau_recomb).expect("write");
    writeln!(txt, "x_e_final = {:.9e}", s.x_e_final).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[microphysics_gate]").expect("write");
    writeln!(txt, "yp_ok = {}", s.yp_ok).expect("write");
    writeln!(txt, "dh_ok = {}", s.dh_ok).expect("write");
    writeln!(txt, "recombination_ok = {}", s.recombination_ok).expect("write");
    writeln!(txt, "opacity_positive_ok = {}", s.opacity_positive_ok).expect("write");
    writeln!(txt, "passes_all = {}", s.passes_all()).expect("write");

    let mut json = File::create(&json_path).expect("create json");
    writeln!(
        json,
        "{{\n  \"inputs\": {{\"h0_km_s_mpc\": {:.9}, \"omega_b0\": {:.12}, \"omega_m0\": {:.12}, \"omega_r0\": {:.12}, \"omega_k0\": {:.12}, \"omega_lambda0\": {:.12}, \"eta10\": {:.9}}},\n  \"windows\": {{\"yp_abs_max\": {:.9}, \"dh_rel_max\": {:.9}, \"z_visibility_min\": {:.9}, \"z_visibility_max\": {:.9}}},\n  \"score\": {{\"yp_network\": {:.9}, \"dh_network\": {:.12e}, \"he3h_network\": {:.12e}, \"bbn_freezeout_seconds\": {:.9}, \"z_visibility_peak\": {:.9}, \"tau_recomb\": {:.9e}, \"x_e_final\": {:.9e}, \"yp_ok\": {}, \"dh_ok\": {}, \"recombination_ok\": {}, \"opacity_positive_ok\": {}, \"passes_all\": {}}}\n}}",
        a.h0_km_s_mpc,
        a.omega_b0,
        a.omega_m0,
        a.omega_r0,
        a.omega_k0,
        a.omega_lambda0,
        a.eta10,
        w.yp_abs_max,
        w.dh_rel_max,
        w.z_visibility_min,
        w.z_visibility_max,
        s.yp_network,
        s.dh_network,
        s.he3h_network,
        s.bbn_freezeout_seconds,
        s.z_visibility_peak,
        s.tau_recomb,
        s.x_e_final,
        s.yp_ok,
        s.dh_ok,
        s.recombination_ok,
        s.opacity_positive_ok,
        s.passes_all(),
    )
    .expect("write json");

    println!("wrote {txt_path}");
    println!("wrote {json_path}");
    println!(
        "Microphysics: Yp={:.5}, D/H={:.3e}, z_vis={:.1}, pass={}",
        s.yp_network,
        s.dh_network,
        s.z_visibility_peak,
        s.passes_all(),
    );
}
