//! GRAND-352 CI gate for explicit BBN/recombination microphysics lane.

use gutoe_physics::{evaluate_microphysics_gate, MicrophysicsAssumptions, MicrophysicsWindows};
use std::fs::{self, File};
use std::io::Write;
use std::process;

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
        std::env::var("GUTOE_MICRO_GATE_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let json_path = format!("{out_dir}/microphysics_ci_gate.json");
    let mut json = File::create(&json_path).expect("create gate json");

    writeln!(
        json,
        "{{\n  \"overall_pass\": {},\n  \"windows\": {{\"yp_abs_max\": {:.9}, \"dh_rel_max\": {:.9}, \"z_visibility_min\": {:.9}, \"z_visibility_max\": {:.9}}},\n  \"score\": {{\"yp_network\": {:.9}, \"dh_network\": {:.12e}, \"z_visibility_peak\": {:.9}, \"tau_recomb\": {:.9e}, \"x_e_final\": {:.9e}, \"yp_ok\": {}, \"dh_ok\": {}, \"recombination_ok\": {}, \"opacity_positive_ok\": {}, \"passes_all\": {}}}\n}}",
        s.passes_all(),
        w.yp_abs_max,
        w.dh_rel_max,
        w.z_visibility_min,
        w.z_visibility_max,
        s.yp_network,
        s.dh_network,
        s.z_visibility_peak,
        s.tau_recomb,
        s.x_e_final,
        s.yp_ok,
        s.dh_ok,
        s.recombination_ok,
        s.opacity_positive_ok,
        s.passes_all(),
    )
    .expect("write gate json");

    println!(
        "Microphysics gate: pass={} (Yp={:.5}, D/H={:.3e}, z_vis={:.1})",
        s.passes_all(),
        s.yp_network,
        s.dh_network,
        s.z_visibility_peak,
    );
    println!("wrote {json_path}");

    if !s.passes_all() {
        eprintln!(
            "FAIL: yp_ok={} dh_ok={} recombination_ok={} opacity_positive_ok={}",
            s.yp_ok, s.dh_ok, s.recombination_ok, s.opacity_positive_ok,
        );
        process::exit(2);
    }
}
