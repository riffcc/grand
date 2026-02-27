//! GRAND-350 CI gate for assembled end-to-end universe simulation.

use gutoe_physics::{evaluate_universe_gate, UniverseAssumptions, UniverseWindows};
use std::fs::{self, File};
use std::io::Write;
use std::process;

fn main() {
    let assumptions = UniverseAssumptions::default();
    let windows = UniverseWindows::default();
    let score = evaluate_universe_gate(assumptions, windows);

    let out_dir =
        std::env::var("GUTOE_UNIVERSE_GATE_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let json_path = format!("{out_dir}/universe_ci_gate.json");
    let mut json = File::create(&json_path).expect("create gate json");

    writeln!(
        json,
        "{{\n  \"overall_pass\": {},\n  \"assumptions\": {{\"omega_r0\": {:.12}, \"omega_k0\": {:.12}, \"h0_ref_km_s_mpc\": {:.9}}},\n  \"windows\": {{\"h0_rel_error_max\": {:.12}, \"age_gyr_min\": {:.9}, \"age_gyr_max\": {:.9}, \"recombination_age_kyr_min\": {:.9}, \"recombination_age_kyr_max\": {:.9}, \"bbn_age_sec_min\": {:.9}, \"bbn_age_sec_max\": {:.9}}},\n  \"score\": {{\"h0_km_s_mpc\": {:.9}, \"h0_rel_error\": {:.12}, \"age_gyr\": {:.9}, \"recombination_age_kyr\": {:.9}, \"bbn_age_seconds\": {:.9}, \"rs_drag_mpc\": {:.9}, \"theta_star_rad\": {:.9e}, \"l_peak1\": {:.9}, \"l_peak2\": {:.9}, \"micro_yp\": {:.9}, \"micro_dh\": {:.12e}, \"micro_z_visibility_peak\": {:.9}, \"inflation_ok\": {}, \"baryogenesis_ok\": {}, \"bbn_ok\": {}, \"microphysics_ok\": {}, \"dark_matter_unified_ok\": {}, \"transfer_ok\": {}, \"h0_ok\": {}, \"age_ok\": {}, \"recombination_ok\": {}, \"bbn_timing_ok\": {}, \"passes_early_universe\": {}, \"passes_late_universe\": {}, \"passes_all\": {}}}\n}}",
        score.passes_all(),
        assumptions.omega_r0,
        assumptions.omega_k0,
        assumptions.h0_ref_km_s_mpc,
        windows.h0_rel_error_max,
        windows.age_gyr_min,
        windows.age_gyr_max,
        windows.recombination_age_kyr_min,
        windows.recombination_age_kyr_max,
        windows.bbn_age_sec_min,
        windows.bbn_age_sec_max,
        score.h0_km_s_mpc,
        score.h0_rel_error,
        score.age_gyr,
        score.recombination_age_kyr,
        score.bbn_age_seconds,
        score.transfer.rs_drag_mpc,
        score.transfer.theta_star_rad,
        score.transfer.l_peak1,
        score.transfer.l_peak2,
        score.microphysics.yp_network,
        score.microphysics.dh_network,
        score.microphysics.z_visibility_peak,
        score.inflation_ok,
        score.baryogenesis_ok,
        score.bbn_ok,
        score.microphysics_ok,
        score.dark_matter_unified_ok,
        score.transfer_ok,
        score.h0_ok,
        score.age_ok,
        score.recombination_ok,
        score.bbn_timing_ok,
        score.passes_early_universe(),
        score.passes_late_universe(),
        score.passes_all(),
    )
    .expect("write gate json");

    println!(
        "Universe gate: pass={} (H0={:.4}, age={:.4} Gyr, recomb={:.2} kyr, BBN={:.2} s, r_s={:.2} Mpc, l1={:.2}, micro Yp={:.4})",
        score.passes_all(),
        score.h0_km_s_mpc,
        score.age_gyr,
        score.recombination_age_kyr,
        score.bbn_age_seconds,
        score.transfer.rs_drag_mpc,
        score.transfer.l_peak1,
        score.microphysics.yp_network,
    );
    println!("wrote {json_path}");

    if !score.passes_all() {
        eprintln!(
            "FAIL: inflation={} baryogenesis={} bbn={} microphysics={} dark_unified={} transfer={} h0={} age={} recomb={} bbn_timing={}",
            score.inflation_ok,
            score.baryogenesis_ok,
            score.bbn_ok,
            score.microphysics_ok,
            score.dark_matter_unified_ok,
            score.transfer_ok,
            score.h0_ok,
            score.age_ok,
            score.recombination_ok,
            score.bbn_timing_ok,
        );
        process::exit(2);
    }
}
