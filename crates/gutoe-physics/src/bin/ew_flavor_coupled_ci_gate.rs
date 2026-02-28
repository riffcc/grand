//! Coupled electroweak + flavor closure gate.
//!
//! This gate tracks two closure targets together:
//! 1) sin²(theta_W) at M_Z
//! 2) neutrino oscillation splitting ratio Δm²32 / Δm²21

use gutoe_em::{
    neutrino_absolute_masses_from_texture, neutrino_hierarchy_prediction, sin2_weinberg,
};
use gutoe_physics::dynamics_map::StandardModelDynamicsMap;
use std::fs::{self, File};
use std::io::Write;
use std::process;

const SIN2_THETA_W_MZ_TARGET: f64 = 0.23122;
const SIN2_THETA_W_MZ_ABS_TOL: f64 = 5.0e-4;

const SOLAR_DM21_TARGET_EV2: f64 = 7.53e-5;
const ATMOSPHERIC_DM32_TARGET_EV2: f64 = 2.453e-3;
const SPLITTING_RATIO_TARGET: f64 = ATMOSPHERIC_DM32_TARGET_EV2 / SOLAR_DM21_TARGET_EV2;
const SPLITTING_RATIO_REL_TOL: f64 = 0.05;

fn rel_err(observed: f64, target: f64) -> f64 {
    if target.abs() < 1.0e-30 {
        0.0
    } else {
        (observed - target) / target
    }
}

fn main() {
    let out_dir = std::env::var("GUTOE_EW_FLAVOR_GATE_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let json_path = format!("{out_dir}/ew_flavor_coupled_ci_gate.json");

    let structural_sin2 = sin2_weinberg();
    let map = StandardModelDynamicsMap::from_clifford_z3();
    let shift = map.sin2_theta_w_mz_shift_structural();
    let mz_sin2 = map.sin2_theta_w_at_mz();
    let mz_abs_err = (mz_sin2 - SIN2_THETA_W_MZ_TARGET).abs();
    let mz_signed_err = mz_sin2 - SIN2_THETA_W_MZ_TARGET;

    let hierarchy = neutrino_hierarchy_prediction();
    let abs = neutrino_absolute_masses_from_texture();
    let m1_ev = abs.m1_ev;
    let m2_ev = abs.m2_ev;
    let m3_ev = abs.m3_ev;
    let dm21_ev2 = abs.dm21_ev2;
    let dm32_ev2 = abs.dm32_ev2;
    let ratio = abs.splitting_ratio_32_over_21;
    let ratio_rel_err = rel_err(ratio, SPLITTING_RATIO_TARGET);
    let clifford_half_dim = map.clifford_dim as f64 / 2.0;

    let ew_ok = mz_abs_err <= SIN2_THETA_W_MZ_ABS_TOL;
    let ratio_ok = ratio_rel_err.abs() <= SPLITTING_RATIO_REL_TOL;
    let hierarchy_ok = hierarchy == "normal";
    let ordering_ok = dm21_ev2 > 0.0 && dm32_ev2 > 0.0;

    let overall_pass = ew_ok && ratio_ok && hierarchy_ok && ordering_ok;

    let mut json = File::create(&json_path).expect("create ew+flavor gate json");
    writeln!(
        json,
        "{{\n  \"overall_pass\": {},\n  \"windows\": {{\"sin2_theta_w_mz_target\": {:.12}, \"sin2_theta_w_mz_abs_tol\": {:.12e}, \"splitting_ratio_target\": {:.12}, \"splitting_ratio_rel_tol\": {:.6}}},\n  \"electroweak\": {{\"sin2_structural\": {:.12}, \"alpha\": {:.12}, \"clifford_half_dim\": {:.12}, \"shift_coeff\": {:.12}, \"delta_sin2\": {:.12e}, \"sin2_mz_bridge\": {:.12}, \"mz_abs_err\": {:.12e}, \"mz_signed_err\": {:.12e}}},\n  \"flavor\": {{\"hierarchy\": \"{}\", \"hierarchy_exponent\": {:.12}, \"m1_ev\": {:.12e}, \"m2_ev\": {:.12e}, \"m3_ev\": {:.12e}, \"dm21_ev2\": {:.12e}, \"dm32_ev2\": {:.12e}, \"splitting_ratio\": {:.12}, \"splitting_ratio_rel_err\": {:.12e}}},\n  \"checks\": {{\"ew_ok\": {}, \"ratio_ok\": {}, \"hierarchy_ok\": {}, \"ordering_ok\": {}}}\n}}",
        if overall_pass { "true" } else { "false" },
        SIN2_THETA_W_MZ_TARGET,
        SIN2_THETA_W_MZ_ABS_TOL,
        SPLITTING_RATIO_TARGET,
        SPLITTING_RATIO_REL_TOL,
        structural_sin2,
        map.alpha_leading_order,
        clifford_half_dim,
        clifford_half_dim,
        shift,
        mz_sin2,
        mz_abs_err,
        mz_signed_err,
        hierarchy,
        abs.hierarchy_exponent,
        m1_ev,
        m2_ev,
        m3_ev,
        dm21_ev2,
        dm32_ev2,
        ratio,
        ratio_rel_err,
        if ew_ok { "true" } else { "false" },
        if ratio_ok { "true" } else { "false" },
        if hierarchy_ok { "true" } else { "false" },
        if ordering_ok { "true" } else { "false" }
    )
    .expect("write ew+flavor gate json");

    println!(
        "ew_flavor_coupled_ci_gate: pass={} sin2_mz={:.9} ratio={:.6} hierarchy={}",
        overall_pass, mz_sin2, ratio, hierarchy
    );
    println!("wrote {json_path}");

    if !overall_pass {
        eprintln!(
            "FAIL: ew_ok={} ratio_ok={} hierarchy_ok={} ordering_ok={} mz_abs_err={:.3e} ratio_rel_err={:.3e}",
            ew_ok, ratio_ok, hierarchy_ok, ordering_ok, mz_abs_err, ratio_rel_err
        );
        process::exit(2);
    }
}
