//! CI gate for multiplicative triangulation stability.
//!
//! This gate checks that the forced-parameter solver remains numerically stable
//! and that structural-vs-required deltas are explicitly surfaced.

use gutoe_em::{
    neutrino_hierarchy_exponent_structural, triangulate_ew_shift_for_target,
    triangulate_neutrino_from_splittings,
};
use std::fs::{self, File};
use std::io::Write;
use std::process;

const SOLAR_DM21_TARGET_EV2: f64 = 7.53e-5;
const ATMOSPHERIC_DM32_TARGET_EV2: f64 = 2.453e-3;
const SIN2_THETA_W_MZ_TARGET: f64 = 0.23122;
const KAPPA_STRUCTURAL: f64 = 60.0 / 11.0;

const RATIO_FIT_REL_MAX: f64 = 1.0e-9;
const ABS_SPLIT_REL_MAX: f64 = 1.0e-9;
const P_STRUCT_REL_MAX: f64 = 5.0e-3; // 0.5%

fn rel_err(observed: f64, target: f64) -> f64 {
    if target.abs() < 1.0e-30 {
        0.0
    } else {
        (observed - target) / target
    }
}

fn main() {
    let out_dir = std::env::var("GUTOE_TRIANGULATION_GATE_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let json_path = format!("{out_dir}/triangulation_ci_gate.json");

    let tri = triangulate_neutrino_from_splittings(SOLAR_DM21_TARGET_EV2, ATMOSPHERIC_DM32_TARGET_EV2);
    let ew = triangulate_ew_shift_for_target(SIN2_THETA_W_MZ_TARGET);
    let p_struct = neutrino_hierarchy_exponent_structural();

    let dm21_rel = rel_err(tri.dm21_ev2, SOLAR_DM21_TARGET_EV2);
    let dm32_rel = rel_err(tri.dm32_ev2, ATMOSPHERIC_DM32_TARGET_EV2);
    let p_struct_rel = rel_err(tri.p_triangulated, p_struct);
    let kappa_vs_struct = rel_err(tri.kappa_geo, KAPPA_STRUCTURAL);

    let ratio_ok = tri.ratio_fit_rel_err.abs() <= RATIO_FIT_REL_MAX;
    let abs_ok = dm21_rel.abs() <= ABS_SPLIT_REL_MAX && dm32_rel.abs() <= ABS_SPLIT_REL_MAX;
    let p_ok = p_struct_rel.abs() <= P_STRUCT_REL_MAX;

    // This gate validates solver integrity, not full structural closure.
    let overall_pass = ratio_ok && abs_ok && p_ok;

    let mut json = File::create(&json_path).expect("create triangulation gate json");
    writeln!(
        json,
        "{{\n  \"overall_pass\": {},\n  \"windows\": {{\"ratio_fit_rel_max\": {:.3e}, \"abs_split_rel_max\": {:.3e}, \"p_struct_rel_max\": {:.6}}},\n  \"triangulated\": {{\"p\": {:.12}, \"ratio_fit\": {:.12}, \"ratio_fit_rel_err\": {:.12e}, \"kappa\": {:.12}, \"dm21_ev2\": {:.12e}, \"dm32_ev2\": {:.12e}}},\n  \"structural\": {{\"p\": {:.12}, \"kappa\": {:.12}, \"ew_coeff\": {:.12}}},\n  \"residuals\": {{\"p_struct_rel\": {:.12e}, \"kappa_vs_struct_rel\": {:.12e}, \"dm21_rel\": {:.12e}, \"dm32_rel\": {:.12e}, \"ew_coeff_rel_delta\": {:.12e}}},\n  \"checks\": {{\"ratio_ok\": {}, \"abs_ok\": {}, \"p_ok\": {}}}\n}}",
        if overall_pass { "true" } else { "false" },
        RATIO_FIT_REL_MAX,
        ABS_SPLIT_REL_MAX,
        P_STRUCT_REL_MAX,
        tri.p_triangulated,
        tri.ratio_fit,
        tri.ratio_fit_rel_err,
        tri.kappa_geo,
        tri.dm21_ev2,
        tri.dm32_ev2,
        p_struct,
        KAPPA_STRUCTURAL,
        ew.coeff_structural,
        p_struct_rel,
        kappa_vs_struct,
        dm21_rel,
        dm32_rel,
        ew.coeff_rel_delta,
        if ratio_ok { "true" } else { "false" },
        if abs_ok { "true" } else { "false" },
        if p_ok { "true" } else { "false" }
    )
    .expect("write triangulation gate json");

    println!(
        "triangulation_ci_gate: pass={} p={:.6} kappa={:.6} ew_coeff_delta={:.3e}",
        overall_pass, tri.p_triangulated, tri.kappa_geo, ew.coeff_rel_delta
    );
    println!("wrote {json_path}");

    if !overall_pass {
        eprintln!(
            "FAIL: ratio_ok={} abs_ok={} p_ok={} ratio_rel={:.3e} dm21_rel={:.3e} dm32_rel={:.3e} p_struct_rel={:.3e}",
            ratio_ok, abs_ok, p_ok, tri.ratio_fit_rel_err, dm21_rel, dm32_rel, p_struct_rel
        );
        process::exit(2);
    }
}
