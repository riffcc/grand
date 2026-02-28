//! Multiplicative triangulation lane for coupled EW + neutrino anchors.
//!
//! Goal:
//! - solve hidden multiplicative factors from independent observables
//! - keep the solution in log-space where possible
//! - report structural-vs-required deltas without per-lane retuning

use gutoe_em::{
    neutrino_absolute_masses_from_texture, neutrino_hierarchy_exponent_structural,
    triangulate_ew_shift_for_target, triangulate_neutrino_from_splittings,
};
use std::fs::{self, File};
use std::io::Write;

const SOLAR_DM21_TARGET_EV2: f64 = 7.53e-5;
const ATMOSPHERIC_DM32_TARGET_EV2: f64 = 2.453e-3;
const SIN2_THETA_W_MZ_TARGET: f64 = 0.23122;
const KAPPA_STRUCTURAL: f64 = 60.0 / 11.0;

fn main() {
    let out_dir = std::env::var("GUTOE_TRIANGULATE_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let txt_path = format!("{out_dir}/triangulate_params_report.txt");
    let json_path = format!("{out_dir}/triangulate_params_report.json");

    let tri = triangulate_neutrino_from_splittings(SOLAR_DM21_TARGET_EV2, ATMOSPHERIC_DM32_TARGET_EV2);
    let ew = triangulate_ew_shift_for_target(SIN2_THETA_W_MZ_TARGET);

    let p_structural = neutrino_hierarchy_exponent_structural();
    let abs_current = neutrino_absolute_masses_from_texture();
    let ratio_struct = abs_current.splitting_ratio_32_over_21;

    let kappa_current = {
        let alpha = abs_current.alpha_physical;
        let alpha4 = alpha.powi(4);
        abs_current.mass_scale_ev / (abs_current.electron_mass_anchor_ev * alpha4)
    };

    let dm21_struct = abs_current.dm21_ev2;
    let dm32_struct = abs_current.dm32_ev2;

    let mut txt = File::create(&txt_path).expect("create triangulation txt");
    writeln!(txt, "triangulation lane: multiplicative/log-space anchors").expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[anchors]").expect("write");
    writeln!(txt, "solar_dm21_target_ev2 = {:.12e}", SOLAR_DM21_TARGET_EV2).expect("write");
    writeln!(
        txt,
        "atmospheric_dm32_target_ev2 = {:.12e}",
        ATMOSPHERIC_DM32_TARGET_EV2
    )
    .expect("write");
    writeln!(txt, "splitting_ratio_target = {:.12}", tri.ratio_target).expect("write");
    writeln!(txt, "sin2_theta_w_mz_target = {:.12}", SIN2_THETA_W_MZ_TARGET).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[texture_eigenvalue_ratios]").expect("write");
    writeln!(txt, "r1 = {:.12e}", tri.r1).expect("write");
    writeln!(txt, "r2 = {:.12e}", tri.r2).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[triangulated_parameters]").expect("write");
    writeln!(txt, "p_ratio = {:.12}", tri.p_triangulated).expect("write");
    writeln!(txt, "kappa_dm21 = {:.12}", tri.kappa_dm21).expect("write");
    writeln!(txt, "kappa_dm32 = {:.12}", tri.kappa_dm32).expect("write");
    writeln!(txt, "kappa_geo = {:.12}", tri.kappa_geo).expect("write");
    writeln!(txt, "ew_coeff_required = {:.12}", ew.coeff_required).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[structural_reference]").expect("write");
    writeln!(txt, "p_structural = {:.12}", p_structural).expect("write");
    writeln!(txt, "kappa_structural = {:.12}", KAPPA_STRUCTURAL).expect("write");
    writeln!(txt, "kappa_current_runtime = {:.12}", kappa_current).expect("write");
    writeln!(txt, "ew_coeff_structural = {:.12}", ew.coeff_structural).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[residuals]").expect("write");
    writeln!(txt, "ratio_fit = {:.12}", tri.ratio_fit).expect("write");
    writeln!(txt, "ratio_fit_rel_err = {:.12e}", tri.ratio_fit_rel_err).expect("write");
    writeln!(
        txt,
        "kappa_consistency_rel = {:.12e}",
        tri.kappa_consistency_rel
    )
    .expect("write");
    writeln!(
        txt,
        "kappa_vs_structural_rel = {:.12e}",
        (tri.kappa_geo - KAPPA_STRUCTURAL) / KAPPA_STRUCTURAL
    )
    .expect("write");
    writeln!(txt, "ew_coeff_delta_rel = {:.12e}", ew.coeff_rel_delta).expect("write");
    writeln!(txt, "sin2_mz_structural = {:.12}", ew.sin2_structural_mz).expect("write");
    writeln!(
        txt,
        "sin2_mz_structural_abs_err = {:.12e}",
        ew.sin2_structural_abs_err
    )
    .expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[reconstructed_from_triangulated]").expect("write");
    writeln!(txt, "m1_ev = {:.12e}", tri.m1_ev).expect("write");
    writeln!(txt, "m2_ev = {:.12e}", tri.m2_ev).expect("write");
    writeln!(txt, "m3_ev = {:.12e}", tri.m3_ev).expect("write");
    writeln!(txt, "dm21_ev2 = {:.12e}", tri.dm21_ev2).expect("write");
    writeln!(txt, "dm32_ev2 = {:.12e}", tri.dm32_ev2).expect("write");
    writeln!(
        txt,
        "ratio = {:.12}",
        tri.dm32_ev2 / tri.dm21_ev2.max(1.0e-30)
    )
    .expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[current_runtime_reference]").expect("write");
    writeln!(txt, "runtime_dm21_ev2 = {:.12e}", dm21_struct).expect("write");
    writeln!(txt, "runtime_dm32_ev2 = {:.12e}", dm32_struct).expect("write");
    writeln!(txt, "runtime_ratio = {:.12}", ratio_struct).expect("write");
    writeln!(
        txt,
        "runtime_ratio_rel_err = {:.12e}",
        (ratio_struct - tri.ratio_target) / tri.ratio_target
    )
    .expect("write");

    let mut json = File::create(&json_path).expect("create triangulation json");
    writeln!(
        json,
        "{{\n  \"anchors\": {{\n    \"dm21_target_ev2\": {:.12e},\n    \"dm32_target_ev2\": {:.12e},\n    \"ratio_target\": {:.12},\n    \"sin2_theta_w_mz_target\": {:.12}\n  }},\n  \"texture\": {{\"r1\": {:.12e}, \"r2\": {:.12e}}},\n  \"triangulated\": {{\n    \"p_ratio\": {:.12},\n    \"kappa_dm21\": {:.12},\n    \"kappa_dm32\": {:.12},\n    \"kappa_geo\": {:.12},\n    \"log_kappa_geo\": {:.12},\n    \"ew_coeff_required\": {:.12},\n    \"log_ew_coeff_required\": {:.12}\n  }},\n  \"structural\": {{\n    \"p\": {:.12},\n    \"kappa_expected\": {:.12},\n    \"kappa_current_runtime\": {:.12},\n    \"ew_coeff\": {:.12},\n    \"sin2_structural\": {:.12},\n    \"sin2_shift_structural\": {:.12e},\n    \"sin2_mz_structural\": {:.12}\n  }},\n  \"residuals\": {{\n    \"ratio_fit\": {:.12},\n    \"ratio_fit_rel_err\": {:.12e},\n    \"kappa_consistency_rel\": {:.12e},\n    \"kappa_vs_structural_rel\": {:.12e},\n    \"ew_coeff_delta_rel\": {:.12e},\n    \"sin2_mz_structural_abs_err\": {:.12e}\n  }},\n  \"reconstructed\": {{\n    \"m1_ev\": {:.12e},\n    \"m2_ev\": {:.12e},\n    \"m3_ev\": {:.12e},\n    \"dm21_ev2\": {:.12e},\n    \"dm32_ev2\": {:.12e},\n    \"ratio\": {:.12}\n  }},\n  \"runtime_reference\": {{\n    \"dm21_ev2\": {:.12e},\n    \"dm32_ev2\": {:.12e},\n    \"ratio\": {:.12}\n  }}\n}}",
        SOLAR_DM21_TARGET_EV2,
        ATMOSPHERIC_DM32_TARGET_EV2,
        tri.ratio_target,
        SIN2_THETA_W_MZ_TARGET,
        tri.r1,
        tri.r2,
        tri.p_triangulated,
        tri.kappa_dm21,
        tri.kappa_dm32,
        tri.kappa_geo,
        tri.kappa_geo.ln(),
        ew.coeff_required,
        ew.coeff_required.ln(),
        p_structural,
        KAPPA_STRUCTURAL,
        kappa_current,
        ew.coeff_structural,
        ew.sin2_structural,
        ew.shift_structural,
        ew.sin2_structural_mz,
        tri.ratio_fit,
        tri.ratio_fit_rel_err,
        tri.kappa_consistency_rel,
        (tri.kappa_geo - KAPPA_STRUCTURAL) / KAPPA_STRUCTURAL,
        ew.coeff_rel_delta,
        ew.sin2_structural_abs_err,
        tri.m1_ev,
        tri.m2_ev,
        tri.m3_ev,
        tri.dm21_ev2,
        tri.dm32_ev2,
        tri.dm32_ev2 / tri.dm21_ev2.max(1.0e-30),
        dm21_struct,
        dm32_struct,
        ratio_struct
    )
    .expect("write triangulation json");

    println!("wrote {txt_path}");
    println!("wrote {json_path}");
    println!(
        "triangulated p={:.6}, kappa={:.6}, ew_coeff={:.6}",
        tri.p_triangulated, tri.kappa_geo, ew.coeff_required
    );
}
