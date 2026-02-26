//! Structural cosmological constant report from Clifford-derived suppression.

use gutoe_physics::constants::{
    lambda_cosmological_full_candidate, lambda_cosmological_signature_candidate,
    lambda_cosmological_structural, lambda_cosmological_suppression, lambda_micro_finite_mode_rescale,
    lambda_micro_mode_count, lambda_micro_mode_count_from_ternary_depth,
    lorentz_signature_factor_from_bivector_split, ALPHA_INV_LEADING_ORDER,
    BIVECTOR_TIMELIKE_SPACELIKE_COUNT, BIVECTOR_TOTAL_COUNT, EWSB_SCALE_FACTOR_STRUCTURAL,
    HIGGS_QUARTIC_STRUCTURAL, LAMBDA_COSMOLOGICAL_OBSERVED, PLANCK_LENGTH, Z3_FIXED_GRADE1_COUNT,
};
use std::f64::consts::SQRT_2;
use std::fs::{self, File};
use std::io::Write;

fn main() {
    let suppression = lambda_cosmological_suppression();
    let lambda_struct = lambda_cosmological_structural();
    let lambda_signature = lambda_cosmological_signature_candidate();
    let lambda_full = lambda_cosmological_full_candidate();
    let micro_mode_count = lambda_micro_mode_count();
    let micro_mode_count_ternary = lambda_micro_mode_count_from_ternary_depth();
    let micro_rescale = lambda_micro_finite_mode_rescale();
    let k_required = lambda_struct / LAMBDA_COSMOLOGICAL_OBSERVED;
    let k_signature = lorentz_signature_factor_from_bivector_split();
    let ratio_struct = lambda_struct / LAMBDA_COSMOLOGICAL_OBSERVED;
    let rel_err_struct = (lambda_struct - LAMBDA_COSMOLOGICAL_OBSERVED).abs() / LAMBDA_COSMOLOGICAL_OBSERVED;
    let ratio_signature = lambda_signature / LAMBDA_COSMOLOGICAL_OBSERVED;
    let rel_err_signature =
        (lambda_signature - LAMBDA_COSMOLOGICAL_OBSERVED).abs() / LAMBDA_COSMOLOGICAL_OBSERVED;
    let ratio_full = lambda_full / LAMBDA_COSMOLOGICAL_OBSERVED;
    let rel_err_full = (lambda_full - LAMBDA_COSMOLOGICAL_OBSERVED).abs() / LAMBDA_COSMOLOGICAL_OBSERVED;

    let out_dir = "/tmp/bh_renders";
    let _ = fs::create_dir_all(out_dir);
    let txt_path = format!("{out_dir}/lambda_cosmological_report.txt");
    let json_path = format!("{out_dir}/lambda_cosmological_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[structural_inputs]").expect("write");
    writeln!(txt, "planck_length_m = {:.12e}", PLANCK_LENGTH).expect("write");
    writeln!(txt, "higgs_quartic = {:.12}", HIGGS_QUARTIC_STRUCTURAL).expect("write");
    writeln!(txt, "alpha_inv_lo = {}", ALPHA_INV_LEADING_ORDER).expect("write");
    writeln!(txt, "suppression = {:.12e}", suppression).expect("write");
    writeln!(txt, "sqrt2 = {:.12}", SQRT_2).expect("write");
    writeln!(txt, "bivector_total = {:.0}", BIVECTOR_TOTAL_COUNT).expect("write");
    writeln!(txt, "bivector_timelike_spacelike = {:.0}", BIVECTOR_TIMELIKE_SPACELIKE_COUNT).expect("write");
    writeln!(txt, "k_signature = {:.12}", k_signature).expect("write");
    writeln!(txt, "k_required = {:.12}", k_required).expect("write");
    writeln!(txt, "k_required_over_k_signature = {:.12}", k_required / k_signature).expect("write");
    writeln!(txt, "ewsb_scale_factor = {:.0}", EWSB_SCALE_FACTOR_STRUCTURAL).expect("write");
    writeln!(txt, "z3_fixed_grade1_count = {:.0}", Z3_FIXED_GRADE1_COUNT).expect("write");
    writeln!(txt, "micro_mode_count = {:.0}", micro_mode_count).expect("write");
    writeln!(txt, "micro_mode_count_ternary = {:.0}", micro_mode_count_ternary).expect("write");
    writeln!(txt, "micro_count_ratio_ewsb_over_ternary = {:.12}", micro_mode_count / micro_mode_count_ternary)
        .expect("write");
    writeln!(txt, "micro_rescale = {:.12}", micro_rescale).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[lambda_cosmological_structural]").expect("write");
    writeln!(txt, "lambda_structural = {:.12e}", lambda_struct).expect("write");
    writeln!(txt, "lambda_observed = {:.12e}", LAMBDA_COSMOLOGICAL_OBSERVED).expect("write");
    writeln!(txt, "ratio_struct_over_obs = {:.12}", ratio_struct).expect("write");
    writeln!(txt, "relative_error = {:.12}", rel_err_struct).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[lambda_cosmological_signature_candidate]").expect("write");
    writeln!(txt, "lambda_signature = {:.12e}", lambda_signature).expect("write");
    writeln!(txt, "ratio_signature_over_obs = {:.12}", ratio_signature).expect("write");
    writeln!(txt, "relative_error = {:.12}", rel_err_signature).expect("write");
    writeln!(txt, "residual_over_sqrt2 = {:.12}", ratio_struct / SQRT_2).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[lambda_cosmological_full_candidate]").expect("write");
    writeln!(txt, "lambda_full = {:.12e}", lambda_full).expect("write");
    writeln!(txt, "ratio_full_over_obs = {:.12}", ratio_full).expect("write");
    writeln!(txt, "relative_error = {:.12}", rel_err_full).expect("write");

    let mut json = File::create(&json_path).expect("create json");
    writeln!(
        json,
        "{{\n  \"planck_length_m\": {:.12e},\n  \"higgs_quartic\": {:.12},\n  \"alpha_inv_lo\": {},\n  \"suppression\": {:.12e},\n  \"sqrt2\": {:.12},\n  \"bivector_total\": {:.0},\n  \"bivector_timelike_spacelike\": {:.0},\n  \"k_signature\": {:.12},\n  \"k_required\": {:.12},\n  \"k_required_over_k_signature\": {:.12},\n  \"ewsb_scale_factor\": {:.0},\n  \"z3_fixed_grade1_count\": {:.0},\n  \"micro_mode_count\": {:.0},\n  \"micro_mode_count_ternary\": {:.0},\n  \"micro_count_ratio_ewsb_over_ternary\": {:.12},\n  \"micro_rescale\": {:.12},\n  \"lambda_structural\": {:.12e},\n  \"lambda_signature_candidate\": {:.12e},\n  \"lambda_full_candidate\": {:.12e},\n  \"lambda_observed\": {:.12e},\n  \"ratio_struct_over_obs\": {:.12},\n  \"ratio_signature_over_obs\": {:.12},\n  \"ratio_full_over_obs\": {:.12},\n  \"residual_over_sqrt2\": {:.12},\n  \"relative_error_structural\": {:.12},\n  \"relative_error_signature_candidate\": {:.12},\n  \"relative_error_full_candidate\": {:.12}\n}}",
        PLANCK_LENGTH,
        HIGGS_QUARTIC_STRUCTURAL,
        ALPHA_INV_LEADING_ORDER,
        suppression,
        SQRT_2,
        BIVECTOR_TOTAL_COUNT,
        BIVECTOR_TIMELIKE_SPACELIKE_COUNT,
        k_signature,
        k_required,
        k_required / k_signature,
        EWSB_SCALE_FACTOR_STRUCTURAL,
        Z3_FIXED_GRADE1_COUNT,
        micro_mode_count,
        micro_mode_count_ternary,
        micro_mode_count / micro_mode_count_ternary,
        micro_rescale,
        lambda_struct,
        lambda_signature,
        lambda_full,
        LAMBDA_COSMOLOGICAL_OBSERVED,
        ratio_struct,
        ratio_signature,
        ratio_full,
        ratio_struct / SQRT_2,
        rel_err_struct,
        rel_err_signature,
        rel_err_full
    )
    .expect("write json");

    println!("wrote {txt_path}");
    println!("wrote {json_path}");
    println!(
        "Λ_struct={:.6e}, Λ_obs={:.6e}, ratio={:.4}, rel_err={:.4}",
        lambda_struct, LAMBDA_COSMOLOGICAL_OBSERVED, ratio_struct, rel_err_struct
    );
    println!(
        "Λ_sig=Λ_struct/sqrt2={:.6e}, ratio={:.4}, rel_err={:.4}",
        lambda_signature, ratio_signature, rel_err_signature
    );
    println!(
        "Λ_full=Λ_sig*(486/485)={:.6e}, ratio={:.6}, rel_err={:.6}",
        lambda_full, ratio_full, rel_err_full
    );
}
