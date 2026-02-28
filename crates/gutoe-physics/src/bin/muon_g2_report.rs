//! Muon g-2 anomaly report from shared GUTOE primitives.
//!
//! This lane targets the unresolved discrepancy `Δa_μ = a_μ(exp) - a_μ(SM)`.
//! Structural candidate:
//!   Δa_μ,cand = α^3 / (N_gauge * N_complement)
//! where
//!   N_gauge = 8 + 3 + 1 = 12
//!   N_complement = 2^4 - 3 = 13

use gutoe_physics::constants::{ALPHA, ALPHA_LEADING_ORDER};
use gutoe_physics::StandardModelDynamicsMap;
use std::f64::consts::PI;
use std::fs::{self, File};
use std::io::Write;

fn env_f64(name: &str, default: f64) -> f64 {
    std::env::var(name)
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(default)
}

fn main() {
    let out_dir =
        std::env::var("GUTOE_MUON_G2_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);

    // Reference values used for this report lane (can be overridden).
    let a_mu_exp = env_f64("GUTOE_AMU_EXP_REF", 0.001_165_920_59);
    let a_mu_sm = env_f64("GUTOE_AMU_SM_REF", 0.001_165_918_10);
    let delta_ref = a_mu_exp - a_mu_sm;

    let sm = StandardModelDynamicsMap::from_clifford_z3();
    let n_gauge = sm.total_gauge_generators as f64; // 12
    let n_complement = (sm.clifford_dim - sm.magnetic_triplet_card) as f64; // 13
    let denom = n_gauge * n_complement; // 156

    // Candidate unresolved-gap term from shared counts.
    let delta_cand_phys_alpha = ALPHA.powi(3) / denom;
    let delta_cand_struct_alpha = ALPHA_LEADING_ORDER.powi(3) / denom;

    // Standard first-order Schwinger term (for scale context).
    let schwinger_phys = ALPHA / (2.0 * PI);
    let schwinger_struct = ALPHA_LEADING_ORDER / (2.0 * PI);

    let a_mu_pred_phys = a_mu_sm + delta_cand_phys_alpha;
    let a_mu_pred_struct = a_mu_sm + delta_cand_struct_alpha;

    let delta_abs_err_phys = (delta_cand_phys_alpha - delta_ref).abs();
    let delta_abs_err_struct = (delta_cand_struct_alpha - delta_ref).abs();
    let delta_rel_err_phys = if delta_ref != 0.0 {
        delta_abs_err_phys / delta_ref.abs()
    } else {
        f64::NAN
    };
    let delta_rel_err_struct = if delta_ref != 0.0 {
        delta_abs_err_struct / delta_ref.abs()
    } else {
        f64::NAN
    };

    let txt_path = format!("{out_dir}/muon_g2_report.txt");
    let json_path = format!("{out_dir}/muon_g2_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "a_mu_exp_ref = {:.12e}", a_mu_exp).expect("write");
    writeln!(txt, "a_mu_sm_ref = {:.12e}", a_mu_sm).expect("write");
    writeln!(txt, "delta_ref = {:.12e}", delta_ref).expect("write");
    writeln!(txt, "n_gauge = {:.0}", n_gauge).expect("write");
    writeln!(txt, "n_complement = {:.0}", n_complement).expect("write");
    writeln!(txt, "denominator = {:.0}", denom).expect("write");
    writeln!(txt, "alpha_physical = {:.15e}", ALPHA).expect("write");
    writeln!(txt, "alpha_structural = {:.15e}", ALPHA_LEADING_ORDER).expect("write");
    writeln!(txt, "schwinger_phys_alpha = {:.12e}", schwinger_phys).expect("write");
    writeln!(txt, "schwinger_struct_alpha = {:.12e}", schwinger_struct).expect("write");
    writeln!(
        txt,
        "delta_candidate_phys_alpha = {:.12e}",
        delta_cand_phys_alpha
    )
    .expect("write");
    writeln!(
        txt,
        "delta_candidate_struct_alpha = {:.12e}",
        delta_cand_struct_alpha
    )
    .expect("write");
    writeln!(txt, "delta_abs_err_phys = {:.12e}", delta_abs_err_phys).expect("write");
    writeln!(
        txt,
        "delta_abs_err_struct = {:.12e}",
        delta_abs_err_struct
    )
    .expect("write");
    writeln!(txt, "delta_rel_err_phys = {:.12e}", delta_rel_err_phys).expect("write");
    writeln!(
        txt,
        "delta_rel_err_struct = {:.12e}",
        delta_rel_err_struct
    )
    .expect("write");
    writeln!(txt, "a_mu_pred_phys = {:.12e}", a_mu_pred_phys).expect("write");
    writeln!(txt, "a_mu_pred_struct = {:.12e}", a_mu_pred_struct).expect("write");

    let mut json = File::create(&json_path).expect("create json");
    writeln!(
        json,
        "{{\n  \"references\": {{\"a_mu_exp\": {:.12e}, \"a_mu_sm\": {:.12e}, \"delta_ref\": {:.12e}}},\n  \"structural_counts\": {{\"n_gauge\": {:.0}, \"n_complement\": {:.0}, \"denominator\": {:.0}}},\n  \"alphas\": {{\"physical\": {:.15e}, \"structural\": {:.15e}}},\n  \"schwinger\": {{\"physical_alpha\": {:.12e}, \"structural_alpha\": {:.12e}}},\n  \"delta_candidate\": {{\"physical_alpha\": {:.12e}, \"structural_alpha\": {:.12e}}},\n  \"delta_error\": {{\"abs_phys\": {:.12e}, \"abs_struct\": {:.12e}, \"rel_phys\": {:.12e}, \"rel_struct\": {:.12e}}},\n  \"a_mu_pred\": {{\"physical_alpha\": {:.12e}, \"structural_alpha\": {:.12e}}}\n}}",
        a_mu_exp,
        a_mu_sm,
        delta_ref,
        n_gauge,
        n_complement,
        denom,
        ALPHA,
        ALPHA_LEADING_ORDER,
        schwinger_phys,
        schwinger_struct,
        delta_cand_phys_alpha,
        delta_cand_struct_alpha,
        delta_abs_err_phys,
        delta_abs_err_struct,
        delta_rel_err_phys,
        delta_rel_err_struct,
        a_mu_pred_phys,
        a_mu_pred_struct
    )
    .expect("write json");

    println!("wrote {txt_path}");
    println!("wrote {json_path}");
    println!(
        "Δa_mu ref={:.6e}, candidate={:.6e} (phys α, rel err {:.3e})",
        delta_ref, delta_cand_phys_alpha, delta_rel_err_phys
    );
}

