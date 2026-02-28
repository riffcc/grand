use gutoe_physics::StandardModelDynamicsMap;
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

// PDG-like central defaults for reporting (explicitly assumptions, not derived).
const MZ_GEV_DEFAULT: f64 = 91.1876;
const ALPHA_S_MZ_DEFAULT: f64 = 0.1181;
const PROTON_MASS_MEV_OBS: f64 = 938.272_088_16;
const ALPHA_INV_STRUCTURAL: f64 = 137.0;
const CLIFFORD_DIM: f64 = 16.0;
const C_INF_STRUCTURAL: f64 = 67.0 / 66.0;

fn env_f64(name: &str, default: f64) -> f64 {
    std::env::var(name)
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(default)
}

fn env_u32(name: &str, default: u32) -> u32 {
    std::env::var(name)
        .ok()
        .and_then(|s| s.parse::<u32>().ok())
        .unwrap_or(default)
}

/// One-loop running in the Lean-compatible form:
/// alpha_s(Q) = 2π / (β0 * ln(Q/Λ)).
fn alpha_s_one_loop(beta0: f64, q_gev: f64, lambda_qcd_gev: f64) -> f64 {
    let x = (q_gev / lambda_qcd_gev).ln();
    (2.0 * std::f64::consts::PI) / (beta0 * x)
}

/// Invert one-loop formula:
/// Λ = Q * exp(-2π/(β0 * alpha_s(Q))).
fn lambda_qcd_from_anchor(beta0: f64, q_gev: f64, alpha_s_q: f64) -> f64 {
    q_gev * (-(2.0 * std::f64::consts::PI) / (beta0 * alpha_s_q)).exp()
}

fn beta0_su3_projected(nf: u32) -> f64 {
    // β0 = (11/3) C_A - (2/3) n_f with C_A = 3 for SU(3).
    11.0 - (2.0 / 3.0) * nf as f64
}

fn beta1_su3_projected(nf: u32) -> f64 {
    // β1 = 102 - (38/3) n_f for SU(3).
    102.0 - (38.0 / 3.0) * nf as f64
}

fn beta2_su3_projected(nf: u32) -> f64 {
    // MS-bar: β2 = 2857/2 - 5033/18 n_f + 325/54 n_f² for SU(3).
    (2857.0 / 2.0) - (5033.0 / 18.0) * nf as f64 + (325.0 / 54.0) * (nf as f64).powi(2)
}

fn alpha_s_two_loop(beta0: f64, beta1: f64, q_gev: f64, lambda_qcd_gev: f64) -> f64 {
    let l = ((q_gev * q_gev) / (lambda_qcd_gev * lambda_qcd_gev)).ln();
    let c = beta1 / (beta0 * beta0);
    let correction = 1.0 - c * l.ln() / l;
    (4.0 * std::f64::consts::PI / (beta0 * l)) * correction
}

fn alpha_s_three_loop(
    beta0: f64,
    beta1: f64,
    beta2: f64,
    q_gev: f64,
    lambda_qcd_gev: f64,
) -> f64 {
    let l = ((q_gev * q_gev) / (lambda_qcd_gev * lambda_qcd_gev)).ln();
    let ln_l = l.ln();
    let c1 = beta1 / (beta0 * beta0);
    let c2 = beta2 / (beta0 * beta0 * beta0);
    let bracket = 1.0 - c1 * ln_l / l + (c1 * c1 * (ln_l * ln_l - ln_l - 1.0) + c2) / (l * l);
    (4.0 * std::f64::consts::PI / (beta0 * l)) * bracket
}

fn lambda_qcd_two_loop_from_anchor(beta0: f64, beta1: f64, q_gev: f64, alpha_s_q: f64) -> f64 {
    // Monotonic bisection in lambda on (0, q): alpha_s rises with lambda.
    let mut lo = 1.0e-6_f64;
    let mut hi = q_gev * 0.999;
    for _ in 0..140 {
        let mid = 0.5 * (lo + hi);
        let a_mid = alpha_s_two_loop(beta0, beta1, q_gev, mid);
        if a_mid.is_nan() || !a_mid.is_finite() || a_mid > alpha_s_q {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    0.5 * (lo + hi)
}

fn lambda_qcd_three_loop_from_anchor(
    beta0: f64,
    beta1: f64,
    beta2: f64,
    q_gev: f64,
    alpha_s_q: f64,
) -> f64 {
    let mut lo = 1.0e-6_f64;
    let mut hi = q_gev * 0.999;
    for _ in 0..160 {
        let mid = 0.5 * (lo + hi);
        let a_mid = alpha_s_three_loop(beta0, beta1, beta2, q_gev, mid);
        if a_mid.is_nan() || !a_mid.is_finite() || a_mid > alpha_s_q {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    0.5 * (lo + hi)
}

fn match_lambda_at_threshold(
    lambda_high: f64,
    nf_high: u32,
    nf_low: u32,
    threshold_gev: f64,
) -> f64 {
    // One-loop continuity: α_s^high(threshold) = α_s^low(threshold).
    let alpha_thr = alpha_s_one_loop(beta0_su3_projected(nf_high), threshold_gev, lambda_high);
    lambda_qcd_from_anchor(beta0_su3_projected(nf_low), threshold_gev, alpha_thr)
}

fn match_lambda_two_loop_at_threshold(
    lambda_high: f64,
    nf_high: u32,
    nf_low: u32,
    threshold_gev: f64,
) -> f64 {
    let beta0_hi = beta0_su3_projected(nf_high);
    let beta1_hi = beta1_su3_projected(nf_high);
    let beta0_lo = beta0_su3_projected(nf_low);
    let beta1_lo = beta1_su3_projected(nf_low);
    let alpha_thr = alpha_s_two_loop(beta0_hi, beta1_hi, threshold_gev, lambda_high);
    lambda_qcd_two_loop_from_anchor(beta0_lo, beta1_lo, threshold_gev, alpha_thr)
}

fn match_lambda_three_loop_at_threshold(
    lambda_high: f64,
    nf_high: u32,
    nf_low: u32,
    threshold_gev: f64,
) -> f64 {
    let beta0_hi = beta0_su3_projected(nf_high);
    let beta1_hi = beta1_su3_projected(nf_high);
    let beta2_hi = beta2_su3_projected(nf_high);
    let beta0_lo = beta0_su3_projected(nf_low);
    let beta1_lo = beta1_su3_projected(nf_low);
    let beta2_lo = beta2_su3_projected(nf_low);
    let alpha_thr = alpha_s_three_loop(beta0_hi, beta1_hi, beta2_hi, threshold_gev, lambda_high);
    lambda_qcd_three_loop_from_anchor(beta0_lo, beta1_lo, beta2_lo, threshold_gev, alpha_thr)
}

fn main() {
    let out_dir = std::env::var("GUTOE_QCD_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/qcd_scale_report".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let sm = StandardModelDynamicsMap::from_clifford_z3();
    let beta0_unified = sm.beta0; // structural 58/3 (Clifford-wide)
    let nf_mz = env_u32("GUTOE_QCD_NF_MZ", 5);
    let beta0_su3 = beta0_su3_projected(nf_mz);
    let beta1_su3 = beta1_su3_projected(nf_mz);
    let beta2_su3 = beta2_su3_projected(nf_mz);
    let q_ref_gev = env_f64("GUTOE_QCD_QREF_GEV", MZ_GEV_DEFAULT);
    let alpha_s_ref = env_f64("GUTOE_QCD_ALPHA_S_REF", ALPHA_S_MZ_DEFAULT);
    let m_b_gev = env_f64("GUTOE_QCD_MB_GEV", 4.18);
    let m_c_gev = env_f64("GUTOE_QCD_MC_GEV", 1.27);

    let lambda_unified_gev = lambda_qcd_from_anchor(beta0_unified, q_ref_gev, alpha_s_ref);
    let alpha_unified_backcheck = alpha_s_one_loop(beta0_unified, q_ref_gev, lambda_unified_gev);
    let alpha_unified_1gev = alpha_s_one_loop(beta0_unified, 1.0, lambda_unified_gev);
    let alpha_unified_2gev = alpha_s_one_loop(beta0_unified, 2.0, lambda_unified_gev);
    let alpha_unified_10gev = alpha_s_one_loop(beta0_unified, 10.0, lambda_unified_gev);

    let lambda_su3_gev = lambda_qcd_from_anchor(beta0_su3, q_ref_gev, alpha_s_ref);
    let alpha_su3_backcheck = alpha_s_one_loop(beta0_su3, q_ref_gev, lambda_su3_gev);
    let alpha_su3_1gev = alpha_s_one_loop(beta0_su3, 1.0, lambda_su3_gev);
    let alpha_su3_2gev = alpha_s_one_loop(beta0_su3, 2.0, lambda_su3_gev);
    let alpha_su3_10gev = alpha_s_one_loop(beta0_su3, 10.0, lambda_su3_gev);
    let lambda_nf4_gev = match_lambda_at_threshold(lambda_su3_gev, 5, 4, m_b_gev);
    let lambda_nf3_gev = match_lambda_at_threshold(lambda_nf4_gev, 4, 3, m_c_gev);
    let alpha_matched_10gev = alpha_s_one_loop(beta0_su3_projected(5), 10.0, lambda_su3_gev);
    let alpha_matched_2gev = alpha_s_one_loop(beta0_su3_projected(4), 2.0, lambda_nf4_gev);
    let alpha_matched_1gev = alpha_s_one_loop(beta0_su3_projected(3), 1.0, lambda_nf3_gev);
    let lambda_su3_2l_gev =
        lambda_qcd_two_loop_from_anchor(beta0_su3, beta1_su3, q_ref_gev, alpha_s_ref);
    let lambda_nf4_2l_gev = match_lambda_two_loop_at_threshold(lambda_su3_2l_gev, 5, 4, m_b_gev);
    let lambda_nf3_2l_gev = match_lambda_two_loop_at_threshold(lambda_nf4_2l_gev, 4, 3, m_c_gev);
    let alpha_2l_matched_10gev = alpha_s_two_loop(
        beta0_su3_projected(5),
        beta1_su3_projected(5),
        10.0,
        lambda_su3_2l_gev,
    );
    let alpha_2l_matched_2gev = alpha_s_two_loop(
        beta0_su3_projected(4),
        beta1_su3_projected(4),
        2.0,
        lambda_nf4_2l_gev,
    );
    let alpha_2l_matched_1gev = alpha_s_two_loop(
        beta0_su3_projected(3),
        beta1_su3_projected(3),
        1.0,
        lambda_nf3_2l_gev,
    );
    let lambda_su3_3l_gev =
        lambda_qcd_three_loop_from_anchor(beta0_su3, beta1_su3, beta2_su3, q_ref_gev, alpha_s_ref);
    let lambda_nf4_3l_gev = match_lambda_three_loop_at_threshold(lambda_su3_3l_gev, 5, 4, m_b_gev);
    let lambda_nf3_3l_gev = match_lambda_three_loop_at_threshold(lambda_nf4_3l_gev, 4, 3, m_c_gev);
    let alpha_3l_matched_10gev = alpha_s_three_loop(
        beta0_su3_projected(5),
        beta1_su3_projected(5),
        beta2_su3_projected(5),
        10.0,
        lambda_su3_3l_gev,
    );
    let alpha_3l_matched_2gev = alpha_s_three_loop(
        beta0_su3_projected(4),
        beta1_su3_projected(4),
        beta2_su3_projected(4),
        2.0,
        lambda_nf4_3l_gev,
    );
    let alpha_3l_matched_1gev = alpha_s_three_loop(
        beta0_su3_projected(3),
        beta1_su3_projected(3),
        beta2_su3_projected(3),
        1.0,
        lambda_nf3_3l_gev,
    );
    let alpha_s_structural_leading = CLIFFORD_DIM / ALPHA_INV_STRUCTURAL; // 16/137
    let alpha_s_structural_corrected = alpha_s_structural_leading * C_INF_STRUCTURAL; // (16/137)*(67/66)
    let lambda_su3_3l_struct_lead_gev = lambda_qcd_three_loop_from_anchor(
        beta0_su3,
        beta1_su3,
        beta2_su3,
        q_ref_gev,
        alpha_s_structural_leading,
    );
    let lambda_nf4_3l_struct_lead_gev =
        match_lambda_three_loop_at_threshold(lambda_su3_3l_struct_lead_gev, 5, 4, m_b_gev);
    let lambda_nf3_3l_struct_lead_gev =
        match_lambda_three_loop_at_threshold(lambda_nf4_3l_struct_lead_gev, 4, 3, m_c_gev);
    let alpha_3l_struct_lead_2gev = alpha_s_three_loop(
        beta0_su3_projected(4),
        beta1_su3_projected(4),
        beta2_su3_projected(4),
        2.0,
        lambda_nf4_3l_struct_lead_gev,
    );
    let lambda_su3_3l_struct_corr_gev = lambda_qcd_three_loop_from_anchor(
        beta0_su3,
        beta1_su3,
        beta2_su3,
        q_ref_gev,
        alpha_s_structural_corrected,
    );
    let lambda_nf4_3l_struct_corr_gev =
        match_lambda_three_loop_at_threshold(lambda_su3_3l_struct_corr_gev, 5, 4, m_b_gev);
    let lambda_nf3_3l_struct_corr_gev =
        match_lambda_three_loop_at_threshold(lambda_nf4_3l_struct_corr_gev, 4, 3, m_c_gev);
    let alpha_3l_struct_corr_2gev = alpha_s_three_loop(
        beta0_su3_projected(4),
        beta1_su3_projected(4),
        beta2_su3_projected(4),
        2.0,
        lambda_nf4_3l_struct_corr_gev,
    );
    let mp_over_lambda_nf3_1l = PROTON_MASS_MEV_OBS / (lambda_nf3_gev * 1000.0);
    let mp_over_lambda_nf3_2l = PROTON_MASS_MEV_OBS / (lambda_nf3_2l_gev * 1000.0);
    let mp_over_lambda_nf3_3l = PROTON_MASS_MEV_OBS / (lambda_nf3_3l_gev * 1000.0);
    let mp_over_lambda_nf3_3l_struct_lead =
        PROTON_MASS_MEV_OBS / (lambda_nf3_3l_struct_lead_gev * 1000.0);
    let mp_over_lambda_nf3_3l_struct_corr =
        PROTON_MASS_MEV_OBS / (lambda_nf3_3l_struct_corr_gev * 1000.0);

    // "Confinement scale in torsion sector" status is currently a gap.
    // We expose this explicitly so downstream reports don't imply a bridge exists.
    let torsion_bridge_status = "missing_explicit_bridge";

    let txt_path = out.join("qcd_scale_report.txt");
    let json_path = out.join("qcd_scale_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[structural_inputs]").expect("write");
    writeln!(txt, "beta0_structural_unified = {:.12}", beta0_unified).expect("write");
    writeln!(txt, "beta0_expected_58_over_3 = {:.12}", 58.0 / 3.0).expect("write");
    writeln!(txt, "beta0_su3_projected = {:.12}", beta0_su3).expect("write");
    writeln!(txt, "beta1_su3_projected = {:.12}", beta1_su3).expect("write");
    writeln!(txt, "beta2_su3_projected = {:.12}", beta2_su3).expect("write");
    writeln!(txt, "nf_mz_projected = {}", nf_mz).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[assumptions]").expect("write");
    writeln!(txt, "q_ref_gev = {:.12}", q_ref_gev).expect("write");
    writeln!(txt, "alpha_s_ref = {:.12}", alpha_s_ref).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[derived_qcd_scale_unified_beta0]").expect("write");
    writeln!(txt, "lambda_qcd_gev = {:.12}", lambda_unified_gev).expect("write");
    writeln!(txt, "lambda_qcd_mev = {:.6}", lambda_unified_gev * 1000.0).expect("write");
    writeln!(txt, "alpha_s_backcheck_ref = {:.12}", alpha_unified_backcheck).expect("write");
    writeln!(txt, "alpha_s_1gev = {:.12}", alpha_unified_1gev).expect("write");
    writeln!(txt, "alpha_s_2gev = {:.12}", alpha_unified_2gev).expect("write");
    writeln!(txt, "alpha_s_10gev = {:.12}", alpha_unified_10gev).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[derived_qcd_scale_su3_projected]").expect("write");
    writeln!(txt, "lambda_qcd_gev = {:.12}", lambda_su3_gev).expect("write");
    writeln!(txt, "lambda_qcd_mev = {:.6}", lambda_su3_gev * 1000.0).expect("write");
    writeln!(txt, "alpha_s_backcheck_ref = {:.12}", alpha_su3_backcheck).expect("write");
    writeln!(txt, "alpha_s_1gev = {:.12}", alpha_su3_1gev).expect("write");
    writeln!(txt, "alpha_s_2gev = {:.12}", alpha_su3_2gev).expect("write");
    writeln!(txt, "alpha_s_10gev = {:.12}", alpha_su3_10gev).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[derived_qcd_scale_su3_threshold_matched]").expect("write");
    writeln!(txt, "m_b_gev = {:.12}", m_b_gev).expect("write");
    writeln!(txt, "m_c_gev = {:.12}", m_c_gev).expect("write");
    writeln!(txt, "lambda_nf5_gev = {:.12}", lambda_su3_gev).expect("write");
    writeln!(txt, "lambda_nf4_gev = {:.12}", lambda_nf4_gev).expect("write");
    writeln!(txt, "lambda_nf3_gev = {:.12}", lambda_nf3_gev).expect("write");
    writeln!(txt, "alpha_s_10gev_matched_nf5 = {:.12}", alpha_matched_10gev).expect("write");
    writeln!(txt, "alpha_s_2gev_matched_nf4 = {:.12}", alpha_matched_2gev).expect("write");
    writeln!(txt, "alpha_s_1gev_matched_nf3 = {:.12}", alpha_matched_1gev).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[derived_qcd_scale_su3_threshold_matched_two_loop]").expect("write");
    writeln!(txt, "lambda_nf5_gev = {:.12}", lambda_su3_2l_gev).expect("write");
    writeln!(txt, "lambda_nf4_gev = {:.12}", lambda_nf4_2l_gev).expect("write");
    writeln!(txt, "lambda_nf3_gev = {:.12}", lambda_nf3_2l_gev).expect("write");
    writeln!(txt, "alpha_s_10gev_matched_nf5 = {:.12}", alpha_2l_matched_10gev).expect("write");
    writeln!(txt, "alpha_s_2gev_matched_nf4 = {:.12}", alpha_2l_matched_2gev).expect("write");
    writeln!(txt, "alpha_s_1gev_matched_nf3 = {:.12}", alpha_2l_matched_1gev).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[derived_qcd_scale_su3_threshold_matched_three_loop]").expect("write");
    writeln!(txt, "lambda_nf5_gev = {:.12}", lambda_su3_3l_gev).expect("write");
    writeln!(txt, "lambda_nf4_gev = {:.12}", lambda_nf4_3l_gev).expect("write");
    writeln!(txt, "lambda_nf3_gev = {:.12}", lambda_nf3_3l_gev).expect("write");
    writeln!(txt, "alpha_s_10gev_matched_nf5 = {:.12}", alpha_3l_matched_10gev).expect("write");
    writeln!(txt, "alpha_s_2gev_matched_nf4 = {:.12}", alpha_3l_matched_2gev).expect("write");
    writeln!(txt, "alpha_s_1gev_matched_nf3 = {:.12}", alpha_3l_matched_1gev).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[alpha_s_mz_structural_candidates]").expect("write");
    writeln!(txt, "alpha_s_mz_structural_leading = {:.12}", alpha_s_structural_leading).expect("write");
    writeln!(txt, "alpha_s_mz_structural_corrected = {:.12}", alpha_s_structural_corrected).expect("write");
    writeln!(txt, "alpha_s_mz_observed_anchor = {:.12}", alpha_s_ref).expect("write");
    writeln!(
        txt,
        "leading_rel_err_percent = {:.6}",
        100.0 * (alpha_s_structural_leading - alpha_s_ref) / alpha_s_ref
    )
    .expect("write");
    writeln!(
        txt,
        "corrected_rel_err_percent = {:.6}",
        100.0 * (alpha_s_structural_corrected - alpha_s_ref) / alpha_s_ref
    )
    .expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[derived_qcd_scale_three_loop_from_structural_alpha_s]").expect("write");
    writeln!(txt, "leading_lambda_nf3_gev = {:.12}", lambda_nf3_3l_struct_lead_gev).expect("write");
    writeln!(txt, "leading_alpha_s_2gev = {:.12}", alpha_3l_struct_lead_2gev).expect("write");
    writeln!(txt, "corrected_lambda_nf3_gev = {:.12}", lambda_nf3_3l_struct_corr_gev).expect("write");
    writeln!(txt, "corrected_alpha_s_2gev = {:.12}", alpha_3l_struct_corr_2gev).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[proton_ratio]").expect("write");
    writeln!(txt, "proton_mass_mev_obs = {:.12}", PROTON_MASS_MEV_OBS).expect("write");
    writeln!(txt, "mp_over_lambda_nf3_one_loop = {:.12}", mp_over_lambda_nf3_1l).expect("write");
    writeln!(txt, "mp_over_lambda_nf3_two_loop = {:.12}", mp_over_lambda_nf3_2l).expect("write");
    writeln!(txt, "mp_over_lambda_nf3_three_loop = {:.12}", mp_over_lambda_nf3_3l).expect("write");
    writeln!(
        txt,
        "mp_over_lambda_nf3_three_loop_structural_leading = {:.12}",
        mp_over_lambda_nf3_3l_struct_lead
    )
    .expect("write");
    writeln!(
        txt,
        "mp_over_lambda_nf3_three_loop_structural_corrected = {:.12}",
        mp_over_lambda_nf3_3l_struct_corr
    )
    .expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[status]").expect("write");
    writeln!(txt, "torsion_confinement_bridge = {}", torsion_bridge_status).expect("write");
    writeln!(
        txt,
        "note = one-loop Λ_QCD derivation currently requires alpha_s(Q_ref) anchor; anchor-free Λ_QCD is open"
    )
    .expect("write");

    let mut json = File::create(&json_path).expect("create json");
    writeln!(
        json,
        "{{\n  \"structural_inputs\": {{\"beta0_unified\": {:.12}, \"beta0_expected_58_over_3\": {:.12}, \"beta0_su3_projected\": {:.12}, \"beta1_su3_projected\": {:.12}, \"beta2_su3_projected\": {:.12}, \"nf_mz_projected\": {}, \"clifford_dim\": {:.0}, \"alpha_inv_structural\": {:.0}, \"c_inf_structural\": {:.12}}},\n  \"assumptions\": {{\"q_ref_gev\": {:.12}, \"alpha_s_ref\": {:.12}, \"m_b_gev\": {:.12}, \"m_c_gev\": {:.12}}},\n  \"alpha_s_mz_structural_candidates\": {{\"leading\": {:.12}, \"corrected\": {:.12}, \"observed_anchor\": {:.12}}},\n  \"derived_qcd_scale_unified_beta0\": {{\"lambda_qcd_gev\": {:.12}, \"lambda_qcd_mev\": {:.6}, \"alpha_s_backcheck_ref\": {:.12}, \"alpha_s_1gev\": {:.12}, \"alpha_s_2gev\": {:.12}, \"alpha_s_10gev\": {:.12}}},\n  \"derived_qcd_scale_su3_projected\": {{\"lambda_qcd_gev\": {:.12}, \"lambda_qcd_mev\": {:.6}, \"alpha_s_backcheck_ref\": {:.12}, \"alpha_s_1gev\": {:.12}, \"alpha_s_2gev\": {:.12}, \"alpha_s_10gev\": {:.12}}},\n  \"derived_qcd_scale_su3_threshold_matched\": {{\"lambda_nf5_gev\": {:.12}, \"lambda_nf4_gev\": {:.12}, \"lambda_nf3_gev\": {:.12}, \"alpha_s_10gev_matched_nf5\": {:.12}, \"alpha_s_2gev_matched_nf4\": {:.12}, \"alpha_s_1gev_matched_nf3\": {:.12}}},\n  \"derived_qcd_scale_su3_threshold_matched_two_loop\": {{\"lambda_nf5_gev\": {:.12}, \"lambda_nf4_gev\": {:.12}, \"lambda_nf3_gev\": {:.12}, \"alpha_s_10gev_matched_nf5\": {:.12}, \"alpha_s_2gev_matched_nf4\": {:.12}, \"alpha_s_1gev_matched_nf3\": {:.12}}},\n  \"derived_qcd_scale_su3_threshold_matched_three_loop\": {{\"lambda_nf5_gev\": {:.12}, \"lambda_nf4_gev\": {:.12}, \"lambda_nf3_gev\": {:.12}, \"alpha_s_10gev_matched_nf5\": {:.12}, \"alpha_s_2gev_matched_nf4\": {:.12}, \"alpha_s_1gev_matched_nf3\": {:.12}}},\n  \"derived_qcd_scale_three_loop_from_structural_alpha_s\": {{\"leading_lambda_nf3_gev\": {:.12}, \"leading_alpha_s_2gev\": {:.12}, \"corrected_lambda_nf3_gev\": {:.12}, \"corrected_alpha_s_2gev\": {:.12}}},\n  \"proton_ratio\": {{\"proton_mass_mev_obs\": {:.12}, \"mp_over_lambda_nf3_one_loop\": {:.12}, \"mp_over_lambda_nf3_two_loop\": {:.12}, \"mp_over_lambda_nf3_three_loop\": {:.12}, \"mp_over_lambda_nf3_three_loop_structural_leading\": {:.12}, \"mp_over_lambda_nf3_three_loop_structural_corrected\": {:.12}}},\n  \"status\": {{\"torsion_confinement_bridge\": \"{}\", \"note\": \"one-loop/two-loop/three-loop Lambda_QCD derivations currently require alpha_s(Q_ref) anchor; structural alpha_s(M_Z) candidates are now reported\"}}\n}}",
        beta0_unified,
        58.0 / 3.0,
        beta0_su3,
        beta1_su3,
        beta2_su3,
        nf_mz,
        CLIFFORD_DIM,
        ALPHA_INV_STRUCTURAL,
        C_INF_STRUCTURAL,
        q_ref_gev,
        alpha_s_ref,
        m_b_gev,
        m_c_gev,
        alpha_s_structural_leading,
        alpha_s_structural_corrected,
        alpha_s_ref,
        lambda_unified_gev,
        lambda_unified_gev * 1000.0,
        alpha_unified_backcheck,
        alpha_unified_1gev,
        alpha_unified_2gev,
        alpha_unified_10gev,
        lambda_su3_gev,
        lambda_su3_gev * 1000.0,
        alpha_su3_backcheck,
        alpha_su3_1gev,
        alpha_su3_2gev,
        alpha_su3_10gev,
        lambda_su3_gev,
        lambda_nf4_gev,
        lambda_nf3_gev,
        alpha_matched_10gev,
        alpha_matched_2gev,
        alpha_matched_1gev,
        lambda_su3_2l_gev,
        lambda_nf4_2l_gev,
        lambda_nf3_2l_gev,
        alpha_2l_matched_10gev,
        alpha_2l_matched_2gev,
        alpha_2l_matched_1gev,
        lambda_su3_3l_gev,
        lambda_nf4_3l_gev,
        lambda_nf3_3l_gev,
        alpha_3l_matched_10gev,
        alpha_3l_matched_2gev,
        alpha_3l_matched_1gev,
        lambda_nf3_3l_struct_lead_gev,
        alpha_3l_struct_lead_2gev,
        lambda_nf3_3l_struct_corr_gev,
        alpha_3l_struct_corr_2gev,
        PROTON_MASS_MEV_OBS,
        mp_over_lambda_nf3_1l,
        mp_over_lambda_nf3_2l,
        mp_over_lambda_nf3_3l,
        mp_over_lambda_nf3_3l_struct_lead,
        mp_over_lambda_nf3_3l_struct_corr,
        torsion_bridge_status
    )
    .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "unified beta0={:.6} -> lambda_qcd={:.6} GeV ({:.1} MeV)",
        beta0_unified,
        lambda_unified_gev,
        lambda_unified_gev * 1000.0
    );
    println!(
        "SU(3)-projected beta0={:.6} (nf={}) -> lambda_qcd={:.6} GeV ({:.1} MeV) from alpha_s({:.3} GeV)={:.4}",
        beta0_su3,
        nf_mz,
        lambda_su3_gev,
        lambda_su3_gev * 1000.0,
        q_ref_gev,
        alpha_s_ref
    );
    println!(
        "threshold-matched one-loop: lambda_nf5={:.6} GeV, lambda_nf4={:.6} GeV, lambda_nf3={:.6} GeV | alpha_s(10 GeV)={:.4}, alpha_s(2 GeV)={:.4}, alpha_s(1 GeV)={:.4}",
        lambda_su3_gev,
        lambda_nf4_gev,
        lambda_nf3_gev,
        alpha_matched_10gev,
        alpha_matched_2gev,
        alpha_matched_1gev
    );
    println!(
        "threshold-matched two-loop: lambda_nf5={:.6} GeV, lambda_nf4={:.6} GeV, lambda_nf3={:.6} GeV | alpha_s(10 GeV)={:.4}, alpha_s(2 GeV)={:.4}, alpha_s(1 GeV)={:.4}",
        lambda_su3_2l_gev,
        lambda_nf4_2l_gev,
        lambda_nf3_2l_gev,
        alpha_2l_matched_10gev,
        alpha_2l_matched_2gev,
        alpha_2l_matched_1gev
    );
    println!(
        "threshold-matched three-loop: lambda_nf5={:.6} GeV, lambda_nf4={:.6} GeV, lambda_nf3={:.6} GeV | alpha_s(10 GeV)={:.4}, alpha_s(2 GeV)={:.4}, alpha_s(1 GeV)={:.4}",
        lambda_su3_3l_gev,
        lambda_nf4_3l_gev,
        lambda_nf3_3l_gev,
        alpha_3l_matched_10gev,
        alpha_3l_matched_2gev,
        alpha_3l_matched_1gev
    );
    println!(
        "structural alpha_s(M_Z) candidates: leading(16/137)={:.6}, corrected((16/137)*(67/66))={:.6}, anchor={:.6}",
        alpha_s_structural_leading,
        alpha_s_structural_corrected,
        alpha_s_ref
    );
    println!(
        "from structural alpha_s (3-loop): leading lambda_nf3={:.6} GeV alpha_s(2GeV)={:.4}; corrected lambda_nf3={:.6} GeV alpha_s(2GeV)={:.4}",
        lambda_nf3_3l_struct_lead_gev,
        alpha_3l_struct_lead_2gev,
        lambda_nf3_3l_struct_corr_gev,
        alpha_3l_struct_corr_2gev
    );
    println!(
        "mp/lambda_nf3: one-loop={:.4}, two-loop={:.4}, three-loop={:.4}, structural-leading={:.4}, structural-corrected={:.4}",
        mp_over_lambda_nf3_1l,
        mp_over_lambda_nf3_2l,
        mp_over_lambda_nf3_3l,
        mp_over_lambda_nf3_3l_struct_lead,
        mp_over_lambda_nf3_3l_struct_corr
    );
}
