use gutoe_em::alpha::{
    electron_mass_from_clifford_improved_with_alpha,
    lepton_masses_from_electron_structural_alpha, lepton_masses_from_electron_with_alpha,
    triangular, ALPHA_INVERSE_PHYSICAL, ALPHA_INVERSE_STRUCTURAL,
};
use gutoe_physics::constants::{ALPHA_LEADING_ORDER, C, G, HBAR};
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

const M_E_OBS: f64 = 0.510_998_950;
const M_MU_OBS: f64 = 105.658_375_5;
const M_TAU_OBS: f64 = 1776.93;
const KG_TO_MEV: f64 = 5.609_588_603e29;
const TWO_TERM_ABS_ERR_MAX: f64 = 1.0e-5;

fn electron_transduction_factor_f() -> f64 {
    let ratio_corrected: f64 = 115.0 / 22.0;
    let c_inf: f64 = 67.0 / 66.0;
    ALPHA_LEADING_ORDER.powi(13) * ratio_corrected.powi(3) * c_inf * 12.0_f64.powi(3)
}

fn g_from_me(me_mev: f64, f: f64) -> (f64, f64) {
    let m_pl_mev = me_mev / f;
    let m_pl_kg = m_pl_mev / KG_TO_MEV;
    let g_pred = HBAR * C / (m_pl_kg * m_pl_kg);
    let rel = (g_pred - G) / G;
    (g_pred, rel)
}

fn main() {
    let out_dir = std::env::var("GUTOE_ALPHA_WEB_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/alpha_web_ci_report".to_string());
    let out = PathBuf::from(out_dir);
    fs::create_dir_all(&out).expect("create output dir");

    let t16 = triangular(1 << 4);
    let alpha_inv_struct = ALPHA_INVERSE_STRUCTURAL;
    let alpha_inv_phys = ALPHA_INVERSE_PHYSICAL;
    let alpha = 1.0 / alpha_inv_phys;
    let delta_target = alpha_inv_phys - alpha_inv_struct;
    let delta_first = 5.0 * alpha;
    let delta_second = 5.0 * alpha - 9.0 * alpha * alpha;
    let delta_first_abs_err = (delta_first - delta_target).abs();
    let delta_second_abs_err = (delta_second - delta_target).abs();

    let masses_struct = lepton_masses_from_electron_structural_alpha(M_E_OBS);
    let masses_phys = lepton_masses_from_electron_with_alpha(M_E_OBS, 1.0 / alpha_inv_phys);
    let [_me_s, mmu_s, mtau_s] = masses_struct;
    let [_me_p, mmu_p, mtau_p] = masses_phys;

    let mu_rel_s = (mmu_s - M_MU_OBS) / M_MU_OBS;
    let tau_rel_s = (mtau_s - M_TAU_OBS) / M_TAU_OBS;
    let mu_rel_p = (mmu_p - M_MU_OBS) / M_MU_OBS;
    let tau_rel_p = (mtau_p - M_TAU_OBS) / M_TAU_OBS;

    let f = electron_transduction_factor_f();
    let me_struct_from_mu_tau = electron_mass_from_clifford_improved_with_alpha(
        M_MU_OBS,
        M_TAU_OBS,
        1.0 / alpha_inv_struct,
    );
    let (g_me_obs, g_me_obs_rel) = g_from_me(M_E_OBS, f);
    let (g_me_structalpha, g_me_structalpha_rel) = g_from_me(me_struct_from_mu_tau, f);

    let alpha_identity_ok = t16 + 1 == 137 && (alpha_inv_struct - 137.0).abs() < 1e-15;
    let structural_lane_sane = mu_rel_s.abs() < 0.01 && tau_rel_s.abs() < 0.01;
    let bridge_lane_sane = g_me_obs_rel.abs() < 0.01;
    let alpha_second_order_improves = delta_second_abs_err < delta_first_abs_err;
    let alpha_second_order_within_band = delta_second_abs_err <= TWO_TERM_ABS_ERR_MAX;
    let passes_all =
        alpha_identity_ok
            && structural_lane_sane
            && bridge_lane_sane
            && alpha_second_order_improves
            && alpha_second_order_within_band;

    let txt_path = out.join("alpha_web_ci_report.txt");
    let json_path = out.join("alpha_web_ci_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[alpha_identity]").expect("write");
    writeln!(txt, "T16 = {}", t16).expect("write");
    writeln!(txt, "T16_plus_1 = {}", t16 + 1).expect("write");
    writeln!(txt, "alpha_inv_structural = {:.12}", alpha_inv_struct).expect("write");
    writeln!(txt, "alpha_inv_physical = {:.12}", alpha_inv_phys).expect("write");
    writeln!(txt, "alpha_inv_rel_offset = {:.12e}", (alpha_inv_phys - alpha_inv_struct) / alpha_inv_struct)
        .expect("write");
    writeln!(txt).expect("write");

    writeln!(txt, "[lepton_lane_from_me_obs]").expect("write");
    writeln!(txt, "structural_alpha_mu_rel = {:.12e}", mu_rel_s).expect("write");
    writeln!(txt, "structural_alpha_tau_rel = {:.12e}", tau_rel_s).expect("write");
    writeln!(txt, "physical_alpha_mu_rel = {:.12e}", mu_rel_p).expect("write");
    writeln!(txt, "physical_alpha_tau_rel = {:.12e}", tau_rel_p).expect("write");
    writeln!(txt).expect("write");

    writeln!(txt, "[alpha_correction_lane]").expect("write");
    writeln!(txt, "delta_target = {:.12e}", delta_target).expect("write");
    writeln!(txt, "delta_first_order_5alpha = {:.12e}", delta_first).expect("write");
    writeln!(
        txt,
        "delta_second_order_5alpha_minus_9alpha2 = {:.12e}",
        delta_second
    )
    .expect("write");
    writeln!(txt, "first_abs_error = {:.12e}", delta_first_abs_err).expect("write");
    writeln!(txt, "second_abs_error = {:.12e}", delta_second_abs_err).expect("write");
    writeln!(
        txt,
        "second_order_abs_error_band_max = {:.12e}",
        TWO_TERM_ABS_ERR_MAX
    )
    .expect("write");
    writeln!(
        txt,
        "second_order_improves = {}",
        alpha_second_order_improves
    )
    .expect("write");
    writeln!(
        txt,
        "second_order_within_band = {}",
        alpha_second_order_within_band
    )
    .expect("write");
    writeln!(txt).expect("write");

    writeln!(txt, "[g_bridge]").expect("write");
    writeln!(txt, "F = {:.15e}", f).expect("write");
    writeln!(txt, "g_pred_me_obs = {:.15e}", g_me_obs).expect("write");
    writeln!(txt, "g_rel_me_obs = {:.15e}", g_me_obs_rel).expect("write");
    writeln!(txt, "g_pred_me_structalpha = {:.15e}", g_me_structalpha).expect("write");
    writeln!(txt, "g_rel_me_structalpha = {:.15e}", g_me_structalpha_rel).expect("write");
    writeln!(txt).expect("write");

    writeln!(txt, "[ci_gate]").expect("write");
    writeln!(txt, "alpha_identity_ok = {}", alpha_identity_ok).expect("write");
    writeln!(txt, "structural_lane_sane = {}", structural_lane_sane).expect("write");
    writeln!(txt, "bridge_lane_sane = {}", bridge_lane_sane).expect("write");
    writeln!(txt, "alpha_second_order_improves = {}", alpha_second_order_improves).expect("write");
    writeln!(
        txt,
        "alpha_second_order_within_band = {}",
        alpha_second_order_within_band
    )
    .expect("write");
    writeln!(txt, "passes_all = {}", passes_all).expect("write");

    let mut json = File::create(&json_path).expect("create json");
    writeln!(
        json,
        "{{\n  \"alpha_identity\": {{\"t16\": {}, \"t16_plus_1\": {}, \"alpha_inv_structural\": {:.12}, \"alpha_inv_physical\": {:.12}, \"alpha_inv_rel_offset\": {:.12e}}},\n  \"lepton_lane_from_me_obs\": {{\"structural_alpha_mu_rel\": {:.12e}, \"structural_alpha_tau_rel\": {:.12e}, \"physical_alpha_mu_rel\": {:.12e}, \"physical_alpha_tau_rel\": {:.12e}}},\n  \"alpha_correction_lane\": {{\"delta_target\": {:.12e}, \"delta_first_order_5alpha\": {:.12e}, \"delta_second_order_5alpha_minus_9alpha2\": {:.12e}, \"first_abs_error\": {:.12e}, \"second_abs_error\": {:.12e}, \"second_abs_error_band_max\": {:.12e}, \"second_order_improves\": {}, \"second_order_within_band\": {}}},\n  \"g_bridge\": {{\"F\": {:.15e}, \"me_struct_from_mu_tau\": {:.12}, \"g_pred_me_obs\": {:.15e}, \"g_rel_me_obs\": {:.15e}, \"g_pred_me_structalpha\": {:.15e}, \"g_rel_me_structalpha\": {:.15e}}},\n  \"ci_gate\": {{\"alpha_identity_ok\": {}, \"structural_lane_sane\": {}, \"bridge_lane_sane\": {}, \"alpha_second_order_improves\": {}, \"alpha_second_order_within_band\": {}, \"passes_all\": {}}}\n}}",
        t16,
        t16 + 1,
        alpha_inv_struct,
        alpha_inv_phys,
        (alpha_inv_phys - alpha_inv_struct) / alpha_inv_struct,
        mu_rel_s,
        tau_rel_s,
        mu_rel_p,
        tau_rel_p,
        delta_target,
        delta_first,
        delta_second,
        delta_first_abs_err,
        delta_second_abs_err,
        TWO_TERM_ABS_ERR_MAX,
        alpha_second_order_improves,
        alpha_second_order_within_band,
        f,
        me_struct_from_mu_tau,
        g_me_obs,
        g_me_obs_rel,
        g_me_structalpha,
        g_me_structalpha_rel,
        alpha_identity_ok,
        structural_lane_sane,
        bridge_lane_sane,
        alpha_second_order_improves,
        alpha_second_order_within_band,
        passes_all
    )
    .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!("ci_gate passes_all={}", passes_all);

    if std::env::var("GUTOE_ALPHA_WEB_STRICT").ok().as_deref() == Some("1") && !passes_all {
        std::process::exit(2);
    }
}
