use gutoe_physics::constants::{ALPHA_LEADING_ORDER, C, G, HBAR};
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

const ELECTRON_MASS_MEV_OBS: f64 = 0.510_998_950;
const MUON_MASS_MEV_OBS: f64 = 105.658_375_5;
const TAU_MASS_MEV_OBS: f64 = 1776.93;
const PROTON_MASS_MEV_OBS: f64 = 938.272_088_16;
const MP_ME_STRUCT: f64 = 1836.0;
const KG_TO_MEV: f64 = 5.609_588_603e29;

/// Electron transduction candidate (dimensionless):
/// F = α^13 * (115/22)^3 * (67/66) * 12^3
fn electron_transduction_factor_f() -> f64 {
    let ratio_corrected: f64 = 115.0 / 22.0;
    let c_inf: f64 = 67.0 / 66.0;
    ALPHA_LEADING_ORDER.powi(13) * ratio_corrected.powi(3) * c_inf * 12.0_f64.powi(3)
}

#[derive(Clone, Copy)]
struct ModeResult {
    m_e_mev: f64,
    m_pl_kg_pred: f64,
    g_pred: f64,
    rel_err_g: f64,
}

fn evaluate_mode(m_e_mev: f64, f: f64) -> ModeResult {
    let m_pl_mev_pred = m_e_mev / f;
    let m_pl_kg_pred = m_pl_mev_pred / KG_TO_MEV;
    let g_pred = HBAR * C / (m_pl_kg_pred * m_pl_kg_pred);
    let rel_err_g = (g_pred - G) / G;
    ModeResult {
        m_e_mev,
        m_pl_kg_pred,
        g_pred,
        rel_err_g,
    }
}

/// Muon-lane electron estimate from the Koide phase with alpha^2 correction:
/// δ = 3π/4 - 5α(13/12) - c2 α², with c2 = 15/16 (current sweep candidate).
fn electron_from_mu_tau_phase_alpha2(m_mu: f64, m_tau: f64) -> f64 {
    use std::f64::consts::PI;
    let alpha = 1.0 / gutoe_em::alpha::ALPHA_INVERSE_PHYSICAL;
    let c2 = 15.0 / 16.0;
    let correction = 5.0 * alpha * (13.0 / 12.0) + c2 * alpha * alpha;
    let delta = 3.0 * PI / 4.0 - correction;
    let two_pi_3 = 2.0 * PI / 3.0;
    let c1 = (delta + two_pi_3).cos();
    let c2c = (delta + 2.0 * two_pi_3).cos();
    let a_mu = m_mu.sqrt();
    let a_tau = m_tau.sqrt();
    let m = (c2c * a_mu - c1 * a_tau) / (c2c - c1);
    let s = (a_mu / m - 1.0) / c1;
    let amp_e = m * (1.0 + s * delta.cos());
    amp_e * amp_e
}

fn main() {
    let out_dir =
        std::env::var("GUTOE_G_BRIDGE_OUT").unwrap_or_else(|_| "/tmp/bh_renders/g_bridge_report".to_string());
    let out = PathBuf::from(out_dir);
    fs::create_dir_all(&out).expect("create output dir");

    let f = electron_transduction_factor_f();

    // Mode A: measured electron mass input.
    let mode_measured = evaluate_mode(ELECTRON_MASS_MEV_OBS, f);

    // Mode B: structural electron lane via proton anchor + mp/me=1836.
    let m_e_structural = PROTON_MASS_MEV_OBS / MP_ME_STRUCT;
    let mode_structural = evaluate_mode(m_e_structural, f);

    // Mode C: muon/Koide phase lane with alpha^2 correction candidate.
    let m_e_muon_lane = electron_from_mu_tau_phase_alpha2(MUON_MASS_MEV_OBS, TAU_MASS_MEV_OBS);
    let mode_muon_lane = evaluate_mode(m_e_muon_lane, f);

    let txt_path = out.join("g_bridge_report.txt");
    let json_path = out.join("g_bridge_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[bridge]").expect("write");
    writeln!(txt, "F_expr = alpha^13 * (115/22)^3 * (67/66) * 12^3").expect("write");
    writeln!(txt, "F = {:.15e}", f).expect("write");
    writeln!(txt, "relation = G_pred = hbar * c * F^2 / m_e^2").expect("write");
    writeln!(txt, "G_codata = {:.15e}", G).expect("write");
    writeln!(txt).expect("write");

    writeln!(txt, "[mode_measured_electron]").expect("write");
    writeln!(txt, "m_e_mev = {:.12}", mode_measured.m_e_mev).expect("write");
    writeln!(txt, "m_pl_kg_pred = {:.15e}", mode_measured.m_pl_kg_pred).expect("write");
    writeln!(txt, "g_pred = {:.15e}", mode_measured.g_pred).expect("write");
    writeln!(txt, "g_rel_error = {:.15e}", mode_measured.rel_err_g).expect("write");
    writeln!(txt).expect("write");

    writeln!(txt, "[mode_structural_electron_from_proton]").expect("write");
    writeln!(txt, "m_e_mev = {:.12}", mode_structural.m_e_mev).expect("write");
    writeln!(txt, "m_pl_kg_pred = {:.15e}", mode_structural.m_pl_kg_pred).expect("write");
    writeln!(txt, "g_pred = {:.15e}", mode_structural.g_pred).expect("write");
    writeln!(txt, "g_rel_error = {:.15e}", mode_structural.rel_err_g).expect("write");
    writeln!(txt).expect("write");

    writeln!(txt, "[mode_muon_phase_alpha2_electron]").expect("write");
    writeln!(txt, "m_e_mev = {:.12}", mode_muon_lane.m_e_mev).expect("write");
    writeln!(txt, "m_pl_kg_pred = {:.15e}", mode_muon_lane.m_pl_kg_pred).expect("write");
    writeln!(txt, "g_pred = {:.15e}", mode_muon_lane.g_pred).expect("write");
    writeln!(txt, "g_rel_error = {:.15e}", mode_muon_lane.rel_err_g).expect("write");

    let mut json = File::create(&json_path).expect("create json");
    writeln!(
        json,
        "{{\n  \"bridge\": {{\"F_expr\": \"alpha^13 * (115/22)^3 * (67/66) * 12^3\", \"relation\": \"G_pred = hbar*c*F^2/m_e^2\", \"F\": {:.15e}, \"G_codata\": {:.15e}}},\n  \"mode_measured_electron\": {{\"m_e_mev\": {:.12}, \"m_pl_kg_pred\": {:.15e}, \"g_pred\": {:.15e}, \"g_rel_error\": {:.15e}}},\n  \"mode_structural_electron_from_proton\": {{\"m_e_mev\": {:.12}, \"m_pl_kg_pred\": {:.15e}, \"g_pred\": {:.15e}, \"g_rel_error\": {:.15e}}},\n  \"mode_muon_phase_alpha2_electron\": {{\"m_e_mev\": {:.12}, \"m_pl_kg_pred\": {:.15e}, \"g_pred\": {:.15e}, \"g_rel_error\": {:.15e}, \"phase_model\": \"delta=3pi/4-5alpha*(13/12)-(15/16)alpha^2\"}}\n}}",
        f,
        G,
        mode_measured.m_e_mev,
        mode_measured.m_pl_kg_pred,
        mode_measured.g_pred,
        mode_measured.rel_err_g,
        mode_structural.m_e_mev,
        mode_structural.m_pl_kg_pred,
        mode_structural.g_pred,
        mode_structural.rel_err_g,
        mode_muon_lane.m_e_mev,
        mode_muon_lane.m_pl_kg_pred,
        mode_muon_lane.g_pred,
        mode_muon_lane.rel_err_g
    )
    .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "measured m_e mode: G_pred={:.6e} (rel={:+.6e})",
        mode_measured.g_pred, mode_measured.rel_err_g
    );
    println!(
        "structural m_e mode: G_pred={:.6e} (rel={:+.6e})",
        mode_structural.g_pred, mode_structural.rel_err_g
    );
    println!(
        "muon-phase m_e mode: G_pred={:.6e} (rel={:+.6e})",
        mode_muon_lane.g_pred, mode_muon_lane.rel_err_g
    );
}
