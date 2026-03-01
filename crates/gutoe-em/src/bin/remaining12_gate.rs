//! Remaining-12 unified CI gate.
//!
//! Combines:
//! 1) Neutrino endgame lane checks (structural + triangulated)
//! 2) Absolute-scale endgame checks (lattice/proton-anchor branch)
//!
//! Exits nonzero if any required check fails.

use gutoe_em::alpha::koide_ratio;
use gutoe_em::{
    electron_mass_from_proton_anchor, electroweak_vev_from_fermi,
    electroweak_vev_from_lattice_order_parameter, higgs_mass_from_vev,
    neutrino_absolute_masses_from_texture, neutrino_hierarchy_prediction,
    triangulate_neutrino_from_splittings, w_mass_from_vev_and_alpha, z_mass_from_vev_and_alpha,
    ALPHA_EW_MZ,
};
use serde::Serialize;
use std::fs::{self, File};
use std::io::Write;
use std::process;

const DM21_TARGET_EV2: f64 = 7.53e-5;
const DM32_TARGET_EV2: f64 = 2.453e-3;
const NEUTRINO_RATIO_TOL: f64 = 0.05;
const NEUTRINO_TRI_TOL: f64 = 1.0e-9;
const G_F: f64 = 1.166_378_7e-5;
const SCALE_TOL: f64 = 0.01;

const M_E_REF_MEV: f64 = 0.510_998_950;
const M_W_REF_GEV: f64 = 80.377;
const M_Z_REF_GEV: f64 = 91.1876;
const M_H_REF_GEV: f64 = 125.25;

#[derive(Debug, Clone, Copy, Serialize)]
struct NeutrinoChecks {
    hierarchy_ok: bool,
    tiny_ok: bool,
    ratio_ok: bool,
    abs_splittings_ok: bool,
    no_fit_pass: bool,
    tri_ratio_ok: bool,
    tri_abs_ok: bool,
    triangulated_pass: bool,
}

#[derive(Debug, Clone, Copy, Serialize)]
struct NeutrinoMetrics {
    m1_ev: f64,
    m2_ev: f64,
    m3_ev: f64,
    sum_ev: f64,
    dm21_rel_err: f64,
    dm32_rel_err: f64,
    ratio_rel_err: f64,
    koide_k: f64,
    koide_s2: f64,
    p_triangulated: f64,
    tri_dm21_rel_err: f64,
    tri_dm32_rel_err: f64,
}

#[derive(Debug, Clone, Copy, Serialize)]
struct AbsoluteChecks {
    electron_ok: bool,
    vev_ok: bool,
    lattice_masses_ok: bool,
    fermi_masses_ok: bool,
    overall_pass: bool,
}

#[derive(Debug, Clone, Copy, Serialize)]
struct AbsoluteMetrics {
    electron_rel_err: f64,
    vev_rel_err: f64,
    m_w_lat_rel_err: f64,
    m_z_lat_rel_err: f64,
    m_h_lat_rel_err: f64,
    m_w_fer_rel_err: f64,
    m_z_fer_rel_err: f64,
    m_h_fer_rel_err: f64,
}

#[derive(Debug, Clone, Serialize)]
struct GateReport {
    overall_pass: bool,
    neutrino: NeutrinoBlock,
    absolute_scale: AbsoluteBlock,
}

#[derive(Debug, Clone, Serialize)]
struct NeutrinoBlock {
    checks: NeutrinoChecks,
    metrics: NeutrinoMetrics,
}

#[derive(Debug, Clone, Serialize)]
struct AbsoluteBlock {
    checks: AbsoluteChecks,
    metrics: AbsoluteMetrics,
}

fn rel_err(obs: f64, target: f64) -> f64 {
    if target.abs() < 1.0e-30 {
        0.0
    } else {
        (obs - target) / target
    }
}

fn pct_ok(rel: f64, tol: f64) -> bool {
    rel.abs() <= tol
}

fn main() {
    let out_dir =
        std::env::var("GUTOE_REMAINING12_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let txt_path = format!("{out_dir}/remaining12_gate.txt");
    let json_path = format!("{out_dir}/remaining12_gate.json");

    // ── Lane 1: neutrino endgame ─────────────────────────────────────────────
    let hierarchy = neutrino_hierarchy_prediction();
    let abs = neutrino_absolute_masses_from_texture();
    let tri = triangulate_neutrino_from_splittings(DM21_TARGET_EV2, DM32_TARGET_EV2);

    let dm21_rel = rel_err(abs.dm21_ev2, DM21_TARGET_EV2);
    let dm32_rel = rel_err(abs.dm32_ev2.abs(), DM32_TARGET_EV2);
    let ratio_target = DM32_TARGET_EV2 / DM21_TARGET_EV2;
    let ratio_rel = rel_err(abs.splitting_ratio_32_over_21, ratio_target);

    let k = koide_ratio([abs.m1_ev, abs.m2_ev, abs.m3_ev]);
    let s2 = 6.0 * k - 2.0;

    let hierarchy_ok = hierarchy == "normal";
    let tiny_ok = abs.m3_ev < 0.8 && abs.sum_ev < 0.12;
    let ratio_ok = pct_ok(ratio_rel, NEUTRINO_RATIO_TOL);
    let abs_splittings_ok = dm21_rel.abs() <= NEUTRINO_RATIO_TOL && dm32_rel.abs() <= NEUTRINO_RATIO_TOL;
    let no_fit_pass = hierarchy_ok && tiny_ok && ratio_ok && abs_splittings_ok;

    let tri_dm21_rel = rel_err(tri.dm21_ev2, DM21_TARGET_EV2);
    let tri_dm32_rel = rel_err(tri.dm32_ev2, DM32_TARGET_EV2);
    let tri_ratio_ok = tri.ratio_fit_rel_err.abs() <= NEUTRINO_TRI_TOL;
    let tri_abs_ok =
        tri_dm21_rel.abs() <= NEUTRINO_TRI_TOL && tri_dm32_rel.abs() <= NEUTRINO_TRI_TOL;
    let triangulated_pass = tri_ratio_ok && tri_abs_ok;

    let n_checks = NeutrinoChecks {
        hierarchy_ok,
        tiny_ok,
        ratio_ok,
        abs_splittings_ok,
        no_fit_pass,
        tri_ratio_ok,
        tri_abs_ok,
        triangulated_pass,
    };
    let n_metrics = NeutrinoMetrics {
        m1_ev: abs.m1_ev,
        m2_ev: abs.m2_ev,
        m3_ev: abs.m3_ev,
        sum_ev: abs.sum_ev,
        dm21_rel_err: dm21_rel,
        dm32_rel_err: dm32_rel,
        ratio_rel_err: ratio_rel,
        koide_k: k,
        koide_s2: s2,
        p_triangulated: tri.p_triangulated,
        tri_dm21_rel_err: tri_dm21_rel,
        tri_dm32_rel_err: tri_dm32_rel,
    };

    // ── Lane 2: absolute-scale endgame ───────────────────────────────────────
    let me_rel = rel_err(electron_mass_from_proton_anchor(), M_E_REF_MEV);
    let v_lattice = electroweak_vev_from_lattice_order_parameter(1.0);
    let v_fermi = electroweak_vev_from_fermi(G_F);
    let v_rel = rel_err(v_lattice, v_fermi);

    let m_w_lat_rel = rel_err(w_mass_from_vev_and_alpha(v_lattice, ALPHA_EW_MZ), M_W_REF_GEV);
    let m_z_lat_rel = rel_err(z_mass_from_vev_and_alpha(v_lattice, ALPHA_EW_MZ), M_Z_REF_GEV);
    let m_h_lat_rel = rel_err(higgs_mass_from_vev(v_lattice), M_H_REF_GEV);

    let m_w_fer_rel = rel_err(w_mass_from_vev_and_alpha(v_fermi, ALPHA_EW_MZ), M_W_REF_GEV);
    let m_z_fer_rel = rel_err(z_mass_from_vev_and_alpha(v_fermi, ALPHA_EW_MZ), M_Z_REF_GEV);
    let m_h_fer_rel = rel_err(higgs_mass_from_vev(v_fermi), M_H_REF_GEV);

    let electron_ok = pct_ok(me_rel, SCALE_TOL);
    let vev_ok = pct_ok(v_rel, SCALE_TOL);
    let lattice_masses_ok = pct_ok(m_w_lat_rel, SCALE_TOL)
        && pct_ok(m_z_lat_rel, SCALE_TOL)
        && pct_ok(m_h_lat_rel, SCALE_TOL);
    let fermi_masses_ok = pct_ok(m_w_fer_rel, SCALE_TOL)
        && pct_ok(m_z_fer_rel, SCALE_TOL)
        && pct_ok(m_h_fer_rel, SCALE_TOL);
    let a_overall = electron_ok && vev_ok && lattice_masses_ok && fermi_masses_ok;

    let a_checks = AbsoluteChecks {
        electron_ok,
        vev_ok,
        lattice_masses_ok,
        fermi_masses_ok,
        overall_pass: a_overall,
    };
    let a_metrics = AbsoluteMetrics {
        electron_rel_err: me_rel,
        vev_rel_err: v_rel,
        m_w_lat_rel_err: m_w_lat_rel,
        m_z_lat_rel_err: m_z_lat_rel,
        m_h_lat_rel_err: m_h_lat_rel,
        m_w_fer_rel_err: m_w_fer_rel,
        m_z_fer_rel_err: m_z_fer_rel,
        m_h_fer_rel_err: m_h_fer_rel,
    };

    let overall_pass = n_checks.no_fit_pass && n_checks.triangulated_pass && a_checks.overall_pass;

    let report = GateReport {
        overall_pass,
        neutrino: NeutrinoBlock {
            checks: n_checks,
            metrics: n_metrics,
        },
        absolute_scale: AbsoluteBlock {
            checks: a_checks,
            metrics: a_metrics,
        },
    };

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[remaining12_gate]").expect("write");
    writeln!(txt, "overall_pass={}", report.overall_pass).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[neutrino_checks]").expect("write");
    writeln!(
        txt,
        "hierarchy_ok={} tiny_ok={} ratio_ok={} abs_splittings_ok={} no_fit_pass={} tri_ratio_ok={} tri_abs_ok={} triangulated_pass={}",
        report.neutrino.checks.hierarchy_ok,
        report.neutrino.checks.tiny_ok,
        report.neutrino.checks.ratio_ok,
        report.neutrino.checks.abs_splittings_ok,
        report.neutrino.checks.no_fit_pass,
        report.neutrino.checks.tri_ratio_ok,
        report.neutrino.checks.tri_abs_ok,
        report.neutrino.checks.triangulated_pass,
    )
    .expect("write");
    writeln!(
        txt,
        "dm21_rel={:+.6e} dm32_rel={:+.6e} ratio_rel={:+.6e} K={:.12e} s2={:.12e} p_tri={:.12e}",
        report.neutrino.metrics.dm21_rel_err,
        report.neutrino.metrics.dm32_rel_err,
        report.neutrino.metrics.ratio_rel_err,
        report.neutrino.metrics.koide_k,
        report.neutrino.metrics.koide_s2,
        report.neutrino.metrics.p_triangulated,
    )
    .expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[absolute_checks]").expect("write");
    writeln!(
        txt,
        "electron_ok={} vev_ok={} lattice_masses_ok={} fermi_masses_ok={} overall={}",
        report.absolute_scale.checks.electron_ok,
        report.absolute_scale.checks.vev_ok,
        report.absolute_scale.checks.lattice_masses_ok,
        report.absolute_scale.checks.fermi_masses_ok,
        report.absolute_scale.checks.overall_pass,
    )
    .expect("write");
    writeln!(
        txt,
        "electron_rel={:+.6e} vev_rel={:+.6e} mW_lat_rel={:+.6e} mZ_lat_rel={:+.6e} mH_lat_rel={:+.6e} mW_fer_rel={:+.6e} mZ_fer_rel={:+.6e} mH_fer_rel={:+.6e}",
        report.absolute_scale.metrics.electron_rel_err,
        report.absolute_scale.metrics.vev_rel_err,
        report.absolute_scale.metrics.m_w_lat_rel_err,
        report.absolute_scale.metrics.m_z_lat_rel_err,
        report.absolute_scale.metrics.m_h_lat_rel_err,
        report.absolute_scale.metrics.m_w_fer_rel_err,
        report.absolute_scale.metrics.m_z_fer_rel_err,
        report.absolute_scale.metrics.m_h_fer_rel_err,
    )
    .expect("write");

    let mut json = File::create(&json_path).expect("create json");
    serde_json::to_writer_pretty(&mut json, &report).expect("write json");
    writeln!(json).expect("newline");

    println!(
        "remaining12_gate: overall_pass={} neutrino(no_fit={},tri={}) absolute={}",
        report.overall_pass,
        report.neutrino.checks.no_fit_pass,
        report.neutrino.checks.triangulated_pass,
        report.absolute_scale.checks.overall_pass
    );
    println!("wrote {txt_path}");
    println!("wrote {json_path}");

    if !report.overall_pass {
        process::exit(2);
    }
}
