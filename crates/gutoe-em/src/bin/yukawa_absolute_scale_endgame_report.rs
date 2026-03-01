//! Endgame absolute-scale closure report:
//! lattice/proton-anchor branch against Fermi/PDG anchors.

use gutoe_em::weak::{
    electron_mass_from_proton_anchor, electroweak_vev_from_fermi,
    electroweak_vev_from_lattice_order_parameter, higgs_mass_from_vev, w_mass_from_vev_and_alpha,
    w_z_mass_ratio, z_mass_from_vev_and_alpha, ALPHA_EW_MZ,
};
use serde::Serialize;
use std::fs::{self, File};
use std::io::Write;

const G_F: f64 = 1.166_378_7e-5;
const M_E_REF_MEV: f64 = 0.510_998_950;
const M_W_REF_GEV: f64 = 80.377;
const M_Z_REF_GEV: f64 = 91.1876;
const M_H_REF_GEV: f64 = 125.25;

#[derive(Debug, Clone, Copy, Serialize)]
struct ScalarComparison {
    predicted: f64,
    reference: f64,
    rel_err: f64,
}

#[derive(Debug, Clone, Serialize)]
struct Report {
    inputs: Inputs,
    electron_anchor: ScalarComparison,
    vev_lattice_vs_fermi: ScalarComparison,
    masses_lattice_branch: Masses,
    masses_fermi_branch: Masses,
    checks: Checks,
}

#[derive(Debug, Clone, Copy, Serialize)]
struct Inputs {
    g_f: f64,
    alpha_mz: f64,
    f0_vac: f64,
    w_over_z_structural: f64,
}

#[derive(Debug, Clone, Copy, Serialize)]
struct Masses {
    m_w_gev: ScalarComparison,
    m_z_gev: ScalarComparison,
    m_h_gev: ScalarComparison,
}

#[derive(Debug, Clone, Copy, Serialize)]
struct Checks {
    electron_anchor_ok_1pct: bool,
    vev_ok_1pct: bool,
    lattice_mass_ok_1pct: bool,
    fermi_mass_ok_1pct: bool,
    overall_pass: bool,
}

fn rel_err(pred: f64, reference: f64) -> f64 {
    if reference.abs() < 1.0e-30 {
        0.0
    } else {
        (pred - reference) / reference
    }
}

fn cmp(pred: f64, reference: f64) -> ScalarComparison {
    ScalarComparison {
        predicted: pred,
        reference,
        rel_err: rel_err(pred, reference),
    }
}

fn pct_ok(x: ScalarComparison) -> bool {
    x.rel_err.abs() <= 0.01
}

fn main() {
    let out_dir = std::env::var("GUTOE_ABSOLUTE_ENDGAME_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let txt_path = format!("{out_dir}/yukawa_absolute_scale_endgame_report.txt");
    let json_path = format!("{out_dir}/yukawa_absolute_scale_endgame_report.json");

    let f0_vac = 1.0;
    let me_anchor = cmp(electron_mass_from_proton_anchor(), M_E_REF_MEV);

    let v_lattice = electroweak_vev_from_lattice_order_parameter(f0_vac);
    let v_fermi = electroweak_vev_from_fermi(G_F);
    let v_cmp = cmp(v_lattice, v_fermi);

    let m_w_lat = cmp(w_mass_from_vev_and_alpha(v_lattice, ALPHA_EW_MZ), M_W_REF_GEV);
    let m_z_lat = cmp(z_mass_from_vev_and_alpha(v_lattice, ALPHA_EW_MZ), M_Z_REF_GEV);
    let m_h_lat = cmp(higgs_mass_from_vev(v_lattice), M_H_REF_GEV);

    let m_w_fer = cmp(w_mass_from_vev_and_alpha(v_fermi, ALPHA_EW_MZ), M_W_REF_GEV);
    let m_z_fer = cmp(z_mass_from_vev_and_alpha(v_fermi, ALPHA_EW_MZ), M_Z_REF_GEV);
    let m_h_fer = cmp(higgs_mass_from_vev(v_fermi), M_H_REF_GEV);

    let lattice_mass_ok = [m_w_lat, m_z_lat, m_h_lat].into_iter().all(pct_ok);
    let fermi_mass_ok = [m_w_fer, m_z_fer, m_h_fer].into_iter().all(pct_ok);

    let checks = Checks {
        electron_anchor_ok_1pct: pct_ok(me_anchor),
        vev_ok_1pct: pct_ok(v_cmp),
        lattice_mass_ok_1pct: lattice_mass_ok,
        fermi_mass_ok_1pct: fermi_mass_ok,
        overall_pass: pct_ok(me_anchor) && pct_ok(v_cmp) && lattice_mass_ok && fermi_mass_ok,
    };

    let report = Report {
        inputs: Inputs {
            g_f: G_F,
            alpha_mz: ALPHA_EW_MZ,
            f0_vac,
            w_over_z_structural: w_z_mass_ratio(),
        },
        electron_anchor: me_anchor,
        vev_lattice_vs_fermi: v_cmp,
        masses_lattice_branch: Masses {
            m_w_gev: m_w_lat,
            m_z_gev: m_z_lat,
            m_h_gev: m_h_lat,
        },
        masses_fermi_branch: Masses {
            m_w_gev: m_w_fer,
            m_z_gev: m_z_fer,
            m_h_gev: m_h_fer,
        },
        checks,
    };

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[absolute_scale_endgame]").expect("write");
    writeln!(
        txt,
        "inputs: G_F={:.12e} alpha_mz={:.12} f0_vac={:.6} w_over_z_structural={:.12}",
        report.inputs.g_f,
        report.inputs.alpha_mz,
        report.inputs.f0_vac,
        report.inputs.w_over_z_structural
    )
    .expect("write");
    writeln!(txt).expect("write");
    writeln!(
        txt,
        "electron_anchor_mev: pred={:.12} ref={:.12} rel={:+.6e}",
        report.electron_anchor.predicted,
        report.electron_anchor.reference,
        report.electron_anchor.rel_err
    )
    .expect("write");
    writeln!(
        txt,
        "vev_lattice_vs_fermi_gev: pred={:.12} ref={:.12} rel={:+.6e}",
        report.vev_lattice_vs_fermi.predicted,
        report.vev_lattice_vs_fermi.reference,
        report.vev_lattice_vs_fermi.rel_err
    )
    .expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[lattice_branch_masses]").expect("write");
    writeln!(
        txt,
        "mW: pred={:.12} ref={:.12} rel={:+.6e}",
        report.masses_lattice_branch.m_w_gev.predicted,
        report.masses_lattice_branch.m_w_gev.reference,
        report.masses_lattice_branch.m_w_gev.rel_err
    )
    .expect("write");
    writeln!(
        txt,
        "mZ: pred={:.12} ref={:.12} rel={:+.6e}",
        report.masses_lattice_branch.m_z_gev.predicted,
        report.masses_lattice_branch.m_z_gev.reference,
        report.masses_lattice_branch.m_z_gev.rel_err
    )
    .expect("write");
    writeln!(
        txt,
        "mH: pred={:.12} ref={:.12} rel={:+.6e}",
        report.masses_lattice_branch.m_h_gev.predicted,
        report.masses_lattice_branch.m_h_gev.reference,
        report.masses_lattice_branch.m_h_gev.rel_err
    )
    .expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[fermi_branch_masses]").expect("write");
    writeln!(
        txt,
        "mW: pred={:.12} ref={:.12} rel={:+.6e}",
        report.masses_fermi_branch.m_w_gev.predicted,
        report.masses_fermi_branch.m_w_gev.reference,
        report.masses_fermi_branch.m_w_gev.rel_err
    )
    .expect("write");
    writeln!(
        txt,
        "mZ: pred={:.12} ref={:.12} rel={:+.6e}",
        report.masses_fermi_branch.m_z_gev.predicted,
        report.masses_fermi_branch.m_z_gev.reference,
        report.masses_fermi_branch.m_z_gev.rel_err
    )
    .expect("write");
    writeln!(
        txt,
        "mH: pred={:.12} ref={:.12} rel={:+.6e}",
        report.masses_fermi_branch.m_h_gev.predicted,
        report.masses_fermi_branch.m_h_gev.reference,
        report.masses_fermi_branch.m_h_gev.rel_err
    )
    .expect("write");
    writeln!(txt).expect("write");
    writeln!(
        txt,
        "checks: electron_1pct={} vev_1pct={} lattice_masses_1pct={} fermi_masses_1pct={} overall={}",
        report.checks.electron_anchor_ok_1pct,
        report.checks.vev_ok_1pct,
        report.checks.lattice_mass_ok_1pct,
        report.checks.fermi_mass_ok_1pct,
        report.checks.overall_pass
    )
    .expect("write");

    let mut json = File::create(&json_path).expect("create json");
    serde_json::to_writer_pretty(&mut json, &report).expect("write json");
    writeln!(json).expect("newline");

    println!(
        "absolute_endgame: overall_pass={} vev_rel={:+.3e} mW_lat_rel={:+.3e} mH_lat_rel={:+.3e}",
        report.checks.overall_pass,
        report.vev_lattice_vs_fermi.rel_err,
        report.masses_lattice_branch.m_w_gev.rel_err,
        report.masses_lattice_branch.m_h_gev.rel_err,
    );
    println!("wrote {txt_path}");
    println!("wrote {json_path}");
}
