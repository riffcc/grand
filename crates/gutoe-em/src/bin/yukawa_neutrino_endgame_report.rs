//! Endgame neutrino lane report:
//! - structural no-fit lane
//! - oscillation-triangulated lane
//! - Koide diagnostics (including K_nu = 1/2 probe)

use gutoe_em::alpha::koide_ratio;
use gutoe_em::{
    neutrino_absolute_masses_from_texture, neutrino_hierarchy_prediction,
    triangulate_neutrino_from_splittings,
};
use serde::Serialize;
use std::fs::{self, File};
use std::io::Write;

const DM21_TARGET_EV2: f64 = 7.53e-5;
const DM32_TARGET_EV2: f64 = 2.453e-3;
const RATIO_TOL: f64 = 0.05;
const ABS_TOL: f64 = 0.05;
const K_NU_TARGET: f64 = 0.5;

#[derive(Debug, Clone, Copy, Serialize)]
struct Lane {
    m1_ev: f64,
    m2_ev: f64,
    m3_ev: f64,
    sum_ev: f64,
    dm21_ev2: f64,
    dm32_ev2: f64,
    dm21_rel_err: f64,
    dm32_rel_err: f64,
    ratio_32_over_21: f64,
    ratio_rel_err: f64,
    koide_k: f64,
    s2: f64,
    k_nu_rel_err: f64,
}

#[derive(Debug, Clone, Copy, Serialize)]
struct Checks {
    hierarchy_ok: bool,
    tiny_ok: bool,
    ratio_ok: bool,
    abs_splittings_ok: bool,
    no_fit_pass: bool,
    triangulated_pass: bool,
}

#[derive(Debug, Clone, Serialize)]
struct Report {
    hierarchy_prediction: String,
    targets: Targets,
    structural: Lane,
    triangulated: Lane,
    triangulated_aux: TriAux,
    checks: Checks,
}

#[derive(Debug, Clone, Copy, Serialize)]
struct Targets {
    dm21_target_ev2: f64,
    dm32_target_ev2: f64,
    ratio_target: f64,
    k_nu_target: f64,
    ratio_tol: f64,
    abs_tol: f64,
}

#[derive(Debug, Clone, Copy, Serialize)]
struct TriAux {
    p_triangulated: f64,
    ratio_fit_rel_err: f64,
    kappa_dm21: f64,
    kappa_dm32: f64,
    kappa_geo: f64,
    kappa_consistency_rel: f64,
}

fn rel_err(obs: f64, target: f64) -> f64 {
    if target.abs() < 1.0e-30 {
        0.0
    } else {
        (obs - target) / target
    }
}

fn lane_from_masses(m1: f64, m2: f64, m3: f64) -> Lane {
    let sum_ev = m1 + m2 + m3;
    let dm21_ev2 = m2 * m2 - m1 * m1;
    let dm32_ev2 = m3 * m3 - m2 * m2;
    let ratio = dm32_ev2.abs() / dm21_ev2.abs().max(1.0e-30);
    let k = koide_ratio([m1, m2, m3]);
    let s2 = 6.0 * k - 2.0;
    Lane {
        m1_ev: m1,
        m2_ev: m2,
        m3_ev: m3,
        sum_ev,
        dm21_ev2,
        dm32_ev2,
        dm21_rel_err: rel_err(dm21_ev2, DM21_TARGET_EV2),
        dm32_rel_err: rel_err(dm32_ev2.abs(), DM32_TARGET_EV2),
        ratio_32_over_21: ratio,
        ratio_rel_err: rel_err(ratio, DM32_TARGET_EV2 / DM21_TARGET_EV2),
        koide_k: k,
        s2,
        k_nu_rel_err: rel_err(k, K_NU_TARGET),
    }
}

fn main() {
    let out_dir = std::env::var("GUTOE_NEUTRINO_ENDGAME_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let txt_path = format!("{out_dir}/yukawa_neutrino_endgame_report.txt");
    let json_path = format!("{out_dir}/yukawa_neutrino_endgame_report.json");

    let hierarchy = neutrino_hierarchy_prediction().to_string();
    let abs = neutrino_absolute_masses_from_texture();
    let tri = triangulate_neutrino_from_splittings(DM21_TARGET_EV2, DM32_TARGET_EV2);

    let structural = lane_from_masses(abs.m1_ev, abs.m2_ev, abs.m3_ev);
    let triangulated = lane_from_masses(tri.m1_ev, tri.m2_ev, tri.m3_ev);

    let hierarchy_ok = hierarchy == "normal";
    let tiny_ok = structural.m3_ev < 0.8 && structural.sum_ev < 0.12;
    let ratio_ok = structural.ratio_rel_err.abs() <= RATIO_TOL;
    let abs_splittings_ok =
        structural.dm21_rel_err.abs() <= ABS_TOL && structural.dm32_rel_err.abs() <= ABS_TOL;
    let no_fit_pass = hierarchy_ok && tiny_ok && ratio_ok;
    let triangulated_pass =
        tri.ratio_fit_rel_err.abs() < 1.0e-9
            && triangulated.dm21_rel_err.abs() < 1.0e-9
            && triangulated.dm32_rel_err.abs() < 1.0e-9;

    let report = Report {
        hierarchy_prediction: hierarchy.clone(),
        targets: Targets {
            dm21_target_ev2: DM21_TARGET_EV2,
            dm32_target_ev2: DM32_TARGET_EV2,
            ratio_target: DM32_TARGET_EV2 / DM21_TARGET_EV2,
            k_nu_target: K_NU_TARGET,
            ratio_tol: RATIO_TOL,
            abs_tol: ABS_TOL,
        },
        structural,
        triangulated,
        triangulated_aux: TriAux {
            p_triangulated: tri.p_triangulated,
            ratio_fit_rel_err: tri.ratio_fit_rel_err,
            kappa_dm21: tri.kappa_dm21,
            kappa_dm32: tri.kappa_dm32,
            kappa_geo: tri.kappa_geo,
            kappa_consistency_rel: tri.kappa_consistency_rel,
        },
        checks: Checks {
            hierarchy_ok,
            tiny_ok,
            ratio_ok,
            abs_splittings_ok,
            no_fit_pass,
            triangulated_pass,
        },
    };

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[neutrino_endgame]").expect("write");
    writeln!(txt, "hierarchy_prediction={}", report.hierarchy_prediction).expect("write");
    writeln!(
        txt,
        "targets: dm21={:.12e} dm32={:.12e} ratio={:.12e} K_nu_target={:.12}",
        report.targets.dm21_target_ev2,
        report.targets.dm32_target_ev2,
        report.targets.ratio_target,
        report.targets.k_nu_target
    )
    .expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[structural_no_fit]").expect("write");
    writeln!(
        txt,
        "masses_ev=({:.12e}, {:.12e}, {:.12e}) sum={:.12e}",
        report.structural.m1_ev,
        report.structural.m2_ev,
        report.structural.m3_ev,
        report.structural.sum_ev
    )
    .expect("write");
    writeln!(
        txt,
        "dm21={:.12e} (rel={:+.6e}) dm32={:.12e} (rel={:+.6e})",
        report.structural.dm21_ev2,
        report.structural.dm21_rel_err,
        report.structural.dm32_ev2,
        report.structural.dm32_rel_err
    )
    .expect("write");
    writeln!(
        txt,
        "ratio32/21={:.12e} (rel={:+.6e}) K={:.12e} s2={:.12e} K_vs_half_rel={:+.6e}",
        report.structural.ratio_32_over_21,
        report.structural.ratio_rel_err,
        report.structural.koide_k,
        report.structural.s2,
        report.structural.k_nu_rel_err
    )
    .expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[triangulated]").expect("write");
    writeln!(
        txt,
        "p={:.12e} kappa_dm21={:.12e} kappa_dm32={:.12e} kappa_geo={:.12e} kappa_consistency_rel={:+.6e}",
        report.triangulated_aux.p_triangulated,
        report.triangulated_aux.kappa_dm21,
        report.triangulated_aux.kappa_dm32,
        report.triangulated_aux.kappa_geo,
        report.triangulated_aux.kappa_consistency_rel
    )
    .expect("write");
    writeln!(
        txt,
        "masses_ev=({:.12e}, {:.12e}, {:.12e}) sum={:.12e}",
        report.triangulated.m1_ev,
        report.triangulated.m2_ev,
        report.triangulated.m3_ev,
        report.triangulated.sum_ev
    )
    .expect("write");
    writeln!(
        txt,
        "dm21={:.12e} (rel={:+.6e}) dm32={:.12e} (rel={:+.6e})",
        report.triangulated.dm21_ev2,
        report.triangulated.dm21_rel_err,
        report.triangulated.dm32_ev2,
        report.triangulated.dm32_rel_err
    )
    .expect("write");
    writeln!(
        txt,
        "ratio32/21={:.12e} (rel={:+.6e}) K={:.12e} s2={:.12e} K_vs_half_rel={:+.6e}",
        report.triangulated.ratio_32_over_21,
        report.triangulated.ratio_rel_err,
        report.triangulated.koide_k,
        report.triangulated.s2,
        report.triangulated.k_nu_rel_err
    )
    .expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[checks]").expect("write");
    writeln!(
        txt,
        "hierarchy_ok={} tiny_ok={} ratio_ok={} abs_splittings_ok={} no_fit_pass={} triangulated_pass={}",
        report.checks.hierarchy_ok,
        report.checks.tiny_ok,
        report.checks.ratio_ok,
        report.checks.abs_splittings_ok,
        report.checks.no_fit_pass,
        report.checks.triangulated_pass
    )
    .expect("write");

    let mut json = File::create(&json_path).expect("create json");
    serde_json::to_writer_pretty(&mut json, &report).expect("write json");
    writeln!(json).expect("newline");

    println!(
        "neutrino_endgame: no_fit_pass={} tri_pass={} K_no_fit={:.6} K_tri={:.6}",
        report.checks.no_fit_pass,
        report.checks.triangulated_pass,
        report.structural.koide_k,
        report.triangulated.koide_k
    );
    println!("wrote {txt_path}");
    println!("wrote {json_path}");
}
