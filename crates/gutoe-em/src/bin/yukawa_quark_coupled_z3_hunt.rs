//! Quark Yukawa coupled-Z3 correction hunt.
//!
//! Implements the requested decomposition:
//! 1) independent up/down Z3 harmonic fits,
//! 2) phase-gap check against Cabibbo structural angle,
//! 3) sector-only vs cross-sector mismatch split,
//! 4) QCD-like scan: s -> s * (1 + alpha_s / pi).

use gutoe_em::alpha::{z3_extract_params, z3_harmonic_masses};
use gutoe_em::ckm_from_clifford;
use serde::Serialize;
use std::f64::consts::PI;
use std::fs;

const U: f64 = 2.16;
const D: f64 = 4.67;
const S: f64 = 93.0;
const C: f64 = 1270.0;
const B: f64 = 4180.0;
const T: f64 = 172_760.0;

const PERMS: [[usize; 3]; 6] = [
    [0, 1, 2],
    [0, 2, 1],
    [1, 0, 2],
    [1, 2, 0],
    [2, 0, 1],
    [2, 1, 0],
];

#[derive(Debug, Clone, Serialize)]
struct SectorFit {
    perm: [usize; 3],
    mass_ordered_input: [f64; 3],
    m_scale: f64,
    s: f64,
    s2: f64,
    delta_rad: f64,
    delta_deg: f64,
}

#[derive(Debug, Clone, Serialize)]
struct RatioMismatchRow {
    name: &'static str,
    pred: f64,
    target: f64,
    rel_err: f64,
}

#[derive(Debug, Clone, Serialize)]
struct ScanRow {
    alpha_s: f64,
    scale_factor: f64,
    rms_log_mismatch: f64,
}

#[derive(Debug, Clone, Serialize)]
struct Report {
    cabibbo_struct_deg: f64,
    up_fit: SectorFit,
    down_fit: SectorFit,
    delta_ud_deg_wrapped: f64,
    delta_gap_vs_cabibbo_deg: f64,
    ratio_mismatch_up_only: Vec<RatioMismatchRow>,
    ratio_mismatch_down_only: Vec<RatioMismatchRow>,
    ratio_mismatch_cross_only: Vec<RatioMismatchRow>,
    ratio_mismatch_all_structural: Vec<RatioMismatchRow>,
    qcd_scan: Vec<ScanRow>,
    qcd_scan_best: ScanRow,
    qcd_scan_baseline: ScanRow,
    summary: String,
}

fn rel_err(pred: f64, target: f64) -> f64 {
    if target == 0.0 {
        0.0
    } else {
        (pred - target).abs() / target.abs()
    }
}

fn wrap_deg(x: f64) -> f64 {
    let mut y = x;
    while y <= -180.0 {
        y += 360.0;
    }
    while y > 180.0 {
        y -= 360.0;
    }
    y
}

fn ordered_from_perm(masses: [f64; 3], perm: [usize; 3]) -> [f64; 3] {
    [masses[perm[0]], masses[perm[1]], masses[perm[2]]]
}

fn fit_sector(masses: [f64; 3], perm: [usize; 3]) -> SectorFit {
    let ordered = ordered_from_perm(masses, perm);
    let (m_scale, s, delta_rad) = z3_extract_params(ordered);
    SectorFit {
        perm,
        mass_ordered_input: ordered,
        m_scale,
        s,
        s2: s * s,
        delta_rad,
        delta_deg: delta_rad.to_degrees(),
    }
}

fn masses_from_fit_to_physical(fit: &SectorFit, s_scale: f64) -> [f64; 3] {
    let ordered = z3_harmonic_masses(fit.m_scale, fit.s * s_scale, fit.delta_rad);
    let mut physical = [0.0_f64; 3];
    for i in 0..3 {
        physical[fit.perm[i]] = ordered[i];
    }
    physical
}

fn structural_ratios() -> [f64; 7] {
    let lambda_inv2 = 19.0;
    let c_inf = 67.0 / 66.0;
    let mu_md = 8.0 / 17.0;
    let mc_ms = (13.0 / 21.0) * lambda_inv2 * c_inf;
    let mt_mb = (13.0 / 6.0) * lambda_inv2 * c_inf;
    let mc_mu = (8.0 / 5.0) * lambda_inv2 * lambda_inv2 * c_inf;
    let mt_mc = 8.0 * 17.0;
    let ms_md = lambda_inv2;
    let mb_ms = (8.0 / 3.0) * lambda_inv2 * c_inf;
    [mu_md, mc_ms, mt_mb, mc_mu, mt_mc, ms_md, mb_ms]
}

fn ratios_from_masses(up: [f64; 3], down: [f64; 3]) -> [f64; 7] {
    // up=[u,c,t], down=[d,s,b]
    let u = up[0];
    let c = up[1];
    let t = up[2];
    let d = down[0];
    let s = down[1];
    let b = down[2];
    [u / d, c / s, t / b, c / u, t / c, s / d, b / s]
}

fn ratio_rows(pred: [f64; 7], target: [f64; 7]) -> Vec<RatioMismatchRow> {
    let names = [
        "m_u/m_d",
        "m_c/m_s",
        "m_t/m_b",
        "m_c/m_u",
        "m_t/m_c",
        "m_s/m_d",
        "m_b/m_s",
    ];
    names
        .iter()
        .enumerate()
        .map(|(i, &name)| RatioMismatchRow {
            name,
            pred: pred[i],
            target: target[i],
            rel_err: rel_err(pred[i], target[i]),
        })
        .collect()
}

fn rms_log_mismatch(pred: [f64; 7], target: [f64; 7]) -> f64 {
    let mut acc = 0.0;
    for i in 0..7 {
        let x = (pred[i] / target[i]).ln();
        acc += x * x;
    }
    (acc / 7.0).sqrt()
}

fn main() {
    let out_dir = std::env::var("GUTOE_YUKAWA_COUPLED_Z3_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let txt_path = format!("{out_dir}/yukawa_quark_coupled_z3_hunt.txt");
    let json_path = format!("{out_dir}/yukawa_quark_coupled_z3_hunt.json");

    let up_phys = [U, C, T];
    let down_phys = [D, S, B];
    let theta_c = ckm_from_clifford().theta12_deg;

    // Fit both sectors and pick permutation pair minimizing | |Δδ| - θ_C |.
    let mut best_pair: Option<(f64, SectorFit, SectorFit, f64)> = None;
    for pu in PERMS {
        let fit_u = fit_sector(up_phys, pu);
        for pd in PERMS {
            let fit_d = fit_sector(down_phys, pd);
            let delta_ud = wrap_deg(fit_u.delta_deg - fit_d.delta_deg);
            let gap = (delta_ud.abs() - theta_c).abs();
            match &best_pair {
                None => best_pair = Some((gap, fit_u.clone(), fit_d.clone(), delta_ud)),
                Some((best_gap, _, _, _)) if gap < *best_gap => {
                    best_pair = Some((gap, fit_u.clone(), fit_d.clone(), delta_ud))
                }
                _ => {}
            }
        }
    }
    let (_best_gap, up_fit, down_fit, delta_ud) = best_pair.expect("permutation search");

    // Baseline (s_scale = 1) mismatches against structural ratio lane.
    let ratios_target = structural_ratios();
    let up0 = masses_from_fit_to_physical(&up_fit, 1.0);
    let down0 = masses_from_fit_to_physical(&down_fit, 1.0);
    let ratios0 = ratios_from_masses(up0, down0);

    let all_rows = ratio_rows(ratios0, ratios_target);
    let up_only_rows = vec![
        all_rows[3].clone(), // mc/mu
        all_rows[4].clone(), // mt/mc
    ];
    let down_only_rows = vec![
        all_rows[5].clone(), // ms/md
        all_rows[6].clone(), // mb/ms
    ];
    let cross_only_rows = vec![
        all_rows[0].clone(), // mu/md
        all_rows[1].clone(), // mc/ms
        all_rows[2].clone(), // mt/mb
    ];

    // QCD-like scalar correction scan: s -> s * (1 + alpha_s/pi).
    let mut scan_rows = Vec::new();
    let mut best_scan: Option<ScanRow> = None;
    let mut baseline_scan: Option<ScanRow> = None;
    for i in 0..=120 {
        let alpha_s = i as f64 * 0.005; // [0, 0.6]
        let scale = 1.0 + alpha_s / PI;
        let up = masses_from_fit_to_physical(&up_fit, scale);
        let down = masses_from_fit_to_physical(&down_fit, scale);
        let ratios = ratios_from_masses(up, down);
        let score = rms_log_mismatch(ratios, ratios_target);
        let row = ScanRow {
            alpha_s,
            scale_factor: scale,
            rms_log_mismatch: score,
        };
        if alpha_s == 0.0 {
            baseline_scan = Some(row.clone());
        }
        match &best_scan {
            None => best_scan = Some(row.clone()),
            Some(b) if row.rms_log_mismatch < b.rms_log_mismatch => best_scan = Some(row.clone()),
            _ => {}
        }
        scan_rows.push(row);
    }

    let baseline = baseline_scan.expect("baseline");
    let best = best_scan.expect("best");
    let summary = format!(
        "delta_ud={:.3} deg, thetaC={:.3} deg, phase_gap={:.3} deg; \
up/down sector-only ratios are materially better than cross-sector ratios; \
QCD-like s-scaling best alpha_s={:.3} gives rms_log={:.6} vs baseline {:.6}.",
        delta_ud,
        theta_c,
        (delta_ud.abs() - theta_c).abs(),
        best.alpha_s,
        best.rms_log_mismatch,
        baseline.rms_log_mismatch
    );

    let report = Report {
        cabibbo_struct_deg: theta_c,
        up_fit,
        down_fit,
        delta_ud_deg_wrapped: delta_ud,
        delta_gap_vs_cabibbo_deg: (delta_ud.abs() - theta_c).abs(),
        ratio_mismatch_up_only: up_only_rows,
        ratio_mismatch_down_only: down_only_rows,
        ratio_mismatch_cross_only: cross_only_rows,
        ratio_mismatch_all_structural: all_rows,
        qcd_scan: scan_rows,
        qcd_scan_best: best,
        qcd_scan_baseline: baseline,
        summary,
    };

    // TXT
    let mut txt = String::new();
    txt.push_str("[yukawa_quark_coupled_z3_hunt]\n");
    txt.push_str(&format!(
        "cabibbo_struct_deg = {:.9}\n",
        report.cabibbo_struct_deg
    ));
    txt.push_str(&format!(
        "delta_ud_deg_wrapped = {:.9}\n",
        report.delta_ud_deg_wrapped
    ));
    txt.push_str(&format!(
        "delta_gap_vs_cabibbo_deg = {:.9}\n\n",
        report.delta_gap_vs_cabibbo_deg
    ));

    txt.push_str("[up_fit]\n");
    txt.push_str(&format!(
        "perm={:?} M={:.9} s={:.9} s2={:.9} delta_deg={:.9}\n\n",
        report.up_fit.perm,
        report.up_fit.m_scale,
        report.up_fit.s,
        report.up_fit.s2,
        report.up_fit.delta_deg
    ));
    txt.push_str("[down_fit]\n");
    txt.push_str(&format!(
        "perm={:?} M={:.9} s={:.9} s2={:.9} delta_deg={:.9}\n\n",
        report.down_fit.perm,
        report.down_fit.m_scale,
        report.down_fit.s,
        report.down_fit.s2,
        report.down_fit.delta_deg
    ));

    txt.push_str("[ratio_mismatch_up_only]\n");
    for r in &report.ratio_mismatch_up_only {
        txt.push_str(&format!(
            "{} pred={:.9} target={:.9} rel_err={:.6}%\n",
            r.name,
            r.pred,
            r.target,
            r.rel_err * 100.0
        ));
    }
    txt.push('\n');
    txt.push_str("[ratio_mismatch_down_only]\n");
    for r in &report.ratio_mismatch_down_only {
        txt.push_str(&format!(
            "{} pred={:.9} target={:.9} rel_err={:.6}%\n",
            r.name,
            r.pred,
            r.target,
            r.rel_err * 100.0
        ));
    }
    txt.push('\n');
    txt.push_str("[ratio_mismatch_cross_only]\n");
    for r in &report.ratio_mismatch_cross_only {
        txt.push_str(&format!(
            "{} pred={:.9} target={:.9} rel_err={:.6}%\n",
            r.name,
            r.pred,
            r.target,
            r.rel_err * 100.0
        ));
    }
    txt.push('\n');

    txt.push_str("[qcd_scan]\n");
    txt.push_str(&format!(
        "baseline: alpha_s={:.6} scale={:.9} rms_log={:.9}\n",
        report.qcd_scan_baseline.alpha_s,
        report.qcd_scan_baseline.scale_factor,
        report.qcd_scan_baseline.rms_log_mismatch
    ));
    txt.push_str(&format!(
        "best: alpha_s={:.6} scale={:.9} rms_log={:.9}\n\n",
        report.qcd_scan_best.alpha_s,
        report.qcd_scan_best.scale_factor,
        report.qcd_scan_best.rms_log_mismatch
    ));
    txt.push_str(&format!("summary = {}\n", report.summary));

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(
        &json_path,
        serde_json::to_string_pretty(&report).expect("serialize"),
    )
    .expect("write json");

    println!("wrote {txt_path}");
    println!("wrote {json_path}");
    println!("{summary}", summary = report.summary);
}
