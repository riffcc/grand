//! Down-sector phase-conjugation test (requested follow-up).
//!
//! Test question:
//!   Does down-sector closure improve if we use a conjugated Z3 harmonic phase
//!   instead of trying CKM amplitude-space rotation?
//!
//! Models tested:
//!   1) baseline observed down masses
//!   2) free-s conjugate phases (delta -> -delta, pi-delta)
//!   3) fixed-sqrt2 conjugate phases (same transforms)
//! where closure is measured against structural targets (ms/md, mb/ms).

use gutoe_em::alpha::{z3_extract_params, z3_harmonic_masses};
use serde::Serialize;
use std::f64::consts::PI;
use std::fs;

const D: f64 = 4.67;
const S: f64 = 93.0;
const B: f64 = 4180.0;

const TARGET_MS_MD: f64 = 19.0;
const TARGET_MB_MS: f64 = (8.0 / 3.0) * 19.0 * (67.0 / 66.0);

const PERMS: [[usize; 3]; 6] = [
    [0, 1, 2],
    [0, 2, 1],
    [1, 0, 2],
    [1, 2, 0],
    [2, 0, 1],
    [2, 1, 0],
];

#[derive(Debug, Clone, Serialize)]
struct ClosureEval {
    name: String,
    source_masses: [f64; 3],
    best_perm: [usize; 3],
    best_masses_reordered: [f64; 3],
    ms_over_md: f64,
    mb_over_ms: f64,
    rms_log_closure: f64,
}

#[derive(Debug, Clone, Serialize)]
struct Report {
    down_fit_m_scale: f64,
    down_fit_s: f64,
    down_fit_s2: f64,
    down_fit_delta_deg: f64,
    target_ms_over_md: f64,
    target_mb_over_ms: f64,
    evaluations: Vec<ClosureEval>,
    best_name: String,
    summary: String,
}

fn closure_score(ms_md: f64, mb_ms: f64) -> f64 {
    let x = (ms_md / TARGET_MS_MD).ln();
    let y = (mb_ms / TARGET_MB_MS).ln();
    ((x * x + y * y) / 2.0).sqrt()
}

fn best_reordered(m: [f64; 3]) -> ([usize; 3], [f64; 3], f64, f64, f64) {
    let mut best = None::<([usize; 3], [f64; 3], f64, f64, f64)>;
    for p in PERMS {
        let r = [m[p[0]], m[p[1]], m[p[2]]];
        if r[0] <= 0.0 || r[1] <= 0.0 || r[2] <= 0.0 {
            continue;
        }
        let ms_md = r[1] / r[0];
        let mb_ms = r[2] / r[1];
        let score = closure_score(ms_md, mb_ms);
        match &best {
            None => best = Some((p, r, ms_md, mb_ms, score)),
            Some((_, _, _, _, s_old)) if score < *s_old => {
                best = Some((p, r, ms_md, mb_ms, score))
            }
            _ => {}
        }
    }
    best.expect("at least one permutation")
}

fn eval(name: &str, masses: [f64; 3]) -> ClosureEval {
    let (perm, reordered, ms_md, mb_ms, score) = best_reordered(masses);
    ClosureEval {
        name: name.to_string(),
        source_masses: masses,
        best_perm: perm,
        best_masses_reordered: reordered,
        ms_over_md: ms_md,
        mb_over_ms: mb_ms,
        rms_log_closure: score,
    }
}

fn main() {
    let out_dir = std::env::var("GUTOE_YUKAWA_DOWN_PHASE_CONJ_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let txt_path = format!("{out_dir}/yukawa_down_phase_conjugation_test.txt");
    let json_path = format!("{out_dir}/yukawa_down_phase_conjugation_test.json");

    let down = [D, S, B];
    let (m_scale, s, delta) = z3_extract_params(down);

    let mut evaluations = Vec::new();
    evaluations.push(eval("baseline_observed", down));

    // Free-s conjugation transforms.
    let free_neg = z3_harmonic_masses(m_scale, s, -delta);
    evaluations.push(eval("free_s_delta_to_neg_delta", free_neg));

    let free_pi_minus = z3_harmonic_masses(m_scale, s, PI - delta);
    evaluations.push(eval("free_s_delta_to_pi_minus_delta", free_pi_minus));

    // Fixed s = sqrt(2) conjugation transforms.
    let s2 = std::f64::consts::SQRT_2;
    let sqrt2_neg = z3_harmonic_masses(m_scale, s2, -delta);
    evaluations.push(eval("sqrt2_delta_to_neg_delta", sqrt2_neg));

    let sqrt2_pi_minus = z3_harmonic_masses(m_scale, s2, PI - delta);
    evaluations.push(eval("sqrt2_delta_to_pi_minus_delta", sqrt2_pi_minus));

    let best = evaluations
        .iter()
        .min_by(|a, b| a.rms_log_closure.total_cmp(&b.rms_log_closure))
        .expect("non-empty");

    let best_name = best.name.clone();
    let best_rms = best.rms_log_closure;
    let baseline_rms = evaluations[0].rms_log_closure;
    let summary = format!(
        "best model={} rms_log={:.6}; baseline rms_log={:.6}.",
        best_name, best_rms, baseline_rms
    );

    let report = Report {
        down_fit_m_scale: m_scale,
        down_fit_s: s,
        down_fit_s2: s * s,
        down_fit_delta_deg: delta.to_degrees(),
        target_ms_over_md: TARGET_MS_MD,
        target_mb_over_ms: TARGET_MB_MS,
        evaluations,
        best_name,
        summary,
    };

    let mut txt = String::new();
    txt.push_str("[yukawa_down_phase_conjugation_test]\n");
    txt.push_str(&format!(
        "down_fit: M={:.9} s={:.9} s2={:.9} delta_deg={:.9}\n",
        report.down_fit_m_scale, report.down_fit_s, report.down_fit_s2, report.down_fit_delta_deg
    ));
    txt.push_str(&format!(
        "targets: ms/md={:.9} mb/ms={:.9}\n\n",
        report.target_ms_over_md, report.target_mb_over_ms
    ));
    for e in &report.evaluations {
        txt.push_str(&format!(
            "[{}]\nperm={:?} reordered={:?}\nms/md={:.9} mb/ms={:.9} rms_log={:.9}\n\n",
            e.name,
            e.best_perm,
            e.best_masses_reordered,
            e.ms_over_md,
            e.mb_over_ms,
            e.rms_log_closure
        ));
    }
    txt.push_str(&format!("best_name={}\nsummary={}\n", report.best_name, report.summary));

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
