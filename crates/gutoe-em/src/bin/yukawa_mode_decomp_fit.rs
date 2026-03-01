//! Decomposition fit for Yukawa full-dynamics scan.
//!
//! Consumes `/tmp/bh_renders/yukawa_full_dynamics_scan.csv` (or env override)
//! and fits:
//!   1) s²_down plateau lane
//!   2) s²_Lg UV lane (gap to 3) with alpha_s and asymptotic forms
//!   3) up/down split lane with accumulated/incremental Yukawa diagnostics
//!
//! Sharp tests added:
//! - gap_lg = 3 - s²_Lg ~ c1 * alpha_s (origin-constrained) and free affine fit
//! - d(Δ_ud)/d(ln μ) vs y_t²(midpoint)
//! - Δ_QCD = s²_down - 2 proximity to 4/9 = C_F/N_c

use serde::Serialize;
use std::collections::HashMap;
use std::f64::consts::{PI, SQRT_2};
use std::fs;

const G_F: f64 = 1.166_378_7e-5;

#[derive(Debug, Clone, Serialize)]
struct InputRow {
    mu_gev: f64,
    alpha_s_mu: f64,
    m_u: f64,
    m_d: f64,
    m_s: f64,
    m_c: f64,
    m_b: f64,
    m_t: f64,
}

#[derive(Debug, Clone, Serialize)]
struct DecompRow {
    mu_gev: f64,
    alpha_s_mu: f64,
    y_t: f64,
    iy_t2_dlnmu: f64,
    k_up: f64,
    k_down: f64,
    k_lg: f64,
    s2_up: f64,
    s2_down: f64,
    s2_lg: f64,
    delta_ud: f64,
    gap_lg_to_3: f64,
}

#[derive(Debug, Clone, Serialize)]
struct SegmentRow {
    mu_lo_gev: f64,
    mu_hi_gev: f64,
    dlnmu: f64,
    y_t2_mid: f64,
    ddelta_ud_dlnmu: f64,
}

#[derive(Debug, Clone, Serialize)]
struct LinearFit {
    intercept: f64,
    slope: f64,
    r2: f64,
    rmse: f64,
}

#[derive(Debug, Clone, Serialize)]
struct OriginFit {
    slope: f64,
    rmse: f64,
}

#[derive(Debug, Clone, Serialize)]
struct ConstantFit {
    mean: f64,
    rmse: f64,
    max_abs_dev: f64,
}

#[derive(Debug, Clone, Serialize)]
struct AsymptoteFit {
    c: f64,
    p: f64,
    r2_log: f64,
    rmse_s2: f64,
}

#[derive(Debug, Clone, Serialize)]
struct CasimirCheck {
    delta_qcd: f64,
    cf_over_nc: f64,
    delta_minus_cf_over_nc: f64,
    rel_percent_vs_cf_over_nc: f64,
    kappa_mt_in_delta_eq_kappa_4as_over_pi: f64,
    kappa_uv_in_delta_eq_kappa_4as_over_pi: f64,
}

#[derive(Debug, Clone, Serialize)]
struct Report {
    source_csv: String,
    output_dir: String,
    vev_gev: f64,
    rows: Vec<DecompRow>,
    segment_rows: Vec<SegmentRow>,
    monotonic_s2_lg_increasing: bool,
    monotonic_gap_lg_decreasing: bool,
    monotonic_delta_ud_increasing: bool,
    down_plateau_constant: ConstantFit,
    down_vs_alpha_s: LinearFit,
    lg_gap_vs_alpha_s: LinearFit,
    lg_gap_vs_alpha_s_origin: OriginFit,
    delta_ud_vs_iyt2: LinearFit,
    delta_ud_vs_yt2: LinearFit,
    delta_ud_vs_yt: LinearFit,
    delta_ud_derivative_vs_yt2_mid: LinearFit,
    lg_gap_asymptote: AsymptoteFit,
    casimir_check: CasimirCheck,
    summary: String,
}

fn electroweak_vev_from_fermi(gf: f64) -> f64 {
    1.0 / (2.0_f64.sqrt() * gf).sqrt()
}

fn koide(vals: [f64; 3]) -> f64 {
    let num = vals[0] + vals[1] + vals[2];
    let den = (vals[0].sqrt() + vals[1].sqrt() + vals[2].sqrt()).powi(2);
    num / den
}

fn s2_from_koide(k: f64) -> f64 {
    6.0 * k - 2.0
}

fn linear_fit(xs: &[f64], ys: &[f64]) -> LinearFit {
    assert!(!xs.is_empty() && xs.len() == ys.len(), "linear_fit inputs");
    let n = xs.len();
    let mx = xs.iter().sum::<f64>() / n as f64;
    let my = ys.iter().sum::<f64>() / n as f64;
    let num = xs
        .iter()
        .zip(ys.iter())
        .map(|(x, y)| (x - mx) * (y - my))
        .sum::<f64>();
    let den = xs.iter().map(|x| (x - mx).powi(2)).sum::<f64>();
    let slope = if den > 0.0 { num / den } else { 0.0 };
    let intercept = my - slope * mx;

    let preds: Vec<f64> = xs.iter().map(|x| intercept + slope * x).collect();
    let mse = ys
        .iter()
        .zip(preds.iter())
        .map(|(y, p)| (y - p).powi(2))
        .sum::<f64>()
        / n as f64;
    let rmse = mse.sqrt();
    let ss_tot = ys.iter().map(|y| (y - my).powi(2)).sum::<f64>();
    let ss_res = ys
        .iter()
        .zip(preds.iter())
        .map(|(y, p)| (y - p).powi(2))
        .sum::<f64>();
    let r2 = if ss_tot > 0.0 { 1.0 - ss_res / ss_tot } else { 1.0 };

    LinearFit {
        intercept,
        slope,
        r2,
        rmse,
    }
}

fn origin_fit(xs: &[f64], ys: &[f64]) -> OriginFit {
    assert!(!xs.is_empty() && xs.len() == ys.len(), "origin_fit inputs");
    let den = xs.iter().map(|x| x * x).sum::<f64>();
    let num = xs.iter().zip(ys.iter()).map(|(x, y)| x * y).sum::<f64>();
    let slope = if den > 0.0 { num / den } else { 0.0 };
    let rmse = (xs
        .iter()
        .zip(ys.iter())
        .map(|(x, y)| (y - slope * x).powi(2))
        .sum::<f64>()
        / xs.len() as f64)
        .sqrt();
    OriginFit { slope, rmse }
}

fn constant_fit(ys: &[f64]) -> ConstantFit {
    assert!(!ys.is_empty(), "constant_fit inputs");
    let n = ys.len() as f64;
    let mean = ys.iter().sum::<f64>() / n;
    let mse = ys.iter().map(|y| (y - mean).powi(2)).sum::<f64>() / n;
    let rmse = mse.sqrt();
    let max_abs_dev = ys
        .iter()
        .map(|y| (y - mean).abs())
        .fold(0.0_f64, |a, b| a.max(b));
    ConstantFit {
        mean,
        rmse,
        max_abs_dev,
    }
}

fn asymptote_fit(mu: &[f64], s2_lg: &[f64]) -> AsymptoteFit {
    // Fit ln(3 - s2) = ln(C) - p ln(ln(mu))
    let mut xs = Vec::new();
    let mut ys = Vec::new();
    for (&m, &s2) in mu.iter().zip(s2_lg.iter()) {
        let gap = 3.0 - s2;
        if m > 1.0 && gap > 0.0 {
            xs.push((m.ln()).ln());
            ys.push(gap.ln());
        }
    }
    let lf = linear_fit(&xs, &ys);
    let p = -lf.slope;
    let c = lf.intercept.exp();

    let mut err2 = 0.0;
    let mut cnt = 0usize;
    for (&m, &s2) in mu.iter().zip(s2_lg.iter()) {
        if m > 1.0 {
            let pred_gap = c * m.ln().powf(-p);
            let pred_s2 = 3.0 - pred_gap;
            err2 += (pred_s2 - s2).powi(2);
            cnt += 1;
        }
    }
    let rmse_s2 = if cnt > 0 {
        (err2 / cnt as f64).sqrt()
    } else {
        0.0
    };

    AsymptoteFit {
        c,
        p,
        r2_log: lf.r2,
        rmse_s2,
    }
}

fn parse_csv(path: &str) -> Vec<InputRow> {
    let raw = fs::read_to_string(path).expect("read source csv");
    let mut lines = raw.lines();
    let header = lines.next().expect("csv header");
    let cols: Vec<&str> = header.split(',').collect();
    let mut idx = HashMap::new();
    for (i, c) in cols.iter().enumerate() {
        idx.insert(*c, i);
    }

    let req = [
        "mu_gev",
        "alpha_s_mu",
        "m_u",
        "m_d",
        "m_s",
        "m_c",
        "m_b",
        "m_t",
    ];
    for k in &req {
        assert!(idx.contains_key(k), "missing csv column: {k}");
    }

    let mut out = Vec::new();
    for ln in lines {
        let line = ln.trim();
        if line.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split(',').collect();
        let g = |k: &str| -> f64 {
            let j = *idx.get(k).expect("column idx");
            parts[j].parse::<f64>().expect("parse f64")
        };
        out.push(InputRow {
            mu_gev: g("mu_gev"),
            alpha_s_mu: g("alpha_s_mu"),
            m_u: g("m_u"),
            m_d: g("m_d"),
            m_s: g("m_s"),
            m_c: g("m_c"),
            m_b: g("m_b"),
            m_t: g("m_t"),
        });
    }
    out.sort_by(|a, b| a.mu_gev.total_cmp(&b.mu_gev));
    out
}

fn main() {
    let source_csv = std::env::var("GUTOE_YUKAWA_FULL_DYN_CSV")
        .unwrap_or_else(|_| "/tmp/bh_renders/yukawa_full_dynamics_scan.csv".to_string());
    let out_dir =
        std::env::var("GUTOE_YUKAWA_DECOMP_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);

    let rows_in = parse_csv(&source_csv);
    let vev = electroweak_vev_from_fermi(G_F);

    let mut rows = Vec::new();
    let mut iy = 0.0_f64;
    let mut prev_mu: Option<f64> = None;
    let mut prev_yt2: Option<f64> = None;

    for r in &rows_in {
        let k_up = koide([r.m_u, r.m_c, r.m_t]);
        let k_down = koide([r.m_d, r.m_s, r.m_b]);
        let l1 = (r.m_u * r.m_d).sqrt();
        let l2 = (r.m_c * r.m_s).sqrt();
        let l3 = (r.m_t * r.m_b).sqrt();
        let k_lg = koide([l1, l2, l3]);

        let s2_up = s2_from_koide(k_up);
        let s2_down = s2_from_koide(k_down);
        let s2_lg = s2_from_koide(k_lg);
        let delta_ud = s2_up - s2_down;
        let gap_lg = 3.0 - s2_lg;
        let y_t = SQRT_2 * (r.m_t / 1000.0) / vev;
        let y_t2 = y_t * y_t;

        if let (Some(m0), Some(y20)) = (prev_mu, prev_yt2) {
            let dt = r.mu_gev.ln() - m0.ln();
            iy += 0.5 * (y20 + y_t2) * dt;
        }
        prev_mu = Some(r.mu_gev);
        prev_yt2 = Some(y_t2);

        rows.push(DecompRow {
            mu_gev: r.mu_gev,
            alpha_s_mu: r.alpha_s_mu,
            y_t,
            iy_t2_dlnmu: iy,
            k_up,
            k_down,
            k_lg,
            s2_up,
            s2_down,
            s2_lg,
            delta_ud,
            gap_lg_to_3: gap_lg,
        });
    }

    let mu: Vec<f64> = rows.iter().map(|r| r.mu_gev).collect();
    let asv: Vec<f64> = rows.iter().map(|r| r.alpha_s_mu).collect();
    let ytv: Vec<f64> = rows.iter().map(|r| r.y_t).collect();
    let yt2v: Vec<f64> = rows.iter().map(|r| r.y_t * r.y_t).collect();
    let iyv: Vec<f64> = rows.iter().map(|r| r.iy_t2_dlnmu).collect();
    let s2d: Vec<f64> = rows.iter().map(|r| r.s2_down).collect();
    let s2lg: Vec<f64> = rows.iter().map(|r| r.s2_lg).collect();
    let dud: Vec<f64> = rows.iter().map(|r| r.delta_ud).collect();
    let gap_lg: Vec<f64> = rows.iter().map(|r| r.gap_lg_to_3).collect();

    let down_plateau = constant_fit(&s2d);
    let down_vs_alpha = linear_fit(&asv, &s2d);
    let lg_gap_vs_alpha = linear_fit(&asv, &gap_lg);
    let lg_gap_vs_alpha_origin = origin_fit(&asv, &gap_lg);
    let delta_vs_iy = linear_fit(&iyv, &dud);
    let delta_vs_yt2 = linear_fit(&yt2v, &dud);
    let delta_vs_yt = linear_fit(&ytv, &dud);
    let lg_asym = asymptote_fit(&mu, &s2lg);

    let mut segment_rows = Vec::new();
    for i in 1..rows.len() {
        let lo = &rows[i - 1];
        let hi = &rows[i];
        let dlnmu = hi.mu_gev.ln() - lo.mu_gev.ln();
        let y_t2_mid = 0.5 * (lo.y_t * lo.y_t + hi.y_t * hi.y_t);
        let ddelta_ud_dlnmu = (hi.delta_ud - lo.delta_ud) / dlnmu;
        segment_rows.push(SegmentRow {
            mu_lo_gev: lo.mu_gev,
            mu_hi_gev: hi.mu_gev,
            dlnmu,
            y_t2_mid,
            ddelta_ud_dlnmu,
        });
    }
    let ddelta_vec: Vec<f64> = segment_rows.iter().map(|s| s.ddelta_ud_dlnmu).collect();
    let yt2_mid_vec: Vec<f64> = segment_rows.iter().map(|s| s.y_t2_mid).collect();
    let delta_deriv_vs_yt2 = linear_fit(&yt2_mid_vec, &ddelta_vec);

    let monotonic_s2_lg_increasing = s2lg.windows(2).all(|w| w[1] >= w[0]);
    let monotonic_gap_lg_decreasing = gap_lg.windows(2).all(|w| w[1] <= w[0]);
    let monotonic_delta_ud_increasing = dud.windows(2).all(|w| w[1] >= w[0]);

    let delta_qcd = down_plateau.mean - 2.0;
    let cf_over_nc = 4.0 / 9.0;
    let delta_minus_cf_over_nc = delta_qcd - cf_over_nc;
    let rel_percent_vs_cf_over_nc = 100.0 * delta_minus_cf_over_nc / cf_over_nc;
    let as_mt = rows.first().map(|r| r.alpha_s_mu).unwrap_or(0.0);
    let as_uv = rows.last().map(|r| r.alpha_s_mu).unwrap_or(0.0);
    let kappa_mt = if as_mt > 0.0 {
        delta_qcd / (4.0 * as_mt / PI)
    } else {
        0.0
    };
    let kappa_uv = if as_uv > 0.0 {
        delta_qcd / (4.0 * as_uv / PI)
    } else {
        0.0
    };
    let casimir_check = CasimirCheck {
        delta_qcd,
        cf_over_nc,
        delta_minus_cf_over_nc,
        rel_percent_vs_cf_over_nc,
        kappa_mt_in_delta_eq_kappa_4as_over_pi: kappa_mt,
        kappa_uv_in_delta_eq_kappa_4as_over_pi: kappa_uv,
    };

    let summary = format!(
        "mode decomposition: s2_down mean={:.9}, s2_lg(1e19)={:.9}, delta_ud(1e19)={:.9}; c1(origin gap~alpha_s)={:.6}; R2(dDelta/dlnmu~yt2mid)={:.6}",
        down_plateau.mean,
        s2lg[s2lg.len() - 1],
        dud[dud.len() - 1],
        lg_gap_vs_alpha_origin.slope,
        delta_deriv_vs_yt2.r2
    );

    let report = Report {
        source_csv: source_csv.clone(),
        output_dir: out_dir.clone(),
        vev_gev: vev,
        rows,
        segment_rows,
        monotonic_s2_lg_increasing,
        monotonic_gap_lg_decreasing,
        monotonic_delta_ud_increasing,
        down_plateau_constant: down_plateau,
        down_vs_alpha_s: down_vs_alpha,
        lg_gap_vs_alpha_s: lg_gap_vs_alpha,
        lg_gap_vs_alpha_s_origin: lg_gap_vs_alpha_origin,
        delta_ud_vs_iyt2: delta_vs_iy,
        delta_ud_vs_yt2: delta_vs_yt2,
        delta_ud_vs_yt: delta_vs_yt,
        delta_ud_derivative_vs_yt2_mid: delta_deriv_vs_yt2,
        lg_gap_asymptote: lg_asym,
        casimir_check,
        summary,
    };

    let txt_path = format!("{out_dir}/yukawa_mode_decomp_fit.txt");
    let csv_path = format!("{out_dir}/yukawa_mode_decomp_fit.csv");
    let json_path = format!("{out_dir}/yukawa_mode_decomp_fit.json");

    let mut txt = String::new();
    txt.push_str("[yukawa_mode_decomp_fit]\n");
    txt.push_str(&format!("source_csv = {}\n", report.source_csv));
    txt.push_str(&format!("vev_gev = {:.9}\n\n", report.vev_gev));
    txt.push_str("rows: mu, alpha_s, y_t, Iy, s2_up, s2_down, s2_lg, delta_ud, gap_lg\n");
    for r in &report.rows {
        txt.push_str(&format!(
            "{:.6e}, {:.9}, {:.9}, {:.9}, {:.9}, {:.9}, {:.9}, {:.9}, {:.9}\n",
            r.mu_gev,
            r.alpha_s_mu,
            r.y_t,
            r.iy_t2_dlnmu,
            r.s2_up,
            r.s2_down,
            r.s2_lg,
            r.delta_ud,
            r.gap_lg_to_3
        ));
    }

    txt.push_str("\n[fits]\n");
    txt.push_str(&format!(
        "s2_down constant: mean={:.9} rmse={:.9} max_abs_dev={:.9}\n",
        report.down_plateau_constant.mean,
        report.down_plateau_constant.rmse,
        report.down_plateau_constant.max_abs_dev
    ));
    txt.push_str(&format!(
        "s2_down vs alpha_s: intercept={:.9} slope={:.9} r2={:.6} rmse={:.9}\n",
        report.down_vs_alpha_s.intercept,
        report.down_vs_alpha_s.slope,
        report.down_vs_alpha_s.r2,
        report.down_vs_alpha_s.rmse
    ));
    txt.push_str(&format!(
        "gap_lg(=3-s2_lg) vs alpha_s: intercept={:.9} slope={:.9} r2={:.6} rmse={:.9}\n",
        report.lg_gap_vs_alpha_s.intercept,
        report.lg_gap_vs_alpha_s.slope,
        report.lg_gap_vs_alpha_s.r2,
        report.lg_gap_vs_alpha_s.rmse
    ));
    txt.push_str(&format!(
        "gap_lg = c1 * alpha_s (origin constrained): c1={:.9} rmse={:.9}\n",
        report.lg_gap_vs_alpha_s_origin.slope,
        report.lg_gap_vs_alpha_s_origin.rmse
    ));
    txt.push_str(&format!(
        "delta_ud vs Iy_t2: intercept={:.9} slope={:.9} r2={:.6} rmse={:.9}\n",
        report.delta_ud_vs_iyt2.intercept,
        report.delta_ud_vs_iyt2.slope,
        report.delta_ud_vs_iyt2.r2,
        report.delta_ud_vs_iyt2.rmse
    ));
    txt.push_str(&format!(
        "delta_ud vs y_t^2: intercept={:.9} slope={:.9} r2={:.6} rmse={:.9}\n",
        report.delta_ud_vs_yt2.intercept,
        report.delta_ud_vs_yt2.slope,
        report.delta_ud_vs_yt2.r2,
        report.delta_ud_vs_yt2.rmse
    ));
    txt.push_str(&format!(
        "delta_ud vs y_t: intercept={:.9} slope={:.9} r2={:.6} rmse={:.9}\n",
        report.delta_ud_vs_yt.intercept,
        report.delta_ud_vs_yt.slope,
        report.delta_ud_vs_yt.r2,
        report.delta_ud_vs_yt.rmse
    ));
    txt.push_str(&format!(
        "d(delta_ud)/dlnmu vs y_t^2(mid): intercept={:.9} slope={:.9} r2={:.6} rmse={:.9}\n",
        report.delta_ud_derivative_vs_yt2_mid.intercept,
        report.delta_ud_derivative_vs_yt2_mid.slope,
        report.delta_ud_derivative_vs_yt2_mid.r2,
        report.delta_ud_derivative_vs_yt2_mid.rmse
    ));
    txt.push_str(&format!(
        "gap_lg=3-s2_lg asymptote: c={:.9} p={:.9} r2_log={:.6} rmse_s2={:.9}\n",
        report.lg_gap_asymptote.c,
        report.lg_gap_asymptote.p,
        report.lg_gap_asymptote.r2_log,
        report.lg_gap_asymptote.rmse_s2
    ));

    txt.push_str("\n[casimir checks]\n");
    txt.push_str(&format!(
        "delta_qcd = s2_down_mean - 2 = {:.9}\n",
        report.casimir_check.delta_qcd
    ));
    txt.push_str(&format!(
        "cf_over_nc = 4/9 = {:.9}\n",
        report.casimir_check.cf_over_nc
    ));
    txt.push_str(&format!(
        "delta_qcd - 4/9 = {:.9} ({:+.3}%)\n",
        report.casimir_check.delta_minus_cf_over_nc,
        report.casimir_check.rel_percent_vs_cf_over_nc
    ));
    txt.push_str(&format!(
        "kappa at mt in delta = kappa*(4*alpha_s/pi): {:.9}\n",
        report.casimir_check.kappa_mt_in_delta_eq_kappa_4as_over_pi
    ));
    txt.push_str(&format!(
        "kappa at 1e19 in delta = kappa*(4*alpha_s/pi): {:.9}\n",
        report.casimir_check.kappa_uv_in_delta_eq_kappa_4as_over_pi
    ));

    txt.push_str("\n[segment derivatives]\n");
    txt.push_str("mu_lo, mu_hi, dlnmu, y_t2_mid, ddelta_ud_dlnmu\n");
    for s in &report.segment_rows {
        txt.push_str(&format!(
            "{:.6e}, {:.6e}, {:.9}, {:.9}, {:.9}\n",
            s.mu_lo_gev, s.mu_hi_gev, s.dlnmu, s.y_t2_mid, s.ddelta_ud_dlnmu
        ));
    }

    txt.push_str("\n[monotonic flags]\n");
    txt.push_str(&format!(
        "s2_lg_increasing = {}\n",
        report.monotonic_s2_lg_increasing
    ));
    txt.push_str(&format!(
        "gap_lg_decreasing = {}\n",
        report.monotonic_gap_lg_decreasing
    ));
    txt.push_str(&format!(
        "delta_ud_increasing = {}\n",
        report.monotonic_delta_ud_increasing
    ));
    txt.push_str(&format!("\nsummary = {}\n", report.summary));
    fs::write(&txt_path, txt).expect("write txt");

    let mut csv = String::new();
    csv.push_str("mu_gev,alpha_s_mu,y_t,iy_t2_dlnmu,s2_up,s2_down,s2_lg,delta_ud,gap_lg_to_3\n");
    for r in &report.rows {
        csv.push_str(&format!(
            "{:.6e},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9}\n",
            r.mu_gev,
            r.alpha_s_mu,
            r.y_t,
            r.iy_t2_dlnmu,
            r.s2_up,
            r.s2_down,
            r.s2_lg,
            r.delta_ud,
            r.gap_lg_to_3
        ));
    }
    fs::write(&csv_path, csv).expect("write csv");

    fs::write(&json_path, serde_json::to_string_pretty(&report).expect("serialize json"))
        .expect("write json");

    println!("wrote {txt_path}");
    println!("wrote {csv_path}");
    println!("wrote {json_path}");
    println!("{}", report.summary);
}
