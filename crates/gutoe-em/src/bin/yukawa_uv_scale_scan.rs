//! UV scale scan for quark Yukawa structure.
//!
//! Goal:
//! 1) Run all quark masses to a common renormalization scale μ (one-loop QCD,
//!    threshold matched at m_c, m_b, m_t).
//! 2) Decompose by generation into:
//!      L_g = sqrt(m_up,g * m_down,g)      (mass ladder mode)
//!      S_g = 0.5 * ln(m_up,g / m_down,g)   (isospin split mode)
//! 3) Evaluate:
//!      - Z3 closure on L_g (free fit and fixed-s^2=2 constrained fit)
//!      - S_g vs generation linearity.
//! 4) Scan μ in {m_t, 1e4, 1e8, 1e12, 1e16} GeV.

use gutoe_em::alpha::{z3_extract_params, z3_harmonic_masses};
use serde::Serialize;
use std::f64::consts::PI;
use std::fs;

// Input masses (MeV) and reference scales (GeV), aligned with project lanes.
const MU_REF_MEV: f64 = 2.16; // at 2 GeV
const MD_REF_MEV: f64 = 4.67; // at 2 GeV
const MS_REF_MEV: f64 = 93.0; // at 2 GeV
const MC_REF_MEV: f64 = 1270.0; // at m_c
const MB_REF_MEV: f64 = 4180.0; // at m_b
const MT_REF_MEV: f64 = 172_760.0; // near m_t (project convention)

const MU_REF_GEV: f64 = 2.0;
const MD_REF_GEV: f64 = 2.0;
const MS_REF_GEV: f64 = 2.0;
const MC_REF_GEV: f64 = 1.27;
const MB_REF_GEV: f64 = 4.18;
const MT_REF_GEV: f64 = 172.76;

const MZ_GEV: f64 = 91.1876;
const ALPHA_S_MZ: f64 = 0.118; // reference anchor

#[derive(Debug, Clone, Serialize)]
struct LambdaSet {
    lambda3: f64,
    lambda4: f64,
    lambda5: f64,
    lambda6: f64,
}

#[derive(Debug, Clone, Serialize)]
struct Z3Fit {
    m_scale: f64,
    s: f64,
    s2: f64,
    delta_deg: f64,
    masses_pred: [f64; 3],
    rms_rel: f64,
}

#[derive(Debug, Clone, Serialize)]
struct FixedSFit {
    s_fixed: f64,
    m_scale: f64,
    delta_deg: f64,
    masses_pred: [f64; 3],
    rms_rel: f64,
}

#[derive(Debug, Clone, Serialize)]
struct LinearFit {
    slope: f64,
    intercept: f64,
    rmse: f64,
    r2: f64,
    y_pred: [f64; 3],
}

#[derive(Debug, Clone, Serialize)]
struct ScaleRow {
    mu_gev: f64,
    masses_mev: [f64; 6], // [u,d,s,c,b,t]
    lg: [f64; 3],         // [g1,g2,g3] = [sqrt(ud), sqrt(cs), sqrt(tb)]
    sg: [f64; 3],         // [0.5 ln(u/d), 0.5 ln(c/s), 0.5 ln(t/b)]
    lg_z3_free: Z3Fit,
    lg_z3_fixed_s2_2: FixedSFit,
    sg_linear: LinearFit,
}

#[derive(Debug, Clone, Serialize)]
struct Report {
    input_masses_mev: [f64; 6],
    input_scales_gev: [f64; 6],
    alpha_s_mz: f64,
    thresholds_gev: [f64; 3],
    lambdas: LambdaSet,
    scan_mu_gev: Vec<f64>,
    rows: Vec<ScaleRow>,
    monotonic: Monotonic,
    summary: String,
}

#[derive(Debug, Clone, Serialize)]
struct Monotonic {
    lg_fixed_rms_nonincreasing: bool,
    sg_linear_rmse_nonincreasing: bool,
    down_split_signs: [i8; 3],
}

fn beta0(nf: i32) -> f64 {
    11.0 - 2.0 * nf as f64 / 3.0
}

fn alpha_s_one_loop(mu: f64, lambda_nf: f64, nf: i32) -> f64 {
    let b0 = beta0(nf);
    let x = (mu * mu / (lambda_nf * lambda_nf)).ln();
    4.0 * PI / (b0 * x)
}

fn infer_lambda_from_alpha(mu: f64, alpha: f64, nf: i32) -> f64 {
    let b0 = beta0(nf);
    mu * (-(2.0 * PI) / (b0 * alpha)).exp()
}

fn build_lambdas() -> LambdaSet {
    // Anchor in n_f=5 at M_Z.
    let lambda5 = infer_lambda_from_alpha(MZ_GEV, ALPHA_S_MZ, 5);

    // Match at m_b between n_f=5 and n_f=4.
    let a_mb_5 = alpha_s_one_loop(MB_REF_GEV, lambda5, 5);
    let lambda4 = infer_lambda_from_alpha(MB_REF_GEV, a_mb_5, 4);

    // Match at m_c between n_f=4 and n_f=3.
    let a_mc_4 = alpha_s_one_loop(MC_REF_GEV, lambda4, 4);
    let lambda3 = infer_lambda_from_alpha(MC_REF_GEV, a_mc_4, 3);

    // Match at m_t between n_f=5 and n_f=6.
    let a_mt_5 = alpha_s_one_loop(MT_REF_GEV, lambda5, 5);
    let lambda6 = infer_lambda_from_alpha(MT_REF_GEV, a_mt_5, 6);

    LambdaSet {
        lambda3,
        lambda4,
        lambda5,
        lambda6,
    }
}

fn alpha_s(mu: f64, l: &LambdaSet) -> f64 {
    if mu >= MT_REF_GEV {
        alpha_s_one_loop(mu, l.lambda6, 6)
    } else if mu >= MB_REF_GEV {
        alpha_s_one_loop(mu, l.lambda5, 5)
    } else if mu >= MC_REF_GEV {
        alpha_s_one_loop(mu, l.lambda4, 4)
    } else {
        alpha_s_one_loop(mu, l.lambda3, 3)
    }
}

fn gamma_exp(nf: i32) -> f64 {
    // One-loop quark-mass exponent: 12 / (33 - 2 n_f)
    12.0 / (33.0 - 2.0 * nf as f64)
}

fn run_mass_segment(m0: f64, mu0: f64, mu1: f64, nf: i32, l: &LambdaSet) -> f64 {
    if (mu1 - mu0).abs() < 1e-15 {
        return m0;
    }
    let a0 = alpha_s_one_loop(mu0, match nf {
        3 => l.lambda3,
        4 => l.lambda4,
        5 => l.lambda5,
        6 => l.lambda6,
        _ => panic!("unsupported nf"),
    }, nf);
    let a1 = alpha_s_one_loop(mu1, match nf {
        3 => l.lambda3,
        4 => l.lambda4,
        5 => l.lambda5,
        6 => l.lambda6,
        _ => panic!("unsupported nf"),
    }, nf);
    m0 * (a1 / a0).powf(gamma_exp(nf))
}

fn run_mass_to(mu0: f64, m0: f64, mu: f64, l: &LambdaSet) -> f64 {
    // Piecewise run with threshold crossings at m_c, m_b, m_t.
    if (mu - mu0).abs() < 1e-15 {
        return m0;
    }

    // Build ascending grid for upward running.
    let mut points = vec![mu0];
    for t in [MC_REF_GEV, MB_REF_GEV, MT_REF_GEV] {
        if t > mu0 && t < mu {
            points.push(t);
        }
    }
    points.push(mu);
    points.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let mut m = m0;
    let mut cur = points[0];
    for &next in points.iter().skip(1) {
        // choose nf by the midpoint scale of segment
        let mid = (cur * next).sqrt();
        let nf = if mid >= MT_REF_GEV {
            6
        } else if mid >= MB_REF_GEV {
            5
        } else if mid >= MC_REF_GEV {
            4
        } else {
            3
        };
        m = run_mass_segment(m, cur, next, nf, l);
        cur = next;
    }
    m
}

fn z3_fit(masses: [f64; 3]) -> Z3Fit {
    let (m_scale, s, delta) = z3_extract_params(masses);
    let pred = z3_harmonic_masses(m_scale, s, delta);
    let rms_rel = {
        let e0 = ((pred[0] - masses[0]) / masses[0]).powi(2);
        let e1 = ((pred[1] - masses[1]) / masses[1]).powi(2);
        let e2 = ((pred[2] - masses[2]) / masses[2]).powi(2);
        ((e0 + e1 + e2) / 3.0).sqrt()
    };
    Z3Fit {
        m_scale,
        s,
        s2: s * s,
        delta_deg: delta.to_degrees(),
        masses_pred: pred,
        rms_rel,
    }
}

fn fit_fixed_s(masses: [f64; 3], s_fixed: f64) -> FixedSFit {
    let a = [masses[0].sqrt(), masses[1].sqrt(), masses[2].sqrt()];
    let mut best_obj = f64::INFINITY;
    let mut best_m = 0.0;
    let mut best_delta = 0.0;

    // Dense deterministic scan in phase.
    let n = 20_000usize;
    for i in 0..n {
        let delta = -PI + (2.0 * PI) * (i as f64) / (n as f64);
        let b0 = 1.0 + s_fixed * (delta).cos();
        let b1 = 1.0 + s_fixed * (delta + 2.0 * PI / 3.0).cos();
        let b2 = 1.0 + s_fixed * (delta + 4.0 * PI / 3.0).cos();
        let denom = b0 * b0 + b1 * b1 + b2 * b2;
        if denom <= 1e-14 {
            continue;
        }
        let m = (a[0] * b0 + a[1] * b1 + a[2] * b2) / denom;
        if m <= 0.0 {
            continue;
        }
        let r0 = a[0] - m * b0;
        let r1 = a[1] - m * b1;
        let r2 = a[2] - m * b2;
        let obj = r0 * r0 + r1 * r1 + r2 * r2;
        if obj < best_obj {
            best_obj = obj;
            best_m = m;
            best_delta = delta;
        }
    }

    let pred = z3_harmonic_masses(best_m, s_fixed, best_delta);
    let rms_rel = {
        let e0 = ((pred[0] - masses[0]) / masses[0]).powi(2);
        let e1 = ((pred[1] - masses[1]) / masses[1]).powi(2);
        let e2 = ((pred[2] - masses[2]) / masses[2]).powi(2);
        ((e0 + e1 + e2) / 3.0).sqrt()
    };

    FixedSFit {
        s_fixed,
        m_scale: best_m,
        delta_deg: best_delta.to_degrees(),
        masses_pred: pred,
        rms_rel,
    }
}

fn linear_fit_sg(sg: [f64; 3]) -> LinearFit {
    let x = [1.0_f64, 2.0, 3.0];
    let y = sg;
    let mx = (x[0] + x[1] + x[2]) / 3.0;
    let my = (y[0] + y[1] + y[2]) / 3.0;
    let num = (x[0] - mx) * (y[0] - my) + (x[1] - mx) * (y[1] - my) + (x[2] - mx) * (y[2] - my);
    let den = (x[0] - mx).powi(2) + (x[1] - mx).powi(2) + (x[2] - mx).powi(2);
    let slope = num / den;
    let intercept = my - slope * mx;
    let y_pred = [
        intercept + slope * x[0],
        intercept + slope * x[1],
        intercept + slope * x[2],
    ];
    let mse = ((y[0] - y_pred[0]).powi(2) + (y[1] - y_pred[1]).powi(2) + (y[2] - y_pred[2]).powi(2)) / 3.0;
    let rmse = mse.sqrt();

    let ss_tot = (y[0] - my).powi(2) + (y[1] - my).powi(2) + (y[2] - my).powi(2);
    let ss_res = (y[0] - y_pred[0]).powi(2) + (y[1] - y_pred[1]).powi(2) + (y[2] - y_pred[2]).powi(2);
    let r2 = if ss_tot > 0.0 { 1.0 - ss_res / ss_tot } else { 1.0 };

    LinearFit {
        slope,
        intercept,
        rmse,
        r2,
        y_pred,
    }
}

fn sign_i8(x: f64) -> i8 {
    if x > 0.0 {
        1
    } else if x < 0.0 {
        -1
    } else {
        0
    }
}

fn main() {
    let out_dir = std::env::var("GUTOE_YUKAWA_UV_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let txt_path = format!("{out_dir}/yukawa_uv_scale_scan.txt");
    let json_path = format!("{out_dir}/yukawa_uv_scale_scan.json");
    let csv_path = format!("{out_dir}/yukawa_uv_scale_scan.csv");

    let lambdas = build_lambdas();
    let scan_mu = vec![MT_REF_GEV, 1e4, 1e8, 1e12, 1e16];

    let mut rows = Vec::new();

    for &mu in &scan_mu {
        let m_u = run_mass_to(MU_REF_GEV, MU_REF_MEV, mu, &lambdas);
        let m_d = run_mass_to(MD_REF_GEV, MD_REF_MEV, mu, &lambdas);
        let m_s = run_mass_to(MS_REF_GEV, MS_REF_MEV, mu, &lambdas);
        let m_c = run_mass_to(MC_REF_GEV, MC_REF_MEV, mu, &lambdas);
        let m_b = run_mass_to(MB_REF_GEV, MB_REF_MEV, mu, &lambdas);
        let m_t = run_mass_to(MT_REF_GEV, MT_REF_MEV, mu, &lambdas);

        // Generation mode decomposition.
        let lg = [
            (m_u * m_d).sqrt(),
            (m_c * m_s).sqrt(),
            (m_t * m_b).sqrt(),
        ];
        let sg = [
            0.5 * (m_u / m_d).ln(),
            0.5 * (m_c / m_s).ln(),
            0.5 * (m_t / m_b).ln(),
        ];

        let lg_fit = z3_fit(lg);
        let lg_fixed = fit_fixed_s(lg, std::f64::consts::SQRT_2);
        let sg_lin = linear_fit_sg(sg);

        rows.push(ScaleRow {
            mu_gev: mu,
            masses_mev: [m_u, m_d, m_s, m_c, m_b, m_t],
            lg,
            sg,
            lg_z3_free: lg_fit,
            lg_z3_fixed_s2_2: lg_fixed,
            sg_linear: sg_lin,
        });
    }

    let lg_monotone = rows
        .windows(2)
        .all(|w| w[1].lg_z3_fixed_s2_2.rms_rel <= w[0].lg_z3_fixed_s2_2.rms_rel + 1e-15);
    let sg_monotone = rows
        .windows(2)
        .all(|w| w[1].sg_linear.rmse <= w[0].sg_linear.rmse + 1e-15);

    let down_split_signs = [
        sign_i8(rows[0].sg[0]),
        sign_i8(rows[0].sg[1]),
        sign_i8(rows[0].sg[2]),
    ];

    let monotonic = Monotonic {
        lg_fixed_rms_nonincreasing: lg_monotone,
        sg_linear_rmse_nonincreasing: sg_monotone,
        down_split_signs,
    };

    let summary = format!(
        "scan complete: lg_fixed_nonincreasing={} sg_linear_nonincreasing={} S_signs_at_mt=[{},{},{}]",
        monotonic.lg_fixed_rms_nonincreasing,
        monotonic.sg_linear_rmse_nonincreasing,
        monotonic.down_split_signs[0],
        monotonic.down_split_signs[1],
        monotonic.down_split_signs[2]
    );

    let report = Report {
        input_masses_mev: [MU_REF_MEV, MD_REF_MEV, MS_REF_MEV, MC_REF_MEV, MB_REF_MEV, MT_REF_MEV],
        input_scales_gev: [MU_REF_GEV, MD_REF_GEV, MS_REF_GEV, MC_REF_GEV, MB_REF_GEV, MT_REF_GEV],
        alpha_s_mz: ALPHA_S_MZ,
        thresholds_gev: [MC_REF_GEV, MB_REF_GEV, MT_REF_GEV],
        lambdas,
        scan_mu_gev: scan_mu,
        rows,
        monotonic,
        summary,
    };

    // TXT
    let mut txt = String::new();
    txt.push_str("[yukawa_uv_scale_scan]\n");
    txt.push_str(&format!("alpha_s_mz = {:.9}\n", report.alpha_s_mz));
    txt.push_str(&format!(
        "lambdas: l3={:.9} l4={:.9} l5={:.9} l6={:.9}\n\n",
        report.lambdas.lambda3,
        report.lambdas.lambda4,
        report.lambdas.lambda5,
        report.lambdas.lambda6
    ));

    for row in &report.rows {
        txt.push_str(&format!("[mu = {:.6e} GeV]\n", row.mu_gev));
        txt.push_str(&format!(
            "masses_mev [u,d,s,c,b,t] = [{:.6}, {:.6}, {:.6}, {:.6}, {:.6}, {:.6}]\n",
            row.masses_mev[0],
            row.masses_mev[1],
            row.masses_mev[2],
            row.masses_mev[3],
            row.masses_mev[4],
            row.masses_mev[5],
        ));
        txt.push_str(&format!(
            "L_g = [{:.6}, {:.6}, {:.6}]\n",
            row.lg[0], row.lg[1], row.lg[2]
        ));
        txt.push_str(&format!(
            "S_g = [{:.6}, {:.6}, {:.6}]\n",
            row.sg[0], row.sg[1], row.sg[2]
        ));
        txt.push_str(&format!(
            "L_free: s2={:.9} delta={:.6} rms={:.3e}\n",
            row.lg_z3_free.s2, row.lg_z3_free.delta_deg, row.lg_z3_free.rms_rel
        ));
        txt.push_str(&format!(
            "L_fixed(s2=2): delta={:.6} rms={:.3e}\n",
            row.lg_z3_fixed_s2_2.delta_deg, row.lg_z3_fixed_s2_2.rms_rel
        ));
        txt.push_str(&format!(
            "S_linear: slope={:.9} intercept={:.9} rmse={:.9} r2={:.9}\n\n",
            row.sg_linear.slope,
            row.sg_linear.intercept,
            row.sg_linear.rmse,
            row.sg_linear.r2
        ));
    }

    txt.push_str(&format!(
        "monotonic: L_fixed_nonincreasing={} S_linear_nonincreasing={}\n",
        report.monotonic.lg_fixed_rms_nonincreasing,
        report.monotonic.sg_linear_rmse_nonincreasing
    ));
    txt.push_str(&format!(
        "S_signs_at_mt = [{}, {}, {}]\n",
        report.monotonic.down_split_signs[0],
        report.monotonic.down_split_signs[1],
        report.monotonic.down_split_signs[2]
    ));
    txt.push_str(&format!("summary = {}\n", report.summary));

    fs::write(&txt_path, txt).expect("write txt");

    // CSV compact scan table
    let mut csv = String::new();
    csv.push_str("mu_gev,alpha_s_mu,m_u,m_d,m_s,m_c,m_b,m_t,L1,L2,L3,S1,S2,S3,L_free_s2,L_fixed_rms,S_linear_rmse,S_linear_r2\n");
    for row in &report.rows {
        let a_mu = alpha_s(row.mu_gev, &report.lambdas);
        csv.push_str(&format!(
            "{:.6e},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9e},{:.9e},{:.9}\n",
            row.mu_gev,
            a_mu,
            row.masses_mev[0], row.masses_mev[1], row.masses_mev[2],
            row.masses_mev[3], row.masses_mev[4], row.masses_mev[5],
            row.lg[0], row.lg[1], row.lg[2],
            row.sg[0], row.sg[1], row.sg[2],
            row.lg_z3_free.s2,
            row.lg_z3_fixed_s2_2.rms_rel,
            row.sg_linear.rmse,
            row.sg_linear.r2
        ));
    }
    fs::write(&csv_path, csv).expect("write csv");

    fs::write(
        &json_path,
        serde_json::to_string_pretty(&report).expect("serialize report"),
    )
    .expect("write json");

    println!("wrote {txt_path}");
    println!("wrote {csv_path}");
    println!("wrote {json_path}");
    println!("{}", report.summary);
}
