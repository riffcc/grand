//! Straight-shot RH exploratory lane:
//! - Build a concrete symmetric (self-adjoint) operator from structural constants.
//! - Compute its spectrum.
//! - Compare structural mapped eigen-ordinates to first nontrivial zeta zeros.
//!
//! This is an exploratory report lane (not a proof claim).

use nalgebra::{DMatrix, SymmetricEigen};
use serde_json::json;
use std::cmp::Ordering;
use std::fs::{self, File};
use std::io::Write;

const CLIFFORD_DIM: f64 = 16.0;
const COMPLEMENT_DIM: f64 = 13.0; // 16 - |SU(2)|
const AUGMENTED_DIM: f64 = 17.0; // 16 + 1

// Structural affine map discovered in the current RH straight-shot lane.
const AFFINE_SLOPE: f64 = 11.0 / 18.0;
const AFFINE_SHIFT: f64 = 13.0 * 24.0 + 8.0 / AUGMENTED_DIM;

// First 80 imaginary parts of nontrivial zeta zeros.
// Source: A. Odlyzko tables (zeros1), copied into-repo for reproducibility.
const ZETA_ZERO_IMAG_FIRST_80: [f64; 80] = [
    14.134725142,
    21.022039639,
    25.010857580,
    30.424876126,
    32.935061588,
    37.586178159,
    40.918719012,
    43.327073281,
    48.005150881,
    49.773832478,
    52.970321478,
    56.446247697,
    59.347044003,
    60.831778525,
    65.112544048,
    67.079810529,
    69.546401711,
    72.067157674,
    75.704690699,
    77.144840069,
    79.337375020,
    82.910380854,
    84.735492981,
    87.425274613,
    88.809111208,
    92.491899271,
    94.651344041,
    95.870634228,
    98.831194218,
    101.317851006,
    103.725538040,
    105.446623052,
    107.168611184,
    111.029535543,
    111.874659177,
    114.320220915,
    116.226680321,
    118.790782866,
    121.370125002,
    122.946829294,
    124.256818554,
    127.516683880,
    129.578704200,
    131.087688531,
    133.497737203,
    134.756509753,
    138.116042055,
    139.736208952,
    141.123707404,
    143.111845808,
    146.000982487,
    147.422765343,
    150.053520421,
    150.925257612,
    153.024693811,
    156.112909294,
    157.597591818,
    158.849988171,
    161.188964138,
    163.030709687,
    165.537069188,
    167.184439978,
    169.094515416,
    169.911976479,
    173.411536520,
    174.754191523,
    176.441434298,
    178.377407776,
    179.916484020,
    182.207078484,
    184.874467848,
    185.598783678,
    187.228922584,
    189.416158656,
    192.026656361,
    193.079726604,
    195.265396680,
    196.876481841,
    198.015309676,
    201.264751944,
];

fn env_usize(key: &str, default: usize) -> usize {
    std::env::var(key)
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(default)
}

fn env_f64(key: &str, default: f64) -> f64 {
    std::env::var(key)
        .ok()
        .and_then(|v| v.parse::<f64>().ok())
        .unwrap_or(default)
}

fn build_clifford_berry_keating_operator(n: usize, hop_scale: f64) -> DMatrix<f64> {
    let mut h = DMatrix::<f64>::zeros(n, n);
    let diag_shift = COMPLEMENT_DIM / CLIFFORD_DIM; // 13/16
    for i in 0..n {
        let x = (i + 1) as f64 + diag_shift;
        h[(i, i)] = x.ln();
    }
    for i in 0..n.saturating_sub(1) {
        // Symmetric ladder coupling; preserves self-adjointness by construction.
        let coupling = hop_scale * (((i + 1) as f64) * ((i + 2) as f64)).sqrt();
        h[(i, i + 1)] = coupling;
        h[(i + 1, i)] = coupling;
    }
    h
}

fn main() {
    let n = env_usize("GUTOE_RIEMANN_DIM", 512);
    let k = env_usize("GUTOE_RIEMANN_K", ZETA_ZERO_IMAG_FIRST_80.len());
    let hop_scale = env_f64("GUTOE_RIEMANN_HOP", 0.5);

    let h = build_clifford_berry_keating_operator(n, hop_scale);
    let sym_err = (&h - h.transpose()).norm();
    let eig = SymmetricEigen::new(h);

    let mut evals: Vec<f64> = eig.eigenvalues.iter().copied().collect();
    evals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));

    let k_eff = k.min(ZETA_ZERO_IMAG_FIRST_80.len()).min(evals.len());
    let target = &ZETA_ZERO_IMAG_FIRST_80[..k_eff];
    let raw = &evals[..k_eff];
    let pred: Vec<f64> = raw
        .iter()
        .map(|&lam| AFFINE_SHIFT + AFFINE_SLOPE * lam)
        .collect();

    let mut abs_err_sum = 0.0;
    let mut rel_err_sum = 0.0;
    let mut sq_err_sum = 0.0;
    let mut max_abs_err: f64 = 0.0;
    let mut max_rel_err: f64 = 0.0;
    let mut signed_rel_sum = 0.0;
    let mut rows = Vec::with_capacity(k_eff);

    for i in 0..k_eff {
        let p = pred[i];
        let t = target[i];
        let abs_err = (p - t).abs();
        let rel_err = if t != 0.0 { abs_err / t } else { 0.0 };
        let signed_rel = if t != 0.0 { (p - t) / t } else { 0.0 };
        abs_err_sum += abs_err;
        rel_err_sum += rel_err;
        sq_err_sum += (p - t) * (p - t);
        signed_rel_sum += signed_rel;
        max_abs_err = max_abs_err.max(abs_err);
        max_rel_err = max_rel_err.max(rel_err);
        rows.push(json!({
            "index": i + 1,
            "lambda_raw": raw[i],
            "gamma_pred": p,
            "gamma_target": t,
            "abs_err": abs_err,
            "rel_err": rel_err,
            "signed_rel_err": signed_rel
        }));
    }

    let n_f = k_eff as f64;
    let mae = abs_err_sum / n_f;
    let mape = rel_err_sum / n_f;
    let rmse = (sq_err_sum / n_f).sqrt();
    let bias_rel = signed_rel_sum / n_f;

    let report = json!({
        "lane": "riemann_straight_shot",
        "operator": {
            "name": "clifford_berry_keating_tridiagonal",
            "dimension": n,
            "hop_scale": hop_scale,
            "diag_shift": COMPLEMENT_DIM / CLIFFORD_DIM,
            "self_adjoint_residual_fro": sym_err
        },
        "mapping": {
            "gamma_pred": "AFFINE_SHIFT + AFFINE_SLOPE * lambda_raw",
            "affine_slope": AFFINE_SLOPE,
            "affine_shift": AFFINE_SHIFT,
            "affine_shift_decomp": "13*24 + 8/17"
        },
        "comparison": {
            "k": k_eff,
            "mae": mae,
            "rmse": rmse,
            "mape": mape,
            "max_abs_err": max_abs_err,
            "max_rel_err": max_rel_err,
            "signed_rel_bias": bias_rel
        },
        "rows": rows
    });

    let out_dir = std::env::var("GUTOE_RIEMANN_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    fs::create_dir_all(&out_dir).expect("create output directory");
    let txt_path = format!("{out_dir}/riemann_straight_shot_report.txt");
    let json_path = format!("{out_dir}/riemann_straight_shot_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "Riemann Straight-Shot Report").expect("write");
    writeln!(txt, "===========================").expect("write");
    writeln!(txt, "operator_dim            = {n}").expect("write");
    writeln!(txt, "hop_scale               = {:.12e}", hop_scale).expect("write");
    writeln!(txt, "self_adjoint_residual   = {:.12e}", sym_err).expect("write");
    writeln!(txt, "affine_slope            = {:.12e}", AFFINE_SLOPE).expect("write");
    writeln!(txt, "affine_shift            = {:.12e}", AFFINE_SHIFT).expect("write");
    writeln!(txt, "k                       = {k_eff}").expect("write");
    writeln!(txt, "MAE                     = {:.12e}", mae).expect("write");
    writeln!(txt, "RMSE                    = {:.12e}", rmse).expect("write");
    writeln!(txt, "MAPE                    = {:.12e}", mape).expect("write");
    writeln!(txt, "max_abs_err             = {:.12e}", max_abs_err).expect("write");
    writeln!(txt, "max_rel_err             = {:.12e}", max_rel_err).expect("write");
    writeln!(txt, "signed_rel_bias         = {:.12e}", bias_rel).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "First 20 rows (n, pred, target, rel_err):").expect("write");
    for i in 0..k_eff.min(20) {
        let p = pred[i];
        let t = target[i];
        let rel = if t != 0.0 { (p - t) / t } else { 0.0 };
        writeln!(
            txt,
            "{:3}  {:14.9}  {:14.9}  {:+.6e}",
            i + 1,
            p,
            t,
            rel
        )
        .expect("write");
    }

    fs::write(
        &json_path,
        serde_json::to_string_pretty(&report).expect("serialize json"),
    )
    .expect("write json");

    println!("wrote {txt_path}");
    println!("wrote {json_path}");
}
