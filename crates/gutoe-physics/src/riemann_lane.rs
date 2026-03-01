//! Shared utilities for exploratory Riemann-operator lanes.
//!
//! This module is intentionally empirical/instrumental:
//! it provides reproducible operator construction + scoring against
//! the first nontrivial zeta-zero ordinates and spacing statistics.

use nalgebra::DMatrix;
use std::cmp::Ordering;
use std::f64::consts::PI;

pub const CLIFFORD_DIM: f64 = 16.0;
pub const COMPLEMENT_DIM: f64 = 13.0; // 16 - |SU(2)|
pub const AUGMENTED_DIM: f64 = 17.0; // 16 + 1

pub const STRUCTURAL_AFFINE_SLOPE: f64 = 11.0 / 18.0;
pub const STRUCTURAL_AFFINE_SHIFT: f64 = 13.0 * 24.0 + 8.0 / AUGMENTED_DIM;

// First 80 imaginary parts of nontrivial zeta zeros (Odlyzko table zeros1).
pub const ZETA_ZERO_IMAG_FIRST_80: [f64; 80] = [
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

#[derive(Debug, Clone, Copy)]
pub struct RiemannOperatorParams {
    pub n: usize,
    pub hop1_scale: f64,
    pub hop2_scale: f64,
    pub potential_scale: f64,
}

impl Default for RiemannOperatorParams {
    fn default() -> Self {
        Self {
            n: 512,
            hop1_scale: 0.5,
            hop2_scale: 0.0,
            potential_scale: 0.0,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct RiemannMapParams {
    pub slope: f64,
    pub shift: f64,
}

impl Default for RiemannMapParams {
    fn default() -> Self {
        Self {
            slope: STRUCTURAL_AFFINE_SLOPE,
            shift: STRUCTURAL_AFFINE_SHIFT,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct RiemannFitStats {
    pub k: usize,
    pub mae: f64,
    pub rmse: f64,
    pub mape: f64,
    pub max_abs_err: f64,
    pub max_rel_err: f64,
    pub signed_rel_bias: f64,
    pub spacing_mape: f64,
    pub spacing_rmse: f64,
    pub spacing_ks: f64,
    pub spacing_var: f64,
    pub spacing_var_abs_err_to_gue: f64,
    pub objective_position: f64,
    pub objective_total: f64,
}

pub fn zeta_zero_reference() -> &'static [f64] {
    &ZETA_ZERO_IMAG_FIRST_80
}

pub fn build_operator(params: RiemannOperatorParams) -> DMatrix<f64> {
    let mut h = DMatrix::<f64>::zeros(params.n, params.n);
    let diag_shift = COMPLEMENT_DIM / CLIFFORD_DIM; // 13/16
    for i in 0..params.n {
        let x = (i + 1) as f64 + diag_shift;
        h[(i, i)] = x.ln() + params.potential_scale / x;
    }
    for i in 0..params.n.saturating_sub(1) {
        let c1 = params.hop1_scale * (((i + 1) as f64) * ((i + 2) as f64)).sqrt();
        h[(i, i + 1)] = c1;
        h[(i + 1, i)] = c1;
    }
    for i in 0..params.n.saturating_sub(2) {
        let c2 = params.hop2_scale * (((i + 1) as f64) * ((i + 3) as f64)).sqrt();
        h[(i, i + 2)] = c2;
        h[(i + 2, i)] = c2;
    }
    h
}

pub fn map_eigenvalues(raw: &[f64], map: RiemannMapParams) -> Vec<f64> {
    raw.iter().map(|&x| map.shift + map.slope * x).collect()
}

fn normalized_spacings(series: &[f64]) -> Vec<f64> {
    if series.len() < 3 {
        return Vec::new();
    }
    let mut s: Vec<f64> = series.windows(2).map(|w| w[1] - w[0]).collect();
    let mean = s.iter().sum::<f64>() / s.len() as f64;
    if mean <= 0.0 {
        return Vec::new();
    }
    for v in &mut s {
        *v /= mean;
    }
    s
}

fn ks_distance(a: &[f64], b: &[f64]) -> f64 {
    if a.is_empty() || b.is_empty() {
        return 1.0;
    }
    let mut xa = a.to_vec();
    let mut xb = b.to_vec();
    xa.sort_by(|x, y| x.partial_cmp(y).unwrap_or(Ordering::Equal));
    xb.sort_by(|x, y| x.partial_cmp(y).unwrap_or(Ordering::Equal));

    let na = xa.len() as f64;
    let nb = xb.len() as f64;
    let mut i = 0usize;
    let mut j = 0usize;
    let mut max_d = 0.0_f64;

    while i < xa.len() || j < xb.len() {
        let next_a = xa.get(i).copied().unwrap_or(f64::INFINITY);
        let next_b = xb.get(j).copied().unwrap_or(f64::INFINITY);
        let x = next_a.min(next_b);
        while i < xa.len() && xa[i] <= x {
            i += 1;
        }
        while j < xb.len() && xb[j] <= x {
            j += 1;
        }
        let fa = i as f64 / na;
        let fb = j as f64 / nb;
        max_d = max_d.max((fa - fb).abs());
    }
    max_d
}

fn variance(v: &[f64]) -> f64 {
    if v.len() < 2 {
        return 0.0;
    }
    let m = v.iter().sum::<f64>() / v.len() as f64;
    v.iter()
        .map(|x| {
            let d = *x - m;
            d * d
        })
        .sum::<f64>()
        / v.len() as f64
}

pub fn fit_against_reference(pred: &[f64], target: &[f64]) -> RiemannFitStats {
    let k = pred.len().min(target.len());
    let mut abs_err_sum = 0.0_f64;
    let mut rel_err_sum = 0.0_f64;
    let mut sq_err_sum = 0.0_f64;
    let mut max_abs_err = 0.0_f64;
    let mut max_rel_err = 0.0_f64;
    let mut signed_rel_sum = 0.0_f64;

    for i in 0..k {
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
    }

    let kf = k as f64;
    let mae = abs_err_sum / kf;
    let rmse = (sq_err_sum / kf).sqrt();
    let mape = rel_err_sum / kf;
    let signed_rel_bias = signed_rel_sum / kf;

    let pred_spacing = normalized_spacings(&pred[..k]);
    let tgt_spacing = normalized_spacings(&target[..k]);
    let s_k = pred_spacing.len().min(tgt_spacing.len()).max(1) as f64;
    let mut spacing_sq_sum = 0.0_f64;
    for i in 0..(s_k as usize).min(pred_spacing.len()).min(tgt_spacing.len()) {
        let de = pred_spacing[i] - tgt_spacing[i];
        spacing_sq_sum += de * de;
    }
    let spacing_mape = if !tgt_spacing.is_empty() {
        pred_spacing
            .iter()
            .zip(tgt_spacing.iter())
            .map(|(p, t)| ((*p - *t).abs()) / t.abs().max(1e-12))
            .sum::<f64>()
            / s_k
    } else {
        1.0
    };
    let spacing_rmse = (spacing_sq_sum / s_k).sqrt();
    let spacing_ks = ks_distance(&pred_spacing, &tgt_spacing);
    let spacing_var = variance(&pred_spacing);
    let gue_var_target = 3.0 * PI / 8.0 - 1.0;
    let spacing_var_abs_err_to_gue = (spacing_var - gue_var_target).abs();

    let objective_position = mape;
    let objective_total = mape
        + 0.50 * spacing_ks
        + 0.25 * spacing_mape
        + 0.25 * spacing_var_abs_err_to_gue;

    RiemannFitStats {
        k,
        mae,
        rmse,
        mape,
        max_abs_err,
        max_rel_err,
        signed_rel_bias,
        spacing_mape,
        spacing_rmse,
        spacing_ks,
        spacing_var,
        spacing_var_abs_err_to_gue,
        objective_position,
        objective_total,
    }
}
