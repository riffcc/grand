//! RH scaling/plateau diagnostic lane:
//! - same branch-locked quadratic map protocol used in prior hardening runs
//! - explicit resolution law fit (error vs n)
//! - explicit floor (plateau) test via err(n) = A * n^{-p} + C
//! - residual-curvature diagnostics on predicted zeros vs references

use gutoe_physics::riemann_lane::{build_operator, zeta_zero_reference, RiemannOperatorParams};
use nalgebra::{DMatrix, DVector, SymmetricEigen};
use serde_json::json;
use std::cmp::Ordering;
use std::fs::{self, File};
use std::io::Write;

const TRAIN_N: usize = 40;
const HOLD_N: usize = 40;
const FREEZE_N: usize = 40;
const CORE_N: usize = TRAIN_N + HOLD_N + FREEZE_N;

#[derive(Clone, Copy, Debug)]
struct Map2 {
    a: f64,
    b: f64,
    c: f64,
}

#[derive(Clone, Debug)]
struct BranchResult {
    n: usize,
    start: usize,
    map: Map2,
    train_mape: f64,
    hold_mape: f64,
    freeze_mape: f64,
    core_hold_plus_freeze: f64,
    long_121_500_mape: Option<f64>,
    long_121_1000_mape: Option<f64>,
    curvature_quad: Option<f64>,
    curvature_r2: Option<f64>,
}

#[derive(Clone, Copy, Debug)]
struct LineFit {
    intercept: f64,
    slope: f64,
    sse: f64,
}

fn parse_ns(default: &[usize]) -> Vec<usize> {
    let raw = std::env::var("GUTOE_RIEMANN_NS").unwrap_or_else(|_| {
        default
            .iter()
            .map(|v| v.to_string())
            .collect::<Vec<_>>()
            .join(",")
    });
    let mut out = Vec::new();
    for part in raw.split(',') {
        if let Ok(v) = part.trim().parse::<usize>() {
            if v > 0 {
                out.push(v);
            }
        }
    }
    if out.is_empty() {
        default.to_vec()
    } else {
        out
    }
}

fn load_reference() -> Vec<f64> {
    if let Ok(path) = std::env::var("GUTOE_RIEMANN_REF_PATH") {
        if let Ok(raw) = fs::read_to_string(&path) {
            let vals: Vec<f64> = raw
                .lines()
                .filter_map(|line| line.trim().parse::<f64>().ok())
                .collect();
            if !vals.is_empty() {
                return vals;
            }
        }
    }
    zeta_zero_reference().to_vec()
}

fn fit_poly2(x: &[f64], y: &[f64]) -> Map2 {
    let m = x.len().min(y.len());
    let mut a = DMatrix::<f64>::zeros(m, 3);
    let mut b = DVector::<f64>::zeros(m);
    for i in 0..m {
        let xi = x[i];
        a[(i, 0)] = 1.0;
        a[(i, 1)] = xi;
        a[(i, 2)] = xi * xi;
        b[i] = y[i];
    }
    let svd = a.svd(true, true);
    let sol = svd.solve(&b, 1.0e-12).expect("svd solve");
    Map2 {
        a: sol[0],
        b: sol[1],
        c: sol[2],
    }
}

fn apply_map(map: Map2, x: &[f64]) -> Vec<f64> {
    x.iter()
        .map(|&v| map.a + map.b * v + map.c * v * v)
        .collect()
}

fn mape(pred: &[f64], target: &[f64]) -> f64 {
    let n = pred.len().min(target.len()).max(1);
    let mut s = 0.0;
    for i in 0..n {
        s += ((pred[i] - target[i]).abs()) / target[i].abs().max(1e-12);
    }
    s / n as f64
}

fn is_monotone_increasing(v: &[f64]) -> bool {
    v.windows(2).all(|w| w[1] > w[0])
}

fn long_mape(pred: &[f64], refs: &[f64], start_idx_1based: usize, end_idx_1based: usize) -> Option<f64> {
    if start_idx_1based == 0 || end_idx_1based < start_idx_1based {
        return None;
    }
    if refs.len() < end_idx_1based || pred.len() < end_idx_1based {
        return None;
    }
    let s = start_idx_1based - 1;
    let e = end_idx_1based;
    Some(mape(&pred[s..e], &refs[s..e]))
}

fn sorted_eigs(n: usize) -> Vec<f64> {
    let op = RiemannOperatorParams {
        n,
        hop1_scale: 0.5,
        hop2_scale: 0.0,
        potential_scale: 0.0,
    };
    let h = build_operator(op);
    let eig = SymmetricEigen::new(h);
    let mut w: Vec<f64> = eig.eigenvalues.iter().copied().collect();
    w.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
    w
}

fn fit_line(xs: &[f64], ys: &[f64]) -> Option<LineFit> {
    let n = xs.len().min(ys.len());
    if n < 2 {
        return None;
    }
    let x_mean = xs.iter().take(n).sum::<f64>() / n as f64;
    let y_mean = ys.iter().take(n).sum::<f64>() / n as f64;
    let mut num = 0.0;
    let mut den = 0.0;
    for i in 0..n {
        let dx = xs[i] - x_mean;
        num += dx * (ys[i] - y_mean);
        den += dx * dx;
    }
    if den <= 0.0 {
        return None;
    }
    let slope = num / den;
    let intercept = y_mean - slope * x_mean;
    let mut sse = 0.0;
    for i in 0..n {
        let r = ys[i] - (intercept + slope * xs[i]);
        sse += r * r;
    }
    Some(LineFit {
        intercept,
        slope,
        sse,
    })
}

fn quadratic_residual_curvature(pred: &[f64], refs: &[f64]) -> Option<(f64, f64)> {
    let n = pred.len().min(refs.len()).min(500);
    if n < 8 {
        return None;
    }
    let mut x = Vec::with_capacity(n);
    let mut y = Vec::with_capacity(n);
    for i in 0..n {
        let t = refs[i];
        if t <= 0.0 {
            continue;
        }
        let l = t.ln();
        let rr = pred[i] / t - 1.0;
        x.push(l);
        y.push(rr);
    }
    let m = x.len();
    if m < 8 {
        return None;
    }
    let mut a = DMatrix::<f64>::zeros(m, 3);
    let mut b = DVector::<f64>::zeros(m);
    for i in 0..m {
        let xi = x[i];
        a[(i, 0)] = 1.0;
        a[(i, 1)] = xi;
        a[(i, 2)] = xi * xi;
        b[i] = y[i];
    }
    let svd = a.svd(true, true);
    let sol = svd.solve(&b, 1.0e-12).ok()?;
    let q0 = sol[0];
    let q1 = sol[1];
    let q2 = sol[2];
    let y_mean = y.iter().sum::<f64>() / m as f64;
    let mut sse = 0.0;
    let mut sst = 0.0;
    for i in 0..m {
        let yi = y[i];
        let yhat = q0 + q1 * x[i] + q2 * x[i] * x[i];
        let r = yi - yhat;
        sse += r * r;
        let d = yi - y_mean;
        sst += d * d;
    }
    let r2 = if sst > 0.0 { 1.0 - sse / sst } else { 1.0 };
    Some((q2, r2))
}

fn evaluate_branch(n: usize, refs: &[f64]) -> BranchResult {
    let data = sorted_eigs(n);
    let max_basis = data.len().min(refs.len());
    assert!(max_basis >= CORE_N, "need >= {} points", CORE_N);

    let max_start = data.len() - CORE_N;
    let mut best: Option<(usize, Map2, f64, f64, f64, f64)> = None;
    for start in 0..=max_start {
        let x120 = &data[start..start + CORE_N];
        let map = fit_poly2(&x120[..TRAIN_N], &refs[..TRAIN_N]);
        let pred120 = apply_map(map, x120);
        if !is_monotone_increasing(&pred120) || pred120[0] <= 0.0 {
            continue;
        }
        let train = mape(&pred120[..TRAIN_N], &refs[..TRAIN_N]);
        let hold = mape(&pred120[TRAIN_N..TRAIN_N + HOLD_N], &refs[TRAIN_N..TRAIN_N + HOLD_N]);
        let freeze = mape(&pred120[TRAIN_N + HOLD_N..CORE_N], &refs[TRAIN_N + HOLD_N..CORE_N]);
        let obj = hold + freeze + 0.05 * train + 0.001 * map.c.abs();
        let keep = best.as_ref().map(|v| obj < v.5).unwrap_or(true);
        if keep {
            best = Some((start, map, train, hold, freeze, obj));
        }
    }
    let (start, map, train, hold, freeze, _obj) = best.expect("no branch candidate");
    let max_eval = max_basis.min(data.len() - start);
    let pred = apply_map(map, &data[start..start + max_eval]);

    let long_500 = long_mape(&pred, refs, 121, 500);
    let long_1000 = long_mape(&pred, refs, 121, 1000);
    let curvature = quadratic_residual_curvature(&pred, refs);
    let (curvature_quad, curvature_r2) = curvature
        .map(|(q2, r2)| (Some(q2), Some(r2)))
        .unwrap_or((None, None));

    BranchResult {
        n,
        start,
        map,
        train_mape: train,
        hold_mape: hold,
        freeze_mape: freeze,
        core_hold_plus_freeze: hold + freeze,
        long_121_500_mape: long_500,
        long_121_1000_mape: long_1000,
        curvature_quad,
        curvature_r2,
    }
}

fn fit_scaling_floor(ns: &[f64], errs: &[f64], c_floor: f64) -> Option<(LineFit, f64)> {
    let mut ys = Vec::with_capacity(errs.len());
    for &e in errs {
        if e <= c_floor {
            return None;
        }
        ys.push((e - c_floor).ln());
    }
    let fit = fit_line(ns, &ys)?;
    Some((fit, -fit.slope))
}

fn main() {
    let refs = load_reference();
    let ns = parse_ns(&[512, 768, 1024, 1536, 2048]);

    let mut rows = Vec::new();
    for &n in &ns {
        let out = evaluate_branch(n, &refs);
        rows.push(out);
    }

    rows.sort_by_key(|r| r.n);

    let log_n: Vec<f64> = rows.iter().map(|r| (r.n as f64).ln()).collect();
    let core_err: Vec<f64> = rows.iter().map(|r| r.core_hold_plus_freeze).collect();
    let pure = fit_scaling_floor(&log_n, &core_err, 0.0).expect("pure scaling fit");
    let pure_fit = pure.0;
    let pure_p = pure.1;

    let min_err = core_err
        .iter()
        .copied()
        .fold(f64::INFINITY, |a, b| a.min(b));
    let mut best_floor_c = 0.0;
    let mut best_floor_fit = pure_fit;
    let mut best_floor_p = pure_p;
    let mut best_sse = pure_fit.sse;

    let grid_n = 2000usize;
    for i in 1..=grid_n {
        let c = min_err * 0.95 * (i as f64 / grid_n as f64);
        if let Some((fit, p)) = fit_scaling_floor(&log_n, &core_err, c) {
            if fit.sse < best_sse {
                best_sse = fit.sse;
                best_floor_c = c;
                best_floor_fit = fit;
                best_floor_p = p;
            }
        }
    }

    let sse_gain_rel = (pure_fit.sse - best_sse) / pure_fit.sse.max(1e-12);
    let report = json!({
        "lane": "riemann_nail_scaling_plateau",
        "reference_len": refs.len(),
        "protocol": {
            "train_n": TRAIN_N,
            "hold_n": HOLD_N,
            "freeze_n": FREEZE_N,
            "core_n": CORE_N
        },
        "rows": rows.iter().map(|r| json!({
            "n": r.n,
            "start": r.start,
            "map2": {"a": r.map.a, "b": r.map.b, "c": r.map.c},
            "train_mape": r.train_mape,
            "hold_mape": r.hold_mape,
            "freeze_mape": r.freeze_mape,
            "core_hold_plus_freeze": r.core_hold_plus_freeze,
            "long_121_500_mape": r.long_121_500_mape,
            "long_121_1000_mape": r.long_121_1000_mape,
            "curvature_quad_coeff": r.curvature_quad,
            "curvature_r2": r.curvature_r2
        })).collect::<Vec<_>>(),
        "scaling_fit_no_floor": {
            "model": "core_err = A * n^{-p}",
            "log_intercept": pure_fit.intercept,
            "log_slope": pure_fit.slope,
            "p": pure_p,
            "sse": pure_fit.sse
        },
        "scaling_fit_with_floor": {
            "model": "core_err = A * n^{-p} + C",
            "C": best_floor_c,
            "log_intercept": best_floor_fit.intercept,
            "log_slope": best_floor_fit.slope,
            "p": best_floor_p,
            "sse": best_sse,
            "sse_gain_rel_vs_no_floor": sse_gain_rel
        }
    });

    let out_dir = "/tmp/bh_renders";
    let _ = fs::create_dir_all(out_dir);
    let json_path = format!("{out_dir}/riemann_nail_scaling_plateau_report.json");
    let txt_path = format!("{out_dir}/riemann_nail_scaling_plateau_report.txt");
    fs::write(
        &json_path,
        serde_json::to_string_pretty(&report).expect("json serialize"),
    )
    .expect("write json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "Riemann Nail Scaling/Plateau Report").ok();
    writeln!(txt, "==================================").ok();
    writeln!(txt, "reference_len: {}", refs.len()).ok();
    writeln!(
        txt,
        "protocol: train={}, hold={}, freeze={}, core={}",
        TRAIN_N, HOLD_N, FREEZE_N, CORE_N
    )
    .ok();
    writeln!(txt).ok();
    for r in &rows {
        writeln!(txt, "n={}", r.n).ok();
        writeln!(
            txt,
            "  start={} map=(a={:.9e}, b={:.9e}, c={:.9e})",
            r.start, r.map.a, r.map.b, r.map.c
        )
        .ok();
        writeln!(
            txt,
            "  mape: train={:.9e} hold={:.9e} freeze={:.9e} core={:.9e}",
            r.train_mape, r.hold_mape, r.freeze_mape, r.core_hold_plus_freeze
        )
        .ok();
        writeln!(
            txt,
            "  long: 121..500={:?} 121..1000={:?}",
            r.long_121_500_mape, r.long_121_1000_mape
        )
        .ok();
        writeln!(
            txt,
            "  curvature: quad={:?} r2={:?}",
            r.curvature_quad, r.curvature_r2
        )
        .ok();
    }
    writeln!(txt).ok();
    writeln!(txt, "Scaling fits").ok();
    writeln!(txt, "------------").ok();
    writeln!(
        txt,
        "no-floor: p={:.6} log_slope={:.6} sse={:.9e}",
        pure_p, pure_fit.slope, pure_fit.sse
    )
    .ok();
    writeln!(
        txt,
        "with-floor: C={:.9e} p={:.6} log_slope={:.6} sse={:.9e} sse_gain_rel={:.6e}",
        best_floor_c, best_floor_p, best_floor_fit.slope, best_sse, sse_gain_rel
    )
    .ok();

    println!("{}", serde_json::to_string_pretty(&report).unwrap());
    println!("wrote: {json_path}");
    println!("wrote: {txt_path}");
}

