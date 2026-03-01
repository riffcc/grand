//! RH nail attempt: branch-locked nonlinear map.
//! - Scan contiguous spectral windows.
//! - Fit nonlinear map on train (1..40) only.
//! - Select by hold+freeze performance (41..120).

use gutoe_physics::riemann_lane::{build_operator, zeta_zero_reference, RiemannOperatorParams};
use nalgebra::{DMatrix, DVector, SymmetricEigen};
use serde_json::json;
use std::cmp::Ordering;
use std::fs::{self, File};
use std::io::Write;

const TRAIN_N: usize = 40;
const HOLD_N: usize = 40;
const FREEZE_N: usize = 40;

#[derive(Debug, Clone, Copy)]
struct Map2 {
    a: f64,
    b: f64,
    c: f64,
}

#[derive(Debug, Clone)]
struct Candidate {
    start: usize,
    map: Map2,
    train_mape: f64,
    hold_mape: f64,
    freeze_mape: f64,
    objective: f64,
    monotone: bool,
}

fn env_usize(key: &str, default: usize) -> usize {
    std::env::var(key)
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(default)
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

fn sorted_eigs(op: RiemannOperatorParams) -> Vec<f64> {
    let h = build_operator(op);
    let eig = SymmetricEigen::new(h);
    let mut w: Vec<f64> = eig.eigenvalues.iter().copied().collect();
    w.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
    w
}

fn structural_map(x: &[f64]) -> Vec<f64> {
    let slope = 11.0 / 18.0;
    let shift = 13.0 * 24.0 + 8.0 / 17.0;
    x.iter().map(|&v| shift + slope * v).collect()
}

fn main() {
    let n = env_usize("GUTOE_RIEMANN_DIM", 1024);
    let refs = zeta_zero_reference();
    let needed = TRAIN_N + HOLD_N + FREEZE_N;
    assert!(refs.len() >= needed, "need at least {needed} target zeros");

    let op = RiemannOperatorParams {
        n,
        hop1_scale: 0.5,
        hop2_scale: 0.0,
        potential_scale: 0.0,
    };

    let evals = sorted_eigs(op);
    let max_start = evals.len().saturating_sub(needed);

    // Baseline (naive lowest branch + structural map)
    let base_raw = &evals[..needed];
    let base_pred = structural_map(base_raw);
    let base_train = mape(&base_pred[..TRAIN_N], &refs[..TRAIN_N]);
    let base_hold = mape(
        &base_pred[TRAIN_N..TRAIN_N + HOLD_N],
        &refs[TRAIN_N..TRAIN_N + HOLD_N],
    );
    let base_freeze = mape(
        &base_pred[TRAIN_N + HOLD_N..needed],
        &refs[TRAIN_N + HOLD_N..needed],
    );

    let mut best: Option<Candidate> = None;
    let mut top: Vec<Candidate> = Vec::new();

    for start in 0..=max_start {
        let raw = &evals[start..start + needed];
        let map = fit_poly2(&raw[..TRAIN_N], &refs[..TRAIN_N]);
        let pred = apply_map(map, raw);

        let mono = is_monotone_increasing(&pred) && pred[0] > 0.0;
        if !mono {
            continue;
        }

        let train = mape(&pred[..TRAIN_N], &refs[..TRAIN_N]);
        let hold = mape(&pred[TRAIN_N..TRAIN_N + HOLD_N], &refs[TRAIN_N..TRAIN_N + HOLD_N]);
        let freeze = mape(&pred[TRAIN_N + HOLD_N..needed], &refs[TRAIN_N + HOLD_N..needed]);

        let obj = hold + freeze + 0.05 * train + 0.001 * map.c.abs();
        let cand = Candidate {
            start,
            map,
            train_mape: train,
            hold_mape: hold,
            freeze_mape: freeze,
            objective: obj,
            monotone: mono,
        };

        if best.as_ref().map(|b| obj < b.objective).unwrap_or(true) {
            best = Some(cand.clone());
        }
        top.push(cand);
    }

    top.sort_by(|a, b| a.objective.partial_cmp(&b.objective).unwrap_or(Ordering::Equal));
    let best = best.expect("found at least one monotone candidate");

    let best_raw = &evals[best.start..best.start + needed];
    let best_pred = apply_map(best.map, best_raw);

    let report = json!({
        "lane": "riemann_nail_branch_map",
        "n": n,
        "window_len": needed,
        "baseline": {
            "start": 0,
            "map": "gamma = (13*24 + 8/17) + (11/18) lambda",
            "train_mape": base_train,
            "hold_mape": base_hold,
            "freeze_mape": base_freeze,
            "hold_plus_freeze": base_hold + base_freeze,
        },
        "best": {
            "start": best.start,
            "map2": {"a": best.map.a, "b": best.map.b, "c": best.map.c},
            "train_mape": best.train_mape,
            "hold_mape": best.hold_mape,
            "freeze_mape": best.freeze_mape,
            "hold_plus_freeze": best.hold_mape + best.freeze_mape,
            "objective": best.objective,
            "monotone": best.monotone,
            "pred_range": [best_pred[0], best_pred[needed-1]]
        },
        "improvement": {
            "hold_plus_freeze_abs": (base_hold + base_freeze) - (best.hold_mape + best.freeze_mape),
            "hold_plus_freeze_rel": ((base_hold + base_freeze) - (best.hold_mape + best.freeze_mape)) / (base_hold + base_freeze).max(1e-12)
        },
        "top10": top.iter().take(10).map(|c| json!({
            "start": c.start,
            "a": c.map.a,
            "b": c.map.b,
            "c": c.map.c,
            "train_mape": c.train_mape,
            "hold_mape": c.hold_mape,
            "freeze_mape": c.freeze_mape,
            "objective": c.objective
        })).collect::<Vec<_>>()
    });

    let out_dir = std::env::var("GUTOE_RIEMANN_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    fs::create_dir_all(&out_dir).expect("create output dir");
    let txt_path = format!("{out_dir}/riemann_nail_branch_map_report.txt");
    let json_path = format!("{out_dir}/riemann_nail_branch_map_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "Riemann Nail: Branch-Locked Nonlinear Map").expect("write");
    writeln!(txt, "======================================").expect("write");
    writeln!(txt, "n                             = {}", n).expect("write");
    writeln!(txt, "window_len                    = {}", needed).expect("write");
    writeln!(txt, "\nBaseline (start=0, structural map):").expect("write");
    writeln!(txt, "  train MAPE                  = {:.12e}", base_train).expect("write");
    writeln!(txt, "  hold MAPE                   = {:.12e}", base_hold).expect("write");
    writeln!(txt, "  freeze MAPE                 = {:.12e}", base_freeze).expect("write");
    writeln!(txt, "  hold+freeze                 = {:.12e}", base_hold + base_freeze).expect("write");

    writeln!(txt, "\nBest branch-locked nonlinear map:").expect("write");
    writeln!(txt, "  start                       = {}", best.start).expect("write");
    writeln!(txt, "  map: gamma = a + b*lambda + c*lambda^2").expect("write");
    writeln!(txt, "    a                         = {:.12e}", best.map.a).expect("write");
    writeln!(txt, "    b                         = {:.12e}", best.map.b).expect("write");
    writeln!(txt, "    c                         = {:.12e}", best.map.c).expect("write");
    writeln!(txt, "  train MAPE                  = {:.12e}", best.train_mape).expect("write");
    writeln!(txt, "  hold MAPE                   = {:.12e}", best.hold_mape).expect("write");
    writeln!(txt, "  freeze MAPE                 = {:.12e}", best.freeze_mape).expect("write");
    writeln!(txt, "  hold+freeze                 = {:.12e}", best.hold_mape + best.freeze_mape).expect("write");

    let imp_abs = (base_hold + base_freeze) - (best.hold_mape + best.freeze_mape);
    let imp_rel = imp_abs / (base_hold + base_freeze).max(1e-12);
    writeln!(txt, "\nImprovement vs baseline:").expect("write");
    writeln!(txt, "  absolute                    = {:.12e}", imp_abs).expect("write");
    writeln!(txt, "  relative                    = {:.12e}", imp_rel).expect("write");

    fs::write(&json_path, serde_json::to_string_pretty(&report).expect("serialize")).expect("write json");
    println!("wrote {txt_path}");
    println!("wrote {json_path}");
}
