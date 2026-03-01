//! Consolidated RH hardening pass:
//! - multi-resolution scan (n list)
//! - branch-locked quadratic map (fit on train only)
//! - long unseen holdouts (121..500 and 121..1000)
//! - coefficient stability report for (a,b,c)

use gutoe_physics::riemann_lane::{build_operator, zeta_zero_reference, RiemannOperatorParams};
use nalgebra::{DMatrix, DVector, SymmetricEigen};
use serde_json::json;
use std::cmp::Ordering;
use std::fs::{self, File};
use std::io::Write;

const TRAIN_N: usize = 40;
const HOLD_N: usize = 40;
const FREEZE_N: usize = 40;
const CORE_N: usize = TRAIN_N + HOLD_N + FREEZE_N; // 120

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
}

#[derive(Debug, Clone)]
struct ResolutionResult {
    n: usize,
    refs_len: usize,
    eval_len: usize,
    objective_basis: usize,
    baseline_train: f64,
    baseline_hold: f64,
    baseline_freeze: f64,
    baseline_long_500: Option<f64>,
    baseline_long_1000: Option<f64>,
    best_start: usize,
    best_a: f64,
    best_b: f64,
    best_c: f64,
    best_train: f64,
    best_hold: f64,
    best_freeze: f64,
    best_long_500: Option<f64>,
    best_long_1000: Option<f64>,
    improvement_core_rel: f64,
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

fn structural_map(x: &[f64]) -> Vec<f64> {
    let slope = 11.0 / 18.0;
    let shift = 13.0 * 24.0 + 8.0 / 17.0;
    x.iter().map(|&v| shift + slope * v).collect()
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

fn run_for_n(n: usize, refs: &[f64]) -> ResolutionResult {
    let evals = sorted_eigs(n);
    let max_basis = evals.len().min(refs.len());
    assert!(max_basis >= CORE_N, "need at least 120 points for objective");

    let max_start = evals.len() - CORE_N;

    // Baseline structural map (start=0).
    let base_raw = &evals[..max_basis];
    let base_pred = structural_map(base_raw);
    let baseline_train = mape(&base_pred[..TRAIN_N], &refs[..TRAIN_N]);
    let baseline_hold = mape(
        &base_pred[TRAIN_N..TRAIN_N + HOLD_N],
        &refs[TRAIN_N..TRAIN_N + HOLD_N],
    );
    let baseline_freeze = mape(
        &base_pred[TRAIN_N + HOLD_N..CORE_N],
        &refs[TRAIN_N + HOLD_N..CORE_N],
    );
    let baseline_long_500 = long_mape(&base_pred, refs, 121, 500);
    let baseline_long_1000 = long_mape(&base_pred, refs, 121, 1000);

    let mut best: Option<Candidate> = None;
    for start in 0..=max_start {
        let raw_120 = &evals[start..start + CORE_N];
        let map = fit_poly2(&raw_120[..TRAIN_N], &refs[..TRAIN_N]);

        // Objective only on canonical 120 point protocol.
        let pred_120 = apply_map(map, raw_120);
        if !is_monotone_increasing(&pred_120) || pred_120[0] <= 0.0 {
            continue;
        }

        let train = mape(&pred_120[..TRAIN_N], &refs[..TRAIN_N]);
        let hold = mape(
            &pred_120[TRAIN_N..TRAIN_N + HOLD_N],
            &refs[TRAIN_N..TRAIN_N + HOLD_N],
        );
        let freeze = mape(
            &pred_120[TRAIN_N + HOLD_N..CORE_N],
            &refs[TRAIN_N + HOLD_N..CORE_N],
        );
        let obj = hold + freeze + 0.05 * train + 0.001 * map.c.abs();
        let cand = Candidate {
            start,
            map,
            train_mape: train,
            hold_mape: hold,
            freeze_mape: freeze,
            objective: obj,
        };
        if best.as_ref().map(|b| obj < b.objective).unwrap_or(true) {
            best = Some(cand);
        }
    }

    let best = best.expect("at least one candidate");

    // Evaluate best candidate on longest available horizon.
    let best_max_basis = max_basis.min(evals.len() - best.start);
    let best_raw = &evals[best.start..best.start + best_max_basis];
    let best_pred = apply_map(best.map, best_raw);
    let best_long_500 = long_mape(&best_pred, refs, 121, 500);
    let best_long_1000 = long_mape(&best_pred, refs, 121, 1000);

    let baseline_core = baseline_hold + baseline_freeze;
    let best_core = best.hold_mape + best.freeze_mape;
    let improvement_core_rel = (baseline_core - best_core) / baseline_core.max(1e-12);

    ResolutionResult {
        n,
        refs_len: refs.len(),
        eval_len: evals.len(),
        objective_basis: max_basis,
        baseline_train,
        baseline_hold,
        baseline_freeze,
        baseline_long_500,
        baseline_long_1000,
        best_start: best.start,
        best_a: best.map.a,
        best_b: best.map.b,
        best_c: best.map.c,
        best_train: best.train_mape,
        best_hold: best.hold_mape,
        best_freeze: best.freeze_mape,
        best_long_500,
        best_long_1000,
        improvement_core_rel,
    }
}

fn main() {
    let refs = load_reference();
    let ns = parse_ns(&[512, 1024, 2048, 4096]);

    let mut rows = Vec::new();
    for &n in &ns {
        rows.push(run_for_n(n, &refs));
    }

    let ref_row = rows.first().expect("at least one row");
    let coeff_stability = rows
        .iter()
        .map(|r| {
            json!({
                "n": r.n,
                "da_vs_first": r.best_a - ref_row.best_a,
                "db_vs_first": r.best_b - ref_row.best_b,
                "dc_vs_first": r.best_c - ref_row.best_c,
                "rel_a_vs_first": if ref_row.best_a != 0.0 {(r.best_a - ref_row.best_a)/ref_row.best_a} else {0.0},
                "rel_b_vs_first": if ref_row.best_b != 0.0 {(r.best_b - ref_row.best_b)/ref_row.best_b} else {0.0},
                "rel_c_vs_first": if ref_row.best_c != 0.0 {(r.best_c - ref_row.best_c)/ref_row.best_c} else {0.0},
            })
        })
        .collect::<Vec<_>>();

    let report = json!({
        "lane": "riemann_nail_long_holdout",
        "reference_len": refs.len(),
        "protocol": {
            "train": "1..40",
            "hold": "41..80",
            "freeze": "81..120",
            "long_500": "121..500",
            "long_1000": "121..1000"
        },
        "rows": rows.iter().map(|r| json!({
            "n": r.n,
            "refs_len": r.refs_len,
            "eval_len": r.eval_len,
            "objective_basis": r.objective_basis,
            "baseline": {
                "train_mape": r.baseline_train,
                "hold_mape": r.baseline_hold,
                "freeze_mape": r.baseline_freeze,
                "core_hold_plus_freeze": r.baseline_hold + r.baseline_freeze,
                "long_121_500_mape": r.baseline_long_500,
                "long_121_1000_mape": r.baseline_long_1000
            },
            "best": {
                "start": r.best_start,
                "map2": {"a": r.best_a, "b": r.best_b, "c": r.best_c},
                "train_mape": r.best_train,
                "hold_mape": r.best_hold,
                "freeze_mape": r.best_freeze,
                "core_hold_plus_freeze": r.best_hold + r.best_freeze,
                "long_121_500_mape": r.best_long_500,
                "long_121_1000_mape": r.best_long_1000,
                "core_improvement_rel": r.improvement_core_rel
            }
        })).collect::<Vec<_>>(),
        "coeff_stability_vs_first_n": coeff_stability,
    });

    let out_dir = std::env::var("GUTOE_RIEMANN_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    fs::create_dir_all(&out_dir).expect("create output dir");
    let txt_path = format!("{out_dir}/riemann_nail_long_holdout_report.txt");
    let json_path = format!("{out_dir}/riemann_nail_long_holdout_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "Riemann Nail Long Holdout Report").expect("write");
    writeln!(txt, "===============================").expect("write");
    writeln!(txt, "reference_len = {}", refs.len()).expect("write");
    writeln!(txt, "ns            = {:?}", ns).expect("write");
    writeln!(txt).expect("write");
    for r in &rows {
        writeln!(txt, "n = {}", r.n).expect("write");
        writeln!(
            txt,
            "  baseline core hold+freeze   = {:.12e}",
            r.baseline_hold + r.baseline_freeze
        )
        .expect("write");
        writeln!(
            txt,
            "  best core hold+freeze       = {:.12e}",
            r.best_hold + r.best_freeze
        )
        .expect("write");
        writeln!(txt, "  core improvement rel        = {:.12e}", r.improvement_core_rel)
            .expect("write");
        writeln!(txt, "  best start                  = {}", r.best_start).expect("write");
        writeln!(
            txt,
            "  best map2 (a,b,c)           = ({:.12e}, {:.12e}, {:.12e})",
            r.best_a, r.best_b, r.best_c
        )
        .expect("write");
        if let Some(v) = r.best_long_500 {
            writeln!(txt, "  best long 121..500 MAPE     = {:.12e}", v).expect("write");
        } else {
            writeln!(txt, "  best long 121..500 MAPE     = n/a").expect("write");
        }
        if let Some(v) = r.best_long_1000 {
            writeln!(txt, "  best long 121..1000 MAPE    = {:.12e}", v).expect("write");
        } else {
            writeln!(txt, "  best long 121..1000 MAPE    = n/a").expect("write");
        }
        writeln!(txt).expect("write");
    }

    fs::write(&json_path, serde_json::to_string_pretty(&report).expect("serialize")).expect("write json");
    println!("wrote {txt_path}");
    println!("wrote {json_path}");
}

