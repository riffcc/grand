//! RH ablation/fairness lane:
//! compare equal-complexity quadratic-map protocol across
//! - true structural spectrum
//! - scrambled spectrum control
//! - linear surrogate control

use gutoe_physics::riemann_lane::{build_operator, zeta_zero_reference, RiemannOperatorParams};
use nalgebra::{DMatrix, DVector, SymmetricEigen};
use serde_json::json;
use std::cmp::Ordering;
use std::fs;

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
struct EvalOut {
    start: usize,
    map: Map2,
    train_mape: f64,
    hold_mape: f64,
    freeze_mape: f64,
    long_500_mape: Option<f64>,
    long_1000_mape: Option<f64>,
    objective: f64,
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

fn evaluate_dataset(data: &[f64], refs: &[f64]) -> EvalOut {
    let max_basis = data.len().min(refs.len());
    assert!(
        max_basis >= CORE_N,
        "need at least {} points in both data and refs",
        CORE_N
    );

    let max_start = data.len() - CORE_N;
    let mut best: Option<EvalOut> = None;

    for start in 0..=max_start {
        let x120 = &data[start..start + CORE_N];
        let map = fit_poly2(&x120[..TRAIN_N], &refs[..TRAIN_N]);
        let pred120 = apply_map(map, x120);
        let train = mape(&pred120[..TRAIN_N], &refs[..TRAIN_N]);
        let hold = mape(&pred120[TRAIN_N..TRAIN_N + HOLD_N], &refs[TRAIN_N..TRAIN_N + HOLD_N]);
        let freeze = mape(&pred120[TRAIN_N + HOLD_N..CORE_N], &refs[TRAIN_N + HOLD_N..CORE_N]);
        let mono_penalty = if is_monotone_increasing(&pred120) && pred120[0] > 0.0 {
            0.0
        } else {
            1.0
        };
        let obj = hold + freeze + 0.05 * train + 0.001 * map.c.abs() + 0.1 * mono_penalty;

        let max_eval = max_basis.min(data.len() - start);
        let pred_full = apply_map(map, &data[start..start + max_eval]);
        let long_500 = long_mape(&pred_full, refs, 121, 500);
        let long_1000 = long_mape(&pred_full, refs, 121, 1000);

        let cand = EvalOut {
            start,
            map,
            train_mape: train,
            hold_mape: hold,
            freeze_mape: freeze,
            long_500_mape: long_500,
            long_1000_mape: long_1000,
            objective: obj,
        };
        if best.as_ref().map(|b| obj < b.objective).unwrap_or(true) {
            best = Some(cand);
        }
    }

    best.expect("at least one valid candidate")
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

fn scrambled_copy(v: &[f64]) -> Vec<f64> {
    let mut idx: Vec<usize> = (0..v.len()).collect();
    idx.sort_by_key(|&i| {
        // Deterministic integer hash for stable, non-learning scramble.
        let x = i as u64;
        x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407)
    });
    idx.into_iter().map(|i| v[i]).collect()
}

fn linear_surrogate(n: usize) -> Vec<f64> {
    (0..n).map(|i| (i + 1) as f64).collect()
}

fn main() {
    let refs = load_reference();
    let ns = parse_ns(&[512, 1024, 2048]);

    let mut rows = Vec::new();
    for &n in &ns {
        let true_spec = sorted_eigs(n);
        let scramble_spec = scrambled_copy(&true_spec);
        let linear_spec = linear_surrogate(true_spec.len());

        let truth = evaluate_dataset(&true_spec, &refs);
        let scramble = evaluate_dataset(&scramble_spec, &refs);
        let linear = evaluate_dataset(&linear_spec, &refs);

        rows.push(json!({
            "n": n,
            "truth": {
                "start": truth.start,
                "map2": {"a": truth.map.a, "b": truth.map.b, "c": truth.map.c},
                "train_mape": truth.train_mape,
                "hold_mape": truth.hold_mape,
                "freeze_mape": truth.freeze_mape,
                "core_hold_plus_freeze": truth.hold_mape + truth.freeze_mape,
                "long_121_500_mape": truth.long_500_mape,
                "long_121_1000_mape": truth.long_1000_mape,
                "objective": truth.objective,
            },
            "scrambled_control": {
                "start": scramble.start,
                "map2": {"a": scramble.map.a, "b": scramble.map.b, "c": scramble.map.c},
                "train_mape": scramble.train_mape,
                "hold_mape": scramble.hold_mape,
                "freeze_mape": scramble.freeze_mape,
                "core_hold_plus_freeze": scramble.hold_mape + scramble.freeze_mape,
                "long_121_500_mape": scramble.long_500_mape,
                "long_121_1000_mape": scramble.long_1000_mape,
                "objective": scramble.objective,
            },
            "linear_control": {
                "start": linear.start,
                "map2": {"a": linear.map.a, "b": linear.map.b, "c": linear.map.c},
                "train_mape": linear.train_mape,
                "hold_mape": linear.hold_mape,
                "freeze_mape": linear.freeze_mape,
                "core_hold_plus_freeze": linear.hold_mape + linear.freeze_mape,
                "long_121_500_mape": linear.long_500_mape,
                "long_121_1000_mape": linear.long_1000_mape,
                "objective": linear.objective,
            },
            "truth_vs_scrambled_core_gain_rel": (scramble.hold_mape + scramble.freeze_mape - (truth.hold_mape + truth.freeze_mape))
                / (scramble.hold_mape + scramble.freeze_mape).max(1e-12),
            "truth_vs_linear_core_gain_rel": (linear.hold_mape + linear.freeze_mape - (truth.hold_mape + truth.freeze_mape))
                / (linear.hold_mape + linear.freeze_mape).max(1e-12),
        }));
    }

    let report = json!({
        "lane": "riemann_nail_ablation",
        "reference_len": refs.len(),
        "protocol": {
            "train": "1..40",
            "hold": "41..80",
            "freeze": "81..120",
            "controls": ["scrambled_eigenvalue_sequence", "linear_index_surrogate"],
            "model": "quadratic_map_same_objective_same_branch_scan"
        },
        "rows": rows,
    });

    let out_dir = std::env::var("GUTOE_RIEMANN_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    fs::create_dir_all(&out_dir).expect("create output dir");
    let txt_path = format!("{out_dir}/riemann_nail_ablation_report.txt");
    let json_path = format!("{out_dir}/riemann_nail_ablation_report.json");

    fs::write(&json_path, serde_json::to_string_pretty(&report).expect("serialize")).expect("write json");
    let mut txt = String::new();
    txt.push_str("Riemann Nail Ablation Report\n");
    txt.push_str("===========================\n");
    txt.push_str(&format!("reference_len = {}\n\n", refs.len()));
    for row in report["rows"].as_array().expect("rows") {
        txt.push_str(&format!("n = {}\n", row["n"]));
        txt.push_str(&format!(
            "  truth core hold+freeze     = {:.12e}\n",
            row["truth"]["core_hold_plus_freeze"].as_f64().unwrap_or(f64::NAN)
        ));
        txt.push_str(&format!(
            "  scrambled core hold+freeze = {:.12e}\n",
            row["scrambled_control"]["core_hold_plus_freeze"]
                .as_f64()
                .unwrap_or(f64::NAN)
        ));
        txt.push_str(&format!(
            "  linear core hold+freeze    = {:.12e}\n",
            row["linear_control"]["core_hold_plus_freeze"]
                .as_f64()
                .unwrap_or(f64::NAN)
        ));
        txt.push_str(&format!(
            "  truth vs scrambled gain    = {:.12e}\n",
            row["truth_vs_scrambled_core_gain_rel"]
                .as_f64()
                .unwrap_or(f64::NAN)
        ));
        txt.push_str(&format!(
            "  truth vs linear gain       = {:.12e}\n",
            row["truth_vs_linear_core_gain_rel"]
                .as_f64()
                .unwrap_or(f64::NAN)
        ));
        txt.push('\n');
    }
    fs::write(&txt_path, txt).expect("write txt");

    println!("wrote {txt_path}");
    println!("wrote {json_path}");
}
