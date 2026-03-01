//! Straight-shot RH exploratory lane:
//! - Build a concrete symmetric (self-adjoint) operator from structural constants.
//! - Compute its spectrum.
//! - Compare mapped eigen-ordinates to first nontrivial zeta zeros.
//! - Include spacing-stat objective so we are not only fitting point locations.

use gutoe_physics::riemann_lane::{
    build_operator, fit_against_reference, map_eigenvalues, zeta_zero_reference, RiemannMapParams,
    RiemannOperatorParams,
};
use nalgebra::SymmetricEigen;
use serde_json::json;
use std::cmp::Ordering;
use std::fs::{self, File};
use std::io::Write;

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

fn main() {
    let op = RiemannOperatorParams {
        n: env_usize("GUTOE_RIEMANN_DIM", 512),
        hop1_scale: env_f64("GUTOE_RIEMANN_HOP", 0.5),
        hop2_scale: env_f64("GUTOE_RIEMANN_HOP2", 0.0),
        potential_scale: env_f64("GUTOE_RIEMANN_POT", 0.0),
    };
    let map = RiemannMapParams {
        slope: env_f64("GUTOE_RIEMANN_SLOPE", 11.0 / 18.0),
        shift: env_f64("GUTOE_RIEMANN_SHIFT", 13.0 * 24.0 + 8.0 / 17.0),
    };

    let h = build_operator(op);
    let sym_err = (&h - h.transpose()).norm();
    let eig = SymmetricEigen::new(h);

    let mut evals: Vec<f64> = eig.eigenvalues.iter().copied().collect();
    evals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));

    let k_req = env_usize("GUTOE_RIEMANN_K", zeta_zero_reference().len());
    let k = k_req.min(zeta_zero_reference().len()).min(evals.len());
    let raw = &evals[..k];
    let pred = map_eigenvalues(raw, map);
    let target = &zeta_zero_reference()[..k];
    let fit = fit_against_reference(&pred, target);

    let rows: Vec<_> = (0..k)
        .map(|i| {
            let p = pred[i];
            let t = target[i];
            let abs_err = (p - t).abs();
            let rel_err = if t != 0.0 { abs_err / t } else { 0.0 };
            let signed_rel = if t != 0.0 { (p - t) / t } else { 0.0 };
            json!({
                "index": i + 1,
                "lambda_raw": raw[i],
                "gamma_pred": p,
                "gamma_target": t,
                "abs_err": abs_err,
                "rel_err": rel_err,
                "signed_rel_err": signed_rel
            })
        })
        .collect();

    let report = json!({
        "lane": "riemann_straight_shot",
        "operator": {
            "name": "clifford_berry_keating_family",
            "dimension": op.n,
            "hop1_scale": op.hop1_scale,
            "hop2_scale": op.hop2_scale,
            "potential_scale": op.potential_scale,
            "self_adjoint_residual_fro": sym_err
        },
        "mapping": {
            "gamma_pred": "shift + slope * lambda_raw",
            "affine_slope": map.slope,
            "affine_shift": map.shift
        },
        "comparison": {
            "k": fit.k,
            "mae": fit.mae,
            "rmse": fit.rmse,
            "mape": fit.mape,
            "max_abs_err": fit.max_abs_err,
            "max_rel_err": fit.max_rel_err,
            "signed_rel_bias": fit.signed_rel_bias,
            "spacing_mape": fit.spacing_mape,
            "spacing_rmse": fit.spacing_rmse,
            "spacing_ks": fit.spacing_ks,
            "spacing_var": fit.spacing_var,
            "spacing_var_abs_err_to_gue": fit.spacing_var_abs_err_to_gue,
            "objective_position": fit.objective_position,
            "objective_total": fit.objective_total
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
    writeln!(txt, "operator_dim            = {}", op.n).expect("write");
    writeln!(txt, "hop1_scale              = {:.12e}", op.hop1_scale).expect("write");
    writeln!(txt, "hop2_scale              = {:.12e}", op.hop2_scale).expect("write");
    writeln!(txt, "potential_scale         = {:.12e}", op.potential_scale).expect("write");
    writeln!(txt, "self_adjoint_residual   = {:.12e}", sym_err).expect("write");
    writeln!(txt, "affine_slope            = {:.12e}", map.slope).expect("write");
    writeln!(txt, "affine_shift            = {:.12e}", map.shift).expect("write");
    writeln!(txt, "k                       = {}", fit.k).expect("write");
    writeln!(txt, "MAE                     = {:.12e}", fit.mae).expect("write");
    writeln!(txt, "RMSE                    = {:.12e}", fit.rmse).expect("write");
    writeln!(txt, "MAPE                    = {:.12e}", fit.mape).expect("write");
    writeln!(txt, "max_abs_err             = {:.12e}", fit.max_abs_err).expect("write");
    writeln!(txt, "max_rel_err             = {:.12e}", fit.max_rel_err).expect("write");
    writeln!(txt, "signed_rel_bias         = {:.12e}", fit.signed_rel_bias).expect("write");
    writeln!(txt, "spacing_mape            = {:.12e}", fit.spacing_mape).expect("write");
    writeln!(txt, "spacing_rmse            = {:.12e}", fit.spacing_rmse).expect("write");
    writeln!(txt, "spacing_ks              = {:.12e}", fit.spacing_ks).expect("write");
    writeln!(txt, "spacing_var             = {:.12e}", fit.spacing_var).expect("write");
    writeln!(txt, "spacing_var_abs_err_gue = {:.12e}", fit.spacing_var_abs_err_to_gue).expect("write");
    writeln!(txt, "objective_position      = {:.12e}", fit.objective_position).expect("write");
    writeln!(txt, "objective_total         = {:.12e}", fit.objective_total).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "First 20 rows (n, pred, target, rel_err):").expect("write");
    for i in 0..k.min(20) {
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
