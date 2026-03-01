//! Phase-2 RH exploratory scan:
//! - Expand operator family beyond baseline tridiagonal.
//! - Score candidates on position fit + spacing statistics.

use gutoe_physics::riemann_lane::{
    build_operator, fit_against_reference, map_eigenvalues, zeta_zero_reference, RiemannMapParams,
    RiemannOperatorParams,
};
use nalgebra::SymmetricEigen;
use serde_json::json;
use std::cmp::Ordering;
use std::fs::{self, File};
use std::io::Write;

#[derive(Debug, Clone)]
struct Candidate {
    family: &'static str,
    op: RiemannOperatorParams,
    fit: gutoe_physics::riemann_lane::RiemannFitStats,
    sym_resid: f64,
}

fn env_usize(key: &str, default: usize) -> usize {
    std::env::var(key)
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(default)
}

fn evaluate_candidate(op: RiemannOperatorParams, map: RiemannMapParams, k: usize) -> Candidate {
    let h = build_operator(op);
    let sym_resid = (&h - h.transpose()).norm();
    let eig = SymmetricEigen::new(h);
    let mut evals: Vec<f64> = eig.eigenvalues.iter().copied().collect();
    evals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
    let k_eff = k.min(evals.len()).min(zeta_zero_reference().len());
    let raw = &evals[..k_eff];
    let pred = map_eigenvalues(raw, map);
    let fit = fit_against_reference(&pred, &zeta_zero_reference()[..k_eff]);
    Candidate {
        family: "",
        op,
        fit,
        sym_resid,
    }
}

fn main() {
    let n = env_usize("GUTOE_RIEMANN_SCAN_DIM", 512);
    let k = env_usize("GUTOE_RIEMANN_SCAN_K", 80);

    let map = RiemannMapParams {
        slope: 11.0 / 18.0,
        shift: 13.0 * 24.0 + 8.0 / 17.0,
    };

    // Tight high-fidelity neighborhood around the baseline lane.
    // Keep candidate count modest for n=512 eigensolves.
    let hop1_vals = vec![0.50, 0.52];
    let hop2_vals = vec![0.0, -0.02];
    let pot_vals = vec![0.0, -0.10];

    let mut all: Vec<Candidate> = Vec::new();

    for &h1 in &hop1_vals {
        let mut c = evaluate_candidate(
            RiemannOperatorParams {
                n,
                hop1_scale: h1,
                hop2_scale: 0.0,
                potential_scale: 0.0,
            },
            map,
            k,
        );
        c.family = "baseline";
        all.push(c);
    }

    for &h1 in &hop1_vals {
        for &h2 in &hop2_vals {
            let mut c = evaluate_candidate(
                RiemannOperatorParams {
                    n,
                    hop1_scale: h1,
                    hop2_scale: h2,
                    potential_scale: 0.0,
                },
                map,
                k,
            );
            c.family = "nnn";
            all.push(c);
        }
    }

    for &h1 in &hop1_vals {
        for &pot in &pot_vals {
            let mut c = evaluate_candidate(
                RiemannOperatorParams {
                    n,
                    hop1_scale: h1,
                    hop2_scale: 0.0,
                    potential_scale: pot,
                },
                map,
                k,
            );
            c.family = "potential";
            all.push(c);
        }
    }

    for &h1 in &hop1_vals {
        for &h2 in &hop2_vals {
            for &pot in &pot_vals {
                let mut c = evaluate_candidate(
                    RiemannOperatorParams {
                        n,
                        hop1_scale: h1,
                        hop2_scale: h2,
                        potential_scale: pot,
                    },
                    map,
                    k,
                );
                c.family = "full";
                all.push(c);
            }
        }
    }

    all.sort_by(|a, b| {
        a.fit
            .objective_total
            .partial_cmp(&b.fit.objective_total)
            .unwrap_or(Ordering::Equal)
    });

    let best = all.first().expect("non-empty scan");
    let mut top = Vec::new();
    for (rank, c) in all.iter().take(20).enumerate() {
        top.push(json!({
            "rank": rank + 1,
            "family": c.family,
            "hop1_scale": c.op.hop1_scale,
            "hop2_scale": c.op.hop2_scale,
            "potential_scale": c.op.potential_scale,
            "sym_resid": c.sym_resid,
            "mape": c.fit.mape,
            "rmse": c.fit.rmse,
            "spacing_ks": c.fit.spacing_ks,
            "spacing_mape": c.fit.spacing_mape,
            "spacing_var": c.fit.spacing_var,
            "objective_total": c.fit.objective_total,
        }));
    }

    let report = json!({
        "lane": "riemann_operator_family_scan",
        "scan": {
            "dimension": n,
            "k": k,
            "families": ["baseline", "nnn", "potential", "full"],
            "hop1_axis": hop1_vals,
            "hop2_axis": hop2_vals,
            "potential_axis": pot_vals,
            "candidate_count": all.len()
        },
        "objective": {
            "formula": "mape + 0.50*spacing_ks + 0.25*spacing_mape + 0.25*spacing_var_abs_err_to_gue"
        },
        "best": {
            "family": best.family,
            "hop1_scale": best.op.hop1_scale,
            "hop2_scale": best.op.hop2_scale,
            "potential_scale": best.op.potential_scale,
            "sym_resid": best.sym_resid,
            "mape": best.fit.mape,
            "rmse": best.fit.rmse,
            "spacing_ks": best.fit.spacing_ks,
            "spacing_mape": best.fit.spacing_mape,
            "spacing_var": best.fit.spacing_var,
            "spacing_var_abs_err_to_gue": best.fit.spacing_var_abs_err_to_gue,
            "objective_total": best.fit.objective_total,
        },
        "top20": top
    });

    let out_dir = std::env::var("GUTOE_RIEMANN_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    fs::create_dir_all(&out_dir).expect("create output directory");
    let txt_path = format!("{out_dir}/riemann_operator_family_scan.txt");
    let json_path = format!("{out_dir}/riemann_operator_family_scan.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "Riemann Operator Family Scan").expect("write");
    writeln!(txt, "===========================").expect("write");
    writeln!(txt, "dimension        = {}", n).expect("write");
    writeln!(txt, "k                = {}", k).expect("write");
    writeln!(txt, "candidate_count  = {}", all.len()).expect("write");
    writeln!(txt, "\nBest candidate:").expect("write");
    writeln!(txt, "family           = {}", best.family).expect("write");
    writeln!(txt, "hop1_scale       = {:.12e}", best.op.hop1_scale).expect("write");
    writeln!(txt, "hop2_scale       = {:.12e}", best.op.hop2_scale).expect("write");
    writeln!(txt, "potential_scale  = {:.12e}", best.op.potential_scale).expect("write");
    writeln!(txt, "mape             = {:.12e}", best.fit.mape).expect("write");
    writeln!(txt, "rmse             = {:.12e}", best.fit.rmse).expect("write");
    writeln!(txt, "spacing_ks       = {:.12e}", best.fit.spacing_ks).expect("write");
    writeln!(txt, "spacing_mape     = {:.12e}", best.fit.spacing_mape).expect("write");
    writeln!(txt, "spacing_var      = {:.12e}", best.fit.spacing_var).expect("write");
    writeln!(txt, "objective_total  = {:.12e}", best.fit.objective_total).expect("write");
    writeln!(txt, "\nTop 10 by objective:").expect("write");
    for (rank, c) in all.iter().take(10).enumerate() {
        writeln!(
            txt,
            "{:2}. {:9} h1={:+.3} h2={:+.3} pot={:+.3}  mape={:.5e}  ks={:.5e}  obj={:.5e}",
            rank + 1,
            c.family,
            c.op.hop1_scale,
            c.op.hop2_scale,
            c.op.potential_scale,
            c.fit.mape,
            c.fit.spacing_ks,
            c.fit.objective_total,
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
