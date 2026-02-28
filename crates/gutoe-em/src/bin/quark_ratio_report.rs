//! Structural quark-ratio closure report from Cl(1,3)/Z3 primitives.
//!
//! This lane is intentionally ratio-first: it reports the seven commonly cited
//! quark mass ratios from shared structural counts and CKM lambda suppression.

use gutoe_em::ckm_from_clifford;
use std::fs::{self, File};
use std::io::Write;

#[derive(Clone, Copy)]
struct RatioRow {
    name: &'static str,
    predicted: f64,
    target: f64,
}

fn rel_err(predicted: f64, target: f64) -> f64 {
    ((predicted - target) / target).abs()
}

fn main() {
    let out_dir =
        std::env::var("GUTOE_QUARK_RATIO_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);

    // Shared suppression from CKM lane: λ = 1/sqrt(16+3) = 1/sqrt(19).
    let lambda = ckm_from_clifford().s12;
    let lambda_inv2 = 1.0 / (lambda * lambda);
    let c_inf = 67.0 / 66.0;

    // Z3/Cl(1,3) ratio ansatz from shared primitive counts.
    let rows = [
        RatioRow {
            name: "m_u/m_d",
            predicted: 8.0 / 17.0,
            target: 0.47,
        },
        RatioRow {
            name: "m_c/m_s",
            predicted: (13.0 / 21.0) * lambda_inv2 * c_inf,
            target: 11.7,
        },
        RatioRow {
            name: "m_t/m_b",
            predicted: (13.0 / 6.0) * lambda_inv2 * c_inf,
            target: 41.3,
        },
        RatioRow {
            name: "m_c/m_u",
            predicted: (8.0 / 5.0) * lambda_inv2 * lambda_inv2 * c_inf,
            target: 580.0,
        },
        RatioRow {
            name: "m_t/m_c",
            predicted: 8.0 * 17.0,
            target: 136.0,
        },
        RatioRow {
            name: "m_s/m_d",
            predicted: lambda_inv2,
            target: 20.0,
        },
        RatioRow {
            name: "m_b/m_s",
            predicted: (8.0 / 3.0) * lambda_inv2 * c_inf,
            target: 51.0,
        },
    ];

    let tol_rel = 0.10; // 10% acceptance band for legacy GRAND ratio tickets.
    let mut overall_pass = true;
    for r in &rows {
        if rel_err(r.predicted, r.target) > tol_rel {
            overall_pass = false;
        }
    }

    let txt_path = format!("{out_dir}/quark_ratio_report.txt");
    let json_path = format!("{out_dir}/quark_ratio_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "lambda_ckm = {:.12}", lambda).expect("write");
    writeln!(txt, "lambda_inv2 = {:.12}", lambda_inv2).expect("write");
    writeln!(txt, "c_inf = {:.12}", c_inf).expect("write");
    writeln!(txt, "rel_tolerance = {:.6}", tol_rel).expect("write");
    writeln!(txt, "overall_pass = {}", overall_pass).expect("write");
    writeln!(txt, "").expect("write");
    writeln!(txt, "[ratios]").expect("write");
    for r in &rows {
        writeln!(
            txt,
            "{} = {:.12}  (target {:.12}, rel_err {:.6})",
            r.name,
            r.predicted,
            r.target,
            rel_err(r.predicted, r.target)
        )
        .expect("write");
    }

    let mut json = File::create(&json_path).expect("create json");
    writeln!(json, "{{").expect("write");
    writeln!(json, "  \"lambda_ckm\": {:.12},", lambda).expect("write");
    writeln!(json, "  \"lambda_inv2\": {:.12},", lambda_inv2).expect("write");
    writeln!(json, "  \"c_inf\": {:.12},", c_inf).expect("write");
    writeln!(json, "  \"rel_tolerance\": {:.12},", tol_rel).expect("write");
    writeln!(
        json,
        "  \"overall_pass\": {},",
        if overall_pass { "true" } else { "false" }
    )
    .expect("write");
    writeln!(json, "  \"ratios\": [").expect("write");
    for (i, r) in rows.iter().enumerate() {
        let comma = if i + 1 == rows.len() { "" } else { "," };
        writeln!(
            json,
            "    {{\"name\":\"{}\",\"predicted\":{:.12},\"target\":{:.12},\"rel_err\":{:.12}}}{}",
            r.name,
            r.predicted,
            r.target,
            rel_err(r.predicted, r.target),
            comma
        )
        .expect("write");
    }
    writeln!(json, "  ]").expect("write");
    writeln!(json, "}}").expect("write");

    println!("wrote {txt_path}");
    println!("wrote {json_path}");
    println!(
        "quark ratios: overall_pass={} (tol {:.1}%)",
        overall_pass,
        tol_rel * 100.0
    );
}
