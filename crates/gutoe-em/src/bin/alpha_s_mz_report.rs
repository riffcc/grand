//! One-loop alpha_s report around the Z pole (GRAND-62 / GRAND-97).
//!
//! This tool:
//! - uses the Clifford-fixed beta coefficient from runtime config (`58/3`),
//! - infers `Lambda_QCD` from a target `alpha_s(M_Z)`,
//! - emits a small table showing asymptotic-freedom trend with energy.

use gutoe_em::config::LatticeConfig;
use std::fs::{self, File};
use std::io::Write;

const MZ_GEV: f64 = 91.1876;
const ALPHA_S_MZ_TARGET: f64 = 0.118;

fn alpha_s_one_loop(q_gev: f64, lambda_qcd_gev: f64, beta0: f64) -> f64 {
    let x = (q_gev / lambda_qcd_gev).ln();
    (2.0 * std::f64::consts::PI) / (beta0 * x)
}

fn infer_lambda_from_alpha_at_q(alpha_q: f64, q_gev: f64, beta0: f64) -> f64 {
    q_gev * (-(2.0 * std::f64::consts::PI) / (beta0 * alpha_q)).exp()
}

fn main() {
    let cfg = LatticeConfig::default();
    let beta0 = cfg.beta_coeff; // expected 58/3 from Clifford grades

    let lambda_qcd = infer_lambda_from_alpha_at_q(ALPHA_S_MZ_TARGET, MZ_GEV, beta0);
    let q_points = [10.0_f64, MZ_GEV, 1_000.0_f64, 10_000.0_f64];

    let out_dir = "/tmp/bh_renders";
    let _ = fs::create_dir_all(out_dir);
    let out_csv = format!("{out_dir}/alpha_s_mz_report.csv");

    let mut f = File::create(&out_csv).expect("create alpha_s report csv");
    writeln!(f, "q_gev,alpha_s_one_loop,beta0,lambda_qcd_gev").expect("header");

    let mut values = Vec::with_capacity(q_points.len());
    for q in q_points {
        let a = alpha_s_one_loop(q, lambda_qcd, beta0);
        values.push((q, a));
        writeln!(f, "{q:.6},{a:.9},{beta0:.9},{lambda_qcd:.9}").expect("row");
    }

    // Basic asymptotic-freedom sanity: alpha_s should decrease with higher Q.
    for w in values.windows(2) {
        let (q1, a1) = w[0];
        let (q2, a2) = w[1];
        assert!(q2 > q1, "Q grid must be increasing");
        assert!(
            a2 < a1,
            "asymptotic-freedom violation: alpha_s({q2})={a2} >= alpha_s({q1})={a1}"
        );
    }

    let alpha_mz = alpha_s_one_loop(MZ_GEV, lambda_qcd, beta0);
    println!("beta0 (Clifford) = {beta0:.9}");
    println!("M_Z = {MZ_GEV:.4} GeV");
    println!("inferred Lambda_QCD = {lambda_qcd:.9} GeV");
    println!("alpha_s(M_Z) = {alpha_mz:.9} (target {ALPHA_S_MZ_TARGET:.6})");
    println!("wrote {out_csv}");
}
