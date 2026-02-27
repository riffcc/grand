//! One-loop coupling unification scan (GRAND-63).
//!
//! Uses standard one-loop RG running for gauge couplings:
//!   alpha_i^{-1}(mu) = alpha_i^{-1}(MZ) - (b_i / 2pi) ln(mu/MZ)
//! and reports the scale where spread among (alpha1, alpha2, alpha3) is minimal.

use std::fs::{self, File};
use std::io::Write;

const MZ_GEV: f64 = 91.1876;

// MS-bar values near M_Z (representative).
const ALPHA1_MZ: f64 = 0.016923; // U(1)_Y in SU(5) normalization
const ALPHA2_MZ: f64 = 0.03374; // SU(2)_L
const ALPHA3_MZ: f64 = 0.1180; // SU(3)_c

// One-loop SM beta coefficients.
const B1: f64 = 41.0 / 10.0;
const B2: f64 = -19.0 / 6.0;
const B3: f64 = -7.0;

fn alpha_inv_running(alpha_inv_mz: f64, b: f64, mu_gev: f64) -> f64 {
    alpha_inv_mz - (b / (2.0 * std::f64::consts::PI)) * (mu_gev / MZ_GEV).ln()
}

fn main() {
    let a1_inv_mz = 1.0 / ALPHA1_MZ;
    let a2_inv_mz = 1.0 / ALPHA2_MZ;
    let a3_inv_mz = 1.0 / ALPHA3_MZ;

    let out_dir = "/tmp/bh_renders";
    let _ = fs::create_dir_all(out_dir);
    let csv_path = format!("{out_dir}/coupling_unification_report.csv");
    let txt_path = format!("{out_dir}/coupling_unification_summary.txt");

    let mut csv = File::create(&csv_path).expect("create unification csv");
    writeln!(csv, "mu_gev,alpha1_inv,alpha2_inv,alpha3_inv,spread_inv").expect("header");

    let mut best_mu = MZ_GEV;
    let mut best_spread = f64::INFINITY;
    let mut best_vals = (a1_inv_mz, a2_inv_mz, a3_inv_mz);

    // Log10 scan from 10^2 to 10^19 GeV.
    for i in 0..=1700 {
        let log10_mu = 2.0 + (17.0 * i as f64) / 1700.0;
        let mu = 10f64.powf(log10_mu);
        let a1 = alpha_inv_running(a1_inv_mz, B1, mu);
        let a2 = alpha_inv_running(a2_inv_mz, B2, mu);
        let a3 = alpha_inv_running(a3_inv_mz, B3, mu);
        let max_v = a1.max(a2).max(a3);
        let min_v = a1.min(a2).min(a3);
        let spread = max_v - min_v;

        writeln!(csv, "{mu:.8e},{a1:.9},{a2:.9},{a3:.9},{spread:.9}").expect("row");

        if spread < best_spread {
            best_spread = spread;
            best_mu = mu;
            best_vals = (a1, a2, a3);
        }
    }

    let mut txt = File::create(&txt_path).expect("create summary txt");
    writeln!(txt, "MZ_GEV={MZ_GEV:.6}").unwrap();
    writeln!(
        txt,
        "alpha1_MZ={ALPHA1_MZ:.9} alpha2_MZ={ALPHA2_MZ:.9} alpha3_MZ={ALPHA3_MZ:.9}"
    )
    .unwrap();
    writeln!(txt, "best_mu_gev={best_mu:.8e}").unwrap();
    writeln!(txt, "best_spread_alpha_inv={best_spread:.9}").unwrap();
    writeln!(
        txt,
        "alpha_inv_at_best=({:.9}, {:.9}, {:.9})",
        best_vals.0, best_vals.1, best_vals.2
    )
    .unwrap();

    println!("best unification-like scale mu = {best_mu:.4e} GeV");
    println!("min spread in alpha^-1 = {best_spread:.6}");
    println!(
        "alpha^-1(mu*): ({:.6}, {:.6}, {:.6})",
        best_vals.0, best_vals.1, best_vals.2
    );
    println!("wrote {csv_path} and {txt_path}");
}
