//! GRAND-354: BBN likelihood contour scan for η10/Ωb sensitivity.

use gutoe_physics::{
    evaluate_bbn_gate, evaluate_microphysics_gate, evaluate_universe_gate, BbnWindows,
    MicrophysicsAssumptions, MicrophysicsWindows, UniverseAssumptions, UniverseWindows,
    DARK_TO_VISIBLE_GEOMETRIC_RATIO,
};
use std::fs::{self, File};
use std::io::Write;

const YP_OBS: f64 = 0.245;
const SIGMA_YP: f64 = 0.003;
const DH_OBS: f64 = 2.547e-5;
const SIGMA_DH: f64 = 0.050e-5;

fn chi2_bbn(yp: f64, dh: f64) -> f64 {
    ((yp - YP_OBS) / SIGMA_YP).powi(2) + ((dh - DH_OBS) / SIGMA_DH).powi(2)
}

fn micro(
    h0_km_s_mpc: f64,
    omega_b0: f64,
    omega_r0: f64,
    omega_k0: f64,
    eta10: f64,
) -> Option<(f64, f64, f64)> {
    let omega_dm0 = omega_b0 * DARK_TO_VISIBLE_GEOMETRIC_RATIO;
    let omega_m0 = omega_b0 + omega_dm0;
    let omega_lambda0 = 1.0 - omega_m0 - omega_r0 - omega_k0;
    if omega_lambda0 <= 0.0 {
        return None;
    }

    let s = evaluate_microphysics_gate(
        MicrophysicsAssumptions {
            h0_km_s_mpc,
            omega_b0,
            omega_m0,
            omega_r0,
            omega_k0,
            omega_lambda0,
            eta10,
        },
        MicrophysicsWindows::default(),
    );
    Some((s.yp_network, s.dh_network, s.z_visibility_peak))
}

fn main() {
    let u = evaluate_universe_gate(UniverseAssumptions::default(), UniverseWindows::default());
    let bbn = evaluate_bbn_gate(BbnWindows::default());

    let h0 = u.h0_km_s_mpc;
    let omega_b0 = u.omega_b0;
    let omega_r0 = u.omega_r0;
    let omega_k0 = u.omega_k0;
    let eta10_0 = bbn.eta10;

    let out_dir =
        std::env::var("GUTOE_BBN_SCAN_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let csv_path = format!("{out_dir}/bbn_likelihood_scan.csv");
    let json_path = format!("{out_dir}/bbn_likelihood_scan.json");

    let mut csv = File::create(&csv_path).expect("create csv");
    writeln!(
        csv,
        "omega_b_factor,eta10_factor,omega_b0,eta10,yp,dh,z_visibility,chi2"
    )
    .expect("write header");

    let n = 81usize;
    let omega_b_min = 0.90;
    let omega_b_max = 1.10;
    let eta_min = 0.80;
    let eta_max = 1.20;

    let mut best_chi2 = f64::INFINITY;
    let mut best = (1.0, 1.0, omega_b0, eta10_0, f64::NAN, f64::NAN, f64::NAN);

    let mut count = 0usize;
    let mut count_1s = 0usize;
    let mut count_2s = 0usize;
    let mut count_3s = 0usize;

    for i in 0..n {
        let fb = omega_b_min + (omega_b_max - omega_b_min) * (i as f64) / ((n - 1) as f64);
        for j in 0..n {
            let fe = eta_min + (eta_max - eta_min) * (j as f64) / ((n - 1) as f64);
            let ob = omega_b0 * fb;
            let et = eta10_0 * fe;

            let Some((yp, dh, zv)) = micro(h0, ob, omega_r0, omega_k0, et) else {
                continue;
            };
            let chi2 = chi2_bbn(yp, dh);
            writeln!(
                csv,
                "{:.9},{:.9},{:.12},{:.12},{:.12},{:.12e},{:.6},{:.12}",
                fb, fe, ob, et, yp, dh, zv, chi2
            )
            .expect("write row");

            count += 1;
            if chi2 <= 2.30 {
                count_1s += 1;
            }
            if chi2 <= 6.18 {
                count_2s += 1;
            }
            if chi2 <= 11.83 {
                count_3s += 1;
            }

            if chi2 < best_chi2 {
                best_chi2 = chi2;
                best = (fb, fe, ob, et, yp, dh, zv);
            }
        }
    }

    let eps = 0.01;
    let base = micro(h0, omega_b0, omega_r0, omega_k0, eta10_0).expect("base micro should work");
    let ob_up = micro(h0, omega_b0 * (1.0 + eps), omega_r0, omega_k0, eta10_0)
        .expect("ob+ micro should work");
    let ob_dn = micro(h0, omega_b0 * (1.0 - eps), omega_r0, omega_k0, eta10_0)
        .expect("ob- micro should work");
    let et_up = micro(h0, omega_b0, omega_r0, omega_k0, eta10_0 * (1.0 + eps))
        .expect("eta+ micro should work");
    let et_dn = micro(h0, omega_b0, omega_r0, omega_k0, eta10_0 * (1.0 - eps))
        .expect("eta- micro should work");

    let dln_dh_dln_ob = ((ob_up.1 / ob_dn.1).ln()) / (((1.0 + eps) / (1.0 - eps)).ln());
    let dln_dh_dln_eta = ((et_up.1 / et_dn.1).ln()) / (((1.0 + eps) / (1.0 - eps)).ln());

    let mut json = File::create(&json_path).expect("create json");
    writeln!(
        json,
        "{{\n  \"baseline\": {{\"omega_b0\": {:.12}, \"eta10\": {:.12}, \"yp\": {:.12}, \"dh\": {:.12e}, \"z_visibility\": {:.6}}},\n  \"scan\": {{\"grid_n\": {}, \"points\": {}, \"omega_b_factor_min\": {:.6}, \"omega_b_factor_max\": {:.6}, \"eta10_factor_min\": {:.6}, \"eta10_factor_max\": {:.6}}},\n  \"best_fit\": {{\"chi2\": {:.12}, \"omega_b_factor\": {:.9}, \"eta10_factor\": {:.9}, \"omega_b0\": {:.12}, \"eta10\": {:.12}, \"yp\": {:.12}, \"dh\": {:.12e}, \"z_visibility\": {:.6}}},\n  \"contours\": {{\"chi2_le_2p30_fraction\": {:.9}, \"chi2_le_6p18_fraction\": {:.9}, \"chi2_le_11p83_fraction\": {:.9}}},\n  \"sensitivity\": {{\"dln_dh_dln_omega_b\": {:.9}, \"dln_dh_dln_eta10\": {:.9}}}\n}}",
        omega_b0,
        eta10_0,
        base.0,
        base.1,
        base.2,
        n,
        count,
        omega_b_min,
        omega_b_max,
        eta_min,
        eta_max,
        best_chi2,
        best.0,
        best.1,
        best.2,
        best.3,
        best.4,
        best.5,
        best.6,
        count_1s as f64 / count.max(1) as f64,
        count_2s as f64 / count.max(1) as f64,
        count_3s as f64 / count.max(1) as f64,
        dln_dh_dln_ob,
        dln_dh_dln_eta,
    )
    .expect("write json");

    println!("wrote {csv_path}");
    println!("wrote {json_path}");
    println!(
        "BBN scan: best chi2={:.4}, fb={:.4}, fe={:.4}, Yp={:.5}, D/H={:.3e}",
        best_chi2, best.0, best.1, best.4, best.5
    );
    println!(
        "D/H sensitivity: dln(D/H)/dln(Ωb)={:.3}, dln(D/H)/dln(η10)={:.3}",
        dln_dh_dln_ob, dln_dh_dln_eta
    );
}
