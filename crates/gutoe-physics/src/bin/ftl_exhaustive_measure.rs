//! FTL exhaustive measurement lane (no hardcoded pass/fail gates).
//!
//! This bin computes margins directly from equations:
//! 1) Casimir negative-energy density sweep vs EW/Higgs restoration scale.
//! 2) Higgs-orientation wall surfing kinematics with front vs rear-face (1/10) tension.
//!
//! Output is numeric evidence, not categorical hand-gating.

use gutoe_physics::constants::{C, HBAR, HIGGS_QUARTIC_STRUCTURAL};
use serde_json::json;
use std::f64::consts::PI;
use std::fs;
use std::path::PathBuf;

const GEV_TO_J: f64 = 1.602_176_634e-10;
const HBARC_GEV_M: f64 = 0.197_326_980_4e-15;
const V_EWSB_GEV: f64 = 245.3;
const FC_VOID: f64 = 3.0 / 16.0;
const REAR_FACE_FACTOR: f64 = 1.0 / 10.0;

fn higgs_restoration_density_j_m3() -> f64 {
    // ΔV = λ v^4 / 4 (natural units GeV^4), then convert to J/m^3.
    let delta_v_gev4 = HIGGS_QUARTIC_STRUCTURAL * V_EWSB_GEV.powi(4) / 4.0;
    let gev4_to_j_m3 = GEV_TO_J / HBARC_GEV_M.powi(3);
    delta_v_gev4 * gev4_to_j_m3
}

fn casimir_density_j_m3(gap_m: f64) -> f64 {
    PI.powi(2) * HBAR * C / (720.0 * gap_m.powi(4))
}

fn wall_tension_front_j_m2(delta_theta: f64, thickness_m: f64) -> f64 {
    // σ = f_c * v^2 * (Δθ)^2 / (2 L) in natural units (GeV^3), then convert to J/m^2.
    let l_nat = thickness_m / HBARC_GEV_M; // GeV^-1
    let sigma_gev3 = FC_VOID * V_EWSB_GEV.powi(2) * delta_theta.powi(2) / (2.0 * l_nat);
    let gev3_to_j_m2 = GEV_TO_J / HBARC_GEV_M.powi(2);
    sigma_gev3 * gev3_to_j_m2
}

fn gamma_from_beta(beta: f64) -> f64 {
    1.0 / (1.0 - beta * beta).sqrt()
}

fn beta_from_areal_drive(e_areal_j_m2: f64, sigma_j_m2: f64) -> f64 {
    if sigma_j_m2 <= 0.0 {
        return f64::NAN;
    }
    let gamma = 1.0 + e_areal_j_m2 / sigma_j_m2;
    if gamma <= 1.0 {
        0.0
    } else {
        (1.0 - 1.0 / (gamma * gamma)).sqrt()
    }
}

fn format_sci(v: f64) -> String {
    format!("{v:.12e}")
}

fn main() {
    let out_dir = std::env::var("GUTOE_FTL_EXHAUSTIVE_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ftl_exhaustive_measure".to_string());
    let out = PathBuf::from(&out_dir);
    let _ = fs::create_dir_all(&out);

    let rho_higgs = higgs_restoration_density_j_m3();

    // Casimir sweep down to extreme gaps.
    let casimir_gaps = [
        1e-6_f64, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15,
    ];
    let mut casimir_rows = Vec::new();
    let mut casimir_max = 0.0_f64;
    let mut casimir_max_gap = 0.0_f64;
    for &g in &casimir_gaps {
        let u = casimir_density_j_m3(g);
        let ratio = u / rho_higgs;
        if u > casimir_max {
            casimir_max = u;
            casimir_max_gap = g;
        }
        casimir_rows.push((g, u, ratio));
    }

    // Wall-surf sweep (front vs rear) for representative thickness/orientation scales.
    let thicknesses_m = [1e-18_f64, 1e-15, 1e-12, 1e-9];
    let dthetas = [0.1_f64, 1.0, PI / 2.0, PI];
    let beta_targets = [0.5_f64, 0.9, 0.99, 0.999, 0.9999];

    let mut wall_rows = Vec::new();
    for &l in &thicknesses_m {
        for &dt in &dthetas {
            let sigma_front = wall_tension_front_j_m2(dt, l);
            let sigma_rear = REAR_FACE_FACTOR * sigma_front;
            for &beta_t in &beta_targets {
                let gamma_t = gamma_from_beta(beta_t);
                let e_front = (gamma_t - 1.0) * sigma_front;
                let e_rear = (gamma_t - 1.0) * sigma_rear;
                wall_rows.push((l, dt, beta_t, sigma_front, sigma_rear, e_front, e_rear));
            }
        }
    }

    // Direct finite-energy sweep for beta ceilings.
    let e_over_sigma_grid = [
        1e-9_f64, 1e-6, 1e-3, 1e-1, 1.0, 10.0, 1e2, 1e3, 1e4, 1e6, 1e9,
    ];
    let mut beta_grid_rows = Vec::new();
    let mut beta_front_max = 0.0_f64;
    let mut beta_rear_max = 0.0_f64;
    let mut x_max = 0.0_f64;
    for &x in &e_over_sigma_grid {
        let beta_f = beta_from_areal_drive(x, 1.0);
        let beta_r = beta_from_areal_drive(x, REAR_FACE_FACTOR);
        beta_front_max = beta_front_max.max(beta_f);
        beta_rear_max = beta_rear_max.max(beta_r);
        x_max = x_max.max(x);
        beta_grid_rows.push((x, beta_f, beta_r));
    }

    let ftl_front_detected = beta_front_max > 1.0;
    let ftl_rear_detected = beta_rear_max > 1.0;
    let front_min_one_minus_beta_sq = 1.0 / (1.0 + x_max).powi(2);
    let rear_min_one_minus_beta_sq = 1.0 / (1.0 + x_max / REAR_FACE_FACTOR).powi(2);

    let txt_path = out.join("ftl_exhaustive_measure.txt");
    let csv_path = out.join("ftl_exhaustive_measure.csv");
    let json_path = out.join("ftl_exhaustive_measure.json");

    let mut txt = String::new();
    txt.push_str("[ftl_exhaustive_measure]\n");
    txt.push_str(&format!("higgs_restoration_density_j_m3 = {}\n", format_sci(rho_higgs)));
    txt.push_str(&format!(
        "casimir_max_density_j_m3 = {} at gap_m={}\n",
        format_sci(casimir_max),
        format_sci(casimir_max_gap)
    ));
    txt.push_str(&format!(
        "casimir_to_higgs_ratio_max = {}\n",
        format_sci(casimir_max / rho_higgs)
    ));
    txt.push_str(&format!(
        "casimir_deficit_orders_of_magnitude = {:.6}\n",
        (rho_higgs / casimir_max).log10()
    ));
    txt.push_str(&format!("rear_face_factor = {:.12}\n", REAR_FACE_FACTOR));
    txt.push_str(&format!("beta_front_max_finite_sweep = {:.12}\n", beta_front_max));
    txt.push_str(&format!("beta_rear_max_finite_sweep = {:.12}\n", beta_rear_max));
    txt.push_str(&format!(
        "front_min_1_minus_beta_sq_from_gamma = {}\n",
        format_sci(front_min_one_minus_beta_sq)
    ));
    txt.push_str(&format!(
        "rear_min_1_minus_beta_sq_from_gamma = {}\n",
        format_sci(rear_min_one_minus_beta_sq)
    ));
    txt.push_str(&format!("front_ftl_detected = {}\n", ftl_front_detected));
    txt.push_str(&format!("rear_ftl_detected = {}\n", ftl_rear_detected));
    txt.push_str("\n[wall_energy_requirements]\n");
    txt.push_str("columns = thickness_m,delta_theta,beta_target,sigma_front_j_m2,sigma_rear_j_m2,e_areal_front_j_m2,e_areal_rear_j_m2\n");
    for (l, dt, bt, sf, sr, ef, er) in &wall_rows {
        txt.push_str(&format!(
            "{},{},{},{},{},{},{}\n",
            format_sci(*l),
            format_sci(*dt),
            format_sci(*bt),
            format_sci(*sf),
            format_sci(*sr),
            format_sci(*ef),
            format_sci(*er)
        ));
    }

    let mut csv = String::from(
        "mode,e_over_sigma,beta_front,beta_rear,thickness_m,delta_theta,beta_target,sigma_front_j_m2,sigma_rear_j_m2,e_areal_front_j_m2,e_areal_rear_j_m2\n",
    );
    for (x, bf, br) in &beta_grid_rows {
        csv.push_str(&format!(
            "beta_sweep,{},{},{} ,,,,,,,\n",
            format_sci(*x),
            format_sci(*bf),
            format_sci(*br)
        ));
    }
    for (l, dt, bt, sf, sr, ef, er) in &wall_rows {
        csv.push_str(&format!(
            "wall_requirements,,,,{},{},{},{},{},{},{}\n",
            format_sci(*l),
            format_sci(*dt),
            format_sci(*bt),
            format_sci(*sf),
            format_sci(*sr),
            format_sci(*ef),
            format_sci(*er)
        ));
    }

    let casimir_json_rows: Vec<_> = casimir_rows
        .iter()
        .map(|(g, u, r)| {
            json!({
                "gap_m": g,
                "casimir_density_j_m3": u,
                "ratio_to_higgs_restoration": r
            })
        })
        .collect();

    let wall_json_rows: Vec<_> = wall_rows
        .iter()
        .map(|(l, dt, bt, sf, sr, ef, er)| {
            json!({
                "thickness_m": l,
                "delta_theta": dt,
                "beta_target": bt,
                "sigma_front_j_m2": sf,
                "sigma_rear_j_m2": sr,
                "e_areal_front_j_m2": ef,
                "e_areal_rear_j_m2": er
            })
        })
        .collect();

    let beta_json_rows: Vec<_> = beta_grid_rows
        .iter()
        .map(|(x, bf, br)| {
            json!({
                "e_over_sigma_front": x,
                "beta_front": bf,
                "beta_rear": br
            })
        })
        .collect();

    let payload = json!({
        "higgs_restoration_density_j_m3": rho_higgs,
        "casimir_max_density_j_m3": casimir_max,
        "casimir_max_gap_m": casimir_max_gap,
        "casimir_to_higgs_ratio_max": casimir_max / rho_higgs,
        "casimir_deficit_orders_of_magnitude": (rho_higgs / casimir_max).log10(),
        "rear_face_factor": REAR_FACE_FACTOR,
        "beta_front_max_finite_sweep": beta_front_max,
        "beta_rear_max_finite_sweep": beta_rear_max,
        "front_min_1_minus_beta_sq_from_gamma": front_min_one_minus_beta_sq,
        "rear_min_1_minus_beta_sq_from_gamma": rear_min_one_minus_beta_sq,
        "front_ftl_detected": ftl_front_detected,
        "rear_ftl_detected": ftl_rear_detected,
        "casimir_sweep": casimir_json_rows,
        "beta_sweep": beta_json_rows,
        "wall_energy_requirements": wall_json_rows
    });

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(&csv_path, csv).expect("write csv");
    fs::write(&json_path, serde_json::to_string_pretty(&payload).expect("json")).expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", csv_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "casimir_ratio_max={:.6e}, beta_front_max={:.9}, beta_rear_max={:.9}",
        casimir_max / rho_higgs,
        beta_front_max,
        beta_rear_max
    );
}
