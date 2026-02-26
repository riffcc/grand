//! GRAND-294: FRW / H(z) phenomenology harness from derived GUTOE Λ terms.

use gutoe_physics::constants::{
    lambda_cosmological_full_candidate, lambda_cosmological_signature_candidate,
    lambda_cosmological_structural, C, LAMBDA_COSMOLOGICAL_OBSERVED,
};
use std::fs::{self, File};
use std::io::Write;

const METER_PER_MPC: f64 = 3.085_677_581_491_367e22;
const DEFAULT_OMEGA_M0: f64 = 0.315;
const DEFAULT_OMEGA_R0: f64 = 9.0e-5;
const DEFAULT_OMEGA_K0: f64 = 0.0;
const PLANCK_H0_KM_S_MPC: f64 = 67.4;
const DISTANCE_LADDER_H0_KM_S_MPC: f64 = 73.0;

#[derive(Debug, Clone, Copy)]
struct FrwAssumptions {
    omega_m0: f64,
    omega_r0: f64,
    omega_k0: f64,
}

fn km_s_mpc_to_s_inv(h0_km_s_mpc: f64) -> f64 {
    (h0_km_s_mpc * 1_000.0) / METER_PER_MPC
}

fn s_inv_to_km_s_mpc(h0_s_inv: f64) -> f64 {
    h0_s_inv * METER_PER_MPC / 1_000.0
}

/// Ω_Λ = Λ c² / (3 H0²).
fn omega_lambda_from_lambda_and_h0(lambda: f64, h0_km_s_mpc: f64) -> f64 {
    let h0 = km_s_mpc_to_s_inv(h0_km_s_mpc);
    lambda * C * C / (3.0 * h0 * h0)
}

/// H0 = c * sqrt(Λ / (3 Ω_Λ)).
fn h0_from_lambda_and_omega_lambda(lambda: f64, omega_lambda: f64) -> Option<f64> {
    if lambda <= 0.0 || omega_lambda <= 0.0 {
        return None;
    }
    let h0_s_inv = C * (lambda / (3.0 * omega_lambda)).sqrt();
    Some(s_inv_to_km_s_mpc(h0_s_inv))
}

/// E(z)^2 = Ω_r(1+z)^4 + Ω_m(1+z)^3 + Ω_k(1+z)^2 + Ω_Λ.
fn e2_of_z(z: f64, omega_r: f64, omega_m: f64, omega_k: f64, omega_lambda: f64) -> f64 {
    let one_plus_z = 1.0 + z;
    omega_r * one_plus_z.powi(4)
        + omega_m * one_plus_z.powi(3)
        + omega_k * one_plus_z.powi(2)
        + omega_lambda
}

fn h_of_z(h0_km_s_mpc: f64, z: f64, omega_r: f64, omega_m: f64, omega_k: f64, omega_lambda: f64) -> f64 {
    let e2 = e2_of_z(z, omega_r, omega_m, omega_k, omega_lambda);
    if e2 <= 0.0 {
        return f64::NAN;
    }
    h0_km_s_mpc * e2.sqrt()
}

fn env_f64(name: &str, default: f64) -> f64 {
    std::env::var(name)
        .ok()
        .and_then(|v| v.parse::<f64>().ok())
        .unwrap_or(default)
}

fn main() {
    let assumptions = FrwAssumptions {
        omega_m0: env_f64("GUTOE_OMEGA_M0", DEFAULT_OMEGA_M0),
        omega_r0: env_f64("GUTOE_OMEGA_R0", DEFAULT_OMEGA_R0),
        omega_k0: env_f64("GUTOE_OMEGA_K0", DEFAULT_OMEGA_K0),
    };
    let omega_lambda_flat = 1.0 - assumptions.omega_m0 - assumptions.omega_r0 - assumptions.omega_k0;

    let lambdas = [
        ("observed", LAMBDA_COSMOLOGICAL_OBSERVED),
        ("structural", lambda_cosmological_structural()),
        ("signature", lambda_cosmological_signature_candidate()),
        ("full", lambda_cosmological_full_candidate()),
    ];

    let mut variants = Vec::new();
    for (name, lambda) in lambdas {
        let h0_flat = h0_from_lambda_and_omega_lambda(lambda, omega_lambda_flat);
        let omega_lambda_planck = omega_lambda_from_lambda_and_h0(lambda, PLANCK_H0_KM_S_MPC);
        let omega_lambda_ladder = omega_lambda_from_lambda_and_h0(lambda, DISTANCE_LADDER_H0_KM_S_MPC);
        let omega_k_planck = 1.0 - assumptions.omega_m0 - assumptions.omega_r0 - omega_lambda_planck;
        let omega_k_ladder = 1.0 - assumptions.omega_m0 - assumptions.omega_r0 - omega_lambda_ladder;
        variants.push((
            name,
            lambda,
            h0_flat,
            omega_lambda_planck,
            omega_lambda_ladder,
            omega_k_planck,
            omega_k_ladder,
        ));
    }

    let lambda_full = lambda_cosmological_full_candidate();
    let h0_full_flat = h0_from_lambda_and_omega_lambda(lambda_full, omega_lambda_flat)
        .expect("positive lambda and omega_lambda_flat required");

    let mut curve_rows = Vec::new();
    for i in 0..=20 {
        let z = i as f64 * 0.25;
        let h = h_of_z(
            h0_full_flat,
            z,
            assumptions.omega_r0,
            assumptions.omega_m0,
            assumptions.omega_k0,
            omega_lambda_flat,
        );
        curve_rows.push((z, h, h / h0_full_flat));
    }

    let out_dir = "/tmp/bh_renders";
    let _ = fs::create_dir_all(out_dir);
    let txt_path = format!("{out_dir}/frw_hz_report.txt");
    let json_path = format!("{out_dir}/frw_hz_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "GRAND-294 FRW / H(z) harness").expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[assumptions]").expect("write");
    writeln!(txt, "omega_m0 = {:.8}", assumptions.omega_m0).expect("write");
    writeln!(txt, "omega_r0 = {:.8}", assumptions.omega_r0).expect("write");
    writeln!(txt, "omega_k0 = {:.8}", assumptions.omega_k0).expect("write");
    writeln!(txt, "omega_lambda_flat = {:.8}", omega_lambda_flat).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[variants]").expect("write");
    for (name, lambda, h0_flat, omega_lambda_planck, omega_lambda_ladder, omega_k_planck, omega_k_ladder) in &variants {
        writeln!(
            txt,
            "{}: lambda={:.12e}, h0_flat={:.6}, omega_lambda_at_planck_h0={:.8}, omega_lambda_at_ladder_h0={:.8}, omega_k_at_planck_h0={:.8}, omega_k_at_ladder_h0={:.8}",
            name,
            lambda,
            h0_flat.unwrap_or(f64::NAN),
            omega_lambda_planck,
            omega_lambda_ladder,
            omega_k_planck,
            omega_k_ladder
        )
        .expect("write");
    }
    writeln!(txt).expect("write");
    writeln!(txt, "[full_candidate_curve]").expect("write");
    writeln!(txt, "h0_full_flat = {:.6}", h0_full_flat).expect("write");
    for (z, h, e) in &curve_rows {
        writeln!(txt, "z={:.2}, H(z)={:.6}, E(z)={:.6}", z, h, e).expect("write");
    }

    let mut json = File::create(&json_path).expect("create json");
    writeln!(
        json,
        "{{\n  \"assumptions\": {{\"omega_m0\": {:.10}, \"omega_r0\": {:.10}, \"omega_k0\": {:.10}, \"omega_lambda_flat\": {:.10}}},",
        assumptions.omega_m0, assumptions.omega_r0, assumptions.omega_k0, omega_lambda_flat
    )
    .expect("write");
    writeln!(
        json,
        "  \"h0_references\": {{\"planck_km_s_mpc\": {:.3}, \"distance_ladder_km_s_mpc\": {:.3}}},",
        PLANCK_H0_KM_S_MPC, DISTANCE_LADDER_H0_KM_S_MPC
    )
    .expect("write");
    writeln!(json, "  \"variants\": [").expect("write");
    for (idx, (name, lambda, h0_flat, omega_lambda_planck, omega_lambda_ladder, omega_k_planck, omega_k_ladder)) in
        variants.iter().enumerate()
    {
        writeln!(
            json,
            "    {{\"name\":\"{}\",\"lambda\":{:.12e},\"h0_flat_km_s_mpc\":{:.9},\"omega_lambda_at_planck_h0\":{:.9},\"omega_lambda_at_ladder_h0\":{:.9},\"omega_k_at_planck_h0\":{:.9},\"omega_k_at_ladder_h0\":{:.9}}}{}",
            name,
            lambda,
            h0_flat.unwrap_or(f64::NAN),
            omega_lambda_planck,
            omega_lambda_ladder,
            omega_k_planck,
            omega_k_ladder,
            if idx + 1 == variants.len() { "" } else { "," }
        )
        .expect("write");
    }
    writeln!(json, "  ],").expect("write");
    writeln!(json, "  \"full_candidate_curve\": [").expect("write");
    for (idx, (z, h, e)) in curve_rows.iter().enumerate() {
        writeln!(
            json,
            "    {{\"z\":{:.6},\"h_km_s_mpc\":{:.9},\"e\":{:.9}}}{}",
            z,
            h,
            e,
            if idx + 1 == curve_rows.len() { "" } else { "," }
        )
        .expect("write");
    }
    writeln!(json, "  ]\n}}").expect("write");

    println!("wrote {txt_path}");
    println!("wrote {json_path}");
    println!(
        "FRW harness: Ω_m0={:.4}, Ω_r0={:.6}, Ω_Λ(flat)={:.6}, H0(full,flat)={:.4} km/s/Mpc",
        assumptions.omega_m0, assumptions.omega_r0, omega_lambda_flat, h0_full_flat
    );
}
