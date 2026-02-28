/*!
 * GUTOE Physics - End-to-End Universe Assembly Harness
 * Copyright (C) 2026  Riff Labs
 *
 * GRAND-350:
 *   Assemble inflation + baryogenesis + BBN + dark-sector + FRW expansion
 *   into one executable simulation lane with explicit falsification gates.
 */

use crate::baryogenesis::{evaluate_baryogenesis_gate, BaryogenesisWindows};
use crate::bbn::{evaluate_bbn_gate, BbnWindows};
use crate::constants::{lambda_cosmological_full_candidate, C, DARK_TO_VISIBLE_GEOMETRIC_RATIO};
use crate::cosmo_transfer::{
    evaluate_transfer_gate, TransferAssumptions, TransferScorecard, TransferWindows,
};
use crate::dark_matter_falsification::{
    evaluate_dark_matter_gate, DarkMatterBranchScorecard, DarkMatterFalsificationWindows,
    OMEGA_BARYON_OBS,
};
use crate::dark_sector::DarkSectorBranch;
use crate::inflation::{evaluate_inflation_gate, InflationWindows};
use crate::microphysics::{
    evaluate_microphysics_gate, MicrophysicsAssumptions, MicrophysicsScorecard, MicrophysicsWindows,
};

/// Unit conversion.
pub const METER_PER_MPC: f64 = 3.085_677_581_491_367e22;
/// CMB temperature today.
const T_CMB0_K: f64 = 2.7255;
/// Seconds per Julian year.
pub const SEC_PER_YEAR: f64 = 31_557_600.0;
/// Seconds per gigayear.
pub const SEC_PER_GYR: f64 = 1.0e9 * SEC_PER_YEAR;
/// Upper redshift cutoff for finite FRW integrals in this harness.
pub const Z_INTEGRAL_MAX: f64 = 1.0e12;

/// Baseline FRW assumptions for the assembled lane.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UniverseAssumptions {
    pub omega_r0: f64,
    pub omega_k0: f64,
    pub h0_ref_km_s_mpc: f64,
}

impl Default for UniverseAssumptions {
    fn default() -> Self {
        Self {
            omega_r0: 9.0e-5,
            omega_k0: 0.0,
            h0_ref_km_s_mpc: 67.4,
        }
    }
}

/// Hard acceptance windows for the assembled universe lane.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UniverseWindows {
    pub h0_rel_error_max: f64,
    pub age_gyr_min: f64,
    pub age_gyr_max: f64,
    pub recombination_age_kyr_min: f64,
    pub recombination_age_kyr_max: f64,
    pub bbn_age_sec_min: f64,
    pub bbn_age_sec_max: f64,
}

impl Default for UniverseWindows {
    fn default() -> Self {
        Self {
            h0_rel_error_max: 0.03,
            age_gyr_min: 13.0,
            age_gyr_max: 14.5,
            recombination_age_kyr_min: 200.0,
            recombination_age_kyr_max: 500.0,
            bbn_age_sec_min: 10.0,
            bbn_age_sec_max: 2_000.0,
        }
    }
}

/// Numerical depth controls for end-to-end universe simulation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UniverseSimulationDepth {
    /// Number of log-redshift history samples between z=0 and `history_z_max`.
    pub history_points: usize,
    /// Maximum redshift represented in the exported history table.
    pub history_z_max: f64,
    /// Maximum redshift used in FRW time integrals (acts as t~0 cutoff).
    pub integral_z_max: f64,
}

impl Default for UniverseSimulationDepth {
    fn default() -> Self {
        Self {
            history_points: 256,
            history_z_max: 1.0e9,
            integral_z_max: Z_INTEGRAL_MAX,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct UniverseEpoch {
    pub name: &'static str,
    pub z: f64,
    pub age_seconds: f64,
    pub temperature_k: f64,
    pub h_km_s_mpc: f64,
    pub omega_r: f64,
    pub omega_m: f64,
    pub omega_lambda: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct UniverseHistoryRow {
    pub z: f64,
    pub age_seconds: f64,
    pub temperature_k: f64,
    pub h_km_s_mpc: f64,
    pub omega_r: f64,
    pub omega_m: f64,
    pub omega_lambda: f64,
}

#[derive(Debug, Clone)]
pub struct UniverseScorecard {
    pub lambda_full: f64,
    pub h0_km_s_mpc: f64,
    pub h0_rel_error: f64,
    pub omega_b0: f64,
    pub omega_dm0: f64,
    pub omega_m0: f64,
    pub omega_r0: f64,
    pub omega_k0: f64,
    pub omega_lambda0: f64,
    pub age_gyr: f64,
    pub recombination_age_kyr: f64,
    pub bbn_age_seconds: f64,
    pub inflation_ok: bool,
    pub baryogenesis_ok: bool,
    pub bbn_ok: bool,
    pub microphysics_ok: bool,
    pub dark_matter_unified_ok: bool,
    pub transfer_ok: bool,
    pub transfer: TransferScorecard,
    pub microphysics: MicrophysicsScorecard,
    pub h0_ok: bool,
    pub age_ok: bool,
    pub recombination_ok: bool,
    pub bbn_timing_ok: bool,
    pub epochs: Vec<UniverseEpoch>,
    pub history: Vec<UniverseHistoryRow>,
}

impl UniverseScorecard {
    pub const fn passes_early_universe(&self) -> bool {
        self.inflation_ok && self.baryogenesis_ok && self.bbn_ok && self.microphysics_ok
    }

    pub const fn passes_late_universe(&self) -> bool {
        self.dark_matter_unified_ok && self.transfer_ok && self.h0_ok && self.age_ok
    }

    pub const fn passes_all(&self) -> bool {
        self.passes_early_universe()
            && self.passes_late_universe()
            && self.recombination_ok
            && self.bbn_timing_ok
    }
}

fn km_s_mpc_to_s_inv(h0_km_s_mpc: f64) -> f64 {
    (h0_km_s_mpc * 1_000.0) / METER_PER_MPC
}

fn s_inv_to_km_s_mpc(h0_s_inv: f64) -> f64 {
    h0_s_inv * METER_PER_MPC / 1_000.0
}

/// `Ω_Λ = Λ c² / (3 H0²)` rearranged to solve for `H0`.
fn h0_from_lambda_and_omega_lambda(lambda: f64, omega_lambda: f64) -> Option<f64> {
    if lambda <= 0.0 || omega_lambda <= 0.0 {
        return None;
    }
    let h0_s_inv = C * (lambda / (3.0 * omega_lambda)).sqrt();
    Some(s_inv_to_km_s_mpc(h0_s_inv))
}

fn e2_of_z(z: f64, omega_r0: f64, omega_m0: f64, omega_k0: f64, omega_lambda0: f64) -> f64 {
    let x = 1.0 + z;
    omega_r0 * x.powi(4) + omega_m0 * x.powi(3) + omega_k0 * x.powi(2) + omega_lambda0
}

fn h_of_z(
    h0_km_s_mpc: f64,
    z: f64,
    omega_r0: f64,
    omega_m0: f64,
    omega_k0: f64,
    omega_lambda0: f64,
) -> f64 {
    let e2 = e2_of_z(z, omega_r0, omega_m0, omega_k0, omega_lambda0);
    if e2 <= 0.0 {
        return f64::NAN;
    }
    h0_km_s_mpc * e2.sqrt()
}

/// Midpoint integration over logarithmic redshift bins for
/// `∫ dz / ((1+z) H(z))`.
fn integrate_time_from_z0_to_z1(
    h0_km_s_mpc: f64,
    z0: f64,
    z1: f64,
    omega_r0: f64,
    omega_m0: f64,
    omega_k0: f64,
    omega_lambda0: f64,
) -> f64 {
    if !(z1 > z0 && z0 >= 0.0) {
        return 0.0;
    }

    let h0_s_inv = km_s_mpc_to_s_inv(h0_km_s_mpc);
    let n = 8_192usize;
    let ln0 = (1.0 + z0).ln();
    let ln1 = (1.0 + z1).ln();
    let dln = (ln1 - ln0) / n as f64;

    let mut acc = 0.0;
    for i in 0..n {
        let ln_mid = ln0 + (i as f64 + 0.5) * dln;
        let x = ln_mid.exp(); // x = 1+z
        let z = x - 1.0;
        let e2 = e2_of_z(z, omega_r0, omega_m0, omega_k0, omega_lambda0);
        if e2 <= 0.0 {
            continue;
        }
        let e = e2.sqrt();
        // dz = x dln; dt = dz / (H0 * x * E) = dln / (H0 * E)
        acc += dln / (h0_s_inv * e);
    }
    acc
}

fn age_of_universe_seconds(
    h0_km_s_mpc: f64,
    omega_r0: f64,
    omega_m0: f64,
    omega_k0: f64,
    omega_lambda0: f64,
    integral_z_max: f64,
) -> f64 {
    let z_cap = integral_z_max.max(1.0);
    integrate_time_from_z0_to_z1(
        h0_km_s_mpc,
        0.0,
        z_cap,
        omega_r0,
        omega_m0,
        omega_k0,
        omega_lambda0,
    )
}

fn age_at_redshift_seconds(
    h0_km_s_mpc: f64,
    z: f64,
    omega_r0: f64,
    omega_m0: f64,
    omega_k0: f64,
    omega_lambda0: f64,
    integral_z_max: f64,
) -> f64 {
    let z_cap = integral_z_max.max(z.max(0.0) + 1.0);
    integrate_time_from_z0_to_z1(
        h0_km_s_mpc,
        z.max(0.0),
        z_cap,
        omega_r0,
        omega_m0,
        omega_k0,
        omega_lambda0,
    )
}

fn omega_components_at_z(
    z: f64,
    omega_r0: f64,
    omega_m0: f64,
    omega_k0: f64,
    omega_lambda0: f64,
) -> (f64, f64, f64) {
    let x = 1.0 + z;
    let er = omega_r0 * x.powi(4);
    let em = omega_m0 * x.powi(3);
    let ek = omega_k0 * x.powi(2);
    let el = omega_lambda0;
    let etot = er + em + ek + el;
    if etot <= 0.0 {
        return (f64::NAN, f64::NAN, f64::NAN);
    }
    (er / etot, em / etot, el / etot)
}

fn build_epoch(
    name: &'static str,
    z: f64,
    h0_km_s_mpc: f64,
    omega_r0: f64,
    omega_m0: f64,
    omega_k0: f64,
    omega_lambda0: f64,
    integral_z_max: f64,
) -> UniverseEpoch {
    let age_seconds = age_at_redshift_seconds(
        h0_km_s_mpc,
        z,
        omega_r0,
        omega_m0,
        omega_k0,
        omega_lambda0,
        integral_z_max,
    );
    let (omega_r, omega_m, omega_lambda) =
        omega_components_at_z(z, omega_r0, omega_m0, omega_k0, omega_lambda0);
    UniverseEpoch {
        name,
        z,
        age_seconds,
        temperature_k: T_CMB0_K * (1.0 + z),
        h_km_s_mpc: h_of_z(h0_km_s_mpc, z, omega_r0, omega_m0, omega_k0, omega_lambda0),
        omega_r,
        omega_m,
        omega_lambda,
    }
}

fn simulate_history(
    n: usize,
    z_max: f64,
    h0_km_s_mpc: f64,
    omega_r0: f64,
    omega_m0: f64,
    omega_k0: f64,
    omega_lambda0: f64,
    integral_z_max: f64,
) -> Vec<UniverseHistoryRow> {
    if n == 0 || z_max <= 0.0 {
        return Vec::new();
    }
    let ln1pz_max = (1.0 + z_max).ln();
    let mut out = Vec::with_capacity(n + 1);
    for i in 0..=n {
        let t = i as f64 / n as f64;
        let ln1pz = t * ln1pz_max;
        let x = ln1pz.exp();
        let z = x - 1.0;
        let age_seconds = age_at_redshift_seconds(
            h0_km_s_mpc,
            z,
            omega_r0,
            omega_m0,
            omega_k0,
            omega_lambda0,
            integral_z_max,
        );
        let (omega_r, omega_m, omega_lambda) =
            omega_components_at_z(z, omega_r0, omega_m0, omega_k0, omega_lambda0);
        out.push(UniverseHistoryRow {
            z,
            age_seconds,
            temperature_k: T_CMB0_K * (1.0 + z),
            h_km_s_mpc: h_of_z(h0_km_s_mpc, z, omega_r0, omega_m0, omega_k0, omega_lambda0),
            omega_r,
            omega_m,
            omega_lambda,
        });
    }
    out
}

pub fn evaluate_universe_gate_with_depth(
    assumptions: UniverseAssumptions,
    windows: UniverseWindows,
    depth: UniverseSimulationDepth,
) -> UniverseScorecard {
    let inflation = evaluate_inflation_gate(InflationWindows::default());
    let baryogenesis = evaluate_baryogenesis_gate(BaryogenesisWindows::default());
    let bbn = evaluate_bbn_gate(BbnWindows::default());
    let dark_scorecards = evaluate_dark_matter_gate(DarkMatterFalsificationWindows::default());
    let unified_dark: DarkMatterBranchScorecard = dark_scorecards
        .into_iter()
        .find(|s| s.branch == DarkSectorBranch::Unified)
        .expect("unified dark matter scorecard should exist");

    let omega_b0 = OMEGA_BARYON_OBS;
    let omega_dm0 = OMEGA_BARYON_OBS * DARK_TO_VISIBLE_GEOMETRIC_RATIO;
    let omega_m0 = omega_b0 + omega_dm0;
    let omega_r0 = assumptions.omega_r0;
    let omega_k0 = assumptions.omega_k0;
    let omega_lambda0 = 1.0 - omega_m0 - omega_r0 - omega_k0;

    let lambda_full = lambda_cosmological_full_candidate();
    let h0_km_s_mpc =
        h0_from_lambda_and_omega_lambda(lambda_full, omega_lambda0).unwrap_or(f64::NAN);
    let h0_rel_error =
        ((h0_km_s_mpc - assumptions.h0_ref_km_s_mpc) / assumptions.h0_ref_km_s_mpc).abs();
    let transfer = evaluate_transfer_gate(
        TransferAssumptions {
            h0_km_s_mpc,
            omega_b0,
            omega_m0,
            omega_r0,
            omega_k0,
            omega_lambda0,
            n_s: inflation.n_s,
            a_s: inflation.a_s,
        },
        TransferWindows::default(),
    );
    let microphysics = evaluate_microphysics_gate(
        MicrophysicsAssumptions {
            h0_km_s_mpc,
            omega_b0,
            omega_m0,
            omega_r0,
            omega_k0,
            omega_lambda0,
            eta10: bbn.eta10,
        },
        MicrophysicsWindows::default(),
    );

    let age_seconds = age_of_universe_seconds(
        h0_km_s_mpc,
        omega_r0,
        omega_m0,
        omega_k0,
        omega_lambda0,
        depth.integral_z_max,
    );
    let age_gyr = age_seconds / SEC_PER_GYR;

    let z_bbn = 1.0e9 / T_CMB0_K - 1.0;
    let z_recomb = 1089.0;
    let z_eq_m_lambda = (omega_lambda0 / omega_m0).powf(1.0 / 3.0) - 1.0;

    let epochs = vec![
        build_epoch(
            "BaryogenesisProxy",
            1.0e12 / T_CMB0_K - 1.0,
            h0_km_s_mpc,
            omega_r0,
            omega_m0,
            omega_k0,
            omega_lambda0,
            depth.integral_z_max,
        ),
        build_epoch(
            "BBN",
            z_bbn,
            h0_km_s_mpc,
            omega_r0,
            omega_m0,
            omega_k0,
            omega_lambda0,
            depth.integral_z_max,
        ),
        build_epoch(
            "Recombination",
            z_recomb,
            h0_km_s_mpc,
            omega_r0,
            omega_m0,
            omega_k0,
            omega_lambda0,
            depth.integral_z_max,
        ),
        build_epoch(
            "MatterLambdaEquality",
            z_eq_m_lambda.max(0.0),
            h0_km_s_mpc,
            omega_r0,
            omega_m0,
            omega_k0,
            omega_lambda0,
            depth.integral_z_max,
        ),
        build_epoch(
            "Today",
            0.0,
            h0_km_s_mpc,
            omega_r0,
            omega_m0,
            omega_k0,
            omega_lambda0,
            depth.integral_z_max,
        ),
    ];
    let history = simulate_history(
        depth.history_points,
        depth.history_z_max,
        h0_km_s_mpc,
        omega_r0,
        omega_m0,
        omega_k0,
        omega_lambda0,
        depth.integral_z_max,
    );

    let bbn_age_seconds = epochs
        .iter()
        .find(|e| e.name == "BBN")
        .map(|e| e.age_seconds)
        .unwrap_or(f64::NAN);
    let recombination_age_kyr = epochs
        .iter()
        .find(|e| e.name == "Recombination")
        .map(|e| e.age_seconds / (1_000.0 * SEC_PER_YEAR))
        .unwrap_or(f64::NAN);

    let inflation_ok = inflation.passes_all();
    let baryogenesis_ok = baryogenesis.passes_all();
    let bbn_ok = bbn.passes_all();
    let microphysics_ok = microphysics.passes_all();
    let dark_matter_unified_ok = unified_dark.passes_all();
    let transfer_ok = transfer.passes_all();

    let h0_ok = h0_rel_error <= windows.h0_rel_error_max;
    let age_ok = age_gyr >= windows.age_gyr_min && age_gyr <= windows.age_gyr_max;
    let recombination_ok = recombination_age_kyr >= windows.recombination_age_kyr_min
        && recombination_age_kyr <= windows.recombination_age_kyr_max;
    let bbn_timing_ok =
        bbn_age_seconds >= windows.bbn_age_sec_min && bbn_age_seconds <= windows.bbn_age_sec_max;

    UniverseScorecard {
        lambda_full,
        h0_km_s_mpc,
        h0_rel_error,
        omega_b0,
        omega_dm0,
        omega_m0,
        omega_r0,
        omega_k0,
        omega_lambda0,
        age_gyr,
        recombination_age_kyr,
        bbn_age_seconds,
        inflation_ok,
        baryogenesis_ok,
        bbn_ok,
        microphysics_ok,
        dark_matter_unified_ok,
        transfer_ok,
        transfer,
        microphysics,
        h0_ok,
        age_ok,
        recombination_ok,
        bbn_timing_ok,
        epochs,
        history,
    }
}

pub fn evaluate_universe_gate(
    assumptions: UniverseAssumptions,
    windows: UniverseWindows,
) -> UniverseScorecard {
    evaluate_universe_gate_with_depth(assumptions, windows, UniverseSimulationDepth::default())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn assembled_universe_gate_passes() {
        let s = evaluate_universe_gate(UniverseAssumptions::default(), UniverseWindows::default());
        assert!(s.passes_all(), "universe gate failed: {s:#?}");
    }

    #[test]
    fn assembled_universe_outputs_are_physical() {
        let s = evaluate_universe_gate(UniverseAssumptions::default(), UniverseWindows::default());
        assert!(s.lambda_full > 0.0);
        assert!(s.h0_km_s_mpc.is_finite() && s.h0_km_s_mpc > 0.0);
        assert!(s.age_gyr.is_finite() && s.age_gyr > 0.0);
        assert!(s.omega_lambda0 > 0.0 && s.omega_lambda0 < 1.0);
        assert!(s.omega_m0 > 0.0 && s.omega_m0 < 1.0);
        assert!(s.transfer_ok);
        assert!(s.microphysics_ok);
        assert!(s.transfer.rs_drag_mpc.is_finite() && s.transfer.rs_drag_mpc > 0.0);
        assert!(s.microphysics.yp_network.is_finite() && s.microphysics.yp_network > 0.0);
        assert_eq!(s.epochs.len(), 5);
        assert!(s.history.len() >= 200);
    }
}
