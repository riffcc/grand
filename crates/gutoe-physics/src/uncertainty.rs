/*!
 * GUTOE Physics - Upstream Uncertainty Propagation
 * Copyright (C) 2026  Riff Labs
 *
 * GRAND-353:
 *   Propagate structured uncertainty across upstream cosmology lanes
 *   (inflation, baryogenesis, BBN, dark sector, transfer, microphysics)
 *   with quantitative distributions, not only pass/fail gates.
 */

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::f64::consts::PI;

use crate::baryogenesis::{evaluate_baryogenesis_gate, BaryogenesisWindows, ETA_B_OBSERVED};
use crate::bbn::{
    primordial_deuterium_ratio, primordial_helium3_ratio, primordial_helium4_mass_fraction,
    primordial_lithium7_ratio, BbnWindows, DH_TARGET, HE3H_TARGET, LI7H_OBSERVED, YP_TARGET,
};
use crate::constants::{lambda_cosmological_full_candidate, C, DARK_TO_VISIBLE_GEOMETRIC_RATIO};
use crate::cosmo_transfer::{evaluate_transfer_gate, TransferAssumptions, TransferWindows};
use crate::dark_matter_falsification::{
    evaluate_dark_matter_gate, DarkMatterFalsificationWindows, DM_FRACTION_OBS, OMEGA_BARYON_OBS,
};
use crate::dark_sector::DarkSectorBranch;
use crate::inflation::{evaluate_inflation_gate, InflationWindows};
use crate::microphysics::{
    evaluate_microphysics_gate, MicrophysicsAssumptions, MicrophysicsWindows,
};
use crate::universe::{UniverseAssumptions, UniverseWindows, SEC_PER_GYR, SEC_PER_YEAR};

const METER_PER_MPC: f64 = 3.085_677_581_491_367e22;
const Z_INTEGRAL_MAX: f64 = 1.0e12;
const T_CMB0_K: f64 = 2.7255;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UncertaintyAssumptions {
    pub samples: usize,
    pub seed: u64,
    pub sigma_scale: f64,
    pub omega_b_rel_sigma: f64,
    pub omega_r_rel_sigma: f64,
}

impl Default for UncertaintyAssumptions {
    fn default() -> Self {
        Self {
            samples: 768,
            seed: 0xA11CE_5EED,
            sigma_scale: 1.0,
            omega_b_rel_sigma: 0.01,
            omega_r_rel_sigma: 0.05,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UncertaintyWindows {
    pub pass_fraction_min: f64,
    pub h0_p95_rel_span_max: f64,
    pub theta_star_p95_rel_span_max: f64,
    pub yp_network_p95_span_max: f64,
}

impl Default for UncertaintyWindows {
    fn default() -> Self {
        Self {
            pass_fraction_min: 0.80,
            h0_p95_rel_span_max: 0.12,
            theta_star_p95_rel_span_max: 0.16,
            yp_network_p95_span_max: 0.03,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DistributionSummary {
    pub mean: f64,
    pub std: f64,
    pub p05: f64,
    pub p50: f64,
    pub p95: f64,
    pub min: f64,
    pub max: f64,
}

impl DistributionSummary {
    fn nan() -> Self {
        Self {
            mean: f64::NAN,
            std: f64::NAN,
            p05: f64::NAN,
            p50: f64::NAN,
            p95: f64::NAN,
            min: f64::NAN,
            max: f64::NAN,
        }
    }

    pub fn rel_span95(self) -> f64 {
        if !self.p50.is_finite() || self.p50.abs() <= f64::EPSILON {
            return f64::NAN;
        }
        (self.p95 - self.p05).abs() / self.p50.abs()
    }

    pub fn abs_span95(self) -> f64 {
        (self.p95 - self.p05).abs()
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct UncertaintyScorecard {
    pub requested_samples: usize,
    pub valid_samples: usize,
    pub pass_fraction: f64,
    pub inflation_pass_fraction: f64,
    pub baryogenesis_pass_fraction: f64,
    pub bbn_pass_fraction: f64,
    pub dark_pass_fraction: f64,
    pub transfer_pass_fraction: f64,
    pub microphysics_pass_fraction: f64,
    pub background_pass_fraction: f64,
    pub n_s: DistributionSummary,
    pub a_s: DistributionSummary,
    pub eta10: DistributionSummary,
    pub dm_fraction: DistributionSummary,
    pub h0_km_s_mpc: DistributionSummary,
    pub age_gyr: DistributionSummary,
    pub rs_drag_mpc: DistributionSummary,
    pub theta_star_rad: DistributionSummary,
    pub l_peak1: DistributionSummary,
    pub l_peak2: DistributionSummary,
    pub yp_network: DistributionSummary,
    pub dh_network: DistributionSummary,
    pub z_visibility_peak: DistributionSummary,
}

#[derive(Debug, Clone, PartialEq)]
pub struct UncertaintyGateScorecard {
    pub summary: UncertaintyScorecard,
    pub pass_fraction_ok: bool,
    pub h0_span_ok: bool,
    pub theta_star_span_ok: bool,
    pub yp_span_ok: bool,
}

impl UncertaintyGateScorecard {
    pub const fn passes_all(&self) -> bool {
        self.pass_fraction_ok && self.h0_span_ok && self.theta_star_span_ok && self.yp_span_ok
    }
}

fn gaussian_sample(rng: &mut StdRng) -> f64 {
    let u1 = (1.0_f64 - rng.gen::<f64>()).clamp(1e-12, 1.0);
    let u2 = rng.gen::<f64>();
    (-2.0_f64 * u1.ln()).sqrt() * (2.0_f64 * PI * u2).cos()
}

fn summarize(xs: &mut [f64]) -> DistributionSummary {
    if xs.is_empty() {
        return DistributionSummary::nan();
    }
    let n = xs.len() as f64;
    let mean = xs.iter().copied().sum::<f64>() / n;
    let var = xs
        .iter()
        .map(|x| {
            let d = *x - mean;
            d * d
        })
        .sum::<f64>()
        / n;
    xs.sort_by(|a, b| a.total_cmp(b));

    let q = |p: f64| -> f64 {
        let idx = ((xs.len() - 1) as f64 * p).round() as usize;
        xs[idx.min(xs.len() - 1)]
    };

    DistributionSummary {
        mean,
        std: var.sqrt(),
        p05: q(0.05),
        p50: q(0.50),
        p95: q(0.95),
        min: xs[0],
        max: xs[xs.len() - 1],
    }
}

fn s_inv_to_km_s_mpc(h0_s_inv: f64) -> f64 {
    h0_s_inv * METER_PER_MPC / 1_000.0
}

fn h0_from_lambda_and_omega_lambda(lambda: f64, omega_lambda: f64) -> Option<f64> {
    if lambda <= 0.0 || omega_lambda <= 0.0 {
        return None;
    }
    let h0_s_inv = C * (lambda / (3.0 * omega_lambda)).sqrt();
    Some(s_inv_to_km_s_mpc(h0_s_inv))
}

fn km_s_mpc_to_s_inv(h0_km_s_mpc: f64) -> f64 {
    (h0_km_s_mpc * 1_000.0) / METER_PER_MPC
}

fn e2_of_z(z: f64, omega_r0: f64, omega_m0: f64, omega_k0: f64, omega_lambda0: f64) -> f64 {
    let x = 1.0 + z;
    omega_r0 * x.powi(4) + omega_m0 * x.powi(3) + omega_k0 * x.powi(2) + omega_lambda0
}

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
    let n = 4_096usize;
    let ln0 = (1.0 + z0).ln();
    let ln1 = (1.0 + z1).ln();
    let dln = (ln1 - ln0) / n as f64;

    let mut acc = 0.0;
    for i in 0..n {
        let ln_mid = ln0 + (i as f64 + 0.5) * dln;
        let x = ln_mid.exp();
        let z = x - 1.0;
        let e2 = e2_of_z(z, omega_r0, omega_m0, omega_k0, omega_lambda0);
        if e2 <= 0.0 {
            continue;
        }
        acc += dln / (h0_s_inv * e2.sqrt());
    }
    acc
}

fn age_of_universe_seconds(
    h0_km_s_mpc: f64,
    omega_r0: f64,
    omega_m0: f64,
    omega_k0: f64,
    omega_lambda0: f64,
) -> f64 {
    integrate_time_from_z0_to_z1(
        h0_km_s_mpc,
        0.0,
        Z_INTEGRAL_MAX,
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
) -> f64 {
    integrate_time_from_z0_to_z1(
        h0_km_s_mpc,
        z.max(0.0),
        Z_INTEGRAL_MAX,
        omega_r0,
        omega_m0,
        omega_k0,
        omega_lambda0,
    )
}

pub fn evaluate_uncertainty(
    universe_assumptions: UniverseAssumptions,
    universe_windows: UniverseWindows,
    uncertainty_assumptions: UncertaintyAssumptions,
) -> UncertaintyScorecard {
    let inflation = evaluate_inflation_gate(InflationWindows::default());
    let inflation_windows = InflationWindows::default();
    let baryogenesis_windows = BaryogenesisWindows::default();
    let bbn_windows = BbnWindows::default();
    let transfer_windows = TransferWindows::default();
    let microphysics_windows = MicrophysicsWindows::default();
    let dark_windows = DarkMatterFalsificationWindows::default();

    let dark_unified = evaluate_dark_matter_gate(dark_windows)
        .into_iter()
        .find(|s| s.branch == DarkSectorBranch::Unified)
        .expect("unified dark scorecard should exist");

    let base_baryo = evaluate_baryogenesis_gate(baryogenesis_windows);
    let base_eta10 = base_baryo.eta_predicted * 1.0e10;

    let lambda = lambda_cosmological_full_candidate();

    let ns_sigma = (inflation_windows.ns_tol / 3.0) * uncertainty_assumptions.sigma_scale;
    let as_rel_sigma = (inflation_windows.as_tol / inflation_windows.as_center / 3.0)
        * uncertainty_assumptions.sigma_scale;
    let eta_rel_sigma =
        (baryogenesis_windows.eta_rel_error_max / 3.0) * uncertainty_assumptions.sigma_scale;
    let dm_ratio_rel_sigma = (dark_windows.dm_fraction_delta_abs_max / DM_FRACTION_OBS / 3.0)
        * uncertainty_assumptions.sigma_scale;

    let mut rng = StdRng::seed_from_u64(uncertainty_assumptions.seed);

    let mut n_s_samples = Vec::with_capacity(uncertainty_assumptions.samples);
    let mut a_s_samples = Vec::with_capacity(uncertainty_assumptions.samples);
    let mut eta10_samples = Vec::with_capacity(uncertainty_assumptions.samples);
    let mut dm_frac_samples = Vec::with_capacity(uncertainty_assumptions.samples);
    let mut h0_samples = Vec::with_capacity(uncertainty_assumptions.samples);
    let mut age_samples = Vec::with_capacity(uncertainty_assumptions.samples);
    let mut rs_samples = Vec::with_capacity(uncertainty_assumptions.samples);
    let mut theta_samples = Vec::with_capacity(uncertainty_assumptions.samples);
    let mut l1_samples = Vec::with_capacity(uncertainty_assumptions.samples);
    let mut l2_samples = Vec::with_capacity(uncertainty_assumptions.samples);
    let mut yp_samples = Vec::with_capacity(uncertainty_assumptions.samples);
    let mut dh_samples = Vec::with_capacity(uncertainty_assumptions.samples);
    let mut zvis_samples = Vec::with_capacity(uncertainty_assumptions.samples);

    let mut inflation_pass_count = 0usize;
    let mut baryogenesis_pass_count = 0usize;
    let mut bbn_pass_count = 0usize;
    let mut dark_pass_count = 0usize;
    let mut transfer_pass_count = 0usize;
    let mut microphysics_pass_count = 0usize;
    let mut background_pass_count = 0usize;
    let mut all_pass_count = 0usize;

    for _ in 0..uncertainty_assumptions.samples {
        let ns = (inflation.n_s + gaussian_sample(&mut rng) * ns_sigma).clamp(0.85, 1.10);
        let as_rel = 1.0 + gaussian_sample(&mut rng) * as_rel_sigma;
        let a_s = (inflation.a_s * as_rel).max(1e-12);
        let eta10 = (base_eta10 * (1.0 + gaussian_sample(&mut rng) * eta_rel_sigma)).max(1e-8);

        let omega_b0 = (OMEGA_BARYON_OBS
            * (1.0 + gaussian_sample(&mut rng) * uncertainty_assumptions.omega_b_rel_sigma))
            .max(1e-8);
        let omega_dm_ratio = (DARK_TO_VISIBLE_GEOMETRIC_RATIO
            * (1.0 + gaussian_sample(&mut rng) * dm_ratio_rel_sigma))
            .max(1e-8);
        let omega_dm0 = omega_b0 * omega_dm_ratio;
        let omega_m0 = omega_b0 + omega_dm0;
        let omega_r0 = (universe_assumptions.omega_r0
            * (1.0 + gaussian_sample(&mut rng) * uncertainty_assumptions.omega_r_rel_sigma))
            .max(1e-10);
        let omega_k0 = universe_assumptions.omega_k0;
        let omega_lambda0 = 1.0 - omega_m0 - omega_r0 - omega_k0;
        if !(omega_lambda0 > 0.0) {
            continue;
        }

        let Some(h0_km_s_mpc) = h0_from_lambda_and_omega_lambda(lambda, omega_lambda0) else {
            continue;
        };

        let transfer = evaluate_transfer_gate(
            TransferAssumptions {
                h0_km_s_mpc,
                omega_b0,
                omega_m0,
                omega_r0,
                omega_k0,
                omega_lambda0,
                n_s: ns,
                a_s,
            },
            transfer_windows,
        );

        let microphysics = evaluate_microphysics_gate(
            MicrophysicsAssumptions {
                h0_km_s_mpc,
                omega_b0,
                omega_m0,
                omega_r0,
                omega_k0,
                omega_lambda0,
                eta10,
            },
            microphysics_windows,
        );

        let age_gyr =
            age_of_universe_seconds(h0_km_s_mpc, omega_r0, omega_m0, omega_k0, omega_lambda0)
                / SEC_PER_GYR;
        let z_recomb = 1089.0;
        let z_bbn = 1.0e9 / T_CMB0_K - 1.0;
        let recomb_age_kyr = age_at_redshift_seconds(
            h0_km_s_mpc,
            z_recomb,
            omega_r0,
            omega_m0,
            omega_k0,
            omega_lambda0,
        ) / (1_000.0 * SEC_PER_YEAR);
        let bbn_age_seconds = age_at_redshift_seconds(
            h0_km_s_mpc,
            z_bbn,
            omega_r0,
            omega_m0,
            omega_k0,
            omega_lambda0,
        );

        // Upstream lane checks under sampled parameters.
        let inflation_ok = (ns - inflation_windows.ns_center).abs() <= inflation_windows.ns_tol
            && (a_s - inflation_windows.as_center).abs() <= inflation_windows.as_tol
            && inflation.r <= inflation_windows.r_max;

        let eta_b = eta10 * 1.0e-10;
        let baryo_rel_error = (eta_b - ETA_B_OBSERVED).abs() / ETA_B_OBSERVED;
        let baryogenesis_ok = baryo_rel_error <= baryogenesis_windows.eta_rel_error_max;

        let yp_pred = primordial_helium4_mass_fraction(eta10);
        let dh_pred = primordial_deuterium_ratio(eta10);
        let he3_pred = primordial_helium3_ratio(eta10);
        let li7_pred = primordial_lithium7_ratio(eta10);

        let bbn_ok = (yp_pred - YP_TARGET).abs() <= bbn_windows.yp_abs_max
            && ((dh_pred - DH_TARGET).abs() / DH_TARGET) <= bbn_windows.dh_rel_max
            && ((he3_pred - HE3H_TARGET).abs() / HE3H_TARGET) <= bbn_windows.he3_rel_max
            && {
                let li_ratio = li7_pred / LI7H_OBSERVED;
                li_ratio >= bbn_windows.li_tension_ratio_min
                    && li_ratio <= bbn_windows.li_tension_ratio_max
            };

        let dm_fraction = omega_dm0 / omega_m0;
        let dark_fraction_ok =
            (dm_fraction - DM_FRACTION_OBS).abs() <= dark_windows.dm_fraction_delta_abs_max;
        let dark_ok = dark_unified.rotation_ok && dark_unified.lensing_ok && dark_fraction_ok;

        let transfer_ok = transfer.passes_all();
        let microphysics_ok = microphysics.passes_all();

        let h0_rel_error = ((h0_km_s_mpc - universe_assumptions.h0_ref_km_s_mpc)
            / universe_assumptions.h0_ref_km_s_mpc)
            .abs();
        let background_ok = h0_rel_error <= universe_windows.h0_rel_error_max
            && age_gyr >= universe_windows.age_gyr_min
            && age_gyr <= universe_windows.age_gyr_max
            && recomb_age_kyr >= universe_windows.recombination_age_kyr_min
            && recomb_age_kyr <= universe_windows.recombination_age_kyr_max
            && bbn_age_seconds >= universe_windows.bbn_age_sec_min
            && bbn_age_seconds <= universe_windows.bbn_age_sec_max;

        let all_ok = inflation_ok
            && baryogenesis_ok
            && bbn_ok
            && dark_ok
            && transfer_ok
            && microphysics_ok
            && background_ok;

        inflation_pass_count += inflation_ok as usize;
        baryogenesis_pass_count += baryogenesis_ok as usize;
        bbn_pass_count += bbn_ok as usize;
        dark_pass_count += dark_ok as usize;
        transfer_pass_count += transfer_ok as usize;
        microphysics_pass_count += microphysics_ok as usize;
        background_pass_count += background_ok as usize;
        all_pass_count += all_ok as usize;

        n_s_samples.push(ns);
        a_s_samples.push(a_s);
        eta10_samples.push(eta10);
        dm_frac_samples.push(dm_fraction);
        h0_samples.push(h0_km_s_mpc);
        age_samples.push(age_gyr);
        rs_samples.push(transfer.rs_drag_mpc);
        theta_samples.push(transfer.theta_star_rad);
        l1_samples.push(transfer.l_peak1);
        l2_samples.push(transfer.l_peak2);
        yp_samples.push(microphysics.yp_network);
        dh_samples.push(microphysics.dh_network);
        zvis_samples.push(microphysics.z_visibility_peak);
    }

    let valid = n_s_samples.len();
    let denom = valid.max(1) as f64;

    UncertaintyScorecard {
        requested_samples: uncertainty_assumptions.samples,
        valid_samples: valid,
        pass_fraction: all_pass_count as f64 / denom,
        inflation_pass_fraction: inflation_pass_count as f64 / denom,
        baryogenesis_pass_fraction: baryogenesis_pass_count as f64 / denom,
        bbn_pass_fraction: bbn_pass_count as f64 / denom,
        dark_pass_fraction: dark_pass_count as f64 / denom,
        transfer_pass_fraction: transfer_pass_count as f64 / denom,
        microphysics_pass_fraction: microphysics_pass_count as f64 / denom,
        background_pass_fraction: background_pass_count as f64 / denom,
        n_s: summarize(&mut n_s_samples),
        a_s: summarize(&mut a_s_samples),
        eta10: summarize(&mut eta10_samples),
        dm_fraction: summarize(&mut dm_frac_samples),
        h0_km_s_mpc: summarize(&mut h0_samples),
        age_gyr: summarize(&mut age_samples),
        rs_drag_mpc: summarize(&mut rs_samples),
        theta_star_rad: summarize(&mut theta_samples),
        l_peak1: summarize(&mut l1_samples),
        l_peak2: summarize(&mut l2_samples),
        yp_network: summarize(&mut yp_samples),
        dh_network: summarize(&mut dh_samples),
        z_visibility_peak: summarize(&mut zvis_samples),
    }
}

pub fn evaluate_uncertainty_gate(
    universe_assumptions: UniverseAssumptions,
    universe_windows: UniverseWindows,
    uncertainty_assumptions: UncertaintyAssumptions,
    uncertainty_windows: UncertaintyWindows,
) -> UncertaintyGateScorecard {
    let summary = evaluate_uncertainty(
        universe_assumptions,
        universe_windows,
        uncertainty_assumptions,
    );

    let pass_fraction_ok = summary.pass_fraction >= uncertainty_windows.pass_fraction_min;
    let h0_span_ok = summary.h0_km_s_mpc.rel_span95() <= uncertainty_windows.h0_p95_rel_span_max;
    let theta_star_span_ok =
        summary.theta_star_rad.rel_span95() <= uncertainty_windows.theta_star_p95_rel_span_max;
    let yp_span_ok = summary.yp_network.abs_span95() <= uncertainty_windows.yp_network_p95_span_max;

    UncertaintyGateScorecard {
        summary,
        pass_fraction_ok,
        h0_span_ok,
        theta_star_span_ok,
        yp_span_ok,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uncertainty_lane_produces_valid_distributions() {
        let s = evaluate_uncertainty(
            UniverseAssumptions::default(),
            UniverseWindows::default(),
            UncertaintyAssumptions {
                samples: 128,
                ..UncertaintyAssumptions::default()
            },
        );
        assert!(s.valid_samples > 0);
        assert!(s.h0_km_s_mpc.mean.is_finite());
        assert!(s.theta_star_rad.p95 >= s.theta_star_rad.p05);
        assert!(s.pass_fraction >= 0.0 && s.pass_fraction <= 1.0);
    }

    #[test]
    fn uncertainty_gate_passes_default() {
        let g = evaluate_uncertainty_gate(
            UniverseAssumptions::default(),
            UniverseWindows::default(),
            UncertaintyAssumptions {
                samples: 128,
                ..UncertaintyAssumptions::default()
            },
            UncertaintyWindows::default(),
        );
        assert!(g.passes_all(), "uncertainty gate failed: {g:#?}");
    }
}
