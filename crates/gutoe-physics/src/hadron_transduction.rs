/*!
 * GUTOE Physics - Hadron Transduction Lane (GRAND-353)
 *
 * Purpose:
 * - Map structural hadron proxies to physical MeV scales with explicit anchors.
 * - Quantify uncertainty bands for p, n, π, K predictions.
 * - Expose deterministic pass/fail thresholds for CI.
 *
 * Scope note:
 * This lane is a reduced-order transduction model. It is a reporting/verification
 * scaffold, not a full lattice-QCD bound-state solver.
 */

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::f64::consts::PI;

use crate::chiral_symmetry_breaking::pion_mass_proxy;
use crate::constants::{
    ALPHA_LEADING_ORDER, CLIFFORD_STATE_COUNT_STRUCTURAL, DARK_STATE_COUNT_STRUCTURAL,
    DARK_TO_VISIBLE_GEOMETRIC_RATIO, VISIBLE_STATE_COUNT_STRUCTURAL,
};
use crate::dynamics_map::StandardModelDynamicsMap;

/// PDG-like anchors and running-scale inputs used by this reduced lane.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HadronReferenceAnchors {
    pub electron_mass_mev: f64,
    pub proton_mass_mev_obs: f64,
    pub neutron_mass_mev_obs: f64,
    pub pion_mass_mev_obs: f64,
    pub kaon_mass_mev_obs: f64,
    pub alpha_s_mz: f64,
    pub q_ref_gev: f64,
    pub mb_gev: f64,
    pub mc_gev: f64,
}

impl Default for HadronReferenceAnchors {
    fn default() -> Self {
        Self {
            electron_mass_mev: 0.510_998_950,
            proton_mass_mev_obs: 938.272_088_16,
            neutron_mass_mev_obs: 939.565_420_52,
            pion_mass_mev_obs: 139.570_39,
            kaon_mass_mev_obs: 493.677,
            alpha_s_mz: 0.1181,
            q_ref_gev: 91.1876,
            mb_gev: 4.18,
            mc_gev: 1.27,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HadronUncertaintyAssumptions {
    pub samples: usize,
    pub seed: u64,
    pub electron_mass_sigma_mev: f64,
    pub alpha_s_mz_sigma: f64,
    pub q_ref_gev_sigma: f64,
    pub mb_gev_sigma: f64,
    pub mc_gev_sigma: f64,
    /// Shared transduction systematic applied to chiral -> hadron outputs.
    pub transduction_rel_sigma: f64,
}

impl Default for HadronUncertaintyAssumptions {
    fn default() -> Self {
        Self {
            samples: 4096,
            seed: 0x353_2026_ABCD,
            electron_mass_sigma_mev: 1.6e-10,
            alpha_s_mz_sigma: 8.0e-4,
            q_ref_gev_sigma: 2.1e-3,
            mb_gev_sigma: 3.0e-2,
            mc_gev_sigma: 2.0e-2,
            transduction_rel_sigma: 1.5e-2,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HadronStructuralFactors {
    /// mp/me = 12 * T(17) = 1836.
    pub mp_me_structural_ratio: f64,
    /// π proxy from GRAND-126.
    pub pion_proxy: f64,
    /// Structural chiral -> hadron scaling factor:
    /// (mp/me) / (dim Cl(1,3) * |grade2|) = 1836/(16*6) = 153/8.
    pub pion_transduction_factor: f64,
    /// Structural damping from the finite visible-sector occupancy:
    /// (11 - 1) / 16 = 5/8.
    pub qcd_visibility_damping_factor: f64,
    /// Δnp factor from pion lane:
    /// α * (mZ/mW)^2 = (1/137) * (13/10) = 13/1370.
    pub delta_np_from_pion_factor: f64,
    /// Corrected dark/visible structural ratio = 115/22.
    pub corrected_dark_to_visible_ratio: f64,
    /// K/π structural scaling:
    /// 1 + (1/2) * corrected_ratio = 1 + 115/44 = 159/44.
    pub kaon_to_pion_factor: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HadronMassPrediction {
    pub qcd_scale_nf3_mev: f64,
    pub qcd_scale_effective_mev: f64,
    pub proton_mev: f64,
    pub neutron_mev: f64,
    pub pion_mev: f64,
    pub kaon_mev: f64,
    pub neutron_proton_split_mev: f64,
}

impl HadronMassPrediction {
    fn nan() -> Self {
        Self {
            qcd_scale_nf3_mev: f64::NAN,
            qcd_scale_effective_mev: f64::NAN,
            proton_mev: f64::NAN,
            neutron_mev: f64::NAN,
            pion_mev: f64::NAN,
            kaon_mev: f64::NAN,
            neutron_proton_split_mev: f64::NAN,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HadronMassResiduals {
    pub proton_rel_error: f64,
    pub neutron_rel_error: f64,
    pub pion_rel_error: f64,
    pub kaon_rel_error: f64,
}

impl HadronMassResiduals {
    fn from_prediction(p: HadronMassPrediction, a: HadronReferenceAnchors) -> Self {
        Self {
            proton_rel_error: (p.proton_mev - a.proton_mass_mev_obs) / a.proton_mass_mev_obs,
            neutron_rel_error: (p.neutron_mev - a.neutron_mass_mev_obs) / a.neutron_mass_mev_obs,
            pion_rel_error: (p.pion_mev - a.pion_mass_mev_obs) / a.pion_mass_mev_obs,
            kaon_rel_error: (p.kaon_mev - a.kaon_mass_mev_obs) / a.kaon_mass_mev_obs,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HadronDistributionSummary {
    pub mean: f64,
    pub std: f64,
    pub p05: f64,
    pub p50: f64,
    pub p95: f64,
    pub min: f64,
    pub max: f64,
}

impl HadronDistributionSummary {
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

    pub fn contains(self, x: f64) -> bool {
        self.p05 <= x && x <= self.p95
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct HadronUncertaintySummary {
    pub requested_samples: usize,
    pub valid_samples: usize,
    pub valid_fraction: f64,
    pub qcd_scale_nf3_mev: HadronDistributionSummary,
    pub proton_mev: HadronDistributionSummary,
    pub neutron_mev: HadronDistributionSummary,
    pub pion_mev: HadronDistributionSummary,
    pub kaon_mev: HadronDistributionSummary,
    pub neutron_proton_split_mev: HadronDistributionSummary,
}

#[derive(Debug, Clone, PartialEq)]
pub struct HadronTransductionScorecard {
    pub anchors: HadronReferenceAnchors,
    pub structural: HadronStructuralFactors,
    pub central: HadronMassPrediction,
    pub residuals: HadronMassResiduals,
    pub uncertainty: HadronUncertaintySummary,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HadronTransductionWindows {
    pub proton_rel_error_abs_max: f64,
    pub neutron_rel_error_abs_max: f64,
    pub pion_rel_error_abs_max: f64,
    pub kaon_rel_error_abs_max: f64,
    pub min_valid_fraction: f64,
    pub pion_rel_span95_max: f64,
    pub kaon_rel_span95_max: f64,
}

impl Default for HadronTransductionWindows {
    fn default() -> Self {
        Self {
            proton_rel_error_abs_max: 5.0e-3,
            neutron_rel_error_abs_max: 5.0e-3,
            pion_rel_error_abs_max: 8.0e-2,
            kaon_rel_error_abs_max: 8.0e-2,
            min_valid_fraction: 0.90,
            pion_rel_span95_max: 0.30,
            kaon_rel_span95_max: 0.30,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct HadronTransductionGateScorecard {
    pub score: HadronTransductionScorecard,
    pub proton_rel_error_ok: bool,
    pub neutron_rel_error_ok: bool,
    pub pion_rel_error_ok: bool,
    pub kaon_rel_error_ok: bool,
    pub valid_fraction_ok: bool,
    pub pion_span95_ok: bool,
    pub kaon_span95_ok: bool,
    pub proton_obs_in_p95: bool,
    pub neutron_obs_in_p95: bool,
    pub pion_obs_in_p95: bool,
    pub kaon_obs_in_p95: bool,
}

impl HadronTransductionGateScorecard {
    pub const fn passes_all(&self) -> bool {
        self.proton_rel_error_ok
            && self.neutron_rel_error_ok
            && self.pion_rel_error_ok
            && self.kaon_rel_error_ok
            && self.valid_fraction_ok
            && self.pion_span95_ok
            && self.kaon_span95_ok
            && self.neutron_obs_in_p95
            && self.pion_obs_in_p95
            && self.kaon_obs_in_p95
    }
}

fn triangular(n: u32) -> u32 {
    n * (n + 1) / 2
}

fn gaussian_sample(rng: &mut StdRng) -> f64 {
    let u1 = (1.0_f64 - rng.gen::<f64>()).clamp(1e-12, 1.0);
    let u2 = rng.gen::<f64>();
    (-2.0_f64 * u1.ln()).sqrt() * (2.0_f64 * PI * u2).cos()
}

fn summarize(xs: &mut [f64]) -> HadronDistributionSummary {
    if xs.is_empty() {
        return HadronDistributionSummary::nan();
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
    HadronDistributionSummary {
        mean,
        std: var.sqrt(),
        p05: q(0.05),
        p50: q(0.50),
        p95: q(0.95),
        min: xs[0],
        max: xs[xs.len() - 1],
    }
}

fn beta0_su3(nf: u32) -> f64 {
    11.0 - (2.0 / 3.0) * nf as f64
}

fn beta1_su3(nf: u32) -> f64 {
    102.0 - (38.0 / 3.0) * nf as f64
}

fn beta2_su3(nf: u32) -> f64 {
    (2857.0 / 2.0) - (5033.0 / 18.0) * nf as f64 + (325.0 / 54.0) * (nf as f64).powi(2)
}

fn alpha_s_three_loop(
    beta0: f64,
    beta1: f64,
    beta2: f64,
    q_gev: f64,
    lambda_qcd_gev: f64,
) -> Option<f64> {
    if !(q_gev > 0.0 && lambda_qcd_gev > 0.0 && q_gev > lambda_qcd_gev) {
        return None;
    }
    let l = ((q_gev * q_gev) / (lambda_qcd_gev * lambda_qcd_gev)).ln();
    if !(l.is_finite() && l > 0.0) {
        return None;
    }
    let ln_l = l.ln();
    let c1 = beta1 / (beta0 * beta0);
    let c2 = beta2 / (beta0 * beta0 * beta0);
    let bracket = 1.0 - c1 * ln_l / l + (c1 * c1 * (ln_l * ln_l - ln_l - 1.0) + c2) / (l * l);
    let a = (4.0 * PI / (beta0 * l)) * bracket;
    if a.is_finite() && a > 0.0 {
        Some(a)
    } else {
        None
    }
}

fn lambda_qcd_three_loop_from_anchor(
    beta0: f64,
    beta1: f64,
    beta2: f64,
    q_gev: f64,
    alpha_s_q: f64,
) -> Option<f64> {
    if !(q_gev > 0.0 && alpha_s_q > 0.0) {
        return None;
    }
    let mut lo = 1.0e-6_f64;
    let mut hi = q_gev * 0.999;
    for _ in 0..180 {
        let mid = 0.5 * (lo + hi);
        let Some(a_mid) = alpha_s_three_loop(beta0, beta1, beta2, q_gev, mid) else {
            hi = mid;
            continue;
        };
        if a_mid > alpha_s_q {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let lambda = 0.5 * (lo + hi);
    if lambda.is_finite() && lambda > 0.0 && lambda < q_gev {
        Some(lambda)
    } else {
        None
    }
}

fn match_lambda_three_loop_at_threshold(
    lambda_high: f64,
    nf_high: u32,
    nf_low: u32,
    threshold_gev: f64,
) -> Option<f64> {
    if !(lambda_high > 0.0 && threshold_gev > lambda_high) {
        return None;
    }
    let b0_hi = beta0_su3(nf_high);
    let b1_hi = beta1_su3(nf_high);
    let b2_hi = beta2_su3(nf_high);
    let b0_lo = beta0_su3(nf_low);
    let b1_lo = beta1_su3(nf_low);
    let b2_lo = beta2_su3(nf_low);
    let alpha_thr = alpha_s_three_loop(b0_hi, b1_hi, b2_hi, threshold_gev, lambda_high)?;
    lambda_qcd_three_loop_from_anchor(b0_lo, b1_lo, b2_lo, threshold_gev, alpha_thr)
}

pub fn qcd_scale_nf3_three_loop_mev(anchors: &HadronReferenceAnchors) -> Option<f64> {
    if !(anchors.alpha_s_mz > 0.0
        && anchors.q_ref_gev > 0.0
        && anchors.mb_gev > 0.0
        && anchors.mc_gev > 0.0
        && anchors.q_ref_gev > anchors.mb_gev
        && anchors.mb_gev > anchors.mc_gev)
    {
        return None;
    }
    let b0_nf5 = beta0_su3(5);
    let b1_nf5 = beta1_su3(5);
    let b2_nf5 = beta2_su3(5);
    let lambda_nf5 = lambda_qcd_three_loop_from_anchor(
        b0_nf5,
        b1_nf5,
        b2_nf5,
        anchors.q_ref_gev,
        anchors.alpha_s_mz,
    )?;
    let lambda_nf4 = match_lambda_three_loop_at_threshold(lambda_nf5, 5, 4, anchors.mb_gev)?;
    let lambda_nf3 = match_lambda_three_loop_at_threshold(lambda_nf4, 4, 3, anchors.mc_gev)?;
    Some(lambda_nf3 * 1000.0)
}

pub fn corrected_dark_to_visible_ratio_structural() -> f64 {
    DARK_TO_VISIBLE_GEOMETRIC_RATIO
        - (DARK_STATE_COUNT_STRUCTURAL / 2.0) / VISIBLE_STATE_COUNT_STRUCTURAL
}

pub fn hadron_structural_factors() -> HadronStructuralFactors {
    let sm = StandardModelDynamicsMap::from_clifford_z3();
    let mp_me = (sm.total_gauge_generators * triangular(sm.clifford_dim + 1)) as f64;
    let grade2_card = (sm.su2_generators + sm.magnetic_triplet_card) as f64;
    let pion_transduction_factor = mp_me / (sm.clifford_dim as f64 * grade2_card);
    let qcd_visibility_damping_factor =
        (VISIBLE_STATE_COUNT_STRUCTURAL - 1.0) / CLIFFORD_STATE_COUNT_STRUCTURAL;
    let corrected_ratio = corrected_dark_to_visible_ratio_structural();
    let kaon_to_pion_factor = 1.0 + 0.5 * corrected_ratio;
    HadronStructuralFactors {
        mp_me_structural_ratio: mp_me,
        pion_proxy: pion_mass_proxy(),
        pion_transduction_factor,
        qcd_visibility_damping_factor,
        delta_np_from_pion_factor: ALPHA_LEADING_ORDER * sm.mz_over_mw_sq,
        corrected_dark_to_visible_ratio: corrected_ratio,
        kaon_to_pion_factor,
    }
}

fn predict_hadron_masses(
    anchors: &HadronReferenceAnchors,
    transduction_scale_multiplier: f64,
) -> Option<HadronMassPrediction> {
    if !(transduction_scale_multiplier.is_finite() && transduction_scale_multiplier > 0.0) {
        return None;
    }
    let qcd_scale_nf3_mev = qcd_scale_nf3_three_loop_mev(anchors)?;
    let f = hadron_structural_factors();
    let qcd_scale_effective_mev = qcd_scale_nf3_mev * f.qcd_visibility_damping_factor;
    let proton_mev = anchors.electron_mass_mev * f.mp_me_structural_ratio;
    let pion_mev = f.pion_proxy
        * qcd_scale_effective_mev
        * f.pion_transduction_factor
        * transduction_scale_multiplier;
    let neutron_proton_split_mev = pion_mev * f.delta_np_from_pion_factor;
    let neutron_mev = proton_mev + neutron_proton_split_mev;
    let kaon_mev = pion_mev * f.kaon_to_pion_factor;
    if !(proton_mev.is_finite()
        && neutron_mev.is_finite()
        && pion_mev.is_finite()
        && kaon_mev.is_finite()
        && proton_mev > 0.0
        && neutron_mev > 0.0
        && pion_mev > 0.0
        && kaon_mev > 0.0)
    {
        return None;
    }
    Some(HadronMassPrediction {
        qcd_scale_nf3_mev,
        qcd_scale_effective_mev,
        proton_mev,
        neutron_mev,
        pion_mev,
        kaon_mev,
        neutron_proton_split_mev,
    })
}

pub fn evaluate_hadron_transduction(
    anchors: HadronReferenceAnchors,
    assumptions: HadronUncertaintyAssumptions,
) -> HadronTransductionScorecard {
    let structural = hadron_structural_factors();
    let central = predict_hadron_masses(&anchors, 1.0).unwrap_or_else(HadronMassPrediction::nan);
    let residuals = HadronMassResiduals::from_prediction(central, anchors);

    let mut rng = StdRng::seed_from_u64(assumptions.seed);
    let mut qcd_scale_samples = Vec::with_capacity(assumptions.samples);
    let mut proton_samples = Vec::with_capacity(assumptions.samples);
    let mut neutron_samples = Vec::with_capacity(assumptions.samples);
    let mut pion_samples = Vec::with_capacity(assumptions.samples);
    let mut kaon_samples = Vec::with_capacity(assumptions.samples);
    let mut split_samples = Vec::with_capacity(assumptions.samples);

    for _ in 0..assumptions.samples {
        let mut a = anchors;
        a.electron_mass_mev = (anchors.electron_mass_mev
            + assumptions.electron_mass_sigma_mev * gaussian_sample(&mut rng))
        .max(1.0e-6);
        a.alpha_s_mz = (anchors.alpha_s_mz
            + assumptions.alpha_s_mz_sigma * gaussian_sample(&mut rng))
        .clamp(1.0e-4, 0.5);
        a.q_ref_gev =
            (anchors.q_ref_gev + assumptions.q_ref_gev_sigma * gaussian_sample(&mut rng)).max(5.0);
        a.mb_gev = (anchors.mb_gev + assumptions.mb_gev_sigma * gaussian_sample(&mut rng))
            .clamp(0.5, a.q_ref_gev - 0.2);
        a.mc_gev = (anchors.mc_gev + assumptions.mc_gev_sigma * gaussian_sample(&mut rng))
            .clamp(0.2, a.mb_gev - 0.1);
        if !(a.q_ref_gev > a.mb_gev && a.mb_gev > a.mc_gev) {
            continue;
        }
        let transduction_scale_multiplier =
            1.0 + assumptions.transduction_rel_sigma * gaussian_sample(&mut rng);
        let Some(p) = predict_hadron_masses(&a, transduction_scale_multiplier.max(0.1)) else {
            continue;
        };
        qcd_scale_samples.push(p.qcd_scale_nf3_mev);
        proton_samples.push(p.proton_mev);
        neutron_samples.push(p.neutron_mev);
        pion_samples.push(p.pion_mev);
        kaon_samples.push(p.kaon_mev);
        split_samples.push(p.neutron_proton_split_mev);
    }

    let valid_samples = qcd_scale_samples.len();
    let valid_fraction = if assumptions.samples > 0 {
        valid_samples as f64 / assumptions.samples as f64
    } else {
        0.0
    };
    let uncertainty = HadronUncertaintySummary {
        requested_samples: assumptions.samples,
        valid_samples,
        valid_fraction,
        qcd_scale_nf3_mev: summarize(&mut qcd_scale_samples),
        proton_mev: summarize(&mut proton_samples),
        neutron_mev: summarize(&mut neutron_samples),
        pion_mev: summarize(&mut pion_samples),
        kaon_mev: summarize(&mut kaon_samples),
        neutron_proton_split_mev: summarize(&mut split_samples),
    };

    HadronTransductionScorecard {
        anchors,
        structural,
        central,
        residuals,
        uncertainty,
    }
}

pub fn evaluate_hadron_transduction_gate(
    anchors: HadronReferenceAnchors,
    assumptions: HadronUncertaintyAssumptions,
    windows: HadronTransductionWindows,
) -> HadronTransductionGateScorecard {
    let score = evaluate_hadron_transduction(anchors, assumptions);

    let proton_rel_error_ok =
        score.residuals.proton_rel_error.abs() <= windows.proton_rel_error_abs_max;
    let neutron_rel_error_ok =
        score.residuals.neutron_rel_error.abs() <= windows.neutron_rel_error_abs_max;
    let pion_rel_error_ok = score.residuals.pion_rel_error.abs() <= windows.pion_rel_error_abs_max;
    let kaon_rel_error_ok = score.residuals.kaon_rel_error.abs() <= windows.kaon_rel_error_abs_max;
    let valid_fraction_ok = score.uncertainty.valid_fraction >= windows.min_valid_fraction;
    let pion_span95_ok = score.uncertainty.pion_mev.rel_span95() <= windows.pion_rel_span95_max;
    let kaon_span95_ok = score.uncertainty.kaon_mev.rel_span95() <= windows.kaon_rel_span95_max;
    let proton_obs_in_p95 = score
        .uncertainty
        .proton_mev
        .contains(score.anchors.proton_mass_mev_obs);
    let neutron_obs_in_p95 = score
        .uncertainty
        .neutron_mev
        .contains(score.anchors.neutron_mass_mev_obs);
    let pion_obs_in_p95 = score
        .uncertainty
        .pion_mev
        .contains(score.anchors.pion_mass_mev_obs);
    let kaon_obs_in_p95 = score
        .uncertainty
        .kaon_mev
        .contains(score.anchors.kaon_mass_mev_obs);

    HadronTransductionGateScorecard {
        score,
        proton_rel_error_ok,
        neutron_rel_error_ok,
        pion_rel_error_ok,
        kaon_rel_error_ok,
        valid_fraction_ok,
        pion_span95_ok,
        kaon_span95_ok,
        proton_obs_in_p95,
        neutron_obs_in_p95,
        pion_obs_in_p95,
        kaon_obs_in_p95,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn structural_factors_match_closed_forms() {
        let f = hadron_structural_factors();
        assert!((f.mp_me_structural_ratio - 1836.0).abs() < 1.0e-12);
        assert!((f.pion_transduction_factor - 153.0 / 8.0).abs() < 1.0e-12);
        assert!((f.qcd_visibility_damping_factor - 5.0 / 8.0).abs() < 1.0e-12);
        assert!((f.delta_np_from_pion_factor - 13.0 / 1370.0).abs() < 1.0e-12);
        assert!((f.corrected_dark_to_visible_ratio - 115.0 / 22.0).abs() < 1.0e-12);
        assert!((f.kaon_to_pion_factor - 159.0 / 44.0).abs() < 1.0e-12);
    }

    #[test]
    fn central_predictions_are_finite_and_positive() {
        let s = evaluate_hadron_transduction(
            HadronReferenceAnchors::default(),
            HadronUncertaintyAssumptions::default(),
        );
        assert!(s.central.qcd_scale_nf3_mev.is_finite() && s.central.qcd_scale_nf3_mev > 0.0);
        assert!(
            s.central.qcd_scale_effective_mev.is_finite()
                && s.central.qcd_scale_effective_mev > 0.0
        );
        assert!(s.central.proton_mev.is_finite() && s.central.proton_mev > 0.0);
        assert!(s.central.neutron_mev.is_finite() && s.central.neutron_mev > 0.0);
        assert!(s.central.pion_mev.is_finite() && s.central.pion_mev > 0.0);
        assert!(s.central.kaon_mev.is_finite() && s.central.kaon_mev > 0.0);
        assert!(s.uncertainty.valid_fraction > 0.95);
    }

    #[test]
    fn hadron_gate_passes_defaults() {
        let g = evaluate_hadron_transduction_gate(
            HadronReferenceAnchors::default(),
            HadronUncertaintyAssumptions::default(),
            HadronTransductionWindows::default(),
        );
        assert!(g.passes_all(), "hadron transduction gate failed: {g:#?}");
    }
}
