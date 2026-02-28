/*!
 * GUTOE Physics - Structural Inflation Harness
 * Copyright (C) 2026  Riff Labs
 *
 * GRAND-347:
 *   Zero-free-parameter inflation lane from shared Cl(1,3) primitives.
 */

use crate::constants::{
    ALPHA_LEADING_ORDER, DARK_FRACTION_TOTAL_STATE_SPLIT, DARK_GEOMETRIC_AMPLIFICATION,
    DARK_STATE_COUNT_STRUCTURAL, DARK_TO_VISIBLE_GEOMETRIC_RATIO, LAMBDA_QG,
};

/// CMB scalar spectral index reference.
pub const NS_OBSERVED: f64 = 0.9649;
/// Conservative first-pass tolerance window for this structural lane.
pub const NS_TOL: f64 = 0.0100;
/// Current observational upper bound on tensor-to-scalar ratio.
pub const R_MAX_OBSERVED: f64 = 0.06;
/// Scalar amplitude reference at CMB pivot scale.
pub const AS_OBSERVED: f64 = 2.10e-9;
/// Scalar-amplitude tolerance (first-pass, structural lane).
pub const AS_TOL: f64 = 0.30e-9;

/// Effective relativistic degrees of freedom at reheating (SM baseline).
pub const G_REH: f64 = 106.75;
/// Reduced Planck mass in GeV.
pub const M_PL_REDUCED_GEV: f64 = 2.435e18;

/// Structural e-fold count from shared Clifford dark-sector counts:
/// N = (geometric amplification) × (dark-state count) = 12 × 5 = 60.
pub fn inflation_efolds_structural() -> f64 {
    DARK_GEOMETRIC_AMPLIFICATION * DARK_STATE_COUNT_STRUCTURAL
}

/// Structural Hubble-scale ratio `H_inf / M_pl` from shared Cl(1,3) primitives.
///
/// Composition:
/// - `α_LO^2` from fine-structure structural suppression
/// - `(60/11)` geometric dark budget ratio
/// - `(1 - λ_QG) = 11/12` lattice survival factor
/// - `(3/6) = 1/2` timelike-spacelike/grade-2 split
/// - `1/sqrt(486)` micro-mode dilution
pub fn inflation_hubble_ratio_structural() -> f64 {
    inflation_hubble_ratio_structural_with_correction(1.0)
}

/// Structural Hubble-scale ratio with an explicit multiplicative correction
/// factor `c_inf` (default lane uses `c_inf = 1`).
pub fn inflation_hubble_ratio_structural_with_correction(c_inf: f64) -> f64 {
    let geometric_budget = DARK_TO_VISIBLE_GEOMETRIC_RATIO;
    let survival = 1.0 - LAMBDA_QG;
    let signature_split = 3.0 / 6.0;
    let micro_dilution = 1.0 / 486.0_f64.sqrt();
    c_inf
        * ALPHA_LEADING_ORDER.powi(2)
        * geometric_budget
        * survival
        * signature_split
        * micro_dilution
}

/// Slow-roll epsilon for plateau-like structural lane.
pub fn slow_roll_epsilon(n_efolds: f64) -> f64 {
    if n_efolds <= 0.0 {
        return f64::NAN;
    }
    3.0 / (4.0 * n_efolds * n_efolds)
}

/// Slow-roll eta for plateau-like structural lane.
pub fn slow_roll_eta(n_efolds: f64) -> f64 {
    if n_efolds <= 0.0 {
        return f64::NAN;
    }
    -1.0 / n_efolds
}

/// Scalar spectral index from first-order slow-roll observables.
pub fn scalar_spectral_index(n_efolds: f64) -> f64 {
    let eps = slow_roll_epsilon(n_efolds);
    let eta = slow_roll_eta(n_efolds);
    1.0 - 6.0 * eps + 2.0 * eta
}

/// Tensor-to-scalar ratio from slow-roll epsilon.
pub fn tensor_to_scalar_ratio(n_efolds: f64) -> f64 {
    16.0 * slow_roll_epsilon(n_efolds)
}

/// Scalar amplitude at pivot from `(H/M_pl)^2 / (8 π² ε)`.
pub fn scalar_amplitude(n_efolds: f64, h_over_mpl: f64) -> f64 {
    let eps = slow_roll_epsilon(n_efolds);
    if !(eps > 0.0) || !(h_over_mpl > 0.0) {
        return f64::NAN;
    }
    h_over_mpl.powi(2) / (8.0 * std::f64::consts::PI.powi(2) * eps)
}

/// Structural reheating equation-of-state proxy from shared dark fraction.
pub fn reheating_w_structural() -> f64 {
    DARK_FRACTION_TOTAL_STATE_SPLIT
}

/// Structural reheating duration in e-folds from shared counts: `N_reh = N / 12 = 5`.
pub fn reheating_efolds_structural() -> f64 {
    inflation_efolds_structural() / DARK_GEOMETRIC_AMPLIFICATION
}

/// End-of-inflation energy density in reduced Planck units (`M_pl^4`).
pub fn rho_end_planck_units(h_over_mpl: f64) -> f64 {
    3.0 * h_over_mpl.powi(2)
}

/// Reheating energy density in reduced Planck units.
pub fn rho_reheat_planck_units(h_over_mpl: f64) -> f64 {
    let rho_end = rho_end_planck_units(h_over_mpl);
    let w = reheating_w_structural();
    let n_reh = reheating_efolds_structural();
    rho_end * (-3.0 * (1.0 + w) * n_reh).exp()
}

/// Reheating temperature in GeV from structural reheating map.
pub fn reheating_temperature_gev(h_over_mpl: f64) -> f64 {
    let rho_reh = rho_reheat_planck_units(h_over_mpl);
    if !(rho_reh > 0.0) {
        return f64::NAN;
    }
    let pref = 30.0 / (std::f64::consts::PI.powi(2) * G_REH);
    let t_over_mpl = (pref * rho_reh).powf(0.25);
    t_over_mpl * M_PL_REDUCED_GEV
}

/// Simple Gaussian CMB proxy χ² over (`n_s`, `A_s`) and one-sided `r`.
pub fn cmb_proxy_chi2(n_s: f64, a_s: f64, r: f64) -> f64 {
    let sigma_ns = 0.0040;
    let sigma_as = 0.15e-9;
    let sigma_r = 0.010;
    let z_ns = (n_s - NS_OBSERVED) / sigma_ns;
    let z_as = (a_s - AS_OBSERVED) / sigma_as;
    let z_r = if r <= R_MAX_OBSERVED {
        0.0
    } else {
        (r - R_MAX_OBSERVED) / sigma_r
    };
    z_ns * z_ns + z_as * z_as + z_r * z_r
}

/// Proxy log-likelihood from χ².
pub fn cmb_proxy_loglike(chi2: f64) -> f64 {
    -0.5 * chi2
}

/// Inflation end condition in this lane: ε >= 1.
/// For ε = 3/(4N²), this gives N_end = sqrt(3)/2.
pub fn structural_n_end() -> f64 {
    3.0_f64.sqrt() / 2.0
}

/// Number of e-folds translated to an expansion factor `exp(N)`.
pub fn inflation_expansion_factor(n_efolds: f64) -> f64 {
    n_efolds.exp()
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InflationWindows {
    pub n_efolds_min: f64,
    pub n_efolds_max: f64,
    pub ns_center: f64,
    pub ns_tol: f64,
    pub r_max: f64,
    pub as_center: f64,
    pub as_tol: f64,
    pub t_reheat_min_gev: f64,
    pub chi2_max: f64,
}

impl Default for InflationWindows {
    fn default() -> Self {
        Self {
            n_efolds_min: 50.0,
            n_efolds_max: 70.0,
            ns_center: NS_OBSERVED,
            ns_tol: NS_TOL,
            r_max: R_MAX_OBSERVED,
            as_center: AS_OBSERVED,
            as_tol: AS_TOL,
            // BBN-safe floor.
            t_reheat_min_gev: 1.0e-3,
            // Conservative first-pass CMB proxy gate.
            chi2_max: 9.0,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InflationScorecard {
    pub n_efolds: f64,
    pub epsilon: f64,
    pub eta: f64,
    pub n_s: f64,
    pub r: f64,
    pub a_s: f64,
    pub h_over_mpl: f64,
    pub n_reheat: f64,
    pub w_reheat: f64,
    pub t_reheat_gev: f64,
    pub cmb_proxy_chi2: f64,
    pub cmb_proxy_loglike: f64,
    pub expansion_factor: f64,
    pub n_end: f64,
    pub n_efolds_ok: bool,
    pub n_s_ok: bool,
    pub r_ok: bool,
    pub a_s_ok: bool,
    pub reheating_ok: bool,
    pub cmb_like_ok: bool,
    pub graceful_exit_ok: bool,
}

impl InflationScorecard {
    pub const fn passes_all(self) -> bool {
        self.n_efolds_ok
            && self.n_s_ok
            && self.r_ok
            && self.a_s_ok
            && self.reheating_ok
            && self.cmb_like_ok
            && self.graceful_exit_ok
    }
}

pub fn evaluate_inflation_gate(w: InflationWindows) -> InflationScorecard {
    let n_efolds = inflation_efolds_structural();
    let epsilon = slow_roll_epsilon(n_efolds);
    let eta = slow_roll_eta(n_efolds);
    let n_s = scalar_spectral_index(n_efolds);
    let r = tensor_to_scalar_ratio(n_efolds);
    let h_over_mpl = inflation_hubble_ratio_structural();
    let a_s = scalar_amplitude(n_efolds, h_over_mpl);
    let n_reheat = reheating_efolds_structural();
    let w_reheat = reheating_w_structural();
    let t_reheat_gev = reheating_temperature_gev(h_over_mpl);
    let chi2 = cmb_proxy_chi2(n_s, a_s, r);
    let loglike = cmb_proxy_loglike(chi2);
    let expansion_factor = inflation_expansion_factor(n_efolds);
    let n_end = structural_n_end();

    let n_efolds_ok = n_efolds >= w.n_efolds_min && n_efolds <= w.n_efolds_max;
    let n_s_ok = (n_s - w.ns_center).abs() <= w.ns_tol;
    let r_ok = r <= w.r_max;
    let a_s_ok = (a_s - w.as_center).abs() <= w.as_tol;
    let reheating_ok = t_reheat_gev >= w.t_reheat_min_gev;
    let cmb_like_ok = chi2 <= w.chi2_max;
    let graceful_exit_ok = slow_roll_epsilon(n_end) >= 1.0 - 1e-12;

    InflationScorecard {
        n_efolds,
        epsilon,
        eta,
        n_s,
        r,
        a_s,
        h_over_mpl,
        n_reheat,
        w_reheat,
        t_reheat_gev,
        cmb_proxy_chi2: chi2,
        cmb_proxy_loglike: loglike,
        expansion_factor,
        n_end,
        n_efolds_ok,
        n_s_ok,
        r_ok,
        a_s_ok,
        reheating_ok,
        cmb_like_ok,
        graceful_exit_ok,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn structural_efolds_are_sixty() {
        assert!((inflation_efolds_structural() - 60.0).abs() < 1e-12);
    }

    #[test]
    fn slow_roll_observables_are_physical() {
        let n = inflation_efolds_structural();
        let eps = slow_roll_epsilon(n);
        let eta = slow_roll_eta(n);
        assert!(eps > 0.0 && eps < 1.0);
        assert!(eta < 0.0);
    }

    #[test]
    fn inflation_gate_passes() {
        let s = evaluate_inflation_gate(InflationWindows::default());
        assert!(s.passes_all(), "inflation gate failed: {:?}", s);
    }

    #[test]
    fn scalar_amplitude_is_near_observed_scale() {
        let s = evaluate_inflation_gate(InflationWindows::default());
        assert!(s.a_s > 1.5e-9 && s.a_s < 3.0e-9);
    }

    #[test]
    fn reheating_is_bbn_safe() {
        let s = evaluate_inflation_gate(InflationWindows::default());
        assert!(s.t_reheat_gev > 1.0e-3);
    }
}
