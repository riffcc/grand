/*!
 * GUTOE Physics - Structural Inflation Harness
 * Copyright (C) 2026  Riff Labs
 *
 * GRAND-347:
 *   Zero-free-parameter inflation lane from shared Cl(1,3) primitives.
 */

use crate::constants::{DARK_GEOMETRIC_AMPLIFICATION, DARK_STATE_COUNT_STRUCTURAL};

/// CMB scalar spectral index reference.
pub const NS_OBSERVED: f64 = 0.9649;
/// Conservative first-pass tolerance window for this structural lane.
pub const NS_TOL: f64 = 0.0100;
/// Current observational upper bound on tensor-to-scalar ratio.
pub const R_MAX_OBSERVED: f64 = 0.06;

/// Structural e-fold count from shared Clifford dark-sector counts:
/// N = (geometric amplification) × (dark-state count) = 12 × 5 = 60.
pub fn inflation_efolds_structural() -> f64 {
    DARK_GEOMETRIC_AMPLIFICATION * DARK_STATE_COUNT_STRUCTURAL
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
}

impl Default for InflationWindows {
    fn default() -> Self {
        Self {
            n_efolds_min: 50.0,
            n_efolds_max: 70.0,
            ns_center: NS_OBSERVED,
            ns_tol: NS_TOL,
            r_max: R_MAX_OBSERVED,
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
    pub expansion_factor: f64,
    pub n_end: f64,
    pub n_efolds_ok: bool,
    pub n_s_ok: bool,
    pub r_ok: bool,
    pub graceful_exit_ok: bool,
}

impl InflationScorecard {
    pub const fn passes_all(self) -> bool {
        self.n_efolds_ok && self.n_s_ok && self.r_ok && self.graceful_exit_ok
    }
}

pub fn evaluate_inflation_gate(w: InflationWindows) -> InflationScorecard {
    let n_efolds = inflation_efolds_structural();
    let epsilon = slow_roll_epsilon(n_efolds);
    let eta = slow_roll_eta(n_efolds);
    let n_s = scalar_spectral_index(n_efolds);
    let r = tensor_to_scalar_ratio(n_efolds);
    let expansion_factor = inflation_expansion_factor(n_efolds);
    let n_end = structural_n_end();

    let n_efolds_ok = n_efolds >= w.n_efolds_min && n_efolds <= w.n_efolds_max;
    let n_s_ok = (n_s - w.ns_center).abs() <= w.ns_tol;
    let r_ok = r <= w.r_max;
    let graceful_exit_ok = slow_roll_epsilon(n_end) >= 1.0 - 1e-12;

    InflationScorecard {
        n_efolds,
        epsilon,
        eta,
        n_s,
        r,
        expansion_factor,
        n_end,
        n_efolds_ok,
        n_s_ok,
        r_ok,
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
}
