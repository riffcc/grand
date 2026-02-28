/*!
 * GUTOE Physics - Chiral Symmetry Breaking Structural Lane
 * Copyright (C) 2026  Riff Labs
 *
 * GRAND-126:
 *   Chiral symmetry breaking from shared Cl(1,3) primitives:
 *   - nonzero quark condensate proxy
 *   - pseudo-Goldstone pion scaling
 *   - confinement-linked positive witness
 */

use crate::constants::{ALPHA_LEADING_ORDER, LAMBDA_QG};
use crate::dynamics_map::StandardModelDynamicsMap;

/// Structural quark condensate proxy from shared primitives:
/// `-(1 - λ_QG) * (|quarkOrbit| / dim Cl(1,3))`.
pub fn quark_condensate_proxy() -> f64 {
    let sm = StandardModelDynamicsMap::from_clifford_z3();
    -((1.0 - LAMBDA_QG) * (sm.generations as f64 / sm.clifford_dim as f64))
}

/// Explicit chiral-symmetry breaking scale used in this structural lane.
pub fn chiral_explicit_breaking_alpha() -> f64 {
    ALPHA_LEADING_ORDER
}

/// Pion mass-squared proxy for the pseudo-Goldstone channel.
///
/// In this reduced structural lane:
/// `m_pi^2 ∝ alpha * (-<qq>)`.
pub fn pion_mass_sq_proxy() -> f64 {
    chiral_explicit_breaking_alpha() * (-quark_condensate_proxy())
}

/// Positive mass proxy derived from the mass-squared lane.
pub fn pion_mass_proxy() -> f64 {
    pion_mass_sq_proxy().max(0.0).sqrt()
}

/// Pseudo-Goldstone ratio:
/// `m_pi^2 / (-<qq>) = alpha`.
pub fn pseudo_goldstone_ratio() -> f64 {
    let condensate_mag = (-quark_condensate_proxy()).max(f64::EPSILON);
    pion_mass_sq_proxy() / condensate_mag
}

/// Chiral-limit map for the pion mass-squared proxy.
pub fn pion_mass_sq_from_explicit_breaking(explicit_breaking: f64) -> f64 {
    explicit_breaking * (-quark_condensate_proxy())
}

/// Positive witness linking chiral condensation and confinement:
/// `beta0 * (-<qq>)`.
pub fn confinement_chiral_link_strength() -> f64 {
    let sm = StandardModelDynamicsMap::from_clifford_z3();
    sm.beta0 * (-quark_condensate_proxy())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quark_condensate_proxy_matches_closed_form() {
        let got = quark_condensate_proxy();
        let expected = -11.0 / 64.0;
        assert!(
            (got - expected).abs() < 1.0e-12,
            "condensate proxy mismatch: got={got:.15e}, expected={expected:.15e}"
        );
        assert!(got < 0.0, "expected negative condensate proxy");
    }

    #[test]
    fn pion_mass_sq_proxy_matches_closed_form() {
        let got = pion_mass_sq_proxy();
        let expected = 11.0 / 8768.0;
        assert!(
            (got - expected).abs() < 1.0e-15,
            "pion mass-squared mismatch: got={got:.15e}, expected={expected:.15e}"
        );
        assert!(got > 0.0, "expected positive pion mass-squared proxy");
    }

    #[test]
    fn pseudo_goldstone_ratio_matches_alpha() {
        let got = pseudo_goldstone_ratio();
        let expected = 1.0 / 137.0;
        assert!(
            (got - expected).abs() < 1.0e-15,
            "pseudo-Goldstone ratio mismatch: got={got:.15e}, expected={expected:.15e}"
        );
    }

    #[test]
    fn chiral_limit_recovers_massless_pion_proxy() {
        let at_zero = pion_mass_sq_from_explicit_breaking(0.0);
        assert!(at_zero.abs() < 1.0e-18);
    }

    #[test]
    fn confinement_chiral_link_is_positive() {
        let link = confinement_chiral_link_strength();
        let expected = 319.0 / 96.0;
        assert!(
            (link - expected).abs() < 1.0e-12,
            "confinement link mismatch: got={link:.15e}, expected={expected:.15e}"
        );
        assert!(link > 0.0);
    }
}
