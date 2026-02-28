/*!
 * GUTOE Physics - Molecular Homochirality (Parity-Violating Energy Difference)
 * Copyright (C) 2026  Riff Labs
 *
 * GRAND-368:
 *   Structural lane for amino-acid enantiomer splitting from Cl(1,3) chirality.
 *
 * Scope:
 *   - Derive a zero-fit, channel-count-based parity proxy from shared primitives.
 *   - Map the proxy to an energy scale using weak + electromagnetic constants.
 *   - Expose signed left/right energy shifts for amino-acid enantiomers.
 *
 * This lane is intentionally explicit about assumptions:
 *   - Chemistry model uses the canonical amino-acid backbone heteroatom motif
 *     (N + 2*O) as the dominant parity-sensitive weak-charge source.
 *   - Magnitude target is order-of-magnitude (PVED scale), not sub-ppm spectroscopy.
 */

use crate::constants::{ALPHA_LEADING_ORDER, CLIFFORD_STATE_COUNT_STRUCTURAL};
use gutoe_em::{electron_mass_from_proton_anchor, sin2_weinberg};

/// Fermi constant in GeV^-2.
pub const FERMI_CONSTANT_GEV_INV2: f64 = 1.166_378_7e-5;

/// SU(2) gauge generator count.
pub const SU2_GENERATOR_COUNT: f64 = 3.0;
/// Total SM gauge-generator count.
pub const TOTAL_GAUGE_GENERATOR_COUNT: f64 = 12.0;

/// Weak-sector share of the full gauge algebra: `3/12 = 1/4`.
pub const WEAK_GAUGE_FRACTION: f64 = SU2_GENERATOR_COUNT / TOTAL_GAUGE_GENERATOR_COUNT;

/// Single-channel chiral projection from the 16-state Cl(1,3) basis: `1/16`.
pub const CHIRAL_PROJECTION_FACTOR: f64 = 1.0 / CLIFFORD_STATE_COUNT_STRUCTURAL;

/// Canonical amino-acid backbone heteroatom composition:
/// one N and two O nuclei near the stereocenter.
pub const BACKBONE_N14_MULTIPLICITY: f64 = 1.0;
pub const BACKBONE_O16_MULTIPLICITY: f64 = 2.0;

pub const NITROGEN_Z: f64 = 7.0;
pub const NITROGEN_N: f64 = 7.0;
pub const OXYGEN_Z: f64 = 8.0;
pub const OXYGEN_N: f64 = 8.0;

/// Handedness label for an enantiomer pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Handedness {
    Left,
    Right,
}

/// Weak nuclear charge in this lane:
/// `Q_W = N - (1 - 4 sin²θ_W) Z`.
///
/// With structural `sin²θ_W = 3/13`, this simplifies to `Q_W = N - Z/13`.
pub fn weak_nuclear_charge(z: f64, n: f64) -> f64 {
    let sin2 = sin2_weinberg();
    n - (1.0 - 4.0 * sin2) * z
}

/// Per-nucleus chiral weak weight used in the parity proxy.
pub fn nucleus_chiral_weight(z: f64, n: f64) -> f64 {
    z.powi(3) * weak_nuclear_charge(z, n)
}

/// Canonical amino-acid backbone weak/chiral source factor:
/// `1*N14 + 2*O16`.
pub fn amino_backbone_nuclear_factor() -> f64 {
    BACKBONE_N14_MULTIPLICITY * nucleus_chiral_weight(NITROGEN_Z, NITROGEN_N)
        + BACKBONE_O16_MULTIPLICITY * nucleus_chiral_weight(OXYGEN_Z, OXYGEN_N)
}

/// Pure structural parity factor before alpha suppression:
/// `(weak gauge share) * (chiral projection) * backbone nuclear factor`.
pub fn amino_backbone_parity_factor() -> f64 {
    WEAK_GAUGE_FRACTION * CHIRAL_PROJECTION_FACTOR * amino_backbone_nuclear_factor()
}

/// Dimensionless parity proxy including electromagnetic suppression:
/// `proxy = parity_factor * α^4`.
pub fn amino_backbone_parity_proxy() -> f64 {
    amino_backbone_parity_factor() * ALPHA_LEADING_ORDER.powi(4)
}

/// Structural Rydberg-scale valence energy from shared alpha + electron anchor:
/// `E_R = (1/2) α^2 m_e c^2`.
///
/// `electron_mass_from_proton_anchor()` returns MeV, so multiply by `1e6` for eV.
pub fn rydberg_energy_structural_ev() -> f64 {
    let me_ev = electron_mass_from_proton_anchor() * 1.0e6;
    0.5 * ALPHA_LEADING_ORDER.powi(2) * me_ev
}

/// Weak/electron dimensionless scale:
/// `G_F * m_e^2` (natural units with energies in GeV).
pub fn weak_electron_scale_dimensionless() -> f64 {
    let me_gev = electron_mass_from_proton_anchor() / 1.0e3;
    FERMI_CONSTANT_GEV_INV2 * me_gev.powi(2)
}

/// Predicted amino-acid enantiomer PVED magnitude in eV.
///
/// This yields an order-of-magnitude lane for biological homochirality bias.
pub fn amino_acid_enantiomer_split_ev() -> f64 {
    rydberg_energy_structural_ev()
        * weak_electron_scale_dimensionless()
        * amino_backbone_parity_proxy()
}

/// Signed energy shift for each enantiomer in the pair.
///
/// Convention used in this lane:
/// - `Left` has lower energy (negative shift),
/// - `Right` has higher energy (positive shift),
/// with total splitting `ΔE = E_R - E_L`.
pub fn handedness_energy_shift_ev(h: Handedness) -> f64 {
    let half = 0.5 * amino_acid_enantiomer_split_ev();
    match h {
        Handedness::Left => -half,
        Handedness::Right => half,
    }
}

/// Preferred handedness in this structural sign convention.
pub fn preferred_amino_handedness() -> Handedness {
    Handedness::Left
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn weak_charge_reduces_to_n_minus_z_over_13() {
        let qn = weak_nuclear_charge(NITROGEN_Z, NITROGEN_N);
        let qo = weak_nuclear_charge(OXYGEN_Z, OXYGEN_N);
        assert!((qn - 84.0 / 13.0).abs() < 1e-12);
        assert!((qo - 96.0 / 13.0).abs() < 1e-12);
    }

    #[test]
    fn structural_parity_proxy_is_positive() {
        assert!(amino_backbone_nuclear_factor() > 0.0);
        assert!(amino_backbone_parity_factor() > 0.0);
        assert!(amino_backbone_parity_proxy() > 0.0);
    }

    #[test]
    fn pved_magnitude_is_in_expected_bio_scale_window() {
        let de = amino_acid_enantiomer_split_ev();
        assert!(
            (1.0e-18..=1.0e-16).contains(&de),
            "homochirality PVED out of expected window: {:.6e} eV",
            de
        );
    }

    #[test]
    fn left_is_lower_than_right() {
        let el = handedness_energy_shift_ev(Handedness::Left);
        let er = handedness_energy_shift_ev(Handedness::Right);
        assert!(el < 0.0 && er > 0.0);
        assert!((er - el - amino_acid_enantiomer_split_ev()).abs() < 1e-30);
        assert_eq!(preferred_amino_handedness(), Handedness::Left);
    }
}

