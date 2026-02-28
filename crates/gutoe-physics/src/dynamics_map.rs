//! Formal-to-runtime bridge map for executable SM dynamics parameters.
//!
//! This module packages theorem-linked constants into a single runtime contract,
//! so simulation kernels can consume one coherent parameter map.

use crate::constants::{ALPHA_LEADING_ORDER, LAMBDA_QG};

#[derive(Debug, Clone, PartialEq)]
pub struct StandardModelDynamicsMap {
    pub clifford_dim: u32,
    pub z3_order: u32,
    pub magnetic_triplet_card: u32,
    pub generations: u32,
    pub sin2_theta_w: f64,
    pub cos2_theta_w: f64,
    pub mz_over_mw_sq: f64,
    pub alpha_leading_order: f64,
    pub lambda_qg: f64,
    pub beta0: f64,
    pub su3_generators: u32,
    pub su2_generators: u32,
    pub u1_generators: u32,
    pub total_gauge_generators: u32,
    pub theta_qcd: f64,
}

impl StandardModelDynamicsMap {
    /// Coupled structural EW bridge coefficient from the triangulated Cl(1,3) lane:
    ///
    /// c_EW = d/2 + |grade₂|/(d-|SU(2)|) - 1/((|grade₂|+1)T(16))
    ///      = 8 + 6/13 - 1/(7*136)
    pub fn ew_shift_coeff_structural(&self) -> f64 {
        let clifford_half_dim = self.clifford_dim as f64 / 2.0; // 16/2 = 8
        let grade2 = 6.0;
        let complement = (self.clifford_dim - self.su2_generators) as f64; // 13
        let lattice_shift = grade2 + 1.0; // 7
        let t16 = 136.0;
        clifford_half_dim + grade2 / complement - 1.0 / (lattice_shift * t16)
    }

    pub fn from_clifford_z3() -> Self {
        let clifford_dim = 16_u32;
        let z3_order = 3_u32;
        let magnetic_triplet_card = 3_u32;
        let sin2_theta_w =
            magnetic_triplet_card as f64 / (clifford_dim - magnetic_triplet_card) as f64; // 3/13
        let cos2_theta_w = 1.0 - sin2_theta_w; // 10/13
        let mz_over_mw_sq = 1.0 / cos2_theta_w; // 13/10
        let su3_generators = z3_order * z3_order - 1; // 8
        let su2_generators = 3_u32;
        let u1_generators = 1_u32;
        Self {
            clifford_dim,
            z3_order,
            magnetic_triplet_card,
            generations: z3_order,
            sin2_theta_w,
            cos2_theta_w,
            mz_over_mw_sq,
            alpha_leading_order: ALPHA_LEADING_ORDER,
            lambda_qg: LAMBDA_QG,
            beta0: 58.0 / 3.0,
            su3_generators,
            su2_generators,
            u1_generators,
            total_gauge_generators: su3_generators + su2_generators + u1_generators,
            // Structural CP-odd QCD phase in current GUTOE map.
            // Runtime bridge keeps this explicit and auditable.
            theta_qcd: 0.0,
        }
    }

    pub fn validate_internal_constraints(&self) -> bool {
        let eps = 1e-12;
        (self.sin2_theta_w - 3.0 / 13.0).abs() < eps
            && (self.cos2_theta_w - 10.0 / 13.0).abs() < eps
            && (self.mz_over_mw_sq - 13.0 / 10.0).abs() < eps
            && (self.alpha_leading_order - 1.0 / 137.0).abs() < eps
            && (self.lambda_qg - 1.0 / 12.0).abs() < eps
            && (self.beta0 - 58.0 / 3.0).abs() < eps
            && self.theta_qcd.abs() < eps
            && self.generations == 3
            && self.total_gauge_generators == 12
    }

    /// Coupled structural EW shift from Cl(1,3) + alpha lane:
    ///
    /// Δsin²θ_W = α² * c_EW
    /// c_EW = 8 + 6/13 - 1/(7*136)
    pub fn sin2_theta_w_mz_shift_structural(&self) -> f64 {
        self.alpha_leading_order.powi(2) * self.ew_shift_coeff_structural()
    }

    /// EW value at M_Z from structural terms only.
    pub fn sin2_theta_w_at_mz(&self) -> f64 {
        self.sin2_theta_w + self.sin2_theta_w_mz_shift_structural()
    }

    /// Minimal neutron EDM bridge from θ_QCD.
    ///
    /// The standard chiral estimate is O(2.4e-16 * θ) e·cm.
    /// Keeping this explicit gives a direct falsifiability hook.
    pub fn neutron_edm_e_cm_from_theta(&self) -> f64 {
        2.4e-16 * self.theta_qcd
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn map_matches_theorem_chain_values() {
        let m = StandardModelDynamicsMap::from_clifford_z3();
        assert!((m.sin2_theta_w - 3.0 / 13.0).abs() < 1e-12);
        assert!((m.mz_over_mw_sq - 13.0 / 10.0).abs() < 1e-12);
        assert_eq!(m.generations, 3);
        assert_eq!(m.total_gauge_generators, 12);
    }

    #[test]
    fn map_internal_constraints_hold() {
        let m = StandardModelDynamicsMap::from_clifford_z3();
        assert!(m.validate_internal_constraints());
    }

    #[test]
    fn ew_bridge_sin2_theta_matches_target_window() {
        let m = StandardModelDynamicsMap::from_clifford_z3();
        let sin2_mz = m.sin2_theta_w_at_mz();
        assert!(
            (0.23100..=0.23140).contains(&sin2_mz),
            "sin²(theta_W) at M_Z out of target window: {sin2_mz:.9}"
        );
    }

    #[test]
    fn ew_bridge_shift_is_positive_and_small() {
        let m = StandardModelDynamicsMap::from_clifford_z3();
        let shift = m.sin2_theta_w_mz_shift_structural();
        assert!(shift > 0.0);
        assert!(shift < 1.0e-3, "unexpectedly large EW shift: {shift:.12e}");
        let expected = (1.0 / 137.0f64).powi(2) * (8.0 + 6.0 / 13.0 - 1.0 / (7.0 * 136.0));
        assert!(
            (shift - expected).abs() < 1.0e-15,
            "structural EW shift mismatch: got {shift:.15e}, expected {expected:.15e}"
        );
    }

    #[test]
    fn ew_bridge_coefficient_matches_coupled_structural_form() {
        let m = StandardModelDynamicsMap::from_clifford_z3();
        let coeff = m.ew_shift_coeff_structural();
        let expected = 8.0 + 6.0 / 13.0 - 1.0 / (7.0 * 136.0);
        assert!(
            (coeff - expected).abs() < 1.0e-15,
            "EW coeff mismatch: got {coeff:.15e}, expected {expected:.15e}"
        );
    }

    #[test]
    fn strong_cp_structural_defaults_are_zero() {
        let m = StandardModelDynamicsMap::from_clifford_z3();
        assert!(m.theta_qcd.abs() < 1e-15);
        assert!(m.neutron_edm_e_cm_from_theta().abs() < 1e-30);
    }
}
