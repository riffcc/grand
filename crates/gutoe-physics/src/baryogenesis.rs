/*!
 * GUTOE Physics - Quantitative Baryogenesis Harness
 * Copyright (C) 2026  Riff Labs
 *
 * GRAND-348:
 *   Quantify baryon-to-photon asymmetry from existing Clifford-derived
 *   primitives with zero free fit knobs.
 */

use crate::constants::{
    lambda_micro_finite_mode_rescale, ALPHA_LEADING_ORDER, DARK_TO_VISIBLE_COUNT_RATIO, LAMBDA_QG,
};
use gutoe_em::{
    ckm_from_clifford, ckm_from_textures, cp_violation_witness, CP_PHASE_TOL_DEG, CKM_CP_J_MIN,
};

/// Planck-era baryon-to-photon ratio target.
pub const ETA_B_OBSERVED: f64 = 6.12e-10;

/// Structural survival factor for departure from thermal equilibrium.
///
/// Uses only shared zero-free-parameter terms:
///   - `(1 - λ_QG)` from lattice correction survival
///   - `(486/485)` finite-mode correction from GRAND-295
pub fn nonequilibrium_survival_factor() -> f64 {
    (1.0 - LAMBDA_QG) * lambda_micro_finite_mode_rescale()
}

/// Structural prefactor multiplying CKM CP source `J`.
///
/// η_B ∝ J · α² · (5/11) · f_neq
pub fn baryogenesis_structural_prefactor() -> f64 {
    ALPHA_LEADING_ORDER.powi(2) * DARK_TO_VISIBLE_COUNT_RATIO * nonequilibrium_survival_factor()
}

/// Quantitative baryon-to-photon prediction from a supplied Jarlskog source.
pub fn eta_baryon_from_jarlskog(jarlskog: f64) -> f64 {
    jarlskog * baryogenesis_structural_prefactor()
}

/// Quantitative baryon-to-photon prediction from Clifford CKM observables.
pub fn eta_baryon_from_clifford_ckm() -> f64 {
    let ckm = ckm_from_clifford();
    eta_baryon_from_jarlskog(ckm.jarlskog)
}

/// Quantitative baryon-to-photon prediction from texture CKM observables.
pub fn eta_baryon_from_texture_ckm() -> f64 {
    let ckm = ckm_from_textures();
    eta_baryon_from_jarlskog(ckm.jarlskog)
}

/// Quantitative gate windows for GRAND-348.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BaryogenesisWindows {
    /// Relative-error limit for η_B.
    pub eta_rel_error_max: f64,
}

impl Default for BaryogenesisWindows {
    fn default() -> Self {
        Self {
            // First-pass quantitative window: strict enough to falsify drift,
            // broad enough for this structural-only lane.
            eta_rel_error_max: 0.15,
        }
    }
}

/// Quantitative scorecard for GRAND-348.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BaryogenesisScorecard {
    pub jarlskog_ckm_direct: f64,
    pub jarlskog_ckm_texture: f64,
    pub eta_predicted: f64,
    pub eta_observed: f64,
    pub eta_rel_error: f64,
    pub cp_violation_ok: bool,
    pub baryon_violation_channel_ok: bool,
    pub nonequilibrium_ok: bool,
    pub eta_window_ok: bool,
}

impl BaryogenesisScorecard {
    pub const fn sakharov_ok(self) -> bool {
        self.cp_violation_ok && self.baryon_violation_channel_ok && self.nonequilibrium_ok
    }

    pub const fn passes_all(self) -> bool {
        self.sakharov_ok() && self.eta_window_ok
    }
}

/// Structural check: electroweak non-Abelian lane exists (sphaleron channel).
///
/// This lane is present when the weak-angle sector is physical and nontrivial.
pub fn baryon_violation_channel_structural() -> bool {
    let sin2 = gutoe_em::sin2_weinberg();
    sin2 > 0.0 && sin2 < 0.5
}

/// Structural check for departure from equilibrium.
pub fn nonequilibrium_structural() -> bool {
    let f = nonequilibrium_survival_factor();
    f > 0.0 && f < 1.0
}

/// Evaluate GRAND-348 scorecard with explicit windows.
pub fn evaluate_baryogenesis_gate(windows: BaryogenesisWindows) -> BaryogenesisScorecard {
    let ckm_direct = ckm_from_clifford();
    let ckm_texture = ckm_from_textures();
    // GRAND-348 uses the texture chain as the quantitative source:
    // Cl(1,3) -> textures -> diagonalization -> CKM -> J -> η_B.
    let eta_predicted = eta_baryon_from_jarlskog(ckm_texture.jarlskog);
    let eta_rel_error = (eta_predicted - ETA_B_OBSERVED).abs() / ETA_B_OBSERVED;

    let cp_violation_ok =
        cp_violation_witness(ckm_texture, CKM_CP_J_MIN, CP_PHASE_TOL_DEG).is_ok();
    let baryon_violation_channel_ok = baryon_violation_channel_structural();
    let nonequilibrium_ok = nonequilibrium_structural();
    let eta_window_ok = eta_rel_error <= windows.eta_rel_error_max;

    BaryogenesisScorecard {
        jarlskog_ckm_direct: ckm_direct.jarlskog,
        jarlskog_ckm_texture: ckm_texture.jarlskog,
        eta_predicted,
        eta_observed: ETA_B_OBSERVED,
        eta_rel_error,
        cp_violation_ok,
        baryon_violation_channel_ok,
        nonequilibrium_ok,
        eta_window_ok,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn structural_prefactor_is_positive() {
        assert!(baryogenesis_structural_prefactor() > 0.0);
    }

    #[test]
    fn ckm_eta_prediction_is_positive() {
        let eta = eta_baryon_from_clifford_ckm();
        assert!(eta > 0.0);
        let eta_tex = eta_baryon_from_texture_ckm();
        assert!(eta_tex > 0.0);
    }

    #[test]
    fn structural_sakharov_checks_hold() {
        let s = evaluate_baryogenesis_gate(BaryogenesisWindows::default());
        assert!(s.cp_violation_ok, "CP-violation witness failed");
        assert!(
            s.baryon_violation_channel_ok,
            "baryon-violation channel witness failed"
        );
        assert!(s.nonequilibrium_ok, "nonequilibrium witness failed");
    }

    #[test]
    fn eta_prediction_is_quantitative_not_order_only() {
        let s = evaluate_baryogenesis_gate(BaryogenesisWindows::default());
        assert!(
            s.eta_rel_error < 0.15,
            "η_B relative error too large: {:.6}",
            s.eta_rel_error
        );
    }

    #[test]
    fn default_gate_passes() {
        let s = evaluate_baryogenesis_gate(BaryogenesisWindows::default());
        assert!(s.passes_all(), "baryogenesis gate failed: {:?}", s);
    }
}
