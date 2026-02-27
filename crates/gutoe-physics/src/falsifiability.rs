//! Quantitative falsifiability gates for core GUTOE observables.
//!
//! These thresholds mirror `findings/009-falsifiable-predictions-catalog.md`.

use crate::{
    constants::{ALPHA, LAMBDA_QG},
    dynamics_map::StandardModelDynamicsMap,
};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GateWindow {
    pub min: f64,
    pub max: f64,
}

impl GateWindow {
    pub const fn contains(self, v: f64) -> bool {
        v >= self.min && v <= self.max
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FalsificationWindows {
    pub sin2_theta_w_ew: GateWindow,
    pub mz_over_mw: GateWindow,
    pub alpha_inverse: GateWindow,
    pub theta_qcd_abs_max: f64,
    pub neutron_edm_abs_max_e_cm: f64,
    pub lambda_qg_abs_tol: f64,
}

impl Default for FalsificationWindows {
    fn default() -> Self {
        Self {
            sin2_theta_w_ew: GateWindow {
                min: 0.23100,
                max: 0.23140,
            },
            mz_over_mw: GateWindow {
                min: 1.1335,
                max: 1.1355,
            },
            alpha_inverse: GateWindow {
                min: 137.034,
                max: 137.038,
            },
            // |d_n| < 1e-26 e*cm implies |theta_qcd| < 1e-26 / (2.4e-16) ~= 4.17e-11.
            // Keep a tiny safety margin in the default gate.
            theta_qcd_abs_max: 4.2e-11,
            neutron_edm_abs_max_e_cm: 1.0e-26,
            lambda_qg_abs_tol: 1e-12,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CorrectedObservables {
    pub sin2_theta_w_ew: f64,
    pub mz_over_mw: f64,
    pub alpha_inverse: f64,
    pub theta_qcd: f64,
    pub neutron_edm_e_cm: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FalsificationScorecard {
    pub structural_ok: bool,
    pub corrected_ok: bool,
    pub gauge_count_ok: bool,
    pub lambda_qg_ok: bool,
}

impl FalsificationScorecard {
    pub const fn passes_all(self) -> bool {
        self.structural_ok && self.corrected_ok && self.gauge_count_ok && self.lambda_qg_ok
    }
}

/// Structural gates that do not require measured-scale correction inputs.
pub fn evaluate_structural() -> FalsificationScorecard {
    let m = StandardModelDynamicsMap::from_clifford_z3();
    let lambda_ok = LAMBDA_QG > 0.0
        && (LAMBDA_QG - 1.0 / 12.0).abs() <= FalsificationWindows::default().lambda_qg_abs_tol;
    FalsificationScorecard {
        structural_ok: m.validate_internal_constraints(),
        corrected_ok: true,
        gauge_count_ok: m.total_gauge_generators == 12,
        lambda_qg_ok: lambda_ok,
    }
}

/// Full gate evaluation against corrected observables.
pub fn evaluate_with_corrected(
    corrected: CorrectedObservables,
    windows: FalsificationWindows,
) -> FalsificationScorecard {
    let structural = evaluate_structural();
    let corrected_ok = windows.sin2_theta_w_ew.contains(corrected.sin2_theta_w_ew)
        && windows.mz_over_mw.contains(corrected.mz_over_mw)
        && windows.alpha_inverse.contains(corrected.alpha_inverse)
        && corrected.theta_qcd.abs() <= windows.theta_qcd_abs_max
        && corrected.neutron_edm_e_cm.abs() <= windows.neutron_edm_abs_max_e_cm;
    FalsificationScorecard {
        structural_ok: structural.structural_ok,
        corrected_ok,
        gauge_count_ok: structural.gauge_count_ok,
        lambda_qg_ok: structural.lambda_qg_ok,
    }
}

/// Default corrected observables currently wired from runtime constants.
///
/// `sin²θ_W` and `mZ/mW` remain provisional until the full no-free-parameter EW
/// correction chain is wired; this function keeps that status explicit.
pub fn provisional_corrected_observables() -> CorrectedObservables {
    let m = StandardModelDynamicsMap::from_clifford_z3();
    CorrectedObservables {
        sin2_theta_w_ew: m.sin2_theta_w_at_mz(),
        mz_over_mw: m.mz_over_mw_sq.sqrt(),
        alpha_inverse: 1.0 / ALPHA,
        theta_qcd: m.theta_qcd,
        neutron_edm_e_cm: m.neutron_edm_e_cm_from_theta(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn structural_gates_hold_in_repo_defaults() {
        let s = evaluate_structural();
        assert!(s.structural_ok);
        assert!(s.gauge_count_ok);
        assert!(s.lambda_qg_ok);
    }

    #[test]
    fn corrected_windows_expose_remaining_bridge_gap() {
        let corr = provisional_corrected_observables();
        let s = evaluate_with_corrected(corr, FalsificationWindows::default());
        assert!(
            !s.corrected_ok,
            "corrected gates unexpectedly pass before full correction bridge is wired: {s:?}"
        );
    }

    #[test]
    fn strong_cp_gate_passes_at_structural_zero_theta() {
        let windows = FalsificationWindows::default();
        let corr = provisional_corrected_observables();
        assert!(corr.theta_qcd.abs() <= windows.theta_qcd_abs_max);
        assert!(corr.neutron_edm_e_cm.abs() <= windows.neutron_edm_abs_max_e_cm);
    }

    #[test]
    fn strong_cp_gate_fails_when_theta_exceeds_bound() {
        let windows = FalsificationWindows::default();
        let corr = CorrectedObservables {
            sin2_theta_w_ew: 0.23122,
            mz_over_mw: 1.1345,
            alpha_inverse: 137.036,
            theta_qcd: windows.theta_qcd_abs_max * 1.05,
            neutron_edm_e_cm: windows.neutron_edm_abs_max_e_cm * 1.05,
        };
        let s = evaluate_with_corrected(corr, windows);
        assert!(
            !s.corrected_ok,
            "theta/EDM overflow must fail corrected gates"
        );
    }
}
