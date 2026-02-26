// GUTOE EM — CKM/PMNS flavor-mixing observables from Cl(1,3) integers
// Copyright (C) 2026 Riff Labs, AGPL-3.0-or-later
//
// This module provides a zero-free-parameter algebraic map for:
//   - CKM angles (theta12, theta23, theta13), CP phase delta, Jarlskog J
//   - PMNS angles (theta12, theta23, theta13), CP phase delta, Jarlskog J
//
// Structural inputs are shared Cl(1,3) counts:
//   clifford_dim = 2^4 = 16
//   su2_dim = 3
//   grade1_dim = 4
//   grade2_dim = 6
//   complement_dim = clifford_dim - su2_dim = 13

use std::f64::consts::PI;

const CLIFFORD_DIM: f64 = 16.0; // 2^4
const SU2_DIM: f64 = 3.0;
const GRADE1_DIM: f64 = 4.0;
const GRADE2_DIM: f64 = 6.0;
const COMPLEMENT_DIM: f64 = CLIFFORD_DIM - SU2_DIM; // 13
const AUGMENTED_DIM: f64 = CLIFFORD_DIM + 1.0; // 17
const LATTICE_SHIFT: f64 = GRADE2_DIM + 1.0; // 7

#[derive(Debug, Clone, Copy)]
pub struct MixingObservables {
    pub s12: f64,
    pub s23: f64,
    pub s13: f64,
    pub delta_rad: f64,
    pub theta12_deg: f64,
    pub theta23_deg: f64,
    pub theta13_deg: f64,
    pub delta_deg: f64,
    pub jarlskog: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct MixingTargets {
    pub theta12_deg: f64,
    pub theta23_deg: f64,
    pub theta13_deg: f64,
    pub delta_deg: f64,
    pub jarlskog: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct MixingResiduals {
    pub d_theta12_deg: f64,
    pub d_theta23_deg: f64,
    pub d_theta13_deg: f64,
    pub d_delta_deg: f64,
    pub d_jarlskog: f64,
}

pub const CKM_TARGET: MixingTargets = MixingTargets {
    theta12_deg: 13.0,
    theta23_deg: 2.4,
    theta13_deg: 0.2,
    delta_deg: 68.0,
    jarlskog: 3.0e-5,
};

pub const PMNS_TARGET: MixingTargets = MixingTargets {
    theta12_deg: 33.4,
    theta23_deg: 49.0,
    theta13_deg: 8.5,
    // Current global fits still carry broad uncertainty on leptonic delta.
    // Use a representative center for diagnostics only.
    delta_deg: 197.0,
    jarlskog: -1.0e-2,
};

fn jarlskog(s12: f64, s23: f64, s13: f64, delta_rad: f64) -> f64 {
    let c12 = (1.0 - s12 * s12).sqrt();
    let c23 = (1.0 - s23 * s23).sqrt();
    let c13 = (1.0 - s13 * s13).sqrt();
    c12 * c23 * c13 * c13 * s12 * s23 * s13 * delta_rad.sin()
}

fn build_observables(s12: f64, s23: f64, s13: f64, delta_rad: f64) -> MixingObservables {
    let theta12_deg = s12.asin().to_degrees();
    let theta23_deg = s23.asin().to_degrees();
    let theta13_deg = s13.asin().to_degrees();
    let delta_deg = delta_rad.to_degrees();
    let j = jarlskog(s12, s23, s13, delta_rad);
    MixingObservables {
        s12,
        s23,
        s13,
        delta_rad,
        theta12_deg,
        theta23_deg,
        theta13_deg,
        delta_deg,
        jarlskog: j,
    }
}

/// CKM observables from Cl(1,3) integer structure.
///
/// Definitions:
///   s12 = 1 / sqrt(16 + 3)         (Cabibbo sector from dim + SU(2) mismatch)
///   s23 = 1 / (4 * 6)              (grade-1 / grade-2 coupling suppression)
///   s13 = 1 / (16 * 17)            (Clifford × augmented-dimension suppression)
///   delta = pi/3 + atan(1 / (6+1)) (Z3 phase with lattice correction)
pub fn ckm_from_clifford() -> MixingObservables {
    let s12 = 1.0 / (CLIFFORD_DIM + SU2_DIM).sqrt();
    let s23 = 1.0 / (GRADE1_DIM * GRADE2_DIM);
    let s13 = 1.0 / (CLIFFORD_DIM * AUGMENTED_DIM);
    let delta = PI / 3.0 + (1.0 / LATTICE_SHIFT).atan();
    build_observables(s12, s23, s13, delta)
}

/// PMNS observables from Cl(1,3) integer structure.
///
/// Definitions:
///   sin²(theta12) = 4/13
///   sin²(theta23) = 4/7
///   sin(theta13)  = 1/7
///   delta         = pi + atan(1/3)
pub fn pmns_from_clifford() -> MixingObservables {
    let s12 = (GRADE1_DIM / COMPLEMENT_DIM).sqrt();
    let s23 = (GRADE1_DIM / LATTICE_SHIFT).sqrt();
    let s13 = 1.0 / LATTICE_SHIFT;
    let delta = PI + (1.0 / SU2_DIM).atan();
    build_observables(s12, s23, s13, delta)
}

pub fn residuals(pred: MixingObservables, target: MixingTargets) -> MixingResiduals {
    MixingResiduals {
        d_theta12_deg: pred.theta12_deg - target.theta12_deg,
        d_theta23_deg: pred.theta23_deg - target.theta23_deg,
        d_theta13_deg: pred.theta13_deg - target.theta13_deg,
        d_delta_deg: pred.delta_deg - target.delta_deg,
        d_jarlskog: pred.jarlskog - target.jarlskog,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ckm_observables_land_in_expected_window() {
        let ckm = ckm_from_clifford();
        let r = residuals(ckm, CKM_TARGET);
        assert!(r.d_theta12_deg.abs() < 1.0, "theta12 drift too large: {}", r.d_theta12_deg);
        assert!(r.d_theta23_deg.abs() < 0.4, "theta23 drift too large: {}", r.d_theta23_deg);
        assert!(r.d_theta13_deg.abs() < 0.08, "theta13 drift too large: {}", r.d_theta13_deg);
        assert!(r.d_delta_deg.abs() < 2.0, "delta drift too large: {}", r.d_delta_deg);
        assert!(r.d_jarlskog.abs() < 5e-6, "J drift too large: {}", r.d_jarlskog);
    }

    #[test]
    fn pmns_observables_capture_large_lepton_mixing() {
        let pmns = pmns_from_clifford();
        assert!(pmns.theta12_deg > 30.0);
        assert!(pmns.theta23_deg > 45.0);
        assert!(pmns.theta13_deg > 7.0);
        // Leptonic mixing should be qualitatively larger than CKM.
        let ckm = ckm_from_clifford();
        assert!(pmns.theta12_deg > ckm.theta12_deg);
        assert!(pmns.theta23_deg > ckm.theta23_deg);
        assert!(pmns.theta13_deg > ckm.theta13_deg);
    }
}
