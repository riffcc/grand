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

use num_complex::Complex64;
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

#[derive(Debug, Clone, Copy)]
pub struct MixingEnvelope {
    pub theta12_min_deg: f64,
    pub theta12_max_deg: f64,
    pub theta23_min_deg: f64,
    pub theta23_max_deg: f64,
    pub theta13_min_deg: f64,
    pub theta13_max_deg: f64,
    pub delta_min_deg: f64,
    pub delta_max_deg: f64,
    pub jarlskog_min: f64,
    pub jarlskog_max: f64,
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

// Broad PDG-style CI envelopes for falsification gates.
pub const CKM_PDG_ENVELOPE: MixingEnvelope = MixingEnvelope {
    theta12_min_deg: 12.7,
    theta12_max_deg: 13.3,
    theta23_min_deg: 2.1,
    theta23_max_deg: 2.7,
    theta13_min_deg: 0.16,
    theta13_max_deg: 0.27,
    delta_min_deg: 60.0,
    delta_max_deg: 76.0,
    jarlskog_min: 2.4e-5,
    jarlskog_max: 3.8e-5,
};

pub const PMNS_PDG_ENVELOPE: MixingEnvelope = MixingEnvelope {
    theta12_min_deg: 30.0,
    theta12_max_deg: 37.0,
    theta23_min_deg: 40.0,
    theta23_max_deg: 56.0,
    theta13_min_deg: 7.0,
    theta13_max_deg: 10.5,
    delta_min_deg: 120.0,
    delta_max_deg: 320.0,
    jarlskog_min: -2.0e-2,
    jarlskog_max: -3.0e-3,
};

pub const CKM_CP_J_MIN: f64 = 1.0e-6;
pub const PMNS_CP_J_MIN: f64 = 1.0e-3;
pub const CP_PHASE_TOL_DEG: f64 = 5.0;

fn jarlskog(s12: f64, s23: f64, s13: f64, delta_rad: f64) -> f64 {
    let c12 = (1.0 - s12 * s12).sqrt();
    let c23 = (1.0 - s23 * s23).sqrt();
    let c13 = (1.0 - s13 * s13).sqrt();
    c12 * c23 * c13 * c13 * s12 * s23 * s13 * delta_rad.sin()
}

type CMat3 = [[Complex64; 3]; 3];

fn c_identity() -> CMat3 {
    [
        [Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0), Complex64::new(0.0, 0.0)],
        [Complex64::new(0.0, 0.0), Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)],
        [Complex64::new(0.0, 0.0), Complex64::new(0.0, 0.0), Complex64::new(1.0, 0.0)],
    ]
}

fn c_mul(a: &CMat3, b: &CMat3) -> CMat3 {
    let mut out = [[Complex64::new(0.0, 0.0); 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let mut s = Complex64::new(0.0, 0.0);
            for k in 0..3 {
                s += a[i][k] * b[k][j];
            }
            out[i][j] = s;
        }
    }
    out
}

fn c_conj_transpose(m: &CMat3) -> CMat3 {
    [
        [m[0][0].conj(), m[1][0].conj(), m[2][0].conj()],
        [m[0][1].conj(), m[1][1].conj(), m[2][1].conj()],
        [m[0][2].conj(), m[1][2].conj(), m[2][2].conj()],
    ]
}

fn clamp_unit(x: f64) -> f64 {
    x.clamp(-1.0, 1.0)
}

fn angles_from_unitary(v: &CMat3) -> (f64, f64, f64) {
    let s13 = clamp_unit(v[0][2].norm());
    let c13 = (1.0 - s13 * s13).sqrt();

    let s12 = clamp_unit(v[0][1].norm() / c13.max(1e-15));
    let s23 = clamp_unit(v[1][2].norm() / c13.max(1e-15));
    (s12, s23, s13)
}

/// Jacobi eigen-decomposition for complex Hermitian 3x3.
fn jacobi_eigen_hermitian(mut a: CMat3) -> ([f64; 3], CMat3) {
    let mut u = c_identity();
    let max_iter = 96usize;
    let tol = 1e-12;

    for _ in 0..max_iter {
        let mut p = 0usize;
        let mut q = 1usize;
        let mut max_off = a[0][1].norm();
        for i in 0..3 {
            for j in (i + 1)..3 {
                let off = a[i][j].norm();
                if off > max_off {
                    max_off = off;
                    p = i;
                    q = j;
                }
            }
        }
        if max_off < tol {
            break;
        }

        let apq = a[p][q];
        let gamma = apq.norm();
        if gamma < tol {
            continue;
        }
        let alpha = a[p][p].re;
        let beta = a[q][q].re;
        let tau = (beta - alpha) / (2.0 * gamma);
        let t = if tau >= 0.0 {
            1.0 / (tau + (1.0 + tau * tau).sqrt())
        } else {
            -1.0 / (-tau + (1.0 + tau * tau).sqrt())
        };
        let c = 1.0 / (1.0 + t * t).sqrt();
        let s_abs = c * t;
        let phase = Complex64::from_polar(1.0, apq.arg());
        let s = phase * s_abs;

        for k in 0..3 {
            if k == p || k == q {
                continue;
            }
            let akp = a[k][p];
            let akq = a[k][q];
            let new_kp = akp * c - akq * s.conj();
            let new_kq = akp * s + akq * c;
            a[k][p] = new_kp;
            a[p][k] = new_kp.conj();
            a[k][q] = new_kq;
            a[q][k] = new_kq.conj();
        }
        a[p][p] = Complex64::new(alpha - t * gamma, 0.0);
        a[q][q] = Complex64::new(beta + t * gamma, 0.0);
        a[p][q] = Complex64::new(0.0, 0.0);
        a[q][p] = Complex64::new(0.0, 0.0);

        for row in &mut u {
            let up = row[p];
            let uq = row[q];
            row[p] = up * c - uq * s.conj();
            row[q] = up * s + uq * c;
        }
    }

    let mut idx = [0usize, 1, 2];
    idx.sort_by(|&i, &j| a[i][i].re.partial_cmp(&a[j][j].re).unwrap_or(std::cmp::Ordering::Equal));

    let evals = [a[idx[0]][idx[0]].re, a[idx[1]][idx[1]].re, a[idx[2]][idx[2]].re];
    let mut evecs = [[Complex64::new(0.0, 0.0); 3]; 3];
    for col_new in 0..3 {
        let col_old = idx[col_new];
        for row in 0..3 {
            evecs[row][col_new] = u[row][col_old];
        }
    }
    (evals, evecs)
}

fn apply_column_phase(m: &mut CMat3, col: usize, phase: f64) {
    let e = Complex64::from_polar(1.0, phase);
    for row in m.iter_mut() {
        row[col] *= e;
    }
}

fn apply_row_phase(m: &mut CMat3, row: usize, phase: f64) {
    let e = Complex64::from_polar(1.0, phase);
    for j in 0..3 {
        m[row][j] *= e;
    }
}

/// Canonical rephasing (PDG-like): fix five removable phases.
fn phase_fix_pdg(mut v: CMat3) -> CMat3 {
    for j in 0..3 {
        let phase = -v[0][j].arg();
        apply_column_phase(&mut v, j, phase);
    }
    for i in 1..3 {
        let phase = -v[i][0].arg();
        apply_row_phase(&mut v, i, phase);
    }
    v
}

fn observables_from_unitary(v: &CMat3) -> MixingObservables {
    let v = phase_fix_pdg(*v);
    let (s12, s23, s13) = angles_from_unitary(&v);
    let c12 = (1.0 - s12 * s12).sqrt();
    let c23 = (1.0 - s23 * s23).sqrt();
    let c13 = (1.0 - s13 * s13).sqrt();

    let j = (v[0][0] * v[1][1] * v[0][1].conj() * v[1][0].conj()).im;
    let denom = (c12 * c23 * c13 * c13 * s12 * s23 * s13).max(1e-18);
    let sin_delta = clamp_unit(j / denom);

    let vtd2 = v[2][0].norm_sqr();
    let top = s12 * s12 * s23 * s23 + c12 * c12 * c23 * c23 * s13 * s13 - vtd2;
    let bot = (2.0 * s12 * s23 * c12 * c23 * s13).max(1e-18);
    let cos_delta = clamp_unit(top / bot);

    let mut delta = sin_delta.atan2(cos_delta);
    if delta < 0.0 {
        delta += 2.0 * PI;
    }
    build_observables(s12, s23, s13, delta)
}

fn ckm_mass_textures_from_clifford() -> (CMat3, CMat3) {
    let lambda = 1.0 / (CLIFFORD_DIM + SU2_DIM).sqrt();
    let eta = 1.0 / (GRADE1_DIM * GRADE2_DIM); // 1/24
    let zeta = 1.0 / (CLIFFORD_DIM * AUGMENTED_DIM); // 1/272
    let phi = PI / 3.0 + (1.0 / LATTICE_SHIFT).atan();
    let phi2 = phi / 2.0;

    // Up-sector texture (hierarchical, weakly mixed) from primitive suppressions.
    let mu12 = Complex64::from_polar(lambda * eta, phi2);
    let mu23 = Complex64::from_polar(lambda.powi(3), -phi2);
    let mu13 = Complex64::new(eta * zeta, 0.0);
    let mu = [
        [Complex64::new(lambda.powi(2), 0.0), mu12, mu13],
        [mu12.conj(), Complex64::new(lambda, 0.0), mu23],
        [mu13.conj(), mu23.conj(), Complex64::new(1.0 + eta, 0.0)],
    ];

    // Down-sector texture (Cabibbo-dominant) with Z3 phase placements.
    let md12 = Complex64::from_polar((GRADE1_DIM / (GRADE1_DIM + 1.0)) * lambda, phi2); // 4/5 λ
    let md23 = Complex64::from_polar(eta, phi);
    let md13 = Complex64::from_polar(zeta, phi2);
    let md = [
        [Complex64::new(lambda, 0.0), md12, md13],
        [md12.conj(), Complex64::new(1.0 / (1.0 + lambda), 0.0), md23],
        [md13.conj(), md23.conj(), Complex64::new(2.0, 0.0)],
    ];
    (mu, md)
}

fn pmns_mass_textures_from_clifford() -> (CMat3, CMat3) {
    let eps = 1.0 / LATTICE_SHIFT;
    let eta = 1.0 / (GRADE1_DIM * GRADE2_DIM); // 1/24
    let s12 = (GRADE1_DIM / COMPLEMENT_DIM).sqrt();
    let s23 = (GRADE1_DIM / LATTICE_SHIFT).sqrt();
    let psi = PI + (1.0 / SU2_DIM).atan();
    let psi2 = psi / 2.0;

    // Charged-lepton texture: near-diagonal with tiny complex couplings.
    let ml12 = Complex64::from_polar(eps * eta, psi2);
    let ml23 = Complex64::from_polar(eps.powi(2) * eta, -(1.0 / SU2_DIM).atan());
    let ml13 = Complex64::from_polar(eps.powi(5), psi2);
    let ml = [
        [Complex64::new(eps.powi(4), 0.0), ml12, ml13],
        [ml12.conj(), Complex64::new(eps.powi(3), 0.0), ml23],
        [ml13.conj(), ml23.conj(), Complex64::new(eps.powi(2), 0.0)],
    ];

    // Neutrino texture: large off-diagonal couplings drive large PMNS angles.
    let mnu12 = Complex64::from_polar((SU2_DIM / GRADE1_DIM) * s12, psi2); // 3/4 * sqrt(4/13)
    let mnu23 = Complex64::from_polar((2.0 / 3.0) * s23, -psi2);
    let mnu13 = Complex64::from_polar(eps, psi2);
    let mnu = [
        [Complex64::new(0.0, 0.0), mnu12, mnu13],
        [mnu12.conj(), Complex64::new(eps.powi(2), 0.0), mnu23],
        [mnu13.conj(), mnu23.conj(), Complex64::new(eps, 0.0)],
    ];
    (ml, mnu)
}

/// Derive CKM by diagonalizing explicit algebraic mass textures.
pub fn ckm_from_textures() -> MixingObservables {
    let (mu, md) = ckm_mass_textures_from_clifford();
    let (_, uu) = jacobi_eigen_hermitian(mu);
    let (_, ud) = jacobi_eigen_hermitian(md);
    let v = c_mul(&c_conj_transpose(&uu), &ud);
    observables_from_unitary(&v)
}

/// Derive PMNS by diagonalizing explicit algebraic mass textures.
pub fn pmns_from_textures() -> MixingObservables {
    let (ml, mnu) = pmns_mass_textures_from_clifford();
    let (_, ul) = jacobi_eigen_hermitian(ml);
    let (_, un) = jacobi_eigen_hermitian(mnu);
    let u = c_mul(&c_conj_transpose(&ul), &un);
    observables_from_unitary(&u)
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

pub fn within_envelope(obs: MixingObservables, env: MixingEnvelope) -> Result<(), String> {
    let checks = [
        (
            "theta12_deg",
            obs.theta12_deg,
            env.theta12_min_deg,
            env.theta12_max_deg,
        ),
        (
            "theta23_deg",
            obs.theta23_deg,
            env.theta23_min_deg,
            env.theta23_max_deg,
        ),
        (
            "theta13_deg",
            obs.theta13_deg,
            env.theta13_min_deg,
            env.theta13_max_deg,
        ),
        (
            "delta_deg",
            obs.delta_deg,
            env.delta_min_deg,
            env.delta_max_deg,
        ),
        (
            "jarlskog",
            obs.jarlskog,
            env.jarlskog_min,
            env.jarlskog_max,
        ),
    ];

    for (label, value, lo, hi) in checks {
        if value < lo || value > hi {
            return Err(format!(
                "{label} out of envelope: value={value:.9}, range=[{lo:.9}, {hi:.9}]"
            ));
        }
    }
    Ok(())
}

fn delta_distance_to_cp_conserving_deg(delta_deg: f64) -> f64 {
    let wrapped = ((delta_deg % 360.0) + 360.0) % 360.0;
    let d0 = wrapped.min((360.0 - wrapped).abs());
    let d180 = (wrapped - 180.0).abs();
    d0.min(d180)
}

pub fn cp_violation_witness(
    obs: MixingObservables,
    j_abs_min: f64,
    phase_tol_deg: f64,
) -> Result<(), String> {
    if obs.jarlskog.abs() <= j_abs_min {
        return Err(format!(
            "Jarlskog too small for CPV witness: |J|={:.12e}, min={:.12e}",
            obs.jarlskog.abs(),
            j_abs_min
        ));
    }

    let dist = delta_distance_to_cp_conserving_deg(obs.delta_deg);
    if dist <= phase_tol_deg {
        return Err(format!(
            "phase too close to CP-conserving branch: delta={:.9} deg, dist={:.9} deg, tol={:.9} deg",
            obs.delta_deg, dist, phase_tol_deg
        ));
    }

    Ok(())
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

    #[test]
    fn texture_diagonalization_recovers_ckm_window() {
        let ckm = ckm_from_textures();
        let r = residuals(ckm, CKM_TARGET);
        assert!(r.d_theta12_deg.abs() < 2.0, "theta12 drift too large: {}", r.d_theta12_deg);
        assert!(r.d_theta23_deg.abs() < 1.0, "theta23 drift too large: {}", r.d_theta23_deg);
        assert!(r.d_theta13_deg.abs() < 0.5, "theta13 drift too large: {}", r.d_theta13_deg);
        assert!(r.d_delta_deg.abs() < 30.0, "delta drift too large: {}", r.d_delta_deg);
        assert!(r.d_jarlskog.abs() < 2e-5, "J drift too large: {}", r.d_jarlskog);
    }

    #[test]
    fn texture_diagonalization_keeps_pmns_large() {
        let pmns = pmns_from_textures();
        assert!(pmns.theta12_deg > 25.0);
        assert!(pmns.theta23_deg > 35.0);
        assert!(pmns.theta13_deg > 4.0);
    }

    #[test]
    fn ckm_direct_within_pdg_envelope() {
        let ckm = ckm_from_clifford();
        within_envelope(ckm, CKM_PDG_ENVELOPE).expect("direct CKM outside PDG envelope");
    }

    #[test]
    fn ckm_texture_within_pdg_envelope() {
        let ckm = ckm_from_textures();
        within_envelope(ckm, CKM_PDG_ENVELOPE).expect("texture CKM outside PDG envelope");
    }

    #[test]
    fn pmns_direct_within_pdg_envelope() {
        let pmns = pmns_from_clifford();
        within_envelope(pmns, PMNS_PDG_ENVELOPE).expect("direct PMNS outside PDG envelope");
    }

    #[test]
    fn pmns_texture_within_pdg_envelope() {
        let pmns = pmns_from_textures();
        within_envelope(pmns, PMNS_PDG_ENVELOPE).expect("texture PMNS outside PDG envelope");
    }

    #[test]
    fn ckm_direct_has_cpv_witness() {
        let ckm = ckm_from_clifford();
        cp_violation_witness(ckm, CKM_CP_J_MIN, CP_PHASE_TOL_DEG)
            .expect("direct CKM missing CPV witness");
    }

    #[test]
    fn ckm_texture_has_cpv_witness() {
        let ckm = ckm_from_textures();
        cp_violation_witness(ckm, CKM_CP_J_MIN, CP_PHASE_TOL_DEG)
            .expect("texture CKM missing CPV witness");
    }

    #[test]
    fn pmns_direct_has_cpv_witness() {
        let pmns = pmns_from_clifford();
        cp_violation_witness(pmns, PMNS_CP_J_MIN, CP_PHASE_TOL_DEG)
            .expect("direct PMNS missing CPV witness");
    }

    #[test]
    fn pmns_texture_has_cpv_witness() {
        let pmns = pmns_from_textures();
        cp_violation_witness(pmns, PMNS_CP_J_MIN, CP_PHASE_TOL_DEG)
            .expect("texture PMNS missing CPV witness");
    }
}
