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

use crate::alpha::{ALPHA_INVERSE_PHYSICAL, ALPHA_INVERSE_STRUCTURAL};
use num_complex::Complex64;
use serde::Serialize;
use std::f64::consts::PI;

const CLIFFORD_DIM: f64 = 16.0; // 2^4
const SU2_DIM: f64 = 3.0;
const GRADE1_DIM: f64 = 4.0;
const GRADE2_DIM: f64 = 6.0;
const COMPLEMENT_DIM: f64 = CLIFFORD_DIM - SU2_DIM; // 13
const AUGMENTED_DIM: f64 = CLIFFORD_DIM + 1.0; // 17
const LATTICE_SHIFT: f64 = GRADE2_DIM + 1.0; // 7
pub const PMNS_THETA23_ALPHA2_COEFF_STRUCTURAL: f64 = ALPHA_INVERSE_STRUCTURAL / 4.0; // 137/4

#[derive(Debug, Clone, Copy, Serialize)]
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

#[derive(Debug, Clone, Copy, Serialize)]
pub struct MixingTargets {
    pub theta12_deg: f64,
    pub theta23_deg: f64,
    pub theta13_deg: f64,
    pub delta_deg: f64,
    pub jarlskog: f64,
}

#[derive(Debug, Clone, Copy, Serialize)]
pub struct MixingResiduals {
    pub d_theta12_deg: f64,
    pub d_theta23_deg: f64,
    pub d_theta13_deg: f64,
    pub d_delta_deg: f64,
    pub d_jarlskog: f64,
}

#[derive(Debug, Clone, Copy, Serialize)]
pub struct NeutrinoAbsoluteMasses {
    pub alpha_physical: f64,
    pub electron_mass_anchor_ev: f64,
    pub mass_scale_ev: f64,
    pub hierarchy_exponent: f64,
    pub m1_ev: f64,
    pub m2_ev: f64,
    pub m3_ev: f64,
    pub sum_ev: f64,
    pub dm21_ev2: f64,
    pub dm32_ev2: f64,
    pub dm31_ev2: f64,
    pub splitting_ratio_32_over_21: f64,
}

#[derive(Debug, Clone, Copy, Serialize)]
pub struct NeutrinoTriangulatedSolution {
    pub r1: f64,
    pub r2: f64,
    pub p_triangulated: f64,
    pub ratio_target: f64,
    pub ratio_fit: f64,
    pub ratio_fit_rel_err: f64,
    pub kappa_dm21: f64,
    pub kappa_dm32: f64,
    pub kappa_geo: f64,
    pub kappa_consistency_rel: f64,
    pub mass_scale_ev: f64,
    pub m1_ev: f64,
    pub m2_ev: f64,
    pub m3_ev: f64,
    pub dm21_ev2: f64,
    pub dm32_ev2: f64,
}

#[derive(Debug, Clone, Copy, Serialize)]
pub struct EwShiftTriangulatedSolution {
    pub alpha: f64,
    pub sin2_structural: f64,
    pub sin2_target_mz: f64,
    pub shift_structural: f64,
    pub coeff_structural: f64,
    pub coeff_required: f64,
    pub coeff_rel_delta: f64,
    pub sin2_structural_mz: f64,
    pub sin2_structural_abs_err: f64,
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
        [
            Complex64::new(1.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
        ],
        [
            Complex64::new(0.0, 0.0),
            Complex64::new(1.0, 0.0),
            Complex64::new(0.0, 0.0),
        ],
        [
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(1.0, 0.0),
        ],
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

/// Oscillation splitting ratio model from texture eigenvalue ratios and exponent.
///
/// Let `r1 = |λ1|/|λ3|`, `r2 = |λ2|/|λ3|`, and `m_i ∝ r_i^p` (with `m3` anchor).
/// Then:
///   Δm²21 ∝ r2^(2p) - r1^(2p)
///   Δm²32 ∝ 1 - r2^(2p)
///   R = Δm²32 / Δm²21 = (1 - r2^(2p)) / (r2^(2p) - r1^(2p))
pub fn neutrino_splitting_ratio_from_exponent(r1: f64, r2: f64, exponent: f64) -> f64 {
    let x = r1.powf(2.0 * exponent);
    let y = r2.powf(2.0 * exponent);
    let denom = y - x;
    if denom.abs() < 1.0e-30 {
        return f64::INFINITY;
    }
    (1.0 - y) / denom
}

/// Solve exponent `p` from target splitting ratio using log-space minimization.
///
/// Returns `(p_best, ratio_fit, signed_rel_err)`.
pub fn solve_neutrino_exponent_for_ratio(
    target_ratio: f64,
    r1: f64,
    r2: f64,
) -> (f64, f64, f64) {
    let mut low = 0.05f64;
    let mut high = 80.0f64;
    let mut best_p = neutrino_hierarchy_exponent_structural();
    let mut best_ratio = neutrino_splitting_ratio_from_exponent(r1, r2, best_p);
    let mut best_err = if best_ratio.is_finite() && best_ratio > 0.0 {
        (best_ratio.ln() - target_ratio.ln()).abs()
    } else {
        f64::INFINITY
    };

    for _ in 0..4 {
        let steps = 40_000usize;
        let span = high - low;
        for i in 0..=steps {
            let p = low + span * (i as f64 / steps as f64);
            let ratio = neutrino_splitting_ratio_from_exponent(r1, r2, p);
            if !(ratio.is_finite() && ratio > 0.0) {
                continue;
            }
            let err = (ratio.ln() - target_ratio.ln()).abs();
            if err < best_err {
                best_err = err;
                best_p = p;
                best_ratio = ratio;
            }
        }

        let window = (span / steps as f64) * 100.0;
        low = (best_p - window).max(1.0e-6);
        high = best_p + window;
    }

    let signed_rel_err = (best_ratio - target_ratio) / target_ratio;
    (best_p, best_ratio, signed_rel_err)
}

/// Triangulate neutrino exponent + mass normalization from oscillation anchors.
///
/// This is a forced-parameter diagnostic lane; it does not claim zero-parameter closure.
pub fn triangulate_neutrino_from_splittings(
    dm21_target_ev2: f64,
    dm32_target_ev2: f64,
) -> NeutrinoTriangulatedSolution {
    let mut raw = neutrino_texture_eigenvalues().map(|x| x.abs());
    raw.sort_by(|a, b| a.total_cmp(b));
    let raw_max = raw[2].max(1.0e-18);
    let r1 = raw[0] / raw_max;
    let r2 = raw[1] / raw_max;

    let ratio_target = dm32_target_ev2 / dm21_target_ev2;
    let (p_triangulated, ratio_fit, ratio_fit_rel_err) =
        solve_neutrino_exponent_for_ratio(ratio_target, r1, r2);

    let y1 = r1.powf(2.0 * p_triangulated);
    let y2 = r2.powf(2.0 * p_triangulated);
    let s21 = (y2 - y1).max(1.0e-30);
    let s32 = (1.0 - y2).max(1.0e-30);

    let alpha = 1.0 / ALPHA_INVERSE_PHYSICAL;
    let electron_mass_anchor_ev = crate::weak::electron_mass_from_proton_anchor() * 1.0e6;
    let alpha4 = alpha.powi(4);

    let mass_scale_dm21 = (dm21_target_ev2 / s21).sqrt();
    let mass_scale_dm32 = (dm32_target_ev2 / s32).sqrt();
    let kappa_dm21 = mass_scale_dm21 / (electron_mass_anchor_ev * alpha4);
    let kappa_dm32 = mass_scale_dm32 / (electron_mass_anchor_ev * alpha4);
    let kappa_geo = (kappa_dm21 * kappa_dm32).sqrt();
    let kappa_consistency_rel = if kappa_geo > 0.0 {
        (kappa_dm32 - kappa_dm21) / kappa_geo
    } else {
        f64::INFINITY
    };

    let mass_scale_ev = electron_mass_anchor_ev * alpha4 * kappa_geo;
    let m1_ev = mass_scale_ev * r1.powf(p_triangulated);
    let m2_ev = mass_scale_ev * r2.powf(p_triangulated);
    let m3_ev = mass_scale_ev;
    let dm21_ev2 = m2_ev * m2_ev - m1_ev * m1_ev;
    let dm32_ev2 = m3_ev * m3_ev - m2_ev * m2_ev;

    NeutrinoTriangulatedSolution {
        r1,
        r2,
        p_triangulated,
        ratio_target,
        ratio_fit,
        ratio_fit_rel_err,
        kappa_dm21,
        kappa_dm32,
        kappa_geo,
        kappa_consistency_rel,
        mass_scale_ev,
        m1_ev,
        m2_ev,
        m3_ev,
        dm21_ev2,
        dm32_ev2,
    }
}

/// Triangulate the EW M_Z bridge coefficient from observed `sin²θ_W(M_Z)`.
pub fn triangulate_ew_shift_for_target(sin2_target_mz: f64) -> EwShiftTriangulatedSolution {
    let alpha = 1.0 / ALPHA_INVERSE_STRUCTURAL;
    let sin2_structural = crate::weak::sin2_weinberg();
    let coeff_structural = CLIFFORD_DIM / 2.0; // d/2 = 8
    let shift_structural = alpha * alpha * coeff_structural;
    let sin2_structural_mz = sin2_structural + shift_structural;
    let coeff_required = (sin2_target_mz - sin2_structural) / (alpha * alpha);
    let coeff_rel_delta = (coeff_required - coeff_structural) / coeff_structural;
    let sin2_structural_abs_err = (sin2_structural_mz - sin2_target_mz).abs();

    EwShiftTriangulatedSolution {
        alpha,
        sin2_structural,
        sin2_target_mz,
        shift_structural,
        coeff_structural,
        coeff_required,
        coeff_rel_delta,
        sin2_structural_mz,
        sin2_structural_abs_err,
    }
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
    idx.sort_by(|&i, &j| {
        a[i][i]
            .re
            .partial_cmp(&a[j][j].re)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let evals = [
        a[idx[0]][idx[0]].re,
        a[idx[1]][idx[1]].re,
        a[idx[2]][idx[2]].re,
    ];
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

/// Sorted neutrino texture eigenvalues (ascending) from the Cl(1,3) PMNS texture lane.
pub fn neutrino_texture_eigenvalues() -> [f64; 3] {
    let (_ml, mnu) = pmns_mass_textures_from_clifford();
    let (evals, _un) = jacobi_eigen_hermitian(mnu);
    evals
}

/// Structural hierarchy exponent from Cl(1,3) counts.
///
/// `p = α^{-1} / (|grade₁| + |grade₂|) = 137 / 10`.
pub fn neutrino_hierarchy_exponent_structural() -> f64 {
    ALPHA_INVERSE_STRUCTURAL / (GRADE1_DIM + GRADE2_DIM)
}

/// Absolute neutrino masses (eV) from texture lane + structural suppression.
///
/// Mass scale:
/// `m_scale = m_e * α^4 * (60/11)`.
///
/// Hierarchy mapping:
/// `m_i = m_scale * (|λ_i|/|λ_max|)^p`, `p = 137/10`, with `m_3 = m_scale`.
pub fn neutrino_absolute_masses_from_texture() -> NeutrinoAbsoluteMasses {
    let mut raw = neutrino_texture_eigenvalues().map(|x| x.abs());
    raw.sort_by(|a, b| a.total_cmp(b));

    let alpha_physical = 1.0 / ALPHA_INVERSE_PHYSICAL;
    let electron_mass_anchor_ev = crate::weak::electron_mass_from_proton_anchor() * 1.0e6;
    let mass_scale_ev = electron_mass_anchor_ev * alpha_physical.powi(4) * (60.0 / 11.0);
    let hierarchy_exponent = neutrino_hierarchy_exponent_structural();
    let raw_max = raw[2].max(1.0e-18);

    let m1_ev = mass_scale_ev * (raw[0] / raw_max).powf(hierarchy_exponent);
    let m2_ev = mass_scale_ev * (raw[1] / raw_max).powf(hierarchy_exponent);
    let m3_ev = mass_scale_ev;
    let sum_ev = m1_ev + m2_ev + m3_ev;

    let dm21_ev2 = m2_ev * m2_ev - m1_ev * m1_ev;
    let dm32_ev2 = m3_ev * m3_ev - m2_ev * m2_ev;
    let dm31_ev2 = m3_ev * m3_ev - m1_ev * m1_ev;
    let splitting_ratio_32_over_21 = dm32_ev2.abs() / dm21_ev2.abs().max(1.0e-30);

    NeutrinoAbsoluteMasses {
        alpha_physical,
        electron_mass_anchor_ev,
        mass_scale_ev,
        hierarchy_exponent,
        m1_ev,
        m2_ev,
        m3_ev,
        sum_ev,
        dm21_ev2,
        dm32_ev2,
        dm31_ev2,
        splitting_ratio_32_over_21,
    }
}

/// Hierarchy prediction from the texture eigenvalue ordering.
///
/// Returns `"normal"` when m1 < m2 < m3 (equivalently Δm^2_31 > 0),
/// `"degenerate"` if eigenvalues are nearly equal, else `"inverted_like"`.
pub fn neutrino_hierarchy_prediction() -> &'static str {
    let m = neutrino_texture_eigenvalues();
    let eps = 1.0e-12;
    if (m[2] - m[0]).abs() < eps {
        return "degenerate";
    }
    if m[0] < m[1] && m[1] < m[2] {
        "normal"
    } else {
        "inverted_like"
    }
}

/// Majorana-symmetry residual for the neutrino texture lane.
///
/// A pure left-handed Majorana mass matrix is complex-symmetric (`M = M^T`).
/// This metric reports the largest entrywise deviation from that symmetry.
pub fn neutrino_majorana_symmetry_residual() -> f64 {
    let (_ml, mnu) = pmns_mass_textures_from_clifford();
    let mut max_resid = 0.0f64;
    for i in 0..3 {
        for j in (i + 1)..3 {
            let r = (mnu[i][j] - mnu[j][i]).norm();
            if r > max_resid {
                max_resid = r;
            }
        }
    }
    max_resid
}

/// Structural neutrino mass-character prediction from the texture lane.
///
/// Returns:
/// - `"majorana_like"` if the texture is symmetric to tolerance
/// - `"dirac"` otherwise (the lane is Hermitian, not symmetric)
pub fn neutrino_dirac_majorana_prediction() -> &'static str {
    let residual = neutrino_majorana_symmetry_residual();
    if residual <= 1.0e-12 {
        "majorana_like"
    } else {
        "dirac"
    }
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

/// Direct structural value used by the PMNS lane:
///   sin²(theta23) = 4/7.
pub fn pmns_theta23_sq_direct() -> f64 {
    GRADE1_DIM / LATTICE_SHIFT
}

/// Corrected structural value used by the PMNS alpha² lane:
///   sin²(theta23) = 4/7 - c_alpha2 * alpha_structural^2.
pub fn pmns_theta23_sq_alpha2_corrected(c_alpha2: f64) -> f64 {
    let alpha_structural = 1.0 / ALPHA_INVERSE_STRUCTURAL;
    pmns_theta23_sq_direct() - c_alpha2 * alpha_structural * alpha_structural
}

/// PMNS observables with an optional second-order correction in the theta23 lane.
///
/// Correction ansatz:
///   sin²(theta23) = 4/7 - c_alpha2 * alpha^2
///
/// where `alpha` is the structural electromagnetic coupling (`1/137` from Lean),
/// and `c_alpha2` is a structural-rational coefficient candidate.
pub fn pmns_from_clifford_theta23_alpha2(c_alpha2: f64) -> MixingObservables {
    let s12 = (GRADE1_DIM / COMPLEMENT_DIM).sqrt();
    let s23_sq = pmns_theta23_sq_alpha2_corrected(c_alpha2);
    let s23 = s23_sq.clamp(0.0, 1.0).sqrt();
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
        ("jarlskog", obs.jarlskog, env.jarlskog_min, env.jarlskog_max),
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
        assert!(
            r.d_theta12_deg.abs() < 1.0,
            "theta12 drift too large: {}",
            r.d_theta12_deg
        );
        assert!(
            r.d_theta23_deg.abs() < 0.4,
            "theta23 drift too large: {}",
            r.d_theta23_deg
        );
        assert!(
            r.d_theta13_deg.abs() < 0.08,
            "theta13 drift too large: {}",
            r.d_theta13_deg
        );
        assert!(
            r.d_delta_deg.abs() < 2.0,
            "delta drift too large: {}",
            r.d_delta_deg
        );
        assert!(
            r.d_jarlskog.abs() < 5e-6,
            "J drift too large: {}",
            r.d_jarlskog
        );
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
        assert!(
            r.d_theta12_deg.abs() < 2.0,
            "theta12 drift too large: {}",
            r.d_theta12_deg
        );
        assert!(
            r.d_theta23_deg.abs() < 1.0,
            "theta23 drift too large: {}",
            r.d_theta23_deg
        );
        assert!(
            r.d_theta13_deg.abs() < 0.5,
            "theta13 drift too large: {}",
            r.d_theta13_deg
        );
        assert!(
            r.d_delta_deg.abs() < 30.0,
            "delta drift too large: {}",
            r.d_delta_deg
        );
        assert!(
            r.d_jarlskog.abs() < 2e-5,
            "J drift too large: {}",
            r.d_jarlskog
        );
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

    #[test]
    fn pmns_theta23_alpha2_default_coeff_matches_structural_lane() {
        let alpha_structural = 1.0 / ALPHA_INVERSE_STRUCTURAL;
        let expected = (GRADE1_DIM / LATTICE_SHIFT)
            - PMNS_THETA23_ALPHA2_COEFF_STRUCTURAL * alpha_structural * alpha_structural;
        let got = pmns_theta23_sq_alpha2_corrected(PMNS_THETA23_ALPHA2_COEFF_STRUCTURAL);
        assert!(
            (got - expected).abs() < 1e-15,
            "structural theta23 alpha2 lane mismatch: got={got:.15e}, expected={expected:.15e}"
        );
    }

    #[test]
    fn pmns_theta23_alpha2_improves_over_direct_by_10x() {
        let direct = pmns_from_clifford();
        let corrected = pmns_from_clifford_theta23_alpha2(PMNS_THETA23_ALPHA2_COEFF_STRUCTURAL);
        let direct_resid = (direct.theta23_deg - PMNS_TARGET.theta23_deg).abs();
        let corr_resid = (corrected.theta23_deg - PMNS_TARGET.theta23_deg).abs();
        assert!(
            corr_resid <= direct_resid / 10.0,
            "theta23 alpha2 correction did not improve by 10x: direct={direct_resid:.9}, corrected={corr_resid:.9}"
        );
        within_envelope(corrected, PMNS_PDG_ENVELOPE)
            .expect("corrected PMNS theta23 alpha2 outside PDG envelope");
    }

    #[test]
    fn neutrino_hierarchy_is_normal_in_texture_lane() {
        assert_eq!(neutrino_hierarchy_prediction(), "normal");
    }

    #[test]
    fn neutrino_texture_prefers_dirac_over_majorana_symmetry() {
        let resid = neutrino_majorana_symmetry_residual();
        assert!(
            resid > 1.0e-6,
            "majorana residual too small for Dirac lane claim: {resid:.12e}"
        );
        assert_eq!(neutrino_dirac_majorana_prediction(), "dirac");
    }

    #[test]
    fn neutrino_hierarchy_exponent_is_structural_137_over_10() {
        let p = neutrino_hierarchy_exponent_structural();
        assert!((p - 13.7).abs() < 1e-12, "unexpected hierarchy exponent: {p:.12}");
    }

    #[test]
    fn neutrino_structural_splitting_ratio_is_close_to_oscillation_target() {
        let abs = neutrino_absolute_masses_from_texture();
        let target = 2.453e-3 / 7.53e-5;
        let rel_err = (abs.splitting_ratio_32_over_21 - target) / target;
        assert!(
            rel_err.abs() < 0.05,
            "structural splitting-ratio drift too large: got {:.9}, target {:.9}, rel_err {:.6}",
            abs.splitting_ratio_32_over_21,
            target,
            rel_err
        );
    }

    #[test]
    fn triangulated_exponent_solves_ratio_to_machine_precision() {
        let tri = triangulate_neutrino_from_splittings(7.53e-5, 2.453e-3);
        assert!(
            tri.ratio_fit_rel_err.abs() < 1e-9,
            "triangulated ratio not closed: fit={} target={} rel_err={:.3e}",
            tri.ratio_fit,
            tri.ratio_target,
            tri.ratio_fit_rel_err
        );
    }

    #[test]
    fn triangulated_mass_scale_reconstructs_absolute_splittings() {
        let tri = triangulate_neutrino_from_splittings(7.53e-5, 2.453e-3);
        let dm21_rel = (tri.dm21_ev2 - 7.53e-5) / 7.53e-5;
        let dm32_rel = (tri.dm32_ev2 - 2.453e-3) / 2.453e-3;
        assert!(
            dm21_rel.abs() < 1e-9 && dm32_rel.abs() < 1e-9,
            "triangulated absolute closure drift: dm21_rel={:.3e} dm32_rel={:.3e}",
            dm21_rel,
            dm32_rel
        );
    }

    #[test]
    fn ew_shift_triangulation_exposes_required_uplift_over_structural() {
        let ew = triangulate_ew_shift_for_target(0.23122);
        assert!(ew.coeff_required > ew.coeff_structural);
        assert!(ew.coeff_rel_delta > 0.0);
        assert!(
            (ew.coeff_required - 8.460487692308).abs() < 1e-9,
            "unexpected EW required coeff: {:.12}",
            ew.coeff_required
        );
    }
}
