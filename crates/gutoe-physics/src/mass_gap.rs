/*!
 * GUTOE Physics - Yang-Mills Mass Gap Scaffolding
 * Copyright (C) 2026  Riff Labs
 *
 * This module provides finite-volume transfer-matrix diagnostics used by
 * GRAND-297/298:
 * - positivity checks (entrywise and conservative Gershgorin SPD bound)
 * - spectral-gap observable m_gap(L, a_t) from λ0, λ1 of transfer matrix T
 * - conservative lower-bound extraction using eigenpair residual intervals
 * - monotone-in-volume checks for finite-volume trend analysis
 */

#[derive(Debug, Clone)]
pub struct DenseSymmetricMatrix {
    dim: usize,
    data: Vec<f64>, // Row-major
}

impl DenseSymmetricMatrix {
    pub fn from_rows(rows: &[Vec<f64>]) -> Option<Self> {
        let dim = rows.len();
        if dim == 0 || rows.iter().any(|r| r.len() != dim) {
            return None;
        }
        let mut data = Vec::with_capacity(dim * dim);
        for r in rows {
            data.extend_from_slice(r);
        }
        Some(Self { dim, data })
    }

    pub fn dim(&self) -> usize {
        self.dim
    }

    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.data[i * self.dim + j]
    }

    pub fn is_symmetric(&self, tol: f64) -> bool {
        for i in 0..self.dim {
            for j in (i + 1)..self.dim {
                if (self.get(i, j) - self.get(j, i)).abs() > tol {
                    return false;
                }
            }
        }
        true
    }

    pub fn is_entrywise_nonnegative(&self, tol: f64) -> bool {
        self.data.iter().all(|x| *x >= -tol)
    }

    /// Conservative lower bound on min eigenvalue via Gershgorin discs.
    /// If this is >0 for a symmetric matrix, matrix is SPD.
    pub fn gershgorin_lower_bound(&self) -> f64 {
        let mut lb = f64::INFINITY;
        for i in 0..self.dim {
            let aii = self.get(i, i);
            let mut radius = 0.0;
            for j in 0..self.dim {
                if i != j {
                    radius += self.get(i, j).abs();
                }
            }
            lb = lb.min(aii - radius);
        }
        lb
    }

    fn mat_vec(&self, v: &[f64], out: &mut [f64]) {
        for (i, out_i) in out.iter_mut().enumerate().take(self.dim) {
            let mut s = 0.0;
            for (j, vj) in v.iter().enumerate().take(self.dim) {
                s += self.get(i, j) * *vj;
            }
            *out_i = s;
        }
    }

    fn dot(a: &[f64], b: &[f64]) -> f64 {
        a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
    }

    fn norm(v: &[f64]) -> f64 {
        Self::dot(v, v).sqrt()
    }

    fn normalize(v: &mut [f64]) -> bool {
        let n = Self::norm(v);
        if n <= 0.0 || !n.is_finite() {
            return false;
        }
        for x in v.iter_mut() {
            *x /= n;
        }
        true
    }

    fn rayleigh(&self, v: &[f64], work: &mut [f64]) -> f64 {
        self.mat_vec(v, work);
        Self::dot(v, work)
    }

    fn residual_norm(&self, v: &[f64], lambda: f64, work: &mut [f64]) -> f64 {
        self.mat_vec(v, work);
        for i in 0..self.dim {
            work[i] -= lambda * v[i];
        }
        Self::norm(work)
    }

    /// Largest eigenpair estimate from power iteration.
    pub fn largest_eigenpair_power(&self, max_iters: usize, tol: f64) -> Option<EigenEstimate> {
        if self.dim == 0 {
            return None;
        }
        let mut v = vec![1.0 / (self.dim as f64).sqrt(); self.dim];
        let mut w = vec![0.0; self.dim];
        let mut lambda_prev = f64::NEG_INFINITY;

        for _ in 0..max_iters {
            self.mat_vec(&v, &mut w);
            if !Self::normalize(&mut w) {
                return None;
            }
            v.clone_from_slice(&w);
            let lambda = self.rayleigh(&v, &mut w);
            let resid = self.residual_norm(&v, lambda, &mut w);
            if (lambda - lambda_prev).abs() < tol && resid < tol {
                return Some(EigenEstimate {
                    value: lambda,
                    vector: v,
                    residual: resid,
                });
            }
            lambda_prev = lambda;
        }

        let lambda = self.rayleigh(&v, &mut w);
        let resid = self.residual_norm(&v, lambda, &mut w);
        Some(EigenEstimate {
            value: lambda,
            vector: v,
            residual: resid,
        })
    }

    /// Second eigenvalue estimate via deflated power iteration (orthogonal to v1).
    pub fn second_eigenvalue_deflated(
        &self,
        v1: &[f64],
        max_iters: usize,
        tol: f64,
    ) -> Option<EigenEstimate> {
        if self.dim < 2 || v1.len() != self.dim {
            return None;
        }
        let mut v = vec![0.0; self.dim];
        v[0] = 1.0;
        v[1] = -1.0;
        let proj0 = Self::dot(&v, v1);
        for i in 0..self.dim {
            v[i] -= proj0 * v1[i];
        }
        if !Self::normalize(&mut v) {
            return None;
        }

        let mut w = vec![0.0; self.dim];
        let mut lambda_prev = f64::NEG_INFINITY;
        for _ in 0..max_iters {
            self.mat_vec(&v, &mut w);
            let proj = Self::dot(&w, v1);
            for i in 0..self.dim {
                w[i] -= proj * v1[i];
            }
            if !Self::normalize(&mut w) {
                return None;
            }
            v.clone_from_slice(&w);
            let lambda = self.rayleigh(&v, &mut w);
            let resid = self.residual_norm(&v, lambda, &mut w);
            if (lambda - lambda_prev).abs() < tol && resid < tol {
                return Some(EigenEstimate {
                    value: lambda,
                    vector: v,
                    residual: resid,
                });
            }
            lambda_prev = lambda;
        }

        let lambda = self.rayleigh(&v, &mut w);
        let resid = self.residual_norm(&v, lambda, &mut w);
        Some(EigenEstimate {
            value: lambda,
            vector: v,
            residual: resid,
        })
    }
}

#[derive(Debug, Clone)]
pub struct EigenEstimate {
    pub value: f64,
    pub vector: Vec<f64>,
    pub residual: f64,
}

#[derive(Debug, Clone)]
pub struct TransferMatrixGapEstimate {
    pub lambda0_est: f64,
    pub lambda1_est: f64,
    pub lambda0_residual: f64,
    pub lambda1_residual: f64,
    pub gap_est: f64,
    pub gap_lower_bound: Option<f64>,
}

#[derive(Debug, Clone)]
pub struct TransferMatrixDiagnostics {
    pub dim: usize,
    pub symmetric: bool,
    pub entrywise_nonnegative: bool,
    pub gershgorin_lower_bound: f64,
    pub gap: Option<TransferMatrixGapEstimate>,
}

pub fn transfer_matrix_diagnostics(
    t: &DenseSymmetricMatrix,
    a_t: f64,
    max_iters: usize,
    tol: f64,
) -> Option<TransferMatrixDiagnostics> {
    if a_t <= 0.0 {
        return None;
    }
    let symmetric = t.is_symmetric(tol);
    let entrywise_nonnegative = t.is_entrywise_nonnegative(tol);
    let gersh_lb = t.gershgorin_lower_bound();

    let e0 = t.largest_eigenpair_power(max_iters, tol)?;
    let e1 = t.second_eigenvalue_deflated(&e0.vector, max_iters, tol)?;

    let lambda0 = e0.value.max(0.0);
    let lambda1 = e1.value.max(0.0);
    let ratio = lambda1 / lambda0;
    let gap_est = if lambda0 > 0.0 && lambda1 > 0.0 && ratio < 1.0 {
        -ratio.ln() / a_t
    } else {
        f64::NAN
    };

    // Conservative lower bound using residual intervals:
    // λ0 ∈ [λ0_est-r0, λ0_est+r0], λ1 ∈ [λ1_est-r1, λ1_est+r1].
    // If λ1_ub < λ0_lb then m_gap ≥ -(1/a_t) ln(λ1_ub/λ0_lb).
    let lambda0_lb = (lambda0 - e0.residual).max(0.0);
    let lambda1_ub = (lambda1 + e1.residual).max(0.0);
    let gap_lower_bound = if lambda0_lb > 0.0 && lambda1_ub > 0.0 && lambda1_ub < lambda0_lb {
        Some(-(lambda1_ub / lambda0_lb).ln() / a_t)
    } else {
        None
    };

    let gap = if gap_est.is_finite() {
        Some(TransferMatrixGapEstimate {
            lambda0_est: lambda0,
            lambda1_est: lambda1,
            lambda0_residual: e0.residual,
            lambda1_residual: e1.residual,
            gap_est,
            gap_lower_bound,
        })
    } else {
        None
    };

    Some(TransferMatrixDiagnostics {
        dim: t.dim(),
        symmetric,
        entrywise_nonnegative,
        gershgorin_lower_bound: gersh_lb,
        gap,
    })
}

#[derive(Debug, Clone, Copy)]
pub struct VolumeGapPoint {
    pub volume_l3: usize,
    pub gap_est: f64,
    pub gap_err: f64,
}

/// Check the expected finite-volume monotone trend for gap estimates:
/// as volume grows, finite-size gaps should not increase materially.
pub fn monotone_nonincreasing_in_volume(points: &[VolumeGapPoint], tol: f64) -> bool {
    if points.len() < 2 {
        return true;
    }
    for w in points.windows(2) {
        if w[1].volume_l3 <= w[0].volume_l3 {
            return false;
        }
        if w[1].gap_est > w[0].gap_est + tol {
            return false;
        }
    }
    true
}

/// Conservative continuum-stability envelope from finite-volume estimates.
pub fn continuum_stability_band(points: &[VolumeGapPoint]) -> Option<(f64, f64)> {
    if points.is_empty() {
        return None;
    }
    let mut lo = f64::NEG_INFINITY;
    let mut hi = f64::INFINITY;
    for p in points {
        lo = lo.max(p.gap_est - p.gap_err);
        hi = hi.min(p.gap_est + p.gap_err);
    }
    if lo.is_finite() && hi.is_finite() && lo <= hi {
        Some((lo.max(0.0), hi.max(0.0)))
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn transfer_matrix_gap_and_lower_bound_on_toy_diagonal() {
        let m_gap: f64 = 0.7;
        let m_excited: f64 = 1.4;
        let rows = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, (-m_gap).exp(), 0.0],
            vec![0.0, 0.0, (-m_excited).exp()],
        ];
        let t = DenseSymmetricMatrix::from_rows(&rows).expect("matrix");
        let d = transfer_matrix_diagnostics(&t, 1.0, 5000, 1e-12).expect("diagnostics");
        assert!(d.symmetric);
        assert!(d.entrywise_nonnegative);
        assert!(d.gershgorin_lower_bound > 0.0);

        let g = d.gap.expect("gap estimate");
        assert!((g.gap_est - m_gap).abs() < 5e-4);
        let lb = g.gap_lower_bound.expect("lower bound");
        assert!(lb > 0.0);
        assert!(lb <= g.gap_est + 1e-6);
    }

    #[test]
    fn volume_monotone_and_band_checks() {
        let pts = [
            VolumeGapPoint {
                volume_l3: 8 * 8 * 8,
                gap_est: 1.20,
                gap_err: 0.20,
            },
            VolumeGapPoint {
                volume_l3: 10 * 10 * 10,
                gap_est: 1.08,
                gap_err: 0.15,
            },
            VolumeGapPoint {
                volume_l3: 12 * 12 * 12,
                gap_est: 1.01,
                gap_err: 0.12,
            },
        ];
        assert!(monotone_nonincreasing_in_volume(&pts, 1e-9));
        let (lo, hi) = continuum_stability_band(&pts).expect("band");
        assert!(lo <= hi);
        assert!(lo >= 0.0);
    }
}
