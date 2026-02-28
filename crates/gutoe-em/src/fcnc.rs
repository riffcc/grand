// GUTOE EM — FCNC suppression and GIM diagnostics from CKM mixing structure
//
// GRAND-132:
//   Why are flavor-changing neutral currents (FCNC) suppressed?
//
// This module quantifies two linked facts from the CKM mixing map:
//   1) Tree-level neutral currents are flavor-diagonal (unitarity cancellation).
//   2) Loop-level FCNC amplitudes keep only flavor-difference pieces
//      (GIM cancellation), yielding strong suppression.

use crate::flavor::{ckm_from_clifford, ckm_from_textures, MixingObservables};
use num_complex::Complex64;
use serde::Serialize;

const DOWN_FLAVORS: [&str; 3] = ["d", "s", "b"];
const UP_FLAVORS: [&str; 3] = ["u", "c", "t"];

/// Exact structural small-angle product from the CKM lane:
/// s23 * s13 = (1/24)*(1/272) = 1/6528.
pub const FCNC_LOOP_PROXY_EXPECTED: f64 = 1.0 / 6528.0;

type CMat3 = [[Complex64; 3]; 3];

#[derive(Debug, Clone, Copy, Serialize)]
pub struct GimChannelMetrics {
    pub from: &'static str,
    pub to: &'static str,
    pub lambda_u_abs: f64,
    pub lambda_c_abs: f64,
    pub lambda_t_abs: f64,
    pub lambda_sum_abs: f64,
    pub degenerate_kernel_abs: f64,
    pub split_kernel_abs: f64,
    pub split_kernel_no_gim_abs: f64,
    pub gim_suppression_ratio: f64,
    pub mass_difference_form_residual_abs: f64,
}

#[derive(Debug, Clone, Copy, Serialize)]
pub struct FcncGimMetrics {
    pub source: &'static str,
    pub ckm_s23: f64,
    pub ckm_s13: f64,
    pub structural_loop_proxy: f64,
    pub neutral_current_offdiag_max_abs: f64,
    pub neutral_current_diag_drift_max_abs: f64,
    pub gim_sum_rule_residual_max_abs: f64,
    pub loop_kernel_u: f64,
    pub loop_kernel_c: f64,
    pub loop_kernel_t: f64,
    pub channels: [GimChannelMetrics; 3],
}

fn ckm_unitary_from_observables(obs: MixingObservables) -> CMat3 {
    let s12 = obs.s12.clamp(0.0, 1.0);
    let s23 = obs.s23.clamp(0.0, 1.0);
    let s13 = obs.s13.clamp(0.0, 1.0);
    let c12 = (1.0 - s12 * s12).max(0.0).sqrt();
    let c23 = (1.0 - s23 * s23).max(0.0).sqrt();
    let c13 = (1.0 - s13 * s13).max(0.0).sqrt();

    let s13_e_pos = Complex64::from_polar(s13, obs.delta_rad);
    let s13_e_neg = s13_e_pos.conj();

    [
        [
            Complex64::new(c12 * c13, 0.0),
            Complex64::new(s12 * c13, 0.0),
            s13_e_neg,
        ],
        [
            Complex64::new(-s12 * c23, 0.0) - s13_e_pos * (c12 * s23),
            Complex64::new(c12 * c23, 0.0) - s13_e_pos * (s12 * s23),
            Complex64::new(s23 * c13, 0.0),
        ],
        [
            Complex64::new(s12 * s23, 0.0) - s13_e_pos * (c12 * c23),
            Complex64::new(-c12 * s23, 0.0) - s13_e_pos * (s12 * c23),
            Complex64::new(c23 * c13, 0.0),
        ],
    ]
}

fn c_conj_transpose(m: &CMat3) -> CMat3 {
    [
        [m[0][0].conj(), m[1][0].conj(), m[2][0].conj()],
        [m[0][1].conj(), m[1][1].conj(), m[2][1].conj()],
        [m[0][2].conj(), m[1][2].conj(), m[2][2].conj()],
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

fn neutral_current_unitarity_metrics(v: &CMat3) -> (f64, f64) {
    // Neutral current kernel in the mass basis: V†V.
    let n = c_mul(&c_conj_transpose(v), v);
    let mut offdiag_max = 0.0f64;
    let mut diag_drift_max = 0.0f64;

    for (i, row) in n.iter().enumerate() {
        for (j, val) in row.iter().enumerate() {
            if i == j {
                diag_drift_max = diag_drift_max.max((*val - Complex64::new(1.0, 0.0)).norm());
            } else {
                offdiag_max = offdiag_max.max(val.norm());
            }
        }
    }

    (offdiag_max, diag_drift_max)
}

fn structural_loop_kernel(obs: MixingObservables) -> [f64; 3] {
    // Dimensionless structural loop map:
    // equal term + small flavor-dependent splittings from CKM suppressions.
    [1.0, 1.0 + obs.s23, 1.0 + obs.s23 + obs.s13]
}

fn down_channel_metrics(v: &CMat3, from: usize, to: usize, kernel: [f64; 3]) -> GimChannelMetrics {
    // λ_i^{jk} = V_{ij}* V_{ik}, i in {u,c,t}, j->k in {d,s,b}.
    let lambda = std::array::from_fn::<_, 3, _>(|i| v[i][from].conj() * v[i][to]);
    let lambda_sum = lambda[0] + lambda[1] + lambda[2];

    // Degenerate loop kernel: exact GIM cancellation target.
    let degenerate_amp = lambda_sum;

    // Split kernel: only flavor differences should survive.
    let mut split_amp = Complex64::new(0.0, 0.0);
    let mut split_no_gim_abs = 0.0f64;
    for i in 0..3 {
        let term = lambda[i] * kernel[i];
        split_amp += term;
        split_no_gim_abs += term.norm();
    }

    // GIM rewrite:
    // Σ λ_i f_i = λ_c(f_c-f_u) + λ_t(f_t-f_u) if Σ λ_i = 0.
    let mass_difference_form =
        lambda[1] * (kernel[1] - kernel[0]) + lambda[2] * (kernel[2] - kernel[0]);

    let gim_suppression_ratio = if split_no_gim_abs > 0.0 {
        split_amp.norm() / split_no_gim_abs
    } else {
        0.0
    };

    GimChannelMetrics {
        from: DOWN_FLAVORS[from],
        to: DOWN_FLAVORS[to],
        lambda_u_abs: lambda[0].norm(),
        lambda_c_abs: lambda[1].norm(),
        lambda_t_abs: lambda[2].norm(),
        lambda_sum_abs: lambda_sum.norm(),
        degenerate_kernel_abs: degenerate_amp.norm(),
        split_kernel_abs: split_amp.norm(),
        split_kernel_no_gim_abs: split_no_gim_abs,
        gim_suppression_ratio,
        mass_difference_form_residual_abs: (split_amp - mass_difference_form).norm(),
    }
}

pub fn fcnc_gim_from_observables(source: &'static str, obs: MixingObservables) -> FcncGimMetrics {
    let v = ckm_unitary_from_observables(obs);
    let (offdiag_max, diag_drift_max) = neutral_current_unitarity_metrics(&v);
    let kernel = structural_loop_kernel(obs);

    let pairs = [(0usize, 1usize), (0usize, 2usize), (1usize, 2usize)];
    let channels = pairs.map(|(from, to)| down_channel_metrics(&v, from, to, kernel));
    let gim_sum_rule_residual_max_abs = channels
        .iter()
        .map(|c| c.lambda_sum_abs)
        .fold(0.0f64, f64::max);

    FcncGimMetrics {
        source,
        ckm_s23: obs.s23,
        ckm_s13: obs.s13,
        structural_loop_proxy: obs.s23 * obs.s13,
        neutral_current_offdiag_max_abs: offdiag_max,
        neutral_current_diag_drift_max_abs: diag_drift_max,
        gim_sum_rule_residual_max_abs,
        loop_kernel_u: kernel[0],
        loop_kernel_c: kernel[1],
        loop_kernel_t: kernel[2],
        channels,
    }
}

pub fn fcnc_gim_from_clifford() -> FcncGimMetrics {
    fcnc_gim_from_observables("ckm_direct_clifford", ckm_from_clifford())
}

pub fn fcnc_gim_from_textures() -> FcncGimMetrics {
    fcnc_gim_from_observables("ckm_texture_diagonalization", ckm_from_textures())
}

pub fn ckm_structural_loop_proxy() -> f64 {
    let ckm = ckm_from_clifford();
    ckm.s23 * ckm.s13
}

pub fn ckm_structural_loop_proxy_matches_expected(tol: f64) -> bool {
    (ckm_structural_loop_proxy() - FCNC_LOOP_PROXY_EXPECTED).abs() <= tol.max(0.0)
}

pub fn channel_label(ch: &GimChannelMetrics) -> &'static str {
    match (ch.from, ch.to) {
        ("d", "s") => "d<->s",
        ("d", "b") => "d<->b",
        ("s", "b") => "s<->b",
        _ => "?",
    }
}

pub fn up_flavors() -> [&'static str; 3] {
    UP_FLAVORS
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn structural_loop_proxy_matches_clifford_fraction() {
        let got = ckm_structural_loop_proxy();
        assert!(
            (got - FCNC_LOOP_PROXY_EXPECTED).abs() < 1.0e-15,
            "loop proxy mismatch: got={got:.15e}, expected={:.15e}",
            FCNC_LOOP_PROXY_EXPECTED
        );
    }

    #[test]
    fn tree_level_neutral_current_is_flavor_diagonal() {
        let direct = fcnc_gim_from_clifford();
        let tex = fcnc_gim_from_textures();
        assert!(direct.neutral_current_offdiag_max_abs < 1.0e-12);
        assert!(tex.neutral_current_offdiag_max_abs < 1.0e-12);
        assert!(direct.neutral_current_diag_drift_max_abs < 1.0e-12);
        assert!(tex.neutral_current_diag_drift_max_abs < 1.0e-12);
    }

    #[test]
    fn gim_sum_rule_cancels_degenerate_kernel() {
        let direct = fcnc_gim_from_clifford();
        let tex = fcnc_gim_from_textures();
        assert!(direct.gim_sum_rule_residual_max_abs < 1.0e-12);
        assert!(tex.gim_sum_rule_residual_max_abs < 1.0e-12);
        for ch in &direct.channels {
            assert!(ch.degenerate_kernel_abs < 1.0e-12);
        }
        for ch in &tex.channels {
            assert!(ch.degenerate_kernel_abs < 1.0e-12);
        }
    }

    #[test]
    fn mass_difference_rewrite_is_exact() {
        let direct = fcnc_gim_from_clifford();
        let tex = fcnc_gim_from_textures();
        for ch in &direct.channels {
            assert!(ch.mass_difference_form_residual_abs < 1.0e-12);
        }
        for ch in &tex.channels {
            assert!(ch.mass_difference_form_residual_abs < 1.0e-12);
        }
    }

    #[test]
    fn split_channels_show_strong_gim_suppression() {
        let direct = fcnc_gim_from_clifford();
        let tex = fcnc_gim_from_textures();
        for ch in &direct.channels {
            assert!(
                ch.gim_suppression_ratio < 0.10,
                "direct {} suppression too weak: {:.6e}",
                channel_label(ch),
                ch.gim_suppression_ratio
            );
        }
        for ch in &tex.channels {
            assert!(
                ch.gim_suppression_ratio < 0.10,
                "texture {} suppression too weak: {:.6e}",
                channel_label(ch),
                ch.gim_suppression_ratio
            );
        }
    }
}
