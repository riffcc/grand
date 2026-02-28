// GUTOE EM — Holonomy diagnostics from lattice parallel transport
// Copyright (C) 2026 Riff Labs, AGPL-3.0-or-later
//
// GRAND-216:
//   - closed-loop parallel transport (holonomy) on the Cl(1,3) lattice
//   - restricted holonomy signature U(1) x SU(2) x SU(3) from structural counts
//   - explicit Wilson-plaquette consistency check for triangular loops
//   - geometric phase (Berry/U(1)) composition diagnostics

use std::collections::HashSet;

use num_complex::Complex64;
use serde::{Deserialize, Serialize};

use crate::config::LatticeConfig;
use crate::geometry::{mesh_neighbours, site_coords};
use crate::su2_gauge::{su2_identity, su2_mul, su2_re_tr, Su2, Su2Links};

/// Structural restricted-holonomy signature from the Clifford/Z3 lane.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct RestrictedHolonomySignature {
    pub z3_order: u32,
    pub u1_generators: u32,
    pub su2_generators: u32,
    pub su3_generators: u32,
    pub total_generators: u32,
}

impl RestrictedHolonomySignature {
    /// U(1) x SU(2) x SU(3) generator counts from the same structural map used
    /// across the physics crate.
    pub fn from_clifford_z3() -> Self {
        let z3_order = 3_u32;
        let u1_generators = 1_u32;
        let su2_generators = 3_u32; // |magneticTriplet|
        let su3_generators = z3_order * z3_order - 1; // 8
        Self {
            z3_order,
            u1_generators,
            su2_generators,
            su3_generators,
            total_generators: u1_generators + su2_generators + su3_generators,
        }
    }

    pub fn recovers_sm(&self) -> bool {
        self.u1_generators == 1
            && self.su2_generators == 3
            && self.su3_generators == 8
            && self.total_generators == 12
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct TriangleHolonomySample {
    pub i: usize,
    pub j: usize,
    pub k: usize,
    pub trace_over_2: f64,
    pub class_angle_rad: f64,
    pub wilson_residual_abs: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct HolonomyDiagnostics {
    pub samples: Vec<TriangleHolonomySample>,
    pub max_wilson_residual_abs: f64,
    pub mean_trace_over_2: f64,
    pub mean_class_angle_rad: f64,
}

fn clamp_unit(x: f64) -> f64 {
    x.clamp(-1.0, 1.0)
}

/// SU(2) conjugacy class angle: trace(U)/2 = cos(theta).
pub fn class_angle_from_trace(trace_over_2: f64) -> f64 {
    clamp_unit(trace_over_2).acos()
}

/// U(1) geometric phase e^{i theta}.
pub fn u1_geometric_phase(theta: f64) -> Complex64 {
    Complex64::new(theta.cos(), theta.sin())
}

/// Composition residual for U(1) geometric phase:
/// exp(i(theta_a + theta_b)) ?= exp(i theta_a) exp(i theta_b).
pub fn u1_phase_composition_residual(theta_a: f64, theta_b: f64) -> f64 {
    let lhs = u1_geometric_phase(theta_a + theta_b);
    let rhs = u1_geometric_phase(theta_a) * u1_geometric_phase(theta_b);
    (lhs - rhs).norm()
}

/// Parallel transport product along an oriented path.
pub fn transport_product(links: &Su2Links, path: &[usize]) -> Option<Su2> {
    if path.len() < 2 {
        return None;
    }
    let mut accum = su2_identity();
    for edge in path.windows(2) {
        accum = su2_mul(&accum, &links.get(edge[0], edge[1]));
    }
    Some(accum)
}

/// Holonomy around a closed loop encoded by vertices with repeated start.
pub fn closed_loop_holonomy(links: &Su2Links, loop_vertices: &[usize]) -> Option<Su2> {
    if loop_vertices.len() < 4 {
        return None;
    }
    if loop_vertices.first() != loop_vertices.last() {
        return None;
    }
    transport_product(links, loop_vertices)
}

/// Triangle-loop holonomy for oriented loop i -> j -> k -> i.
pub fn triangle_loop_holonomy(links: &Su2Links, i: usize, j: usize, k: usize) -> Su2 {
    closed_loop_holonomy(links, &[i, j, k, i]).unwrap_or_else(su2_identity)
}

/// Normalized trace for triangle holonomy.
pub fn triangle_loop_trace_over_2(links: &Su2Links, i: usize, j: usize, k: usize) -> f64 {
    su2_re_tr(&triangle_loop_holonomy(links, i, j, k)) / 2.0
}

/// Consistency residual against the direct Wilson plaquette implementation.
pub fn triangle_wilson_residual_abs(links: &Su2Links, i: usize, j: usize, k: usize) -> f64 {
    let via_transport = triangle_loop_trace_over_2(links, i, j, k);
    let via_plaquette = links.plaquette_triangle(i, j, k);
    (via_transport - via_plaquette).abs()
}

/// Enumerate unique triangles (i < j < k) in the lattice mesh.
pub fn enumerate_triangles(cfg: &LatticeConfig) -> Vec<[usize; 3]> {
    let n = cfg.n_sites();
    let mut triangles = Vec::new();
    for i in 0..n {
        let (ri, ci, zi) = site_coords(i, cfg);
        let nbrs_i = mesh_neighbours(ri, ci, zi, cfg);
        for &j in &nbrs_i {
            if j <= i {
                continue;
            }
            let (rj, cj, zj) = site_coords(j, cfg);
            let nbrs_j = mesh_neighbours(rj, cj, zj, cfg);
            let nbrs_i_set: HashSet<usize> = nbrs_i.iter().copied().collect();
            for &k in &nbrs_j {
                if k <= j {
                    continue;
                }
                if nbrs_i_set.contains(&k) {
                    triangles.push([i, j, k]);
                }
            }
        }
    }
    triangles
}

/// Sample triangle-loop holonomy diagnostics.
pub fn sample_holonomy_diagnostics(
    links: &Su2Links,
    cfg: &LatticeConfig,
    max_samples: usize,
) -> HolonomyDiagnostics {
    let triangles = enumerate_triangles(cfg);
    if triangles.is_empty() {
        return HolonomyDiagnostics {
            samples: Vec::new(),
            max_wilson_residual_abs: 0.0,
            mean_trace_over_2: 0.0,
            mean_class_angle_rad: 0.0,
        };
    }

    let n_take = triangles.len().min(max_samples.max(1));
    let mut samples = Vec::with_capacity(n_take);
    let mut trace_sum = 0.0;
    let mut angle_sum = 0.0;
    let mut max_residual: f64 = 0.0;

    for &[i, j, k] in triangles.iter().take(n_take) {
        let trace_over_2 = triangle_loop_trace_over_2(links, i, j, k);
        let class_angle_rad = class_angle_from_trace(trace_over_2);
        let wilson_residual_abs = triangle_wilson_residual_abs(links, i, j, k);
        max_residual = max_residual.max(wilson_residual_abs);
        trace_sum += trace_over_2;
        angle_sum += class_angle_rad;
        samples.push(TriangleHolonomySample {
            i,
            j,
            k,
            trace_over_2,
            class_angle_rad,
            wilson_residual_abs,
        });
    }

    let denom = n_take as f64;
    HolonomyDiagnostics {
        samples,
        max_wilson_residual_abs: max_residual,
        mean_trace_over_2: trace_sum / denom,
        mean_class_angle_rad: angle_sum / denom,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::LatticeConfig;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    fn small_cfg() -> LatticeConfig {
        LatticeConfig {
            hex_rows: 8,
            hex_cols: 8,
            layers: 1,
            ..Default::default()
        }
    }

    #[test]
    fn restricted_holonomy_signature_matches_sm_counts() {
        let sig = RestrictedHolonomySignature::from_clifford_z3();
        assert_eq!(sig.u1_generators, 1);
        assert_eq!(sig.su2_generators, 3);
        assert_eq!(sig.su3_generators, 8);
        assert_eq!(sig.total_generators, 12);
        assert!(sig.recovers_sm());
    }

    #[test]
    fn triangle_holonomy_matches_wilson_plaquette() {
        let cfg = small_cfg();
        let mut rng = StdRng::seed_from_u64(137);
        let links = Su2Links::hot_start(&mut rng, &cfg);
        let tris = enumerate_triangles(&cfg);
        assert!(!tris.is_empty(), "expected at least one triangle");
        let [i, j, k] = tris[0];

        let residual = triangle_wilson_residual_abs(&links, i, j, k);
        assert!(
            residual < 1e-12,
            "triangle holonomy trace must match Wilson plaquette, residual={residual:.3e}"
        );
    }

    #[test]
    fn u1_geometric_phase_is_unitary_and_composes() {
        let theta_a = std::f64::consts::PI / 7.0;
        let theta_b = -std::f64::consts::PI / 11.0;
        let phase = u1_geometric_phase(theta_a);
        assert!((phase.norm() - 1.0).abs() < 1e-12);
        let residual = u1_phase_composition_residual(theta_a, theta_b);
        assert!(
            residual < 1e-12,
            "phase composition residual too large: {residual:.3e}"
        );
    }

    #[test]
    fn closed_loop_requires_explicit_closure() {
        let cfg = small_cfg();
        let links = Su2Links::cold_start(&cfg);
        let open_path = [0usize, 1usize, 2usize];
        assert!(closed_loop_holonomy(&links, &open_path).is_none());
    }
}
