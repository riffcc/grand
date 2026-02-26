/*!
 * GUTOE Physics - Dark Sector Candidate Harness
 * Copyright (C) 2026  Riff Labs
 *
 * Structural GRAND-346 lane:
 *   - particle branch density from Lean-derived count ratio 5/11
 *   - geometric branch density with curvature proxy
 *   - simple rotation/lensing proxies for falsifiable reports
 */

use crate::constants::{
    C, G, DARK_TO_VISIBLE_COUNT_RATIO, DARK_TO_VISIBLE_GEOMETRIC_RATIO, LAMBDA_QG,
};

/// Dark-sector modeling branch.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DarkSectorBranch {
    Particle,
    Geometric,
}

/// Structural particle-branch dark density: ρ_dark = (5/11) * ρ_visible.
pub fn dark_density_particle(visible_density: f64) -> f64 {
    DARK_TO_VISIBLE_COUNT_RATIO * visible_density
}

/// Curvature proxy from the lattice correction scale.
///
/// κ(r) = 1 + λ_QG (r_core / r)^2, clamped to nonnegative inputs.
pub fn curvature_proxy(r: f64, r_core: f64) -> f64 {
    if r <= 0.0 || r_core <= 0.0 {
        return 1.0;
    }
    1.0 + LAMBDA_QG * (r_core / r).powi(2)
}

/// Geometric-branch effective dark density from structural amplification
/// (60/11 from Lean parity) and curvature proxy.
pub fn dark_density_geometric(visible_density: f64, curvature_factor: f64) -> f64 {
    let k = curvature_factor.max(0.0);
    DARK_TO_VISIBLE_GEOMETRIC_RATIO * k * visible_density
}

/// Branch-dispatched effective dark density.
pub fn dark_density(branch: DarkSectorBranch, visible_density: f64, curvature_factor: f64) -> f64 {
    match branch {
        DarkSectorBranch::Particle => dark_density_particle(visible_density),
        DarkSectorBranch::Geometric => dark_density_geometric(visible_density, curvature_factor),
    }
}

/// Total density = visible + dark.
pub fn total_density(branch: DarkSectorBranch, visible_density: f64, curvature_factor: f64) -> f64 {
    visible_density + dark_density(branch, visible_density, curvature_factor)
}

/// Spherical constant-density enclosed mass.
pub fn enclosed_mass_constant_density(density: f64, radius: f64) -> f64 {
    if density <= 0.0 || radius <= 0.0 {
        return 0.0;
    }
    (4.0 / 3.0) * std::f64::consts::PI * radius.powi(3) * density
}

/// Circular velocity from Newtonian balance v² = G M / r.
pub fn circular_velocity(enclosed_mass: f64, radius: f64) -> Option<f64> {
    if enclosed_mass <= 0.0 || radius <= 0.0 {
        return None;
    }
    Some((G * enclosed_mass / radius).sqrt())
}

/// Weak-field lensing deflection proxy α = 4GM/(c² b).
pub fn lensing_deflection(enclosed_mass: f64, impact_parameter: f64) -> Option<f64> {
    if enclosed_mass <= 0.0 || impact_parameter <= 0.0 {
        return None;
    }
    Some(4.0 * G * enclosed_mass / (C * C * impact_parameter))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::{
        DARK_FRACTION_GEOMETRIC_STRUCTURAL, DARK_FRACTION_TOTAL_STATE_SPLIT,
        DARK_TO_VISIBLE_COUNT_RATIO, DARK_TO_VISIBLE_GEOMETRIC_RATIO,
    };

    #[test]
    fn particle_branch_ratio_is_structural() {
        let rho_v = 2.4;
        let rho_d = dark_density_particle(rho_v);
        assert!((rho_d / rho_v - DARK_TO_VISIBLE_COUNT_RATIO).abs() < 1e-15);
    }

    #[test]
    fn total_density_is_nonnegative_for_nonnegative_inputs() {
        let rho_v = 1.2;
        let rho_t = total_density(DarkSectorBranch::Particle, rho_v, 1.0);
        assert!(rho_t >= 0.0);
    }

    #[test]
    fn structural_state_split_fraction_matches_lean_lane() {
        // Lean theorem: dark fraction = 5/16.
        assert!((DARK_FRACTION_TOTAL_STATE_SPLIT - 5.0 / 16.0).abs() < 1e-15);
    }

    #[test]
    fn geometric_branch_ratio_matches_structural_lane() {
        let rho_v = 2.4;
        let rho_d = dark_density_geometric(rho_v, 1.0);
        assert!((rho_d / rho_v - DARK_TO_VISIBLE_GEOMETRIC_RATIO).abs() < 1e-15);
        assert!((DARK_FRACTION_GEOMETRIC_STRUCTURAL - 60.0 / 71.0).abs() < 1e-15);
    }

    #[test]
    fn lensing_proxy_positive_for_positive_mass_and_impact() {
        let alpha = lensing_deflection(1.0e40, 1.0e20).expect("positive inputs");
        assert!(alpha > 0.0);
    }
}
