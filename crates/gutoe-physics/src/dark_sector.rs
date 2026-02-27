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
    C, DARK_TO_VISIBLE_COUNT_RATIO, DARK_TO_VISIBLE_GEOMETRIC_RATIO, G, LAMBDA_COSMOLOGICAL,
    LAMBDA_QG, PLANCK_LENGTH,
};

/// Unit conversion.
pub const METER_PER_KPC: f64 = 3.085_677_581_491_367e19;

/// Dark-sector modeling branch.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DarkSectorBranch {
    Particle,
    Geometric,
    /// Unified branch: local clustering from particle-like lane, total budget
    /// from geometric lane.
    Unified,
}

/// Structural particle-branch dark density: ρ_dark = (5/11) * ρ_visible.
pub fn dark_density_particle(visible_density: f64) -> f64 {
    DARK_TO_VISIBLE_COUNT_RATIO * visible_density
}

/// Curvature proxy from the lattice correction scale.
///
/// ρ_Λ = Λ c² / (8πG).
pub fn vacuum_energy_density_from_lambda(lambda: f64) -> Option<f64> {
    if lambda < 0.0 {
        return None;
    }
    Some(lambda * C * C / (8.0 * std::f64::consts::PI * G))
}

/// Structural vacuum-energy density from the derived cosmological term.
pub fn vacuum_energy_density_structural() -> f64 {
    vacuum_energy_density_from_lambda(LAMBDA_COSMOLOGICAL).unwrap_or(0.0)
}

/// Baryonic density estimate from circular velocity and radius:
/// M(r) = v²r/G and ρ = 3M/(4πr³) = 3v²/(4πGr²).
pub fn baryon_density_from_rotation(v_baryon_kms: f64, radius_kpc: f64) -> Option<f64> {
    if v_baryon_kms <= 0.0 || radius_kpc <= 0.0 {
        return None;
    }
    let v = v_baryon_kms * 1.0e3;
    let r = radius_kpc * METER_PER_KPC;
    Some(3.0 * v * v / (4.0 * std::f64::consts::PI * G * r * r))
}

/// Curvature amplification derived from modified Einstein + cosmology terms:
///
/// κ(r) = (1 + λ_QG (l_P/r)^2) * (1 + ρ_Λ / ρ_visible(r))
///
/// where ρ_Λ = Λ c² /(8πG). This replaces the report-level ad hoc proxy.
pub fn curvature_factor_from_einstein_cosmology(rho_visible: f64, radius_m: f64) -> f64 {
    if rho_visible <= 0.0 || radius_m <= 0.0 {
        return 1.0;
    }
    let rho_lambda = vacuum_energy_density_structural();
    let uv = 1.0 + LAMBDA_QG * (PLANCK_LENGTH / radius_m).powi(2);
    let source = 1.0 + rho_lambda / rho_visible;
    uv * source
}

/// Row-wise κ(r) from baryonic rotation decomposition.
pub fn curvature_factor_from_rotation(v_baryon_kms: f64, radius_kpc: f64) -> f64 {
    let rho_visible = match baryon_density_from_rotation(v_baryon_kms, radius_kpc) {
        Some(v) if v > 0.0 => v,
        _ => return 1.0,
    };
    let radius_m = radius_kpc * METER_PER_KPC;
    curvature_factor_from_einstein_cosmology(rho_visible, radius_m)
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
        DarkSectorBranch::Unified => {
            // Local clustering behavior follows the particle lane, modulated by
            // the derived κ(r) profile.
            let k = curvature_factor.max(0.0);
            DARK_TO_VISIBLE_COUNT_RATIO * k * visible_density
        }
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
