/*!
 * GRAND-128: Singularity resolution lane (black-hole core + Big-Bang bounce).
 *
 * Black-hole side:
 *   r_eff = sqrt(r^2 + r_core^2),  r_core = sqrt(C_inf) * l_P
 *   This keeps curvature scalars finite at r = 0.
 *
 * Cosmology side:
 *   H^2 = (8πG/3) * ρ * (1 - ρ/ρ_crit)
 *   with lattice critical density ρ_crit = C_inf * ρ_P.
 *   This replaces the Big-Bang singularity with a finite-density bounce.
 */

use crate::constants::{C, G, PLANCK_LENGTH, PLANCK_MASS};

/// Lattice Richardson constant from the gravity metric lane.
pub const C_INF: f64 = 5466.0 / 10000.0;

/// Lattice core radius (m): r_core = sqrt(C_inf) * l_P.
pub fn lattice_core_radius_m(l_p: f64) -> f64 {
    C_INF.sqrt() * l_p.max(0.0)
}

/// Effective areal radius with UV floor.
pub fn effective_areal_radius_m(r_m: f64, l_p: f64) -> f64 {
    let r_core = lattice_core_radius_m(l_p);
    (r_m * r_m + r_core * r_core).sqrt()
}

/// Regularized Schwarzschild g_tt using r_eff floor.
pub fn regularized_g_tt(r_m: f64, r_s_m: f64, l_p: f64) -> f64 {
    -(1.0 - r_s_m / effective_areal_radius_m(r_m, l_p).max(1.0e-300))
}

/// Classical Schwarzschild radius from mass.
pub fn schwarzschild_radius_m(mass_kg: f64) -> f64 {
    if mass_kg <= 0.0 {
        return f64::NAN;
    }
    2.0 * G * mass_kg / (C * C)
}

/// Classical Kretschmann scalar for Schwarzschild: K = 12 r_s^2 / r^6.
pub fn kretschmann_classical_m4(r_m: f64, r_s_m: f64) -> Option<f64> {
    if r_m <= 0.0 || r_s_m <= 0.0 {
        return None;
    }
    Some(12.0 * r_s_m * r_s_m / r_m.powi(6))
}

/// Regularized Kretschmann scalar with r_eff floor.
pub fn kretschmann_regularized_m4(r_m: f64, r_s_m: f64, l_p: f64) -> f64 {
    let r_eff = effective_areal_radius_m(r_m, l_p).max(1.0e-300);
    12.0 * r_s_m * r_s_m / r_eff.powi(6)
}

/// Planck density (kg/m^3): ρ_P = m_P / l_P^3.
pub fn planck_density_kg_m3() -> f64 {
    PLANCK_MASS / PLANCK_LENGTH.powi(3)
}

/// Lattice critical density for bounce onset.
pub fn lattice_critical_density_kg_m3() -> f64 {
    C_INF * planck_density_kg_m3()
}

/// Bounce kernel κ(ρ) = ρ(1 - ρ/ρ_crit).
pub fn bounce_kernel(rho_kg_m3: f64, rho_crit_kg_m3: f64) -> f64 {
    if rho_crit_kg_m3 <= 0.0 {
        return 0.0;
    }
    rho_kg_m3 * (1.0 - rho_kg_m3 / rho_crit_kg_m3)
}

/// Lattice-regularized Friedmann RHS.
pub fn hubble_sq_bounce_si(rho_kg_m3: f64, rho_crit_kg_m3: f64) -> f64 {
    if rho_kg_m3 <= 0.0 || rho_crit_kg_m3 <= 0.0 {
        return 0.0;
    }
    ((8.0 * std::f64::consts::PI * G) / 3.0 * bounce_kernel(rho_kg_m3, rho_crit_kg_m3)).max(0.0)
}

/// For ρ(a) = ρ(a=1) * a^{-3(1+w)}, compute bounce scale factor where ρ = ρ_crit.
pub fn bounce_scale_factor(rho_a1_kg_m3: f64, w: f64, rho_crit_kg_m3: f64) -> Option<f64> {
    if rho_a1_kg_m3 <= 0.0 || rho_crit_kg_m3 <= 0.0 {
        return None;
    }
    let exponent = 3.0 * (1.0 + w);
    if exponent <= 0.0 {
        return None;
    }
    Some((rho_a1_kg_m3 / rho_crit_kg_m3).powf(1.0 / exponent))
}

/// Minimal comoving volume fraction at bounce (relative to a=1 reference).
pub fn bounce_volume_fraction(rho_a1_kg_m3: f64, w: f64, rho_crit_kg_m3: f64) -> Option<f64> {
    bounce_scale_factor(rho_a1_kg_m3, w, rho_crit_kg_m3).map(|a_b| a_b.powi(3))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn black_hole_origin_is_finite_under_regularization() {
        let l_p = PLANCK_LENGTH;
        let r_s = 10_000.0;
        let g0 = regularized_g_tt(0.0, r_s, l_p);
        let k0 = kretschmann_regularized_m4(0.0, r_s, l_p);
        assert!(g0.is_finite());
        assert!(k0.is_finite());
        assert!(k0 > 0.0);
    }

    #[test]
    fn classical_kretschmann_grows_faster_than_regularized_near_origin() {
        let r_s = 1000.0;
        let l_p = PLANCK_LENGTH;
        let r = 1.0e-20;
        let k_class = kretschmann_classical_m4(r, r_s).expect("classical K");
        let k_reg = kretschmann_regularized_m4(r, r_s, l_p);
        assert!(k_class >= k_reg);
    }

    #[test]
    fn bounce_kernel_zeroes_at_critical_density() {
        let rho_c = lattice_critical_density_kg_m3();
        let k0 = bounce_kernel(0.0, rho_c);
        let kc = bounce_kernel(rho_c, rho_c);
        let khalf = bounce_kernel(0.5 * rho_c, rho_c);
        assert_eq!(k0, 0.0);
        assert!(kc.abs() < 1.0e-12 * rho_c);
        assert!(khalf > 0.0);
    }

    #[test]
    fn hubble_square_vanishes_at_bounce() {
        let rho_c = lattice_critical_density_kg_m3();
        let h2 = hubble_sq_bounce_si(rho_c, rho_c);
        assert!(h2.abs() < 1.0e-18);
    }

    #[test]
    fn bounce_scale_factor_is_finite_and_positive() {
        let rho_c = lattice_critical_density_kg_m3();
        let rho_a1 = 1.0e-12 * rho_c;
        let a_b = bounce_scale_factor(rho_a1, 1.0 / 3.0, rho_c).expect("a_b");
        let v_b = bounce_volume_fraction(rho_a1, 1.0 / 3.0, rho_c).expect("v_b");
        assert!(a_b > 0.0 && a_b < 1.0);
        assert!(v_b > 0.0 && v_b < 1.0);
    }
}
