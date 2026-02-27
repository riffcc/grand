// GUTOE Kerr Foundations — rotating black hole geometry helpers
// Copyright (C) 2026 Riff Labs, AGPL-3.0-or-later
//
// This module provides a physically-auditable Kerr baseline in Boyer-Lindquist
// variables as groundwork for GRAND-159. It intentionally does NOT yet include
// a lattice spin correction; it is the GR rotating reference layer we can test
// against before extending to a GUTOE-Kerr model.

/// Dimensionless spin `a_* = a / M` with physical Kerr bound |a_*| <= 1.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct KerrSpin {
    pub a_star: f64,
}

impl KerrSpin {
    pub fn new(a_star: f64) -> Option<Self> {
        if a_star.abs() <= 1.0 {
            Some(Self { a_star })
        } else {
            None
        }
    }
}

/// Kerr baseline parameterization using Schwarzschild radius `r_s = 2M`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct KerrMetric {
    pub r_s: f64,
    pub a_star: f64,
}

impl KerrMetric {
    pub fn new(r_s: f64, a_star: f64) -> Option<Self> {
        if r_s > 0.0 && a_star.abs() <= 1.0 {
            Some(Self { r_s, a_star })
        } else {
            None
        }
    }

    /// Mass parameter in geometric units (G=c=1): M = r_s / 2.
    pub fn mass(&self) -> f64 {
        0.5 * self.r_s
    }

    /// Kerr spin length parameter `a = a_* M`.
    pub fn a(&self) -> f64 {
        self.a_star * self.mass()
    }

    /// Returns the dimensionless spin `a/M` (same as `a_star`).
    pub fn a_over_m(&self) -> f64 {
        self.a_star
    }

    /// Σ = r² + a² cos²θ
    pub fn sigma(&self, r: f64, theta: f64) -> f64 {
        r * r + self.a() * self.a() * theta.cos().powi(2)
    }

    /// Δ = r² - 2Mr + a² = r² - r_s r + a²
    pub fn delta(&self, r: f64) -> f64 {
        r * r - self.r_s * r + self.a() * self.a()
    }

    /// Outer/inner Kerr horizons r_± = M ± sqrt(M² - a²).
    pub fn horizons(&self) -> (f64, f64) {
        let m = self.mass();
        let rad = (m * m - self.a() * self.a()).sqrt();
        (m + rad, m - rad)
    }

    /// Horizon angular velocity Ω_H = a / (r_+² + a²).
    pub fn horizon_angular_velocity(&self) -> f64 {
        let (r_plus, _) = self.horizons();
        self.a() / (r_plus * r_plus + self.a() * self.a())
    }

    /// Static limit (ergosurface) radius r_erg(θ) = M + sqrt(M² - a² cos²θ).
    pub fn ergosphere_radius(&self, theta: f64) -> f64 {
        let m = self.mass();
        let rad = (m * m - self.a() * self.a() * theta.cos().powi(2)).sqrt();
        m + rad
    }

    /// Equatorial prograde/retrograde circular photon orbit radii (GR Kerr):
    /// r_ph± = 2M * (1 + cos((2/3) arccos(∓a/M))).
    /// `prograde=true` gives the smaller radius branch for a>0.
    pub fn equatorial_photon_orbit_radius(&self, prograde: bool) -> f64 {
        let m = self.mass();
        let x = self.a() / m;
        let arg = if prograde { -x } else { x }.clamp(-1.0, 1.0);
        2.0 * m * (1.0 + ((2.0 / 3.0) * arg.acos()).cos())
    }

    /// Frame-dragging angular velocity for zero-angular-momentum observers:
    /// ω = -g_{tφ}/g_{φφ} = 2Mar / A
    /// with A = (r²+a²)² - a²Δ sin²θ.
    pub fn frame_dragging_omega(&self, r: f64, theta: f64) -> f64 {
        let a = self.a();
        let s2 = theta.sin().powi(2);
        let delta = self.delta(r);
        let a_cap = (r * r + a * a).powi(2) - a * a * delta * s2;
        if a_cap.abs() < 1e-15 {
            0.0
        } else {
            2.0 * self.mass() * a * r / a_cap
        }
    }

    /// Convert observer image-plane coordinates (alpha, beta) at inclination `theta_obs`
    /// to Kerr impact parameters (xi=Lz/E, eta=Q/E²) for null geodesics with E=1.
    ///
    /// Conventions follow the standard Bardeen/Carter setup:
    ///   xi  = -alpha * sin(theta_obs)
    ///   eta = beta² + (alpha² - a²) cos²(theta_obs)
    pub fn image_to_constants(&self, alpha: f64, beta: f64, theta_obs: f64) -> (f64, f64) {
        let s = theta_obs.sin();
        let c = theta_obs.cos();
        let a = self.a();
        let xi = -alpha * s;
        let eta = beta * beta + (alpha * alpha - a * a) * c * c;
        (xi, eta)
    }

    /// Kerr radial potential for null geodesics (Carter form, E=1):
    ///
    /// R(r) = [(r² + a²) - a*xi]² - Δ * [(xi - a)² + eta]
    ///
    /// Real radial motion requires R(r) >= 0.
    pub fn radial_potential(&self, r: f64, xi: f64, eta: f64) -> f64 {
        let a = self.a();
        let delta = self.delta(r);
        let t = (r * r + a * a) - a * xi;
        t * t - delta * ((xi - a) * (xi - a) + eta)
    }

    /// Polar potential for null geodesics (E=1):
    ///
    /// Θ(θ) = eta + a² cos²θ - xi² cot²θ
    ///
    /// Real polar motion requires Θ(θ) >= 0.
    pub fn polar_potential(&self, theta: f64, xi: f64, eta: f64) -> f64 {
        let c = theta.cos();
        let s = theta.sin();
        let a = self.a();
        let cot2 = if s.abs() < 1e-15 {
            f64::INFINITY
        } else {
            (c / s).powi(2)
        };
        eta + a * a * c * c - xi * xi * cot2
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, PI};

    #[test]
    fn schwarzschild_limit_from_zero_spin() {
        let k = KerrMetric::new(2.0, 0.0).expect("valid");
        let (r_plus, r_minus) = k.horizons();
        assert!((r_plus - 2.0).abs() < 1e-12);
        assert!(r_minus.abs() < 1e-12);
        assert!((k.ergosphere_radius(FRAC_PI_2) - 2.0).abs() < 1e-12);
        assert!((k.equatorial_photon_orbit_radius(true) - 3.0).abs() < 1e-12);
        assert!((k.equatorial_photon_orbit_radius(false) - 3.0).abs() < 1e-12);
    }

    #[test]
    fn extremal_kerr_horizons_coincide() {
        let k = KerrMetric::new(2.0, 1.0).expect("valid");
        let (r_plus, r_minus) = k.horizons();
        assert!((r_plus - 1.0).abs() < 1e-12);
        assert!((r_minus - 1.0).abs() < 1e-12);
    }

    #[test]
    fn ergosphere_touches_horizon_at_poles_and_bulges_at_equator() {
        let k = KerrMetric::new(2.0, 0.9).expect("valid");
        let (r_plus, _) = k.horizons();
        let r_pole = k.ergosphere_radius(0.0);
        let r_equator = k.ergosphere_radius(PI / 2.0);
        assert!((r_pole - r_plus).abs() < 1e-12);
        assert!(r_equator > r_plus);
    }

    #[test]
    fn prograde_photon_orbit_is_smaller_than_retrograde_for_positive_spin() {
        let k = KerrMetric::new(2.0, 0.9).expect("valid");
        let r_pro = k.equatorial_photon_orbit_radius(true);
        let r_ret = k.equatorial_photon_orbit_radius(false);
        assert!(r_pro < r_ret);
    }

    #[test]
    fn extremal_photon_orbits_match_known_values() {
        // a*=1, M=1 for r_s=2: prograde photon orbit at r=M, retrograde at r=4M.
        let k = KerrMetric::new(2.0, 1.0).expect("valid");
        let m = k.mass();
        let r_pro = k.equatorial_photon_orbit_radius(true);
        let r_ret = k.equatorial_photon_orbit_radius(false);
        assert!((r_pro - m).abs() < 1e-12, "r_pro={r_pro}, M={m}");
        assert!(
            (r_ret - 4.0 * m).abs() < 1e-12,
            "r_ret={r_ret}, 4M={}",
            4.0 * m
        );
    }

    #[test]
    fn horizon_angular_velocity_zero_for_non_spinning_case() {
        let k = KerrMetric::new(2.0, 0.0).expect("valid");
        assert!(k.horizon_angular_velocity().abs() < 1e-12);
        assert!(k.frame_dragging_omega(4.0, PI / 2.0).abs() < 1e-12);
    }

    #[test]
    fn image_plane_to_constants_schwarzschild_limit() {
        let k = KerrMetric::new(2.0, 0.0).expect("valid");
        let theta_obs = PI / 3.0;
        let (xi, eta) = k.image_to_constants(2.0, 1.5, theta_obs);
        assert!((xi + 2.0 * theta_obs.sin()).abs() < 1e-12);
        // a=0 => eta = beta^2 + alpha^2 cos^2(theta_obs)
        let expected = 1.5_f64.powi(2) + 2.0_f64.powi(2) * theta_obs.cos().powi(2);
        assert!((eta - expected).abs() < 1e-12);
    }

    #[test]
    fn image_plane_xi_is_beta_invariant() {
        // Parity with Lean: Gutoe.Geodesic3DProjection.kerrXi_beta_invariant
        let k = KerrMetric::new(2.0, 0.7).expect("valid");
        let alpha = 1.25;
        let theta_obs = PI / 5.0;
        let (xi1, _) = k.image_to_constants(alpha, -0.8, theta_obs);
        let (xi2, _) = k.image_to_constants(alpha, 2.2, theta_obs);
        assert!((xi1 - xi2).abs() < 1e-12, "xi1={xi1}, xi2={xi2}");
    }

    #[test]
    fn image_plane_eta_equatorial_is_beta_squared() {
        // Parity with Lean: Gutoe.Geodesic3DProjection.kerrEta_equatorial_from_ray
        let k = KerrMetric::new(2.0, 0.9).expect("valid");
        let alpha = -1.7;
        let beta = 0.6;
        let (_, eta) = k.image_to_constants(alpha, beta, FRAC_PI_2);
        assert!(
            (eta - beta * beta).abs() < 1e-12,
            "eta={eta}, beta^2={}",
            beta * beta
        );
    }

    #[test]
    fn radial_potential_matches_schwarzschild_form_when_a_zero() {
        let k = KerrMetric::new(2.0, 0.0).expect("valid");
        let r = 8.0;
        let xi = 1.25;
        let eta = 3.5;
        let rpot = k.radial_potential(r, xi, eta);
        // a=0 => R = r^4 - Δ*(xi^2 + eta), with Δ=r^2-r_s*r
        let delta = r * r - k.r_s * r;
        let expected = r.powi(4) - delta * (xi * xi + eta);
        assert!(
            (rpot - expected).abs() < 1e-10,
            "rpot={rpot}, expected={expected}"
        );
    }
}
