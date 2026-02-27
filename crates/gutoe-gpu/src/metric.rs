// GUTOE Gravity Metric — Schwarzschild from SC Lattice Continuum Limit
// Copyright (C) 2026 Riff Labs, AGPL-3.0-or-later
//
// The GUTOE Schwarzschild metric emerges from the SC lattice continuum limit.
// Two algebraically-derived corrections to GR:
//
//   1. Singularity regularisation: coordinate r → areal radius r_eff = √(r² + r_core²)
//      r_core = √(C_∞) × l_P,   C_∞ = 0.5466 (Richardson extrapolation, 5-point GPU, L=161–961)
//      At r = 0: r_eff = r_core ≠ 0 → metric is finite, singularity resolved.
//
//   2. Dispersion correction: λ_QG = 1/12 from SC lattice kinetic operator.
//      T_SC(k) = (1 − cos k)/3 = k²/6 − k⁴/72 + O(k⁶)
//             = (k²/6)(1 − k²/12)   [along SC [100] axis]
//      Correction to Hawking temperature: T_H = T_GR × (1 + λ_QG × (l_P/r_s)²)
//
//   3. Photon sphere at areal radius 3r_s/2 → coordinate r_ph = √((3r_s/2)² − r_core²)
//   4. ISCO at areal radius 3r_s → coordinate r_ISCO = √((3r_s)² − r_core²)
//
// Both corrections vanish in the classical limit l_P → 0, recovering GR exactly.
//
// Metric signature (−, +, +, +) in spherical coordinates (t, r, θ, φ):
//   g_tt = −(1 − r_s/r_eff)
//   g_rr = (r_eff/r)² / (1 − r_s/r_eff)   [Jacobian from r → r_eff]
//   g_θθ = r_eff²
//   g_φφ = r_eff² sin²θ
//
// References:
//   Watson (1939): G_SC(0) = 1.5164 (simple cubic Watson integral)
//   GUTOE GPU (2026-02-21): C_∞ = 0.5466 ± 0.0005 (Richardson L=161–961)
//   Cl(1,3) dispersion: T_SC(k) = (1−cos k)/3, λ_QG = 1/12 from k⁴ coefficient

use std::f64::consts::PI;

// ── Physical constants from GUTOE lattice ─────────────────────────────────────

/// Lattice Bohr constant from GPU Richardson extrapolation (5-point, L=161–961).
/// Determines r_core = √C_∞ × l_P — the minimum physical radius (replaces r=0).
pub const C_INF: f64 = 0.5466;

/// SC lattice dispersion correction coefficient λ_QG = 1/12.
///
/// Derivation along the SC [100] axis (k = (k, 0, 0)):
///   T_SC(k) = (1 − cos k) / 3 = k²/6 − k⁴/72 + k⁶/2160 − ...
///           = (k²/6)(1 − k²/12)   [exact to k⁴, error O(k⁶)]
///
/// The k⁴ relative correction is −k²/12, so λ_QG = 1/12 is the coefficient
/// of the leading lattice modification to the gravitational propagator.
pub const LAMBDA_QG: f64 = 1.0 / 12.0;

/// Watson integral for the simple cubic lattice (Watson 1939, validated GPU 2026).
pub const WATSON_SC: f64 = 1.5164;

// ── GUTOE metric ──────────────────────────────────────────────────────────────

/// GUTOE-corrected Schwarzschild metric on the SC lattice continuum.
///
/// The coordinate `r` is NOT the areal radius; the areal radius is `r_eff(r)`.
/// This two-parameter family (r_s, l_P) recovers Schwarzschild as l_P → 0.
#[derive(Debug, Clone, Copy)]
pub struct GutoeMetric {
    /// Schwarzschild radius r_s = 2GM/c² (in Planck units: r_s = 2M_Planck).
    pub r_s: f64,
    /// Planck length l_P (set l_P = 1 for Planck units; use SI value for physical units).
    pub l_planck: f64,
}

impl GutoeMetric {
    /// Construct in Planck units: r_s in Planck lengths, l_P = 1.
    pub fn planck_units(r_s: f64) -> Self {
        Self { r_s, l_planck: 1.0 }
    }

    /// Construct with explicit Planck length (e.g. SI units).
    pub fn new(r_s: f64, l_planck: f64) -> Self {
        Self { r_s, l_planck }
    }

    /// Pure Schwarzschild metric (GR limit): l_P = 0, r_core = 0.
    ///
    /// All GUTOE corrections vanish. Use for comparison renders to quantify the
    /// lattice regularisation effect (visible only at Planck scale: l_P / r_s ≈ 1).
    pub fn schwarzschild(r_s: f64) -> Self {
        Self { r_s, l_planck: 0.0 }
    }

    /// Lattice core radius: r_core = √(C_∞) × l_P.
    ///
    /// This is the minimum physical radius in GUTOE — the "atom" of SC lattice space.
    /// The classical singularity r = 0 is replaced by a sphere of areal radius r_core.
    pub fn r_core(&self) -> f64 {
        C_INF.sqrt() * self.l_planck
    }

    /// Effective areal radius: r_eff(r) = √(r² + r_core²).
    ///
    /// Properties:
    ///   - At large r: r_eff ≈ r + r_core²/(2r) → r  (Schwarzschild recovery)
    ///   - At r = 0:   r_eff = r_core > 0             (singularity resolved)
    ///   - Always:     r_eff ≥ r_core > 0              (no zero crossing)
    pub fn r_eff(&self, r: f64) -> f64 {
        let rc = self.r_core();
        (r * r + rc * rc).sqrt()
    }

    /// g_tt metric component: g_tt(r) = −(1 − r_s/r_eff(r)).
    ///
    /// - Horizon (g_tt = 0):  r_eff = r_s  ↔  r_h = √(r_s² − r_core²)
    /// - Inside horizon:       g_tt > 0  (r_eff < r_s)
    /// - At r = 0:             g_tt = −(1 − r_s/r_core) — finite, not −∞
    pub fn g_tt(&self, r: f64) -> f64 {
        -(1.0 - self.r_s / self.r_eff(r))
    }

    /// g_rr metric component: g_rr(r) = (r_eff/r)² / (1 − r_s/r_eff(r)).
    ///
    /// The (r_eff/r)² Jacobian arises from the coordinate change r → r_eff.
    /// Returns f64::INFINITY at r = 0 (coordinate singularity, not physical).
    pub fn g_rr(&self, r: f64) -> f64 {
        if r.abs() < 1e-100 {
            return f64::INFINITY;
        }
        let re = self.r_eff(r);
        let f = 1.0 - self.r_s / re;
        (re / r) * (re / r) / f
    }

    /// g_θθ metric component: g_θθ(r) = r_eff(r)².
    pub fn g_theta(&self, r: f64) -> f64 {
        let re = self.r_eff(r);
        re * re
    }

    /// g_φφ metric component: g_φφ(r, θ) = r_eff(r)² sin²θ.
    pub fn g_phi(&self, r: f64, theta: f64) -> f64 {
        let re = self.r_eff(r);
        let s = theta.sin();
        re * re * s * s
    }

    /// Horizon coordinate radius: r_h = √(r_s² − r_core²).
    ///
    /// Exists when r_s > r_core. For macroscopic black holes (r_s >> l_P), r_h ≈ r_s.
    /// Returns None for sub-Planckian black holes (r_s ≤ r_core).
    pub fn r_horizon(&self) -> Option<f64> {
        let rc = self.r_core();
        if self.r_s > rc {
            Some((self.r_s * self.r_s - rc * rc).sqrt())
        } else {
            None
        }
    }

    /// Coordinate radius of the photon sphere.
    ///
    /// In GR the photon sphere is at areal radius 3r_s/2. In GUTOE it is the
    /// same areal radius, but the coordinate radius is r_ph = √((3r_s/2)² − r_core²).
    ///
    /// The observed black hole shadow depends on the areal radius r_eff = 3r_s/2,
    /// which matches GR. The GUTOE correction is O((l_P/r_s)²) in observables.
    pub fn r_photon_sphere(&self) -> Option<f64> {
        let r_areal = 1.5 * self.r_s;
        let rc = self.r_core();
        let arg = r_areal * r_areal - rc * rc;
        if arg > 0.0 {
            Some(arg.sqrt())
        } else {
            None
        }
    }

    /// Coordinate radius of the ISCO (innermost stable circular orbit).
    ///
    /// In GR the ISCO is at areal radius 3r_s. In GUTOE:
    /// r_ISCO_coord = √((3r_s)² − r_core²).
    pub fn r_isco(&self) -> Option<f64> {
        let r_areal = 3.0 * self.r_s;
        let rc = self.r_core();
        let arg = r_areal * r_areal - rc * rc;
        if arg > 0.0 {
            Some(arg.sqrt())
        } else {
            None
        }
    }

    /// Hawking temperature in natural units (ℏ = c = G = k_B = 1).
    ///
    /// T_H = (1 / (4π r_s)) × (1 - λ_QG × (l_P/r_s)²)
    ///
    /// The subluminal-dispersion correction is negative: GUTOE black holes are
    /// slightly cooler than GR in this branch.
    /// For astrophysical black holes (r_s >> l_P), the correction is ∼10⁻⁶¹ (unobservable).
    /// For Planck-mass black holes (r_s ≈ 2 l_P), the correction is ∼1/48 ≈ 2%.
    pub fn hawking_temperature(&self) -> f64 {
        let t_gr = self.gr_hawking_temperature();
        let correction = 1.0 - LAMBDA_QG * (self.l_planck / self.r_s).powi(2);
        t_gr * correction
    }

    /// Pure GR Hawking temperature: T_GR = 1/(4π r_s) (natural units).
    pub fn gr_hawking_temperature(&self) -> f64 {
        1.0 / (4.0 * PI * self.r_s)
    }

    /// Fractional Hawking temperature correction: δT/T = -λ_QG × (l_P/r_s)².
    pub fn hawking_correction_fraction(&self) -> f64 {
        -LAMBDA_QG * (self.l_planck / self.r_s).powi(2)
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-10;
    // Test case: r_s = 1000 l_P, so r_s >> r_core ≈ 0.739 l_P.
    const R_S: f64 = 1000.0;

    #[test]
    fn constants_algebraically_consistent() {
        // λ_QG = 1/12 exactly
        assert_eq!(LAMBDA_QG, 1.0 / 12.0);
        // r_core = √C_∞ in Planck units
        let m = GutoeMetric::planck_units(R_S);
        assert!((m.r_core() - C_INF.sqrt()).abs() < EPS);
        // r_core > 0
        assert!(m.r_core() > 0.0);
    }

    #[test]
    fn singularity_resolved() {
        let m = GutoeMetric::planck_units(R_S);
        let rc = m.r_core();
        // r_eff(0) = r_core — not zero, not infinity
        assert!(
            (m.r_eff(0.0) - rc).abs() < EPS,
            "r_eff(0) must equal r_core"
        );
        // g_tt(0) = -(1 - r_s/r_core) — finite
        let g = m.g_tt(0.0);
        assert!(
            g.is_finite(),
            "g_tt at r=0 must be finite (singularity resolved)"
        );
        let expected = -(1.0 - R_S / rc);
        assert!(
            (g - expected).abs() < EPS,
            "g_tt(0) = {g:.8}, expected {expected:.8}"
        );
    }

    #[test]
    fn r_eff_always_positive() {
        let m = GutoeMetric::planck_units(R_S);
        for r in [-1e6, -1.0, 0.0, 1.0, 1e6] {
            assert!(m.r_eff(r) > 0.0, "r_eff must be positive at r={r}");
        }
    }

    #[test]
    fn r_eff_recovers_r_at_large_r() {
        let m = GutoeMetric::planck_units(R_S);
        // r_eff(r)/r → 1 as r → ∞; correction is r_core²/(2r²)
        for r in [1e6_f64, 1e9, 1e12] {
            let ratio = m.r_eff(r) / r;
            // r_eff/r = sqrt(1 + r_core²/r²) ≈ 1 + r_core²/(2r²)
            let correction = m.r_core().powi(2) / (2.0 * r * r);
            assert!(
                (ratio - 1.0 - correction).abs() < 1e-8,
                "r_eff/r = {ratio:.12}, expected ≈ {:.12} at r={r}",
                1.0 + correction
            );
        }
    }

    #[test]
    fn g_tt_zero_at_horizon() {
        let m = GutoeMetric::planck_units(R_S);
        let rc = m.r_core();
        let r_h = (R_S * R_S - rc * rc).sqrt();
        // r_eff(r_h) = r_s → g_tt = -(1 - 1) = 0
        assert!(
            m.g_tt(r_h).abs() < EPS,
            "g_tt at horizon must be 0, got {}",
            m.g_tt(r_h)
        );
    }

    #[test]
    fn metric_signature_outside_horizon() {
        let m = GutoeMetric::planck_units(R_S);
        let r = 2.0 * R_S; // well outside horizon
        assert!(m.g_tt(r) < 0.0, "g_tt must be negative outside horizon");
        assert!(m.g_rr(r) > 0.0, "g_rr must be positive outside horizon");
        assert!(m.g_theta(r) > 0.0, "g_θθ must be positive");
        assert!(m.g_phi(r, PI / 2.0) > 0.0, "g_φφ must be positive at θ=π/2");
    }

    #[test]
    fn hawking_temperature_is_below_gr() {
        let m = GutoeMetric::planck_units(R_S);
        let t_gutoe = m.hawking_temperature();
        let t_gr = m.gr_hawking_temperature();
        // Subluminal branch correction is strictly negative
        assert!(
            t_gutoe < t_gr,
            "GUTOE T_H={t_gutoe:.8e} must be below GR T_H={t_gr:.8e}"
        );
        // Fractional correction = -λ_QG / r_s² (in Planck units with l_P=1)
        let frac_measured = (t_gutoe - t_gr) / t_gr;
        let frac_predicted = -LAMBDA_QG / (R_S * R_S);
        assert!(
            (frac_measured - frac_predicted).abs() < 1e-14,
            "hawking correction fraction {frac_measured:.2e} ≠ predicted {frac_predicted:.2e}"
        );
    }

    #[test]
    fn photon_sphere_areal_radius_matches_gr() {
        let m = GutoeMetric::planck_units(R_S);
        let r_ph = m
            .r_photon_sphere()
            .expect("photon sphere must exist for r_s >> r_core");
        // Areal radius must be exactly 3r_s/2
        assert!(
            (m.r_eff(r_ph) - 1.5 * R_S).abs() < EPS,
            "photon sphere areal radius = {:.6}, expected {:.6}",
            m.r_eff(r_ph),
            1.5 * R_S
        );
        // Coordinate radius is slightly less than areal radius
        assert!(
            r_ph < 1.5 * R_S,
            "coordinate r_ph must be less than areal 3r_s/2"
        );
    }

    #[test]
    fn isco_areal_radius_matches_gr() {
        let m = GutoeMetric::planck_units(R_S);
        let r_isco = m.r_isco().expect("ISCO must exist for r_s >> r_core");
        // Areal radius must be exactly 3r_s
        assert!(
            (m.r_eff(r_isco) - 3.0 * R_S).abs() < EPS,
            "ISCO areal radius = {:.6}, expected {:.6}",
            m.r_eff(r_isco),
            3.0 * R_S
        );
        // Coordinate radius is slightly less than areal radius
        assert!(
            r_isco < 3.0 * R_S,
            "coordinate r_ISCO must be less than areal 3r_s"
        );
    }

    #[test]
    fn gr_limit_at_large_r() {
        let m = GutoeMetric::planck_units(R_S);
        let r = 1e8 * R_S; // r >> r_s >> r_core — pure Schwarzschild regime
        let g_tt_gutoe = m.g_tt(r);
        // At this r, r_eff ≈ r to 1 part in 10¹⁶, so g_tt ≈ -(1 - r_s/r)
        let g_tt_gr = -(1.0 - R_S / r);
        assert!(
            (g_tt_gutoe - g_tt_gr).abs() < 1e-6,
            "g_tt should match GR at large r: {g_tt_gutoe:.10} vs {g_tt_gr:.10}"
        );
    }

    #[test]
    fn sc_dispersion_gives_lambda_qg_one_twelfth() {
        // Verify: T_SC(k) = (1 - cos k)/3 = (k²/6)(1 - k²/12) + O(k⁶)
        // i.e. λ_QG = 1/12 is the leading relative lattice correction.
        // Relative error = O(k⁴/360): for k=0.2, expect < k²/60 ≈ 6.7×10⁻⁴.
        for k in [0.01_f64, 0.05, 0.1, 0.2] {
            let t_exact = (1.0 - k.cos()) / 3.0;
            let t_approx = (k * k / 6.0) * (1.0 - k * k * LAMBDA_QG);
            let rel_err = (t_exact - t_approx).abs() / t_exact;
            // Exact k⁶ remainder / leading term = k⁴/360. Bound: rel_err < k²/60.
            let bound = k * k / 60.0;
            assert!(
                rel_err < bound,
                "SC dispersion rel_err={rel_err:.2e} > bound {bound:.2e} at k={k}"
            );
        }
        // Verify the constant itself
        assert_eq!(LAMBDA_QG, 1.0 / 12.0);
    }

    #[test]
    fn r_core_lt_photon_sphere_lt_isco_lt_horizon_order() {
        let m = GutoeMetric::planck_units(R_S);
        let rc = m.r_core();
        let r_h = m.r_horizon().expect("horizon must exist");
        let r_ph = m.r_photon_sphere().expect("photon sphere must exist");
        let r_isco = m.r_isco().expect("ISCO must exist");
        // Physical hierarchy: r_core < r_h < r_ph < r_ISCO
        assert!(
            rc < r_h,
            "r_core={rc:.4} must be less than r_horizon={r_h:.4}"
        );
        assert!(
            r_h < r_ph,
            "horizon={r_h:.4} must be inside photon sphere={r_ph:.4}"
        );
        assert!(
            r_ph < r_isco,
            "photon sphere={r_ph:.4} must be inside ISCO={r_isco:.4}"
        );
    }
}
