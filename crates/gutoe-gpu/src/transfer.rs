// GUTOE covariant synchrotron transfer scaffold (Lean parity)
// Copyright (C) 2026 Riff Labs, AGPL-3.0-or-later

use crate::synchrotron::{thermal_synchrotron_absorption, thermal_synchrotron_emissivity};

/// GR invariant scaling: j_nu / nu^2 is invariant, so j_obs = g^2 * j_em.
#[inline]
pub fn covariant_emissivity(j_local: f64, g: f64) -> f64 {
    let g_safe = g.max(1e-12);
    j_local.max(0.0) * g_safe * g_safe
}

/// GR invariant scaling: alpha_nu * nu is invariant, so alpha_obs = alpha_em / g.
#[inline]
pub fn covariant_absorption(alpha_local: f64, g: f64) -> f64 {
    let g_safe = g.max(1e-12);
    alpha_local.max(0.0) / g_safe
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SynchrotronTransferCoefficients {
    pub j_em: f64,
    pub alpha_em: f64,
    pub source_em: f64,
    pub j_obs: f64,
    pub alpha_obs: f64,
    pub source_obs: f64,
}

/// Covariant thermal synchrotron coefficients from a literature-fit local model.
///
/// Local `j_nu` uses a Mahadevan-style relativistic fit and `alpha_nu` comes
/// from Kirchhoff (`alpha_nu = j_nu / B_nu`). Then redshift scaling is applied
/// with GR invariants to obtain observer-frame coefficients.
#[inline]
pub fn covariant_synchrotron_coefficients(
    n_e_m3: f64,
    b_tesla: f64,
    te_kelvin: f64,
    nu_obs_hz: f64,
    g: f64,
    sin_pitch: f64,
) -> SynchrotronTransferCoefficients {
    let g_safe = g.max(1e-12);
    let nu_em_hz = (nu_obs_hz / g_safe).max(0.0);
    let j_em = thermal_synchrotron_emissivity(n_e_m3, b_tesla, te_kelvin, nu_em_hz, sin_pitch);
    let alpha_em = thermal_synchrotron_absorption(n_e_m3, b_tesla, te_kelvin, nu_em_hz, sin_pitch);
    let source_em = if alpha_em > 1e-40 { j_em / alpha_em } else { j_em };
    let j_obs = covariant_emissivity(j_em, g_safe);
    let alpha_obs = covariant_absorption(alpha_em, g_safe);
    let source_obs = if alpha_obs > 1e-40 {
        j_obs / alpha_obs
    } else {
        j_obs
    };
    SynchrotronTransferCoefficients {
        j_em,
        alpha_em,
        source_em,
        j_obs,
        alpha_obs,
        source_obs,
    }
}

/// One-step transfer map:
/// I_out = I_in * exp(-tau) + S * (1 - exp(-tau)).
#[inline]
pub fn transfer_step(i_in: f64, source: f64, tau: f64) -> f64 {
    let e = (-tau).exp();
    i_in * e + source * (1.0 - e)
}

/// Stokes vector (I,Q,U,V) transport scaffold for polarized synchrotron.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Stokes {
    pub i: f64,
    pub q: f64,
    pub u: f64,
    pub v: f64,
}

/// Faraday rotation of linear polarization by angle `psi` (radians).
#[inline]
pub fn faraday_rotate(s: Stokes, psi: f64) -> Stokes {
    let c = (2.0 * psi).cos();
    let ss = (2.0 * psi).sin();
    Stokes {
        i: s.i,
        q: s.q * c - s.u * ss,
        u: s.q * ss + s.u * c,
        v: s.v,
    }
}

/// One-step polarized transfer with shared optical depth scalar.
#[inline]
pub fn transfer_step_stokes(i_in: Stokes, source: Stokes, tau: f64) -> Stokes {
    let e = (-tau).exp();
    Stokes {
        i: i_in.i * e + source.i * (1.0 - e),
        q: i_in.q * e + source.q * (1.0 - e),
        u: i_in.u * e + source.u * (1.0 - e),
        v: i_in.v * e + source.v * (1.0 - e),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn covariant_terms_nonnegative_for_nonnegative_inputs() {
        let j = covariant_emissivity(2.0, 0.7);
        let a = covariant_absorption(0.3, 0.7);
        assert!(j >= 0.0);
        assert!(a >= 0.0);
    }

    #[test]
    fn transfer_step_stays_between_input_and_source() {
        let i_in = 0.2;
        let source = 0.9;
        for tau in [0.0, 0.05, 0.5, 1.0, 3.0, 8.0] {
            let out = transfer_step(i_in, source, tau);
            assert!(out >= i_in - 1e-12, "out={out} i_in={i_in} tau={tau}");
            assert!(out <= source + 1e-12, "out={out} source={source} tau={tau}");
        }
    }

    #[test]
    fn transfer_limits_match_expected() {
        let i_in = 0.33;
        let source = 0.81;
        let near_zero = transfer_step(i_in, source, 1e-9);
        let large_tau = transfer_step(i_in, source, 40.0);
        assert!((near_zero - i_in).abs() < 1e-8);
        assert!((large_tau - source).abs() < 1e-8);
    }

    #[test]
    fn faraday_rotation_preserves_linear_polarization_norm() {
        let s = Stokes {
            i: 1.0,
            q: 0.3,
            u: -0.4,
            v: 0.02,
        };
        let p0 = (s.q * s.q + s.u * s.u).sqrt();
        let r = faraday_rotate(s, 0.37);
        let p1 = (r.q * r.q + r.u * r.u).sqrt();
        assert!((p0 - p1).abs() < 1e-10);
        assert_eq!(r.i, s.i);
        assert_eq!(r.v, s.v);
    }

    #[test]
    fn covariant_scaling_matches_invariant_rules() {
        let j = 2.5;
        let a = 0.4;
        let g = 0.5;
        let j_obs = covariant_emissivity(j, g);
        let a_obs = covariant_absorption(a, g);
        assert!((j_obs - j * g * g).abs() < 1e-12);
        assert!((a_obs - a / g).abs() < 1e-12);
    }

    #[test]
    fn thermal_covariant_coefficients_are_finite_and_nonnegative() {
        let c = covariant_synchrotron_coefficients(1.0e11, 30.0, 7.0e10, 230.0e9, 0.72, 0.7);
        assert!(c.j_em.is_finite() && c.j_em >= 0.0);
        assert!(c.alpha_em.is_finite() && c.alpha_em >= 0.0);
        assert!(c.j_obs.is_finite() && c.j_obs >= 0.0);
        assert!(c.alpha_obs.is_finite() && c.alpha_obs >= 0.0);
        assert!(c.source_obs.is_finite() && c.source_obs >= 0.0);
    }
}
