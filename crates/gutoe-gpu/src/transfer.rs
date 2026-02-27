// GUTOE covariant synchrotron transfer scaffold (Lean parity)
// Copyright (C) 2026 Riff Labs, AGPL-3.0-or-later

/// Covariant emissivity proxy: j_obs = j_local * g^3.
#[inline]
pub fn covariant_emissivity(j_local: f64, g: f64) -> f64 {
    j_local * g.powi(3)
}

/// Covariant absorption proxy: alpha_obs = alpha_local * g.
#[inline]
pub fn covariant_absorption(alpha_local: f64, g: f64) -> f64 {
    alpha_local * g
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
}
