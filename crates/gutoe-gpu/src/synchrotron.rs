// GUTOE GRMHD Synchrotron scaffold (Lean parity layer)
// Copyright (C) 2026 Riff Labs, AGPL-3.0-or-later
//
// Mirrors lean/Gutoe/SynchrotronGRMHD.lean:
//   - magnetic sector multiplicity from Z3 triplet: 3
//   - lattice UV factor from lambda_qg = 1/12
//   - Kerr boost factor (1 + |omega|) using frame dragging

use crate::kerr::KerrMetric;
use crate::metric::LAMBDA_QG;
use std::f64::consts::PI;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum RenderSpectrum {
    Bolometric = 0,
    Radio = 1,
    Millimeter = 2,
    Infrared = 3,
    Optical = 4,
    Ultraviolet = 5,
    Xray = 6,
    Gamma = 7,
}

impl RenderSpectrum {
    pub fn parse(raw: &str) -> Option<Self> {
        match raw.to_ascii_lowercase().as_str() {
            "bolometric" | "bolo" | "all" => Some(Self::Bolometric),
            "radio" => Some(Self::Radio),
            "mm" | "millimeter" | "millimetre" | "230ghz" => Some(Self::Millimeter),
            "ir" | "infrared" => Some(Self::Infrared),
            "optical" | "visible" => Some(Self::Optical),
            "uv" | "ultraviolet" => Some(Self::Ultraviolet),
            "xray" | "x-ray" => Some(Self::Xray),
            "gamma" | "gamma-ray" | "gammaray" => Some(Self::Gamma),
            _ => None,
        }
    }

    pub fn from_env() -> Self {
        std::env::var("BH_SPECTRUM")
            .ok()
            .and_then(|s| Self::parse(&s))
            .unwrap_or(Self::Bolometric)
    }

    pub fn as_label(self) -> &'static str {
        match self {
            Self::Bolometric => "bolometric",
            Self::Radio => "radio",
            Self::Millimeter => "millimeter",
            Self::Infrared => "infrared",
            Self::Optical => "optical",
            Self::Ultraviolet => "ultraviolet",
            Self::Xray => "xray",
            Self::Gamma => "gamma",
        }
    }

    pub fn as_i32(self) -> i32 {
        self as i32
    }
}

#[inline]
pub fn magnetic_sector_weight() -> f64 {
    3.0
}

/// Speed of light (m/s).
pub const C_LIGHT: f64 = 299_792_458.0;
/// Planck constant (J*s).
pub const H_PLANCK: f64 = 6.626_070_15e-34;
/// Boltzmann constant (J/K).
pub const K_BOLTZMANN: f64 = 1.380_649e-23;
/// Electron rest mass (kg).
pub const M_ELECTRON: f64 = 9.109_383_701_5e-31;
/// Elementary charge (C).
pub const E_CHARGE: f64 = 1.602_176_634e-19;
/// Electron rest-energy temperature m_e c^2 / k_B in kelvin.
pub const ELECTRON_REST_TEMP_K: f64 = 5.930_867_40e9;

/// Representative observing frequency per rendered band.
#[inline]
pub fn band_frequency_hz(band: RenderSpectrum) -> f64 {
    match band {
        RenderSpectrum::Bolometric => 230.0e9,
        RenderSpectrum::Radio => 43.0e9,
        RenderSpectrum::Millimeter => 230.0e9,
        RenderSpectrum::Infrared => 3.0e14,
        RenderSpectrum::Optical => 5.0e14,
        RenderSpectrum::Ultraviolet => 1.2e15,
        RenderSpectrum::Xray => 5.0e17,
        RenderSpectrum::Gamma => 1.0e20,
    }
}

#[inline]
pub fn electron_theta_e(te_kelvin: f64) -> f64 {
    (te_kelvin / ELECTRON_REST_TEMP_K).max(1e-8)
}

#[inline]
pub fn cyclotron_frequency_hz(b_tesla: f64) -> f64 {
    (E_CHARGE * b_tesla.abs()) / (2.0 * PI * M_ELECTRON)
}

#[inline]
fn bessel_k2_inv_theta_approx(theta_e: f64) -> f64 {
    let theta = theta_e.max(1e-8);
    let x = 1.0 / theta;
    if x > 2.0 {
        // Large-x asymptotic K₂(x) ≈ sqrt(pi/(2x)) e^{-x} (1 + 15/(8x) + 105/(128x²)).
        let pref = (PI / (2.0 * x)).sqrt() * (-x).exp();
        pref * (1.0 + 15.0 / (8.0 * x) + 105.0 / (128.0 * x * x))
    } else {
        // Small-x relativistic limit K₂(x) ≈ 2/x².
        2.0 * theta * theta
    }
}

#[inline]
pub fn mahadevan_iprime(x_m: f64) -> f64 {
    let x = x_m.max(1e-12);
    let x13 = x.powf(1.0 / 3.0);
    let x14 = x.powf(1.0 / 4.0);
    let x16 = x.powf(1.0 / 6.0);
    let xsq = x.sqrt();
    let fit = 4.0505 / x16 * (1.0 + 0.40 / x14 + 0.5316 / xsq) * (-1.8899 * x13).exp();
    fit.max(0.0)
}

/// Planck source function B_nu(T) in SI units.
#[inline]
pub fn planck_b_nu(te_kelvin: f64, nu_hz: f64) -> f64 {
    let t = te_kelvin.max(1.0);
    let nu = nu_hz.max(0.0);
    if nu <= 0.0 {
        return 0.0;
    }
    let x = H_PLANCK * nu / (K_BOLTZMANN * t);
    if x < 1e-3 {
        // Rayleigh-Jeans limit to avoid catastrophic cancellation.
        2.0 * K_BOLTZMANN * t * nu * nu / (C_LIGHT * C_LIGHT)
    } else if x > 700.0 {
        0.0
    } else {
        2.0 * H_PLANCK * nu.powi(3) / (C_LIGHT * C_LIGHT * (x.exp() - 1.0))
    }
}

#[inline]
pub fn one_plus_lambda_qg_sq(l_planck: f64, r_s: f64) -> f64 {
    1.0 + LAMBDA_QG * (l_planck / r_s).powi(2)
}

/// Thermal synchrotron emissivity j_nu (Mahadevan-style fit, relativistic electrons).
///
/// Inputs:
/// - `n_e_m3`: electron number density [m^-3]
/// - `b_tesla`: magnetic field [T]
/// - `te_kelvin`: electron temperature [K]
/// - `nu_hz`: emitted-frame frequency [Hz]
/// - `sin_pitch`: angle-averaged pitch factor in [0,1]
#[inline]
pub fn thermal_synchrotron_emissivity(
    n_e_m3: f64,
    b_tesla: f64,
    te_kelvin: f64,
    nu_hz: f64,
    sin_pitch: f64,
) -> f64 {
    let n_e = n_e_m3.max(0.0);
    let nu = nu_hz.max(0.0);
    if n_e <= 0.0 || nu <= 0.0 {
        return 0.0;
    }

    let theta_e = electron_theta_e(te_kelvin);
    let nu_c = cyclotron_frequency_hz(b_tesla).max(1e-30);
    let sin_p = sin_pitch.clamp(1e-3, 1.0);
    let nu_s = ((2.0 / 9.0) * nu_c * theta_e * theta_e * sin_p).max(1e-30);
    let x_m = (nu / nu_s).max(1e-12);
    let k2 = bessel_k2_inv_theta_approx(theta_e).max(1e-30);
    let pref = (2.0_f64.sqrt() * PI * E_CHARGE * E_CHARGE * n_e * nu_s) / (3.0 * C_LIGHT * k2);
    (pref * mahadevan_iprime(x_m)).max(0.0)
}

/// Thermal synchrotron absorptivity alpha_nu from Kirchhoff's law.
#[inline]
pub fn thermal_synchrotron_absorption(
    n_e_m3: f64,
    b_tesla: f64,
    te_kelvin: f64,
    nu_hz: f64,
    sin_pitch: f64,
) -> f64 {
    let j_nu = thermal_synchrotron_emissivity(n_e_m3, b_tesla, te_kelvin, nu_hz, sin_pitch);
    let b_nu = planck_b_nu(te_kelvin, nu_hz);
    if b_nu > 1e-300 {
        (j_nu / b_nu).max(0.0)
    } else {
        0.0
    }
}

/// Backward-compatible wrapper used by older call-sites.
///
/// `nu` is interpreted as a normalized frequency relative to 230 GHz.
#[inline]
pub fn synchrotron_emissivity(b: f64, nu: f64, r_s: f64, l_planck: f64) -> f64 {
    let nu_hz = nu.max(0.0) * 230.0e9;
    let ne_fid = 1.0e12;
    let te_fid = 6.0e10;
    magnetic_sector_weight()
        * one_plus_lambda_qg_sq(l_planck, r_s.max(1e-12))
        * thermal_synchrotron_emissivity(ne_fid, b, te_fid, nu_hz, 0.785_398_163_39)
}

#[inline]
pub fn synchrotron_emissivity_kerr(
    b: f64,
    nu: f64,
    r_s: f64,
    l_planck: f64,
    kerr: &KerrMetric,
    r: f64,
    theta: f64,
) -> f64 {
    let base = synchrotron_emissivity(b, nu, r_s, l_planck);
    base * (1.0 + kerr.frame_dragging_omega(r, theta).abs())
}

/// Band-limited weight from a local blackbody patch.
///
/// `t_rel` is local disk temperature relative to ISCO.
#[inline]
pub fn band_weight_with_exposure(
    band: RenderSpectrum,
    t_rel: f64,
    fixed_exposure: Option<f64>,
) -> f64 {
    let t = t_rel.max(1e-6);
    let (x0, mut exposure) = match band {
        RenderSpectrum::Bolometric => return 1.0,
        RenderSpectrum::Radio => (0.02, 5.0),
        RenderSpectrum::Millimeter => (0.08, 3.5),
        RenderSpectrum::Infrared => (0.40, 2.0),
        RenderSpectrum::Optical => (1.00, 1.4),
        RenderSpectrum::Ultraviolet => (2.00, 1.1),
        RenderSpectrum::Xray => (6.00, 2.4),
        RenderSpectrum::Gamma => (20.0, 4.0),
    };
    if let Some(v) = fixed_exposure {
        exposure = v.max(0.0);
    }
    let x = x0 / t;
    let planck = if x > 80.0 {
        0.0
    } else {
        x.powi(3) / (x.exp() - 1.0)
    };
    let planck_ref = x0.powi(3) / (x0.exp() - 1.0);
    (exposure * planck / planck_ref.max(1e-12)).clamp(0.0, 64.0)
}

#[inline]
pub fn band_weight(band: RenderSpectrum, t_rel: f64) -> f64 {
    band_weight_with_exposure(band, t_rel, None)
}

#[inline]
pub fn band_tint(band: RenderSpectrum) -> [f64; 3] {
    match band {
        RenderSpectrum::Bolometric => [1.0, 1.0, 1.0],
        RenderSpectrum::Radio => [1.1, 0.45, 0.20],
        RenderSpectrum::Millimeter => [1.2, 0.70, 0.30],
        RenderSpectrum::Infrared => [1.2, 0.35, 0.25],
        RenderSpectrum::Optical => [1.0, 1.0, 1.0],
        RenderSpectrum::Ultraviolet => [0.55, 0.80, 1.20],
        RenderSpectrum::Xray => [0.35, 0.90, 1.25],
        RenderSpectrum::Gamma => [0.90, 0.55, 1.20],
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parser_aliases() {
        assert_eq!(
            RenderSpectrum::parse("mm"),
            Some(RenderSpectrum::Millimeter)
        );
        assert_eq!(RenderSpectrum::parse("x-ray"), Some(RenderSpectrum::Xray));
        assert_eq!(
            RenderSpectrum::parse("all"),
            Some(RenderSpectrum::Bolometric)
        );
        assert_eq!(RenderSpectrum::parse("???"), None);
    }

    #[test]
    fn emissivity_nonnegative_for_nonnegative_frequency() {
        let j = synchrotron_emissivity(2.0, 0.3, 1.0, 1.0);
        assert!(j >= 0.0);
    }

    #[test]
    fn thermal_emissivity_and_absorption_are_nonnegative() {
        let j = thermal_synchrotron_emissivity(1.0e11, 40.0, 8.0e10, 230.0e9, 0.7);
        let a = thermal_synchrotron_absorption(1.0e11, 40.0, 8.0e10, 230.0e9, 0.7);
        assert!(j >= 0.0, "j_nu should be nonnegative");
        assert!(a >= 0.0, "alpha_nu should be nonnegative");
    }

    #[test]
    fn emissivity_scales_roughly_with_density() {
        let j0 = thermal_synchrotron_emissivity(1.0e10, 30.0, 6.0e10, 230.0e9, 0.7);
        let j1 = thermal_synchrotron_emissivity(2.0e10, 30.0, 6.0e10, 230.0e9, 0.7);
        assert!(j1 > j0, "higher density should increase j_nu");
    }

    #[test]
    fn kerr_boost_not_smaller_than_base() {
        let k = KerrMetric::new(1.0, 0.8).expect("valid");
        let base = synchrotron_emissivity(1.2, 0.5, 1.0, 1.0);
        let boosted = synchrotron_emissivity_kerr(1.2, 0.5, 1.0, 1.0, &k, 2.0, 1.1);
        assert!(boosted >= base);
    }

    #[test]
    fn band_frequency_mm_matches_eht_anchor() {
        assert_eq!(band_frequency_hz(RenderSpectrum::Millimeter), 230.0e9);
    }

    #[test]
    fn fixed_exposure_override_changes_band_weight() {
        let t_rel = 1.0;
        let baseline = band_weight(RenderSpectrum::Optical, t_rel);
        let dimmer = band_weight_with_exposure(RenderSpectrum::Optical, t_rel, Some(0.7));
        let brighter = band_weight_with_exposure(RenderSpectrum::Optical, t_rel, Some(2.1));
        assert!(dimmer < baseline, "dimmer={dimmer} baseline={baseline}");
        assert!(
            brighter > baseline,
            "brighter={brighter} baseline={baseline}"
        );
    }
}
