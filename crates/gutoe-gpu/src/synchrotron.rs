// GUTOE GRMHD Synchrotron scaffold (Lean parity layer)
// Copyright (C) 2026 Riff Labs, AGPL-3.0-or-later
//
// Mirrors lean/Gutoe/SynchrotronGRMHD.lean:
//   - magnetic sector multiplicity from Z3 triplet: 3
//   - lattice UV factor from lambda_qg = 1/12
//   - Kerr boost factor (1 + |omega|) using frame dragging

use crate::kerr::KerrMetric;
use crate::metric::LAMBDA_QG;

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

#[inline]
pub fn one_plus_lambda_qg_sq(l_planck: f64, r_s: f64) -> f64 {
    1.0 + LAMBDA_QG * (l_planck / r_s).powi(2)
}

#[inline]
pub fn synchrotron_emissivity(b: f64, nu: f64, r_s: f64, l_planck: f64) -> f64 {
    magnetic_sector_weight()
        * LAMBDA_QG
        * b.powi(2)
        * nu
        * (-(nu * r_s)).exp()
        * one_plus_lambda_qg_sq(l_planck, r_s)
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
pub fn band_weight_with_exposure(band: RenderSpectrum, t_rel: f64, fixed_exposure: Option<f64>) -> f64 {
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
        assert_eq!(RenderSpectrum::parse("mm"), Some(RenderSpectrum::Millimeter));
        assert_eq!(RenderSpectrum::parse("x-ray"), Some(RenderSpectrum::Xray));
        assert_eq!(RenderSpectrum::parse("all"), Some(RenderSpectrum::Bolometric));
        assert_eq!(RenderSpectrum::parse("???"), None);
    }

    #[test]
    fn emissivity_nonnegative_for_nonnegative_frequency() {
        let j = synchrotron_emissivity(2.0, 0.3, 1.0, 1.0);
        assert!(j >= 0.0);
    }

    #[test]
    fn kerr_boost_not_smaller_than_base() {
        let k = KerrMetric::new(1.0, 0.8).expect("valid");
        let base = synchrotron_emissivity(1.2, 0.5, 1.0, 1.0);
        let boosted = synchrotron_emissivity_kerr(1.2, 0.5, 1.0, 1.0, &k, 2.0, 1.1);
        assert!(boosted >= base);
    }

    #[test]
    fn fixed_exposure_override_changes_band_weight() {
        let t_rel = 1.0;
        let baseline = band_weight(RenderSpectrum::Optical, t_rel);
        let dimmer = band_weight_with_exposure(RenderSpectrum::Optical, t_rel, Some(0.7));
        let brighter = band_weight_with_exposure(RenderSpectrum::Optical, t_rel, Some(2.1));
        assert!(dimmer < baseline, "dimmer={dimmer} baseline={baseline}");
        assert!(brighter > baseline, "brighter={brighter} baseline={baseline}");
    }
}
