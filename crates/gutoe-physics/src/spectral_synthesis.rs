//! Composition-driven spectral synthesis (continuum + line proxies).

use std::collections::HashMap;

use crate::Species;

#[derive(Debug, Clone, PartialEq)]
pub struct SpectrumSample {
    pub wavelength_nm: f64,
    pub intensity: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SpectralLine {
    pub name: &'static str,
    pub wavelength_nm: f64,
    pub strength: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SynthesizedSpectrum {
    pub continuum: Vec<SpectrumSample>,
    pub lines: Vec<SpectralLine>,
}

pub fn synthesize_spectrum(
    abund: &HashMap<Species, f64>,
    temperature_k: f64,
    samples: usize,
) -> SynthesizedSpectrum {
    let t = temperature_k.max(1.0);
    let n = samples.max(16);
    let wl_min = 100.0;
    let wl_max = 2000.0;
    let mut continuum = Vec::with_capacity(n);
    for i in 0..n {
        let x = i as f64 / (n.saturating_sub(1).max(1) as f64);
        let wl = wl_min + x * (wl_max - wl_min);
        let bb = planck_proxy(wl, t);
        continuum.push(SpectrumSample {
            wavelength_nm: wl,
            intensity: bb,
        });
    }

    let h = abund.get(&Species::P1).copied().unwrap_or(0.0);
    let he = abund.get(&Species::He4).copied().unwrap_or(0.0);
    let c = abund.get(&Species::C12).copied().unwrap_or(0.0);
    let n14 = abund.get(&Species::N14).copied().unwrap_or(0.0);
    let o15 = abund.get(&Species::O15).copied().unwrap_or(0.0);

    let lines = vec![
        SpectralLine {
            name: "H-alpha",
            wavelength_nm: 656.28,
            strength: h * line_excitation(t, 10_000.0),
        },
        SpectralLine {
            name: "He I 587.6",
            wavelength_nm: 587.6,
            strength: he * line_excitation(t, 20_000.0),
        },
        SpectralLine {
            name: "[C II] 133.5",
            wavelength_nm: 133.5,
            strength: c * line_excitation(t, 8_000.0),
        },
        SpectralLine {
            name: "[N II] 658.4",
            wavelength_nm: 658.4,
            strength: n14 * line_excitation(t, 9_000.0),
        },
        SpectralLine {
            name: "[O III] 500.7",
            wavelength_nm: 500.7,
            strength: o15 * line_excitation(t, 12_000.0),
        },
    ];

    SynthesizedSpectrum { continuum, lines }
}

fn planck_proxy(wavelength_nm: f64, temperature_k: f64) -> f64 {
    let wl_m = wavelength_nm * 1e-9;
    let c2 = 1.4387769e-2;
    let x = (c2 / (wl_m * temperature_k)).clamp(1e-6, 200.0);
    let expx = x.exp();
    (1.0 / wl_m.powi(5)) / (expx - 1.0 + 1e-12)
}

fn line_excitation(t: f64, t0: f64) -> f64 {
    let x = (t / t0).max(1e-9);
    (x / (1.0 + x)).clamp(0.0, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    #[test]
    fn synthesis_outputs_continuum_and_lines() {
        let mut abund = HashMap::new();
        abund.insert(Species::P1, 0.7);
        abund.insert(Species::He4, 0.28);
        let s = synthesize_spectrum(&abund, 5800.0, 64);
        assert_eq!(s.continuum.len(), 64);
        assert!(!s.lines.is_empty());
        assert!(s.continuum.iter().all(|p| p.intensity.is_finite()));
    }

    #[test]
    fn hotter_star_boosts_blue_continuum_proxy() {
        let abund = HashMap::new();
        let cool = synthesize_spectrum(&abund, 4000.0, 128);
        let hot = synthesize_spectrum(&abund, 12000.0, 128);
        let blue_idx = 20usize;
        assert!(hot.continuum[blue_idx].intensity > cool.continuum[blue_idx].intensity);
    }
}
