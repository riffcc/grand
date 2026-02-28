/*!
 * GUTOE Physics - Differential CMB Envelope Operator
 * Copyright (C) 2026  Riff Labs
 *
 * GRAND-344 / GRAND-355:
 *   Apply a differential (not absolute) high-l envelope correction built from
 *   microphysics-derived diffusion scale relative to CLASS baseline damping.
 */

use crate::cmb_class::ClassTtPoint;
use crate::cosmo_transfer::L_PEAK1_OBS;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DifferentialEnvelope {
    pub ell_diff_struct: f64,
    pub ell_diff_class: f64,
    pub ell_transition: f64,
}

/// Canonical projection factor from 3D diffusion geometry into multipole-space
/// effective damping in the current GRAND-344 lane.
pub fn structural_projection_factor_sqrt3_over_2() -> f64 {
    (3.0_f64).sqrt() / 2.0
}

/// Project raw structural damping scale into the effective multipole-space
/// scale used by the differential envelope operator.
pub fn projected_structural_ell_diff(raw_ell_diff: f64) -> f64 {
    raw_ell_diff * structural_projection_factor_sqrt3_over_2()
}

/// Estimate an effective CLASS damping scale from the high-l TT envelope.
///
/// Fit model in log-space: `ln D_l ~= a - l^2 / ell_diff^2` over a smoothed tail.
pub fn estimate_class_ell_diff(
    tt: &[ClassTtPoint],
    ell_min: u32,
    ell_max: u32,
) -> Result<f64, String> {
    if tt.len() < 50 {
        return Err("TT spectrum too short to estimate damping scale".to_string());
    }

    let mut tail: Vec<(f64, f64)> = tt
        .iter()
        .filter(|p| p.ell >= ell_min && p.ell <= ell_max && p.d_ell_tt_uk2 > 0.0)
        .map(|p| (p.ell as f64, p.d_ell_tt_uk2))
        .collect();
    if tail.len() < 100 {
        return Err(format!(
            "insufficient TT points in damping-fit window [{ell_min}, {ell_max}]"
        ));
    }
    tail.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    // Local moving-average smoothing to suppress acoustic oscillations.
    // Use a wide window so acoustic peak/trough oscillations do not bias the
    // damping-tail slope estimate.
    let w = 121usize;
    let mut smooth: Vec<(f64, f64)> = Vec::with_capacity(tail.len());
    for i in 0..tail.len() {
        let lo = i.saturating_sub(w);
        let hi = (i + w + 1).min(tail.len());
        let mut sum = 0.0;
        let mut n = 0usize;
        for (_, y) in tail.iter().take(hi).skip(lo) {
            sum += *y;
            n += 1;
        }
        let y = if n > 0 { sum / n as f64 } else { tail[i].1 };
        smooth.push((tail[i].0, y.max(1e-12)));
    }

    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut sxx = 0.0;
    let mut sxy = 0.0;
    let mut n = 0.0;
    for (ell, dl) in &smooth {
        let x = ell * ell;
        let y = dl.ln();
        sx += x;
        sy += y;
        sxx += x * x;
        sxy += x * y;
        n += 1.0;
    }

    let denom = n * sxx - sx * sx;
    if denom.abs() < 1e-12 {
        return Err("degenerate damping fit (denominator ~ 0)".to_string());
    }
    let slope = (n * sxy - sx * sy) / denom;
    if !(slope < 0.0) {
        return Err(format!(
            "non-negative damping slope from CLASS tail fit: {slope:.6e}"
        ));
    }

    let ell_diff = (-1.0 / slope).sqrt();
    if !ell_diff.is_finite() || ell_diff <= 0.0 {
        return Err(format!("invalid ell_diff from CLASS tail fit: {ell_diff}"));
    }
    Ok(ell_diff)
}

/// Apply a differential high-l envelope correction:
///
/// `F(l) = exp( -l^2 * (1/ell_struct^2 - 1/ell_class^2) * gate(l) )`
/// where `gate(l)` suppresses correction below the first acoustic peak.
pub fn apply_differential_envelope(
    spectrum: &[ClassTtPoint],
    env: DifferentialEnvelope,
) -> Vec<ClassTtPoint> {
    spectrum
        .iter()
        .map(|p| {
            let ell = p.ell as f64;
            let gate = 1.0 - (-(ell / env.ell_transition).powi(2)).exp();
            let inv_s2 = 1.0 / (env.ell_diff_struct * env.ell_diff_struct);
            let inv_c2 = 1.0 / (env.ell_diff_class * env.ell_diff_class);
            let expo = -(ell * ell) * (inv_s2 - inv_c2) * gate;
            let factor = expo.exp();
            let f = if factor.is_finite() { factor } else { 1.0 };
            ClassTtPoint {
                ell: p.ell,
                d_ell_tt_uk2: p.d_ell_tt_uk2 * f,
            }
        })
        .collect()
}

pub fn default_transition_scale() -> f64 {
    L_PEAK1_OBS
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn class_tail_fit_is_reasonable_on_synthetic_data() {
        let ell_diff_true = 1400.0;
        let mut tt = Vec::new();
        for ell in 2..=2500u32 {
            let l = ell as f64;
            let osc = 1.0 + 0.08 * (l / 45.0).sin();
            let dl = 4000.0 * (-(l * l) / (ell_diff_true * ell_diff_true)).exp() * osc.max(0.2);
            tt.push(ClassTtPoint {
                ell,
                d_ell_tt_uk2: dl,
            });
        }
        let est = estimate_class_ell_diff(&tt, 1200, 2200).expect("estimate");
        assert!((est - ell_diff_true).abs() / ell_diff_true < 0.12, "est={est}");
    }

    #[test]
    fn differential_operator_is_identity_when_scales_match() {
        let src = vec![
            ClassTtPoint {
                ell: 100,
                d_ell_tt_uk2: 10.0,
            },
            ClassTtPoint {
                ell: 1000,
                d_ell_tt_uk2: 2.0,
            },
        ];
        let env = DifferentialEnvelope {
            ell_diff_struct: 1500.0,
            ell_diff_class: 1500.0,
            ell_transition: default_transition_scale(),
        };
        let out = apply_differential_envelope(&src, env);
        assert!((out[0].d_ell_tt_uk2 - src[0].d_ell_tt_uk2).abs() < 1e-12);
        assert!((out[1].d_ell_tt_uk2 - src[1].d_ell_tt_uk2).abs() < 1e-12);
    }

    #[test]
    fn projection_factor_matches_sqrt3_over_2() {
        let f = structural_projection_factor_sqrt3_over_2();
        assert!((f - 0.866_025_403_784_438_6).abs() < 1e-12);
        let raw = 1_598.6;
        let proj = projected_structural_ell_diff(raw);
        assert!((proj - 1_384.479_219_291_401_8).abs() < 1e-9);
    }
}
