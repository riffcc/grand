/*!
 * Cyclosporine PK/PD bridge and safety proxy utilities for the MS lane.
 *
 * Scope:
 * - Bridge modeled effect-site concentration targets (nM) to a measurable
 *   whole-blood concentration proxy (ng/mL) through an uncertainty factor.
 * - Quantify exposure uncertainty and safety-window exceedance probabilities.
 *
 * This is a reduced-order translational scaffold, not clinical dosing guidance.
 */

use rand::{rngs::StdRng, Rng, SeedableRng};
use std::f64::consts::PI;

#[derive(Clone, Copy, Debug)]
pub struct CyclosporinePkBridgeInput {
    /// Modeled effect-site concentration target from mechanistic lane (nM).
    pub site_target_nanomolar: f64,
    /// Cyclosporine molecular weight (g/mol), default ~1202.61.
    pub molecular_weight_g_mol: f64,
    /// Median multiplier from modeled site concentration to measured blood concentration.
    pub blood_to_site_gain_median: f64,
    /// Geometric SD of the gain multiplier (lognormal uncertainty).
    pub blood_to_site_gain_gsd: f64,
    pub samples: usize,
    pub seed: u64,
}

#[derive(Clone, Debug)]
pub struct CyclosporinePkBridgeEnsemble {
    pub input: CyclosporinePkBridgeInput,
    pub blood_concentration_nanomolar: Vec<f64>,
    pub blood_concentration_ng_ml: Vec<f64>,
}

#[derive(Clone, Copy, Debug)]
pub struct CyclosporinePkBridgeSummary {
    pub p05_nanomolar: f64,
    pub p25_nanomolar: f64,
    pub p50_nanomolar: f64,
    pub p75_nanomolar: f64,
    pub p95_nanomolar: f64,
    pub p05_ng_ml: f64,
    pub p25_ng_ml: f64,
    pub p50_ng_ml: f64,
    pub p75_ng_ml: f64,
    pub p95_ng_ml: f64,
    pub mean_ng_ml: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct CyclosporineSafetyWindows {
    pub target_zone_low_ng_ml: f64,
    pub target_zone_high_ng_ml: f64,
    pub renal_caution_ng_ml: f64,
    pub renal_high_ng_ml: f64,
    pub neuro_caution_ng_ml: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct CyclosporineSafetyGateInput {
    pub windows: CyclosporineSafetyWindows,
    pub max_prob_above_renal_caution: f64,
    pub max_prob_above_renal_high: f64,
    pub max_prob_above_neuro_caution: f64,
    pub min_prob_in_target_zone: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct CyclosporineSafetyGateScore {
    pub prob_in_target_zone: f64,
    pub prob_above_target_zone: f64,
    pub prob_above_renal_caution: f64,
    pub prob_above_renal_high: f64,
    pub prob_above_neuro_caution: f64,
    pub target_zone_ok: bool,
    pub renal_caution_ok: bool,
    pub renal_high_ok: bool,
    pub neuro_caution_ok: bool,
    pub overall_pass: bool,
}

pub fn default_cyclosporine_pk_bridge_input() -> CyclosporinePkBridgeInput {
    CyclosporinePkBridgeInput {
        site_target_nanomolar: 20.0,
        molecular_weight_g_mol: 1202.61,
        blood_to_site_gain_median: 8.0,
        blood_to_site_gain_gsd: 1.6,
        samples: 50_000,
        seed: 1337,
    }
}

pub fn default_cyclosporine_safety_windows() -> CyclosporineSafetyWindows {
    CyclosporineSafetyWindows {
        target_zone_low_ng_ml: 80.0,
        target_zone_high_ng_ml: 300.0,
        renal_caution_ng_ml: 350.0,
        renal_high_ng_ml: 500.0,
        neuro_caution_ng_ml: 450.0,
    }
}

pub fn default_cyclosporine_safety_gate_input() -> CyclosporineSafetyGateInput {
    CyclosporineSafetyGateInput {
        windows: default_cyclosporine_safety_windows(),
        max_prob_above_renal_caution: 0.15,
        max_prob_above_renal_high: 0.05,
        max_prob_above_neuro_caution: 0.08,
        min_prob_in_target_zone: 0.50,
    }
}

pub fn nanomolar_to_ng_ml(conc_nanomolar: f64, molecular_weight_g_mol: f64) -> f64 {
    conc_nanomolar.max(0.0) * molecular_weight_g_mol.max(0.0) / 1000.0
}

fn standard_normal(rng: &mut StdRng) -> f64 {
    let u1 = rng.gen::<f64>().clamp(1.0e-12, 1.0 - 1.0e-12);
    let u2 = rng.gen::<f64>().clamp(1.0e-12, 1.0 - 1.0e-12);
    (-2.0_f64 * u1.ln()).sqrt() * (2.0_f64 * PI * u2).cos()
}

fn quantile(sorted: &[f64], q: f64) -> f64 {
    if sorted.is_empty() {
        return f64::NAN;
    }
    let qq = q.clamp(0.0, 1.0);
    let idx = ((sorted.len() - 1) as f64 * qq).round() as usize;
    sorted[idx.min(sorted.len() - 1)]
}

pub fn simulate_cyclosporine_pk_bridge(
    input: CyclosporinePkBridgeInput,
) -> CyclosporinePkBridgeEnsemble {
    let mut rng = StdRng::seed_from_u64(input.seed);
    let n = input.samples.max(64);
    let median = input.blood_to_site_gain_median.max(1.0e-6);
    let gsd = input.blood_to_site_gain_gsd.max(1.0 + 1.0e-6);
    let mu = median.ln();
    let sigma = gsd.ln();

    let mut blood_n_m = Vec::with_capacity(n);
    let mut blood_ng_ml = Vec::with_capacity(n);
    for _ in 0..n {
        let z = standard_normal(&mut rng);
        let gain = (mu + sigma * z).exp();
        let c_nm = input.site_target_nanomolar.max(0.0) * gain;
        let c_ng_ml = nanomolar_to_ng_ml(c_nm, input.molecular_weight_g_mol);
        blood_n_m.push(c_nm);
        blood_ng_ml.push(c_ng_ml);
    }

    CyclosporinePkBridgeEnsemble {
        input,
        blood_concentration_nanomolar: blood_n_m,
        blood_concentration_ng_ml: blood_ng_ml,
    }
}

pub fn summarize_cyclosporine_pk_bridge(
    ensemble: &CyclosporinePkBridgeEnsemble,
) -> CyclosporinePkBridgeSummary {
    let mut nms = ensemble.blood_concentration_nanomolar.clone();
    let mut ngs = ensemble.blood_concentration_ng_ml.clone();
    nms.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    ngs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let mean_ng_ml = if ngs.is_empty() {
        f64::NAN
    } else {
        ngs.iter().sum::<f64>() / ngs.len() as f64
    };

    CyclosporinePkBridgeSummary {
        p05_nanomolar: quantile(&nms, 0.05),
        p25_nanomolar: quantile(&nms, 0.25),
        p50_nanomolar: quantile(&nms, 0.50),
        p75_nanomolar: quantile(&nms, 0.75),
        p95_nanomolar: quantile(&nms, 0.95),
        p05_ng_ml: quantile(&ngs, 0.05),
        p25_ng_ml: quantile(&ngs, 0.25),
        p50_ng_ml: quantile(&ngs, 0.50),
        p75_ng_ml: quantile(&ngs, 0.75),
        p95_ng_ml: quantile(&ngs, 0.95),
        mean_ng_ml,
    }
}

pub fn probability_above_ng_ml(
    ensemble: &CyclosporinePkBridgeEnsemble,
    threshold_ng_ml: f64,
) -> f64 {
    if ensemble.blood_concentration_ng_ml.is_empty() {
        return f64::NAN;
    }
    let t = threshold_ng_ml;
    let hits = ensemble
        .blood_concentration_ng_ml
        .iter()
        .filter(|&&x| x >= t)
        .count();
    hits as f64 / ensemble.blood_concentration_ng_ml.len() as f64
}

pub fn probability_between_ng_ml(
    ensemble: &CyclosporinePkBridgeEnsemble,
    low_ng_ml: f64,
    high_ng_ml: f64,
) -> f64 {
    if ensemble.blood_concentration_ng_ml.is_empty() {
        return f64::NAN;
    }
    let lo = low_ng_ml.min(high_ng_ml);
    let hi = low_ng_ml.max(high_ng_ml);
    let hits = ensemble
        .blood_concentration_ng_ml
        .iter()
        .filter(|&&x| x >= lo && x <= hi)
        .count();
    hits as f64 / ensemble.blood_concentration_ng_ml.len() as f64
}

pub fn evaluate_cyclosporine_safety_gate(
    ensemble: &CyclosporinePkBridgeEnsemble,
    gate: CyclosporineSafetyGateInput,
) -> CyclosporineSafetyGateScore {
    let p_in_zone = probability_between_ng_ml(
        ensemble,
        gate.windows.target_zone_low_ng_ml,
        gate.windows.target_zone_high_ng_ml,
    );
    let p_above_zone = probability_above_ng_ml(ensemble, gate.windows.target_zone_high_ng_ml);
    let p_renal_caution = probability_above_ng_ml(ensemble, gate.windows.renal_caution_ng_ml);
    let p_renal_high = probability_above_ng_ml(ensemble, gate.windows.renal_high_ng_ml);
    let p_neuro_caution = probability_above_ng_ml(ensemble, gate.windows.neuro_caution_ng_ml);

    let target_ok = p_in_zone >= gate.min_prob_in_target_zone;
    let renal_caution_ok = p_renal_caution <= gate.max_prob_above_renal_caution;
    let renal_high_ok = p_renal_high <= gate.max_prob_above_renal_high;
    let neuro_caution_ok = p_neuro_caution <= gate.max_prob_above_neuro_caution;
    let overall = target_ok && renal_caution_ok && renal_high_ok && neuro_caution_ok;

    CyclosporineSafetyGateScore {
        prob_in_target_zone: p_in_zone,
        prob_above_target_zone: p_above_zone,
        prob_above_renal_caution: p_renal_caution,
        prob_above_renal_high: p_renal_high,
        prob_above_neuro_caution: p_neuro_caution,
        target_zone_ok: target_ok,
        renal_caution_ok,
        renal_high_ok,
        neuro_caution_ok,
        overall_pass: overall,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn conversion_is_consistent() {
        // 1000 nM = 1 uM => 1e-6 mol/L * MW g/mol => MW mg/L => MW ng/mL.
        let mw = 1202.61;
        let c = nanomolar_to_ng_ml(1000.0, mw);
        assert!((c - mw).abs() < 1.0e-9);
    }

    #[test]
    fn higher_gain_raises_median_exposure() {
        let mut low = default_cyclosporine_pk_bridge_input();
        low.blood_to_site_gain_median = 6.0;
        low.seed = 11;

        let mut high = low;
        high.blood_to_site_gain_median = 10.0;
        high.seed = 11;

        let s_low = summarize_cyclosporine_pk_bridge(&simulate_cyclosporine_pk_bridge(low));
        let s_high = summarize_cyclosporine_pk_bridge(&simulate_cyclosporine_pk_bridge(high));
        assert!(s_high.p50_ng_ml > s_low.p50_ng_ml);
    }

    #[test]
    fn safety_gate_flags_excessive_exposure() {
        let mut input = default_cyclosporine_pk_bridge_input();
        input.blood_to_site_gain_median = 16.0;
        input.seed = 7;
        let ens = simulate_cyclosporine_pk_bridge(input);
        let gate = default_cyclosporine_safety_gate_input();
        let score = evaluate_cyclosporine_safety_gate(&ens, gate);
        assert!(score.prob_above_renal_caution > 0.15 || score.prob_above_renal_high > 0.05);
    }
}
