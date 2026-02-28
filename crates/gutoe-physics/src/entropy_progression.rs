/*!
 * GUTOE Physics - Entropy Progression Lane
 * Copyright (C) 2026  Riff Labs
 *
 * Thermodynamic progression lane:
 * - Computes a universe-scale dissipative proxy over cosmic time.
 * - Activates abiotic/biotic/intelligent channels at derived epoch fractions.
 * - Gates strict stage progression and intelligence-era entropy jump.
 */

use crate::abiogenesis::{evaluate_abiogenesis_gate, AbiogenesisWindows};
use crate::constants::{
    ALPHA_LEADING_ORDER, C, DARK_TO_VISIBLE_COUNT_RATIO, DARK_TO_VISIBLE_GEOMETRIC_RATIO,
};
use crate::universe::{
    evaluate_universe_gate_with_depth, UniverseAssumptions, UniverseSimulationDepth,
    UniverseWindows,
};
use std::f64::consts::PI;

/// Effective stellar source temperature anchor (solar-like photosphere).
pub const STELLAR_SOURCE_TEMPERATURE_K: f64 = 5_772.0;
/// Planetary thermal sink floor used for habitable-surface dissipation.
pub const PLANETARY_SINK_FLOOR_K: f64 = 255.0;
/// Absorbed-flux anchor (Earth-like) used to scale per-area entropy production.
pub const ABSORBED_FLUX_ANCHOR_W_M2: f64 = 239.0;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum DissipativeStage {
    BareRock,
    PrebioticChemistry,
    AutocatalyticLife,
    PhotosyntheticBiosphere,
    MulticellularEcosystem,
    TechnologicalIntelligence,
}

impl DissipativeStage {
    pub const fn as_str(self) -> &'static str {
        match self {
            DissipativeStage::BareRock => "bare_rock",
            DissipativeStage::PrebioticChemistry => "prebiotic_chemistry",
            DissipativeStage::AutocatalyticLife => "autocatalytic_life",
            DissipativeStage::PhotosyntheticBiosphere => "photosynthetic_biosphere",
            DissipativeStage::MulticellularEcosystem => "multicellular_ecosystem",
            DissipativeStage::TechnologicalIntelligence => "technological_intelligence",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EntropyProgressionWindows {
    /// Minimum relative rise required between stage plateaus.
    pub stage_monotone_rel_min: f64,
    /// Intelligence-era step must exceed this factor over previous step.
    pub intelligence_step_factor_min: f64,
    /// Require at least this many local maxima/minima in total history.
    pub extrema_count_min: usize,
}

impl Default for EntropyProgressionWindows {
    fn default() -> Self {
        Self {
            stage_monotone_rel_min: 0.02,
            intelligence_step_factor_min: 2.0,
            extrema_count_min: 1,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StageActivation {
    pub stage: DissipativeStage,
    pub activation_age_gyr: f64,
    /// Incremental gain relative to the bare-rock dissipative baseline.
    pub incremental_gain: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EntropySample {
    pub age_gyr: f64,
    pub z: f64,
    pub star_formation_proxy: f64,
    pub baseline_per_area_w_m2_k: f64,
    pub prebiotic_per_area_w_m2_k: f64,
    pub autocatalytic_per_area_w_m2_k: f64,
    pub photosynthetic_per_area_w_m2_k: f64,
    pub multicellular_per_area_w_m2_k: f64,
    pub intelligence_per_area_w_m2_k: f64,
    pub total_per_area_w_m2_k: f64,
    pub total_universe_w_k: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct StagePlateauSummary {
    pub stage: DissipativeStage,
    pub age_start_gyr: f64,
    pub age_end_gyr: f64,
    pub mean_total_per_area_w_m2_k: f64,
    pub mean_effective_multiplier: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct EntropyProgressionScorecard {
    pub universe_age_gyr: f64,
    pub h0_km_s_mpc: f64,
    pub hubble_radius_m: f64,
    pub hubble_surface_area_m2: f64,
    pub stage_activations: Vec<StageActivation>,
    pub samples: Vec<EntropySample>,
    pub stage_plateaus: Vec<StagePlateauSummary>,
    pub local_maxima_count: usize,
    pub local_minima_count: usize,
    pub max_positive_step_age_gyr: f64,
    pub max_positive_step_w_m2_k: f64,
    pub monotone_stage_plateaus: bool,
    pub intelligence_step_dominant: bool,
    pub extrema_present: bool,
}

impl EntropyProgressionScorecard {
    pub const fn passes_all(&self) -> bool {
        self.monotone_stage_plateaus && self.intelligence_step_dominant && self.extrema_present
    }
}

fn km_s_mpc_to_s_inv(h0_km_s_mpc: f64) -> f64 {
    // 1 Mpc in meters.
    let meter_per_mpc = 3.085_677_581_491_367e22;
    (h0_km_s_mpc * 1_000.0) / meter_per_mpc
}

fn horizon_radius_from_h0_m(h0_km_s_mpc: f64) -> f64 {
    let h0_s_inv = km_s_mpc_to_s_inv(h0_km_s_mpc).max(1.0e-30);
    C / h0_s_inv
}

fn star_formation_proxy_from_z(z: f64) -> f64 {
    // Madau-Lilly style shape with structuralized exponents:
    // ψ(z) ∝ (1+z)^(11/4) / (1 + ((1+z)/3)^6)
    let x = (1.0 + z).max(1.0e-12);
    let numerator = x.powf(11.0 / 4.0);
    let denominator = 1.0 + (x / 3.0).powi(6);
    numerator / denominator
}

fn entropy_efficiency_per_joule(t_sink_k: f64) -> f64 {
    let t_sink = t_sink_k.max(1.0);
    (1.0 / t_sink - 1.0 / STELLAR_SOURCE_TEMPERATURE_K).max(0.0)
}

fn stage_epoch_fractions() -> [(DissipativeStage, f64); 5] {
    // Derived from finite Cl(1,3) channel fractions.
    [
        (DissipativeStage::PrebioticChemistry, 1.0 / 16.0),
        (DissipativeStage::AutocatalyticLife, 5.0 / 16.0),
        (DissipativeStage::PhotosyntheticBiosphere, 11.0 / 16.0),
        (DissipativeStage::MulticellularEcosystem, 13.0 / 16.0),
        (DissipativeStage::TechnologicalIntelligence, 15.0 / 16.0),
    ]
}

pub fn stage_incremental_gains_from_structure() -> [StageActivation; 5] {
    let abiogenesis = evaluate_abiogenesis_gate(AbiogenesisWindows::default(), 298.15);
    let closure_margin = abiogenesis.closure.closure_excess.max(0.0);
    [
        StageActivation {
            stage: DissipativeStage::PrebioticChemistry,
            activation_age_gyr: f64::NAN, // filled later
            incremental_gain: 11.0 * ALPHA_LEADING_ORDER,
        },
        StageActivation {
            stage: DissipativeStage::AutocatalyticLife,
            activation_age_gyr: f64::NAN, // filled later
            incremental_gain: closure_margin,
        },
        StageActivation {
            stage: DissipativeStage::PhotosyntheticBiosphere,
            activation_age_gyr: f64::NAN, // filled later
            incremental_gain: DARK_TO_VISIBLE_GEOMETRIC_RATIO
                / (1.0 + DARK_TO_VISIBLE_GEOMETRIC_RATIO),
        },
        StageActivation {
            stage: DissipativeStage::MulticellularEcosystem,
            activation_age_gyr: f64::NAN, // filled later
            incremental_gain: DARK_TO_VISIBLE_COUNT_RATIO,
        },
        StageActivation {
            stage: DissipativeStage::TechnologicalIntelligence,
            activation_age_gyr: f64::NAN, // filled later
            incremental_gain: DARK_TO_VISIBLE_GEOMETRIC_RATIO,
        },
    ]
}

fn build_stage_activations(universe_age_gyr: f64) -> Vec<StageActivation> {
    let mut gains = stage_incremental_gains_from_structure();
    for (i, (_, frac)) in stage_epoch_fractions().iter().enumerate() {
        gains[i].activation_age_gyr = universe_age_gyr * frac;
    }
    gains.to_vec()
}

fn stage_active(age_gyr: f64, activation_age_gyr: f64) -> bool {
    age_gyr >= activation_age_gyr
}

fn summarize_stage_plateaus(
    samples: &[EntropySample],
    activations: &[StageActivation],
) -> Vec<StagePlateauSummary> {
    let mut out = Vec::new();
    let mut starts = vec![(DissipativeStage::BareRock, 0.0)];
    for a in activations {
        starts.push((a.stage, a.activation_age_gyr));
    }
    starts.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    for i in 0..starts.len() {
        let (stage, age_start) = starts[i];
        let age_end = if i + 1 < starts.len() {
            starts[i + 1].1
        } else {
            f64::INFINITY
        };
        let mut sum = 0.0;
        let mut sum_mul = 0.0;
        let mut count = 0usize;
        for s in samples {
            if s.age_gyr >= age_start && s.age_gyr < age_end {
                sum += s.total_per_area_w_m2_k;
                let mul = if s.baseline_per_area_w_m2_k > 0.0 {
                    s.total_per_area_w_m2_k / s.baseline_per_area_w_m2_k
                } else {
                    1.0
                };
                sum_mul += mul;
                count += 1;
            }
        }
        let mean = if count > 0 {
            sum / count as f64
        } else {
            f64::NAN
        };
        let mean_mul = if count > 0 {
            sum_mul / count as f64
        } else {
            f64::NAN
        };
        out.push(StagePlateauSummary {
            stage,
            age_start_gyr: age_start,
            age_end_gyr: age_end,
            mean_total_per_area_w_m2_k: mean,
            mean_effective_multiplier: mean_mul,
        });
    }
    out
}

fn detect_local_extrema(samples: &[EntropySample]) -> (usize, usize) {
    if samples.len() < 3 {
        return (0, 0);
    }
    let mut nmax = 0usize;
    let mut nmin = 0usize;
    for i in 1..samples.len() - 1 {
        let a = samples[i - 1].total_per_area_w_m2_k;
        let b = samples[i].total_per_area_w_m2_k;
        let c = samples[i + 1].total_per_area_w_m2_k;
        if b > a && b > c {
            nmax += 1;
        }
        if b < a && b < c {
            nmin += 1;
        }
    }
    (nmax, nmin)
}

pub fn evaluate_entropy_progression_gate(
    assumptions: UniverseAssumptions,
    universe_windows: UniverseWindows,
    depth: UniverseSimulationDepth,
    windows: EntropyProgressionWindows,
) -> EntropyProgressionScorecard {
    let universe = evaluate_universe_gate_with_depth(assumptions, universe_windows, depth);
    let universe_age_gyr = universe.age_gyr.max(1.0e-12);

    let hubble_radius_m = horizon_radius_from_h0_m(universe.h0_km_s_mpc);
    let hubble_surface_area_m2 = 4.0 * PI * hubble_radius_m * hubble_radius_m;

    let activations = build_stage_activations(universe_age_gyr);

    let sf_max = universe
        .history
        .iter()
        .map(|h| star_formation_proxy_from_z(h.z))
        .fold(0.0_f64, f64::max)
        .max(1.0e-30);

    let mut samples = Vec::with_capacity(universe.history.len());
    for row in &universe.history {
        let age_gyr = row.age_seconds / (365.25 * 86_400.0 * 1.0e9);
        let sf_proxy = (star_formation_proxy_from_z(row.z) / sf_max).clamp(0.0, 1.0);

        let t_sink = row.temperature_k.max(PLANETARY_SINK_FLOOR_K);
        let entropy_eff = entropy_efficiency_per_joule(t_sink);
        let baseline = ABSORBED_FLUX_ANCHOR_W_M2 * sf_proxy * entropy_eff;

        let mut prebiotic = 0.0;
        let mut autocatalytic = 0.0;
        let mut photosynthetic = 0.0;
        let mut multicellular = 0.0;
        let mut intelligence = 0.0;

        for a in &activations {
            if !stage_active(age_gyr, a.activation_age_gyr) {
                continue;
            }
            let channel = baseline * a.incremental_gain;
            match a.stage {
                DissipativeStage::PrebioticChemistry => prebiotic = channel,
                DissipativeStage::AutocatalyticLife => autocatalytic = channel,
                DissipativeStage::PhotosyntheticBiosphere => photosynthetic = channel,
                DissipativeStage::MulticellularEcosystem => multicellular = channel,
                DissipativeStage::TechnologicalIntelligence => intelligence = channel,
                DissipativeStage::BareRock => {}
            }
        }

        let total_per_area =
            baseline + prebiotic + autocatalytic + photosynthetic + multicellular + intelligence;
        let total_universe = total_per_area * hubble_surface_area_m2;

        samples.push(EntropySample {
            age_gyr,
            z: row.z,
            star_formation_proxy: sf_proxy,
            baseline_per_area_w_m2_k: baseline,
            prebiotic_per_area_w_m2_k: prebiotic,
            autocatalytic_per_area_w_m2_k: autocatalytic,
            photosynthetic_per_area_w_m2_k: photosynthetic,
            multicellular_per_area_w_m2_k: multicellular,
            intelligence_per_area_w_m2_k: intelligence,
            total_per_area_w_m2_k: total_per_area,
            total_universe_w_k: total_universe,
        });
    }

    // Universe history is sampled in increasing z; progression gates require
    // chronological order from early epoch -> late epoch.
    samples.sort_by(|a, b| {
        a.age_gyr
            .partial_cmp(&b.age_gyr)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let stage_plateaus = summarize_stage_plateaus(&samples, &activations);

    let monotone_stage_plateaus = stage_plateaus.windows(2).all(|w2| {
        let a = w2[0].mean_effective_multiplier;
        let b = w2[1].mean_effective_multiplier;
        a.is_finite() && b.is_finite() && b > a * (1.0 + windows.stage_monotone_rel_min)
    });

    let mut max_step = f64::NEG_INFINITY;
    let mut pre_intelligence_max_step = 0.0_f64;
    let mut max_step_age = 0.0;
    let mut intelligence_step = f64::NAN;
    let intelligence_activation_age = activations
        .iter()
        .find(|a| a.stage == DissipativeStage::TechnologicalIntelligence)
        .map(|a| a.activation_age_gyr)
        .unwrap_or(f64::INFINITY);

    for i in 1..samples.len() {
        let dy = samples[i].total_per_area_w_m2_k - samples[i - 1].total_per_area_w_m2_k;
        if dy > max_step {
            max_step = dy;
            max_step_age = samples[i].age_gyr;
        }
        if samples[i].age_gyr < intelligence_activation_age {
            pre_intelligence_max_step = pre_intelligence_max_step.max(dy.max(0.0));
        }
        if samples[i - 1].age_gyr < intelligence_activation_age
            && samples[i].age_gyr >= intelligence_activation_age
        {
            intelligence_step = dy;
        }
    }
    if !intelligence_step.is_finite() {
        intelligence_step = 0.0;
    }

    let intelligence_step_dominant = intelligence_step
        >= windows.intelligence_step_factor_min * pre_intelligence_max_step.max(1.0e-30);

    let (local_maxima_count, local_minima_count) = detect_local_extrema(&samples);
    let extrema_present = local_maxima_count >= windows.extrema_count_min
        && local_minima_count >= windows.extrema_count_min;

    EntropyProgressionScorecard {
        universe_age_gyr,
        h0_km_s_mpc: universe.h0_km_s_mpc,
        hubble_radius_m,
        hubble_surface_area_m2,
        stage_activations: activations,
        samples,
        stage_plateaus,
        local_maxima_count,
        local_minima_count,
        max_positive_step_age_gyr: max_step_age,
        max_positive_step_w_m2_k: max_step.max(0.0),
        monotone_stage_plateaus,
        intelligence_step_dominant,
        extrema_present,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stage_epoch_fractions_are_strictly_increasing() {
        let fs = stage_epoch_fractions();
        for i in 1..fs.len() {
            assert!(fs[i].1 > fs[i - 1].1);
        }
    }

    #[test]
    fn progression_gate_passes_default() {
        let s = evaluate_entropy_progression_gate(
            UniverseAssumptions::default(),
            UniverseWindows::default(),
            UniverseSimulationDepth::default(),
            EntropyProgressionWindows::default(),
        );
        assert!(s.passes_all(), "entropy progression gate failed: {s:#?}");
    }

    #[test]
    fn intelligence_channel_is_nonzero_by_present_epoch() {
        let s = evaluate_entropy_progression_gate(
            UniverseAssumptions::default(),
            UniverseWindows::default(),
            UniverseSimulationDepth::default(),
            EntropyProgressionWindows::default(),
        );
        let last = s.samples.last().expect("sample");
        assert!(
            last.intelligence_per_area_w_m2_k > 0.0,
            "intelligence channel should activate by late epoch"
        );
    }
}
