/*!
 * GUTOE Physics - Galactic Life Map Lane
 * Copyright (C) 2026  Riff Labs
 *
 * Synthetic Milky Way life-map lane driven by entropy-stage progression.
 * Produces present-epoch stage tags plus forward Kardashev extrapolation.
 */

use crate::constants::DARK_TO_VISIBLE_GEOMETRIC_RATIO;
use crate::entropy_progression::{
    evaluate_entropy_progression_gate, stage_incremental_gains_from_structure, DissipativeStage,
    EntropyProgressionWindows,
};
use crate::universe::{UniverseAssumptions, UniverseSimulationDepth, UniverseWindows};
use crate::Z_INTEGRAL_MAX;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::f64::consts::PI;

pub const MILKY_WAY_STELLAR_COUNT_ESTIMATE: f64 = 2.0e11;
pub const DISK_RADIUS_LY: f64 = 45_000.0;
pub const DISK_SCALE_LENGTH_LY: f64 = 9_800.0;
pub const DISK_SCALE_HEIGHT_LY: f64 = 980.0;
pub const BULGE_SCALE_RADIUS_LY: f64 = 5_500.0;
pub const HALO_INNER_RADIUS_LY: f64 = 8_000.0;
pub const HALO_OUTER_RADIUS_LY: f64 = 80_000.0;
pub const SOLAR_GALACTIC_RADIUS_LY: f64 = 26_000.0;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum GalacticComponent {
    Disk,
    Bulge,
    Halo,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CivilizationStage {
    BareRock,
    PrebioticChemistry,
    AutocatalyticLife,
    PhotosyntheticBiosphere,
    MulticellularEcosystem,
    TechnologicalIntelligence,
    KardashevTypeI,
    KardashevTypeII,
    KardashevTypeIII,
}

impl CivilizationStage {
    pub const fn as_str(self) -> &'static str {
        match self {
            CivilizationStage::BareRock => "bare_rock",
            CivilizationStage::PrebioticChemistry => "prebiotic_chemistry",
            CivilizationStage::AutocatalyticLife => "autocatalytic_life",
            CivilizationStage::PhotosyntheticBiosphere => "photosynthetic_biosphere",
            CivilizationStage::MulticellularEcosystem => "multicellular_ecosystem",
            CivilizationStage::TechnologicalIntelligence => "technological_intelligence",
            CivilizationStage::KardashevTypeI => "kardashev_type_i",
            CivilizationStage::KardashevTypeII => "kardashev_type_ii",
            CivilizationStage::KardashevTypeIII => "kardashev_type_iii",
        }
    }

    pub const fn rank(self) -> u8 {
        match self {
            CivilizationStage::BareRock => 0,
            CivilizationStage::PrebioticChemistry => 1,
            CivilizationStage::AutocatalyticLife => 2,
            CivilizationStage::PhotosyntheticBiosphere => 3,
            CivilizationStage::MulticellularEcosystem => 4,
            CivilizationStage::TechnologicalIntelligence => 5,
            CivilizationStage::KardashevTypeI => 6,
            CivilizationStage::KardashevTypeII => 7,
            CivilizationStage::KardashevTypeIII => 8,
        }
    }

    pub const fn is_signal(self) -> bool {
        self.rank() >= CivilizationStage::TechnologicalIntelligence.rank()
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GalacticLifeMapConfig {
    pub sample_count: usize,
    pub rng_seed: u64,
    /// Forecast offsets in Gyr beyond the present age.
    pub forecast_offsets_gyr: [f64; 4],
}

impl Default for GalacticLifeMapConfig {
    fn default() -> Self {
        Self {
            sample_count: 120_000,
            rng_seed: 7_031_337,
            forecast_offsets_gyr: [0.0, 0.5, 1.0, 2.0],
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LifeThresholds {
    pub prebiotic_age_gyr: f64,
    pub autocatalytic_age_gyr: f64,
    pub photosynthetic_age_gyr: f64,
    pub multicellular_age_gyr: f64,
    pub intelligence_age_gyr: f64,
    pub kardashev_i_age_gyr: f64,
    pub kardashev_ii_age_gyr: f64,
    pub kardashev_iii_age_gyr: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StageEntropyMultipliers {
    pub bare_rock: f64,
    pub prebiotic: f64,
    pub autocatalytic: f64,
    pub photosynthetic: f64,
    pub multicellular: f64,
    pub intelligence: f64,
    pub kardashev_i: f64,
    pub kardashev_ii: f64,
    pub kardashev_iii: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GalacticLifeSeed {
    pub id: u64,
    pub component: GalacticComponent,
    pub x_ly: f64,
    pub y_ly: f64,
    pub z_ly: f64,
    pub galactic_radius_ly: f64,
    pub mass_solar: f64,
    pub age_gyr: f64,
    pub metallicity: f64,
    pub main_sequence_lifetime_gyr: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GalacticLifePoint {
    pub seed: GalacticLifeSeed,
    pub habitable: bool,
    pub habitability_score: f64,
    pub stage: CivilizationStage,
    pub entropy_multiplier: f64,
    pub signal: bool,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ForecastSnapshot {
    pub epoch_age_gyr: f64,
    pub delta_gyr: f64,
    pub habitable_count: usize,
    pub signal_count: usize,
    pub type_i_count: usize,
    pub type_ii_count: usize,
    pub type_iii_count: usize,
    pub mean_entropy_multiplier_habitable: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct GalacticLifeMapScorecard {
    pub config: GalacticLifeMapConfig,
    pub universe_age_gyr: f64,
    pub thresholds: LifeThresholds,
    pub multipliers: StageEntropyMultipliers,
    pub points: Vec<GalacticLifePoint>,
    pub stage_counts_present: Vec<(CivilizationStage, usize)>,
    pub habitable_count_present: usize,
    pub signal_count_present: usize,
    pub predicted_signal_count_milky_way_present: f64,
    pub forecasts: Vec<ForecastSnapshot>,
}

impl GalacticLifeMapScorecard {
    pub const fn present_signal_fraction(&self) -> f64 {
        if self.points.is_empty() {
            0.0
        } else {
            self.signal_count_present as f64 / self.points.len() as f64
        }
    }
}

fn exponential_truncated(rng: &mut StdRng, scale: f64, max: f64) -> f64 {
    let u = rng.gen::<f64>().clamp(1.0e-12, 1.0 - 1.0e-12);
    let norm = 1.0 - (-max / scale).exp();
    -scale * (1.0 - u * norm).ln()
}

fn laplace(rng: &mut StdRng, scale: f64) -> f64 {
    let u = rng.gen::<f64>().clamp(1.0e-12, 1.0 - 1.0e-12) - 0.5;
    let sign = if u >= 0.0 { 1.0 } else { -1.0 };
    -sign * scale * (1.0 - 2.0 * u.abs()).ln()
}

fn sample_component(rng: &mut StdRng) -> GalacticComponent {
    let u = rng.gen::<f64>();
    if u < 0.84 {
        GalacticComponent::Disk
    } else if u < 0.96 {
        GalacticComponent::Bulge
    } else {
        GalacticComponent::Halo
    }
}

fn sample_position(component: GalacticComponent, rng: &mut StdRng) -> (f64, f64, f64) {
    match component {
        GalacticComponent::Disk => {
            let r = exponential_truncated(rng, DISK_SCALE_LENGTH_LY, DISK_RADIUS_LY);
            let theta = 2.0 * PI * rng.gen::<f64>();
            let z = laplace(rng, DISK_SCALE_HEIGHT_LY).clamp(-2_000.0, 2_000.0);
            (r * theta.cos(), r * theta.sin(), z)
        }
        GalacticComponent::Bulge => {
            let rr = BULGE_SCALE_RADIUS_LY * (-rng.gen::<f64>().clamp(1.0e-12, 1.0).ln()).sqrt();
            let phi = 2.0 * PI * rng.gen::<f64>();
            let cos_t = 2.0 * rng.gen::<f64>() - 1.0;
            let sin_t = (1.0 - cos_t * cos_t).sqrt();
            let x = rr * sin_t * phi.cos();
            let y = rr * sin_t * phi.sin();
            let z = rr * cos_t * 0.55; // bulge flattening
            (x, y, z)
        }
        GalacticComponent::Halo => {
            // p(r) ~ 1/r^2 between inner and outer radii
            let u = rng.gen::<f64>().clamp(1.0e-12, 1.0 - 1.0e-12);
            let rinv = 1.0 / HALO_INNER_RADIUS_LY;
            let routv = 1.0 / HALO_OUTER_RADIUS_LY;
            let r = 1.0 / (rinv - u * (rinv - routv));
            let phi = 2.0 * PI * rng.gen::<f64>();
            let cos_t = 2.0 * rng.gen::<f64>() - 1.0;
            let sin_t = (1.0 - cos_t * cos_t).sqrt();
            let x = r * sin_t * phi.cos();
            let y = r * sin_t * phi.sin();
            let z = r * cos_t;
            (x, y, z)
        }
    }
}

fn sample_imf_powerlaw(rng: &mut StdRng, m_min: f64, m_max: f64, alpha: f64) -> f64 {
    let u = rng.gen::<f64>().clamp(1.0e-12, 1.0 - 1.0e-12);
    let p = 1.0 - alpha;
    let a = m_min.powf(p);
    let b = m_max.powf(p);
    (a + u * (b - a)).powf(1.0 / p)
}

pub fn main_sequence_lifetime_gyr(mass_solar: f64) -> f64 {
    let m = mass_solar.max(0.08);
    let t = if m < 0.43 {
        10.0 * m.powf(-2.0)
    } else {
        10.0 * m.powf(-2.5)
    };
    t.clamp(0.05, 180.0)
}

fn sample_age_gyr(
    component: GalacticComponent,
    radius_ly: f64,
    universe_age_gyr: f64,
    lifetime_gyr: f64,
    rng: &mut StdRng,
) -> f64 {
    let oldness_shape = match component {
        GalacticComponent::Disk => {
            let radial_oldness = (1.0 - radius_ly / DISK_RADIUS_LY).clamp(0.0, 1.0);
            1.2 + 1.6 * radial_oldness
        }
        GalacticComponent::Bulge => 3.2,
        GalacticComponent::Halo => 2.6,
    };
    let u = rng.gen::<f64>().clamp(1.0e-12, 1.0 - 1.0e-12);
    let age = universe_age_gyr * u.powf(1.0 / oldness_shape);
    age.min(0.98 * lifetime_gyr).max(0.01)
}

fn sample_metallicity(
    component: GalacticComponent,
    radius_ly: f64,
    age_gyr: f64,
    universe_age_gyr: f64,
    rng: &mut StdRng,
) -> f64 {
    let radial_gradient = (-(radius_ly - SOLAR_GALACTIC_RADIUS_LY) / 14_500.0).exp();
    let component_factor = match component {
        GalacticComponent::Disk => 1.0,
        GalacticComponent::Bulge => 1.25,
        GalacticComponent::Halo => 0.25,
    };
    let enrichment = (1.15 - 0.6 * (age_gyr / universe_age_gyr).clamp(0.0, 1.0)).clamp(0.2, 1.2);
    let jitter = 0.75 + 0.55 * rng.gen::<f64>();
    (0.0142 * radial_gradient * component_factor * enrichment * jitter).clamp(1.0e-4, 0.04)
}

pub fn derive_thresholds_and_multipliers() -> (f64, LifeThresholds, StageEntropyMultipliers) {
    let assumptions = UniverseAssumptions::default();
    let universe_windows = UniverseWindows::default();
    let depth = UniverseSimulationDepth {
        history_points: 512,
        history_z_max: 1.0e9,
        integral_z_max: Z_INTEGRAL_MAX,
    };
    let entropy = evaluate_entropy_progression_gate(
        assumptions,
        universe_windows,
        depth,
        EntropyProgressionWindows::default(),
    );
    let universe_age_gyr = entropy.universe_age_gyr;

    let mut prebiotic = f64::NAN;
    let mut auto = f64::NAN;
    let mut photo = f64::NAN;
    let mut multi = f64::NAN;
    let mut intel = f64::NAN;
    for a in entropy.stage_activations {
        match a.stage {
            DissipativeStage::PrebioticChemistry => prebiotic = a.activation_age_gyr,
            DissipativeStage::AutocatalyticLife => auto = a.activation_age_gyr,
            DissipativeStage::PhotosyntheticBiosphere => photo = a.activation_age_gyr,
            DissipativeStage::MulticellularEcosystem => multi = a.activation_age_gyr,
            DissipativeStage::TechnologicalIntelligence => intel = a.activation_age_gyr,
            DissipativeStage::BareRock => {}
        }
    }

    let dt = prebiotic.max(1.0e-9); // = universe_age/16 in current schedule
    let k1_age = intel + dt;
    let k2_age = k1_age + 0.5 * dt;
    let k3_age = k2_age + 0.25 * dt;

    let mut m = 1.0;
    let increments = stage_incremental_gains_from_structure();
    let prebiotic_mul = {
        m += increments[0].incremental_gain;
        m
    };
    let auto_mul = {
        m += increments[1].incremental_gain;
        m
    };
    let photo_mul = {
        m += increments[2].incremental_gain;
        m
    };
    let multi_mul = {
        m += increments[3].incremental_gain;
        m
    };
    let intel_mul = {
        m += increments[4].incremental_gain;
        m
    };

    let r = DARK_TO_VISIBLE_GEOMETRIC_RATIO;
    let k1_mul = intel_mul * (1.0 + r);
    let k2_mul = k1_mul * (1.0 + r * r);
    let k3_mul = k2_mul * (1.0 + r * r * r);

    (
        universe_age_gyr,
        LifeThresholds {
            prebiotic_age_gyr: prebiotic,
            autocatalytic_age_gyr: auto,
            photosynthetic_age_gyr: photo,
            multicellular_age_gyr: multi,
            intelligence_age_gyr: intel,
            kardashev_i_age_gyr: k1_age,
            kardashev_ii_age_gyr: k2_age,
            kardashev_iii_age_gyr: k3_age,
        },
        StageEntropyMultipliers {
            bare_rock: 1.0,
            prebiotic: prebiotic_mul,
            autocatalytic: auto_mul,
            photosynthetic: photo_mul,
            multicellular: multi_mul,
            intelligence: intel_mul,
            kardashev_i: k1_mul,
            kardashev_ii: k2_mul,
            kardashev_iii: k3_mul,
        },
    )
}

pub fn classify_stage(age_gyr: f64, th: LifeThresholds) -> CivilizationStage {
    if age_gyr < th.prebiotic_age_gyr {
        CivilizationStage::BareRock
    } else if age_gyr < th.autocatalytic_age_gyr {
        CivilizationStage::PrebioticChemistry
    } else if age_gyr < th.photosynthetic_age_gyr {
        CivilizationStage::AutocatalyticLife
    } else if age_gyr < th.multicellular_age_gyr {
        CivilizationStage::PhotosyntheticBiosphere
    } else if age_gyr < th.intelligence_age_gyr {
        CivilizationStage::MulticellularEcosystem
    } else if age_gyr < th.kardashev_i_age_gyr {
        CivilizationStage::TechnologicalIntelligence
    } else if age_gyr < th.kardashev_ii_age_gyr {
        CivilizationStage::KardashevTypeI
    } else if age_gyr < th.kardashev_iii_age_gyr {
        CivilizationStage::KardashevTypeII
    } else {
        CivilizationStage::KardashevTypeIII
    }
}

pub fn stage_entropy_multiplier(stage: CivilizationStage, m: StageEntropyMultipliers) -> f64 {
    match stage {
        CivilizationStage::BareRock => m.bare_rock,
        CivilizationStage::PrebioticChemistry => m.prebiotic,
        CivilizationStage::AutocatalyticLife => m.autocatalytic,
        CivilizationStage::PhotosyntheticBiosphere => m.photosynthetic,
        CivilizationStage::MulticellularEcosystem => m.multicellular,
        CivilizationStage::TechnologicalIntelligence => m.intelligence,
        CivilizationStage::KardashevTypeI => m.kardashev_i,
        CivilizationStage::KardashevTypeII => m.kardashev_ii,
        CivilizationStage::KardashevTypeIII => m.kardashev_iii,
    }
}

pub fn habitability_score(seed: GalacticLifeSeed) -> f64 {
    let mass_w = (-((seed.mass_solar - 0.95) / 0.45).powi(2)).exp();
    let metal_low = (seed.metallicity / 0.008).clamp(0.0, 1.2);
    let metal_high = (1.0 - ((seed.metallicity - 0.03).max(0.0) / 0.03)).clamp(0.2, 1.0);
    let metal_w = (metal_low * metal_high).clamp(0.0, 1.2);
    let radial_band = (-((seed.galactic_radius_ly - SOLAR_GALACTIC_RADIUS_LY) / 12_500.0).powi(2))
        .exp();
    let inner_hazard_suppression =
        1.0 / (1.0 + ((14_000.0 - seed.galactic_radius_ly) / 2_300.0).exp());
    let age_w = ((seed.age_gyr - 0.4) / 8.0).clamp(0.0, 1.0);
    let remaining = ((seed.main_sequence_lifetime_gyr - seed.age_gyr)
        / seed.main_sequence_lifetime_gyr)
        .clamp(0.0, 1.0);

    0.24 * mass_w + 0.22 * metal_w + 0.26 * radial_band * inner_hazard_suppression + 0.16 * age_w
        + 0.12 * remaining
}

pub fn is_habitable(seed: GalacticLifeSeed, score: f64) -> bool {
    score >= 0.50
        && seed.mass_solar >= 0.55
        && seed.mass_solar <= 1.35
        && seed.metallicity >= 0.0025
        && seed.age_gyr >= 0.4
        && seed.age_gyr <= 0.98 * seed.main_sequence_lifetime_gyr
}

fn stage_counts(points: &[GalacticLifePoint]) -> Vec<(CivilizationStage, usize)> {
    let mut out = Vec::new();
    for stage in [
        CivilizationStage::BareRock,
        CivilizationStage::PrebioticChemistry,
        CivilizationStage::AutocatalyticLife,
        CivilizationStage::PhotosyntheticBiosphere,
        CivilizationStage::MulticellularEcosystem,
        CivilizationStage::TechnologicalIntelligence,
        CivilizationStage::KardashevTypeI,
        CivilizationStage::KardashevTypeII,
        CivilizationStage::KardashevTypeIII,
    ] {
        let n = points.iter().filter(|p| p.stage == stage).count();
        out.push((stage, n));
    }
    out
}

pub fn infer_component_from_position(x_ly: f64, y_ly: f64, z_ly: f64) -> GalacticComponent {
    let r = (x_ly * x_ly + y_ly * y_ly).sqrt();
    if r <= 8_000.0 && z_ly.abs() <= 4_000.0 {
        GalacticComponent::Bulge
    } else if z_ly.abs() <= 3_000.0 && r <= 55_000.0 {
        GalacticComponent::Disk
    } else {
        GalacticComponent::Halo
    }
}

pub fn evaluate_galactic_life_map(config: GalacticLifeMapConfig) -> GalacticLifeMapScorecard {
    let (universe_age_gyr, thresholds, multipliers) = derive_thresholds_and_multipliers();
    let mut rng = StdRng::seed_from_u64(config.rng_seed);
    let mut points = Vec::with_capacity(config.sample_count);

    for id in 0..config.sample_count {
        let component = sample_component(&mut rng);
        let (x, y, z) = sample_position(component, &mut rng);
        let radius = (x * x + y * y).sqrt();
        let mass = sample_imf_powerlaw(&mut rng, 0.08, 60.0, 2.35);
        let t_ms = main_sequence_lifetime_gyr(mass);
        let age = sample_age_gyr(component, radius, universe_age_gyr, t_ms, &mut rng);
        let metallicity = sample_metallicity(component, radius, age, universe_age_gyr, &mut rng);
        let seed = GalacticLifeSeed {
            id: id as u64,
            component,
            x_ly: x,
            y_ly: y,
            z_ly: z,
            galactic_radius_ly: radius,
            mass_solar: mass,
            age_gyr: age,
            metallicity,
            main_sequence_lifetime_gyr: t_ms,
        };
        let h_score = habitability_score(seed);
        let habitable = is_habitable(seed, h_score);
        let stage = classify_stage(age, thresholds);
        let entropy = stage_entropy_multiplier(stage, multipliers);
        let signal = habitable && stage.is_signal();
        points.push(GalacticLifePoint {
            seed,
            habitable,
            habitability_score: h_score,
            stage,
            entropy_multiplier: entropy,
            signal,
        });
    }

    let stage_counts_present = stage_counts(&points);
    let habitable_count_present = points.iter().filter(|p| p.habitable).count();
    let signal_count_present = points.iter().filter(|p| p.signal).count();
    let signal_fraction = if config.sample_count > 0 {
        signal_count_present as f64 / config.sample_count as f64
    } else {
        0.0
    };
    let predicted_signal_count_milky_way_present = signal_fraction * MILKY_WAY_STELLAR_COUNT_ESTIMATE;

    let mut forecasts = Vec::with_capacity(config.forecast_offsets_gyr.len());
    for dt in config.forecast_offsets_gyr {
        let epoch_age = universe_age_gyr + dt.max(0.0);
        let mut habitable_count = 0usize;
        let mut signal_count = 0usize;
        let mut type_i_count = 0usize;
        let mut type_ii_count = 0usize;
        let mut type_iii_count = 0usize;
        let mut entropy_sum = 0.0;

        for p in &points {
            let aged = p.seed.age_gyr + dt;
            let alive = aged <= 0.98 * p.seed.main_sequence_lifetime_gyr;
            if !alive {
                continue;
            }
            let aged_seed = GalacticLifeSeed {
                age_gyr: aged,
                ..p.seed
            };
            let h_score = habitability_score(aged_seed);
            if !is_habitable(aged_seed, h_score) {
                continue;
            }
            habitable_count += 1;
            let stage = classify_stage(aged, thresholds);
            let ent = stage_entropy_multiplier(stage, multipliers);
            entropy_sum += ent;
            if stage.is_signal() {
                signal_count += 1;
            }
            if stage.rank() >= CivilizationStage::KardashevTypeI.rank() {
                type_i_count += 1;
            }
            if stage.rank() >= CivilizationStage::KardashevTypeII.rank() {
                type_ii_count += 1;
            }
            if stage.rank() >= CivilizationStage::KardashevTypeIII.rank() {
                type_iii_count += 1;
            }
        }

        forecasts.push(ForecastSnapshot {
            epoch_age_gyr: epoch_age,
            delta_gyr: dt,
            habitable_count,
            signal_count,
            type_i_count,
            type_ii_count,
            type_iii_count,
            mean_entropy_multiplier_habitable: if habitable_count > 0 {
                entropy_sum / habitable_count as f64
            } else {
                0.0
            },
        });
    }

    GalacticLifeMapScorecard {
        config,
        universe_age_gyr,
        thresholds,
        multipliers,
        points,
        stage_counts_present,
        habitable_count_present,
        signal_count_present,
        predicted_signal_count_milky_way_present,
        forecasts,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn life_map_is_seed_deterministic() {
        let cfg = GalacticLifeMapConfig {
            sample_count: 4096,
            rng_seed: 2026,
            forecast_offsets_gyr: [0.0, 1.0, 2.0, 3.0],
        };
        let a = evaluate_galactic_life_map(cfg);
        let b = evaluate_galactic_life_map(cfg);
        assert_eq!(a.habitable_count_present, b.habitable_count_present);
        assert_eq!(a.signal_count_present, b.signal_count_present);
        assert_eq!(a.stage_counts_present, b.stage_counts_present);
    }

    #[test]
    fn threshold_order_is_strict() {
        let (_, th, _) = derive_thresholds_and_multipliers();
        assert!(th.prebiotic_age_gyr < th.autocatalytic_age_gyr);
        assert!(th.autocatalytic_age_gyr < th.photosynthetic_age_gyr);
        assert!(th.photosynthetic_age_gyr < th.multicellular_age_gyr);
        assert!(th.multicellular_age_gyr < th.intelligence_age_gyr);
        assert!(th.intelligence_age_gyr < th.kardashev_i_age_gyr);
        assert!(th.kardashev_i_age_gyr < th.kardashev_ii_age_gyr);
        assert!(th.kardashev_ii_age_gyr < th.kardashev_iii_age_gyr);
    }

    #[test]
    fn forward_signal_count_is_non_decreasing_over_short_horizon() {
        let cfg = GalacticLifeMapConfig {
            sample_count: 8192,
            rng_seed: 4242,
            forecast_offsets_gyr: [0.0, 0.5, 1.0, 2.0],
        };
        let s = evaluate_galactic_life_map(cfg);
        for w in s.forecasts.windows(2) {
            assert!(w[1].signal_count >= w[0].signal_count.saturating_sub(16));
        }
    }
}
