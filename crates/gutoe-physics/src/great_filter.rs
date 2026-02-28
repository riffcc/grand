/*!
 * GUTOE Physics - Great Filter Lane
 * Copyright (C) 2026  Riff Labs
 *
 * Monte Carlo survival model for technological civilizations using
 * derived stellar/thermodynamic inputs from existing lanes.
 */

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::f64::consts::PI;

pub const Z_SOLAR_METALLICITY: f64 = 0.0142;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GreatFilterCivilizationInput {
    /// Civilization stage rank from galactic lane (`TechnologicalIntelligence` >= 5).
    pub stage_rank: u8,
    /// Host-star mass in solar masses.
    pub mass_solar: f64,
    /// Host metallicity mass fraction Z.
    pub metallicity_z: f64,
    /// Stellar age in Gyr.
    pub age_gyr: f64,
    /// Habitability score in [0,1]-like range from the Gaia lane.
    pub habitability_score: f64,
    /// Entropy multiplier from progression lane.
    pub entropy_multiplier: f64,
    /// Local Kauffman closure control `N*p`.
    pub local_n_times_p: f64,
    /// Cylindrical galactic radius from center (ly).
    pub galactic_radius_ly: f64,
    /// Height from galactic plane (ly).
    pub galactic_z_ly: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GreatFilterWindows {
    pub trials_per_civilization: usize,
    pub transition_sigma: f64,
    pub atmosphere_sigma: f64,
    pub governance_sigma: f64,
    pub pressure_sigma: f64,
    pub threshold_conflict: f64,
    pub self_guard_multiplier: f64,
    pub environment_guard_multiplier: f64,
}

impl Default for GreatFilterWindows {
    fn default() -> Self {
        Self {
            trials_per_civilization: 512,
            transition_sigma: 0.35,
            atmosphere_sigma: 0.30,
            governance_sigma: 0.26,
            pressure_sigma: 0.22,
            threshold_conflict: 1.15,
            self_guard_multiplier: 1.00,
            environment_guard_multiplier: 1.00,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GreatFilterScorecard {
    pub trials: usize,
    pub mean_transition_years: f64,
    pub energy_pass_fraction: f64,
    pub conflict_pass_fraction: f64,
    pub self_destruction_pass_fraction: f64,
    pub environment_pass_fraction: f64,
    pub strict_pass_fraction: f64,
    pub stellar_stability_likelihood: f64,
    pub orbital_architecture_likelihood: f64,
    pub metallicity_band_likelihood: f64,
    pub galactic_environment_likelihood: f64,
    pub transition_likelihood: f64,
    pub self_destruction_likelihood: f64,
    pub survival_fraction: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct DerivedRiskBaselines {
    transition_years: f64,
    atmosphere_window_years: f64,
    governance_window_years: f64,
    resource_pressure: f64,
    weapons_pressure: f64,
    governance_capacity: f64,
    climate_forcing: f64,
    resilience_capacity: f64,
    stellar_stability: f64,
    orbital_architecture: f64,
    metallicity_band: f64,
    galactic_environment: f64,
}

fn standard_normal(rng: &mut StdRng) -> f64 {
    // Box-Muller (single sample).
    let u1 = rng.gen::<f64>().clamp(1.0e-12, 1.0 - 1.0e-12);
    let u2 = rng.gen::<f64>().clamp(1.0e-12, 1.0 - 1.0e-12);
    (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos()
}

fn sample_lognormal(base: f64, sigma: f64, rng: &mut StdRng) -> f64 {
    let b = base.max(1.0e-9);
    b * (sigma.max(0.0) * standard_normal(rng)).exp()
}

fn logistic(x: f64) -> f64 {
    if x >= 0.0 {
        1.0 / (1.0 + (-x).exp())
    } else {
        let e = x.exp();
        e / (1.0 + e)
    }
}

fn derive_risk_baselines(input: GreatFilterCivilizationInput) -> DerivedRiskBaselines {
    let metal = (input.metallicity_z / Z_SOLAR_METALLICITY).clamp(0.15, 4.0);
    let age_norm = (input.age_gyr / 13.8).clamp(0.05, 1.4);
    let hab = input.habitability_score.clamp(0.05, 1.25);
    let stage_adv = (input.stage_rank as f64 - 5.0).max(0.0);
    let entropy = input.entropy_multiplier.max(1.0).ln_1p();
    let closure = input.local_n_times_p.max(1.0e-6).ln_1p();
    let radius = input.galactic_radius_ly.max(10.0);
    let zabs = input.galactic_z_ly.abs();

    // Energy transition pressure: richer chemistry + mature entropy infrastructure
    // tends to compress pre-fusion bottleneck timescales.
    let transition_years =
        1_550.0 / (metal.powf(0.72) * (1.0 + 0.20 * entropy) * (1.0 + 0.08 * stage_adv));

    // Atmosphere/climate damage window before irreversible lock-in.
    let atmosphere_window_years = 1_350.0
        * hab.powf(0.60)
        * (1.0 + 0.23 * metal.sqrt())
        / input.mass_solar.max(0.2).powf(0.34);

    // Governance response window against capability growth.
    let governance_window_years = 420.0 * (1.0 + 0.56 * age_norm + 0.18 * hab)
        / (1.0 + 0.14 * stage_adv + 0.12 * entropy);

    // Conflict pressure from scarcity and extraction burden.
    let resource_pressure = 1.05 / metal.powf(0.40)
        + 0.23 * (input.mass_solar - 1.0).abs()
        + 0.43 * (1.0 - hab);

    // Capability pressure from tech stack and energy-density access.
    let weapons_pressure = 0.82 + 0.24 * stage_adv + 0.17 * closure + 0.11 * entropy;

    // Institutional capacity from maturity + stable environment.
    let governance_capacity = 1.00 + 0.59 * age_norm + 0.36 * hab + 0.23 * metal.sqrt();

    // Environmental forcing pressure from industrial entropy throughput.
    let climate_forcing = 0.88 + 0.26 * entropy + 0.14 * (input.mass_solar - 1.0).abs();

    // Biospheric/planetary resilience capacity.
    let resilience_capacity = 0.95 + 0.58 * hab + 0.21 * metal.sqrt();

    // Stellar stability proxy: M-dwarf flare suppression and activity aging.
    let mut flare_hazard = if input.mass_solar < 0.60 {
        1.45 * (0.60 / input.mass_solar.max(0.12)).powf(1.2)
    } else if input.mass_solar < 0.80 {
        0.85
    } else if input.mass_solar < 1.20 {
        0.42
    } else {
        0.50 + 0.30 * (input.mass_solar - 1.20)
    };
    flare_hazard *= (1.30 - 0.42 * age_norm).clamp(0.60, 1.40);
    let stellar_stability = (-0.62 * flare_hazard).exp().clamp(0.02, 1.0);

    // Metallicity sweet-spot: low Z lacks rocky chemistry; very high Z overproduces
    // giant-planet migration/dynamical disruption.
    let sweet_center = 1.10;
    let sweet_sigma = 0.55;
    let metal_dev = ((metal / sweet_center).ln() / sweet_sigma).powi(2);
    let metallicity_band = (-0.5 * metal_dev).exp().clamp(0.02, 1.0);

    // Orbital architecture proxy with gas-giant shielding vs migration penalty.
    let bodyguard_prob = 1.0 - (-0.95 * metal.powf(0.80)).exp();
    let migration_penalty = ((metal / 2.35).powi(2) / (1.0 + (metal / 2.35).powi(2))).clamp(0.0, 1.0);
    let orbital_architecture =
        ((0.36 + 0.64 * bodyguard_prob) * (1.0 - 0.60 * migration_penalty)).clamp(0.02, 1.0);

    // Galactic environment hazard: center + plane exposure.
    let center_hazard = (12_000.0 / (radius + 2_000.0)).powf(1.2).clamp(0.05, 3.0);
    let plane_hazard = (-(zabs / 850.0)).exp().clamp(0.0, 1.0);
    let galactic_environment = (-(0.38 * center_hazard + 0.24 * plane_hazard)).exp().clamp(0.02, 1.0);

    DerivedRiskBaselines {
        transition_years: transition_years.clamp(80.0, 12_000.0),
        atmosphere_window_years: atmosphere_window_years.clamp(120.0, 20_000.0),
        governance_window_years: governance_window_years.clamp(40.0, 8_000.0),
        resource_pressure: resource_pressure.clamp(0.10, 4.00),
        weapons_pressure: weapons_pressure.clamp(0.10, 4.00),
        governance_capacity: governance_capacity.clamp(0.10, 4.50),
        climate_forcing: climate_forcing.clamp(0.10, 4.50),
        resilience_capacity: resilience_capacity.clamp(0.10, 4.50),
        stellar_stability,
        orbital_architecture,
        metallicity_band,
        galactic_environment,
    }
}

pub fn evaluate_great_filter(
    input: GreatFilterCivilizationInput,
    windows: GreatFilterWindows,
    seed: u64,
) -> GreatFilterScorecard {
    let mut rng = StdRng::seed_from_u64(seed);
    let base = derive_risk_baselines(input);
    let trials = windows.trials_per_civilization.max(1);

    let mut mean_transition = 0.0;
    let mut energy_passes = 0usize;
    let mut conflict_passes = 0usize;
    let mut self_passes = 0usize;
    let mut env_passes = 0usize;
    let mut all_passes = 0usize;
    let mut stellar_like_sum = 0.0;
    let mut orbital_like_sum = 0.0;
    let mut metal_like_sum = 0.0;
    let mut galactic_like_sum = 0.0;
    let mut transition_like_sum = 0.0;
    let mut self_like_sum = 0.0;
    let mut survival_like_sum = 0.0;

    for _ in 0..trials {
        let dt_transition = sample_lognormal(base.transition_years, windows.transition_sigma, &mut rng);
        let dt_atmosphere =
            sample_lognormal(base.atmosphere_window_years, windows.atmosphere_sigma, &mut rng);
        let dt_governance =
            sample_lognormal(base.governance_window_years, windows.governance_sigma, &mut rng);

        let resource_pressure = (base.resource_pressure
            + windows.pressure_sigma * standard_normal(&mut rng))
        .clamp(0.10, 4.00);
        let weapons_pressure = (base.weapons_pressure
            + windows.pressure_sigma * standard_normal(&mut rng))
        .clamp(0.10, 4.00);
        let governance_capacity = (base.governance_capacity
            + windows.pressure_sigma * standard_normal(&mut rng))
        .clamp(0.10, 4.50);
        let climate_forcing = (base.climate_forcing
            + windows.pressure_sigma * standard_normal(&mut rng))
        .clamp(0.10, 4.50);
        let resilience_capacity = (base.resilience_capacity
            + windows.pressure_sigma * standard_normal(&mut rng))
        .clamp(0.10, 4.50);

        let energy_pass = dt_transition <= dt_atmosphere;
        let conflict_metric = resource_pressure * (dt_transition / dt_governance.max(1.0));
        let conflict_pass = conflict_metric <= windows.threshold_conflict;
        let self_pass = weapons_pressure <= governance_capacity * windows.self_guard_multiplier;
        let env_metric = climate_forcing * (dt_transition / dt_atmosphere.max(1.0)).clamp(0.6, 1.8);
        let env_pass = env_metric <= resilience_capacity * windows.environment_guard_multiplier;
        let all_pass = energy_pass && conflict_pass && self_pass && env_pass;

        // Continuous likelihood channels.
        let stellar_like =
            (base.stellar_stability + 0.10 * standard_normal(&mut rng)).clamp(0.0, 1.0);
        let orbital_like =
            (base.orbital_architecture + 0.12 * standard_normal(&mut rng)).clamp(0.0, 1.0);
        let metal_like = (base.metallicity_band + 0.10 * standard_normal(&mut rng)).clamp(0.0, 1.0);
        let galactic_like =
            (base.galactic_environment + 0.12 * standard_normal(&mut rng)).clamp(0.0, 1.0);
        let transition_like = logistic(
            2.4 * (dt_atmosphere / dt_transition.max(1.0)).ln()
                + 1.7 * (windows.threshold_conflict - conflict_metric)
                + 1.5 * (windows.environment_guard_multiplier - env_metric / resilience_capacity.max(1.0e-6)),
        )
        .clamp(0.0, 1.0);
        let self_like = logistic(
            2.8 * (governance_capacity * windows.self_guard_multiplier - weapons_pressure),
        )
        .clamp(0.0, 1.0);
        let survival_like =
            stellar_like * orbital_like * metal_like * galactic_like * transition_like * self_like;

        mean_transition += dt_transition;
        energy_passes += usize::from(energy_pass);
        conflict_passes += usize::from(conflict_pass);
        self_passes += usize::from(self_pass);
        env_passes += usize::from(env_pass);
        all_passes += usize::from(all_pass);
        stellar_like_sum += stellar_like;
        orbital_like_sum += orbital_like;
        metal_like_sum += metal_like;
        galactic_like_sum += galactic_like;
        transition_like_sum += transition_like;
        self_like_sum += self_like;
        survival_like_sum += survival_like;
    }

    let t = trials as f64;
    GreatFilterScorecard {
        trials,
        mean_transition_years: mean_transition / t,
        energy_pass_fraction: energy_passes as f64 / t,
        conflict_pass_fraction: conflict_passes as f64 / t,
        self_destruction_pass_fraction: self_passes as f64 / t,
        environment_pass_fraction: env_passes as f64 / t,
        strict_pass_fraction: all_passes as f64 / t,
        stellar_stability_likelihood: stellar_like_sum / t,
        orbital_architecture_likelihood: orbital_like_sum / t,
        metallicity_band_likelihood: metal_like_sum / t,
        galactic_environment_likelihood: galactic_like_sum / t,
        transition_likelihood: transition_like_sum / t,
        self_destruction_likelihood: self_like_sum / t,
        survival_fraction: survival_like_sum / t,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_input(metallicity_z: f64) -> GreatFilterCivilizationInput {
        GreatFilterCivilizationInput {
            stage_rank: 5,
            mass_solar: 0.95,
            metallicity_z,
            age_gyr: 8.0,
            habitability_score: 0.7,
            entropy_multiplier: 8.19,
            local_n_times_p: 1.35,
            galactic_radius_ly: 26_000.0,
            galactic_z_ly: 150.0,
        }
    }

    #[test]
    fn scorecard_fractions_are_bounded() {
        let s = evaluate_great_filter(sample_input(0.0142), GreatFilterWindows::default(), 42);
        assert!((0.0..=1.0).contains(&s.energy_pass_fraction));
        assert!((0.0..=1.0).contains(&s.conflict_pass_fraction));
        assert!((0.0..=1.0).contains(&s.self_destruction_pass_fraction));
        assert!((0.0..=1.0).contains(&s.environment_pass_fraction));
        assert!((0.0..=1.0).contains(&s.strict_pass_fraction));
        assert!((0.0..=1.0).contains(&s.stellar_stability_likelihood));
        assert!((0.0..=1.0).contains(&s.orbital_architecture_likelihood));
        assert!((0.0..=1.0).contains(&s.metallicity_band_likelihood));
        assert!((0.0..=1.0).contains(&s.galactic_environment_likelihood));
        assert!((0.0..=1.0).contains(&s.transition_likelihood));
        assert!((0.0..=1.0).contains(&s.self_destruction_likelihood));
        assert!((0.0..=1.0).contains(&s.survival_fraction));
    }

    #[test]
    fn higher_metallicity_improves_survival_in_default_lane() {
        let low = evaluate_great_filter(sample_input(0.004), GreatFilterWindows::default(), 7);
        let high = evaluate_great_filter(sample_input(0.020), GreatFilterWindows::default(), 7);
        assert!(
            high.survival_fraction >= low.survival_fraction,
            "expected higher metallicity to improve survival: low={low:?}, high={high:?}"
        );
    }
}
