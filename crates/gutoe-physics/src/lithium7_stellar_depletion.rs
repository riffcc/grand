//! Pop-II stellar lithium depletion transduction for Li-7 closure.

use crate::{
    eta10_from_baryogenesis, lithium7_residual_post_bbn_depletion_factor, RateEngine,
    BIVECTOR_TOTAL_COUNT, CLIFFORD_STATE_COUNT_STRUCTURAL, DARK_GEOMETRIC_AMPLIFICATION,
    DARK_STATE_COUNT_STRUCTURAL, HELIUM3_ETA_EXP,
};

pub const SOLAR_CORE_TEMPERATURE_K: f64 = 1.57e7;
pub const SOLAR_METALLICITY_Z: f64 = 0.0134;
pub const LITHIUM_BURN_ONSET_TEMPERATURE_K: f64 = 2.5e6;
pub const PRE_MAIN_SEQUENCE_YEARS: f64 = 1.0e7;
pub const INTEGRATION_STEP_YEARS: f64 = 1.0e5;
pub const ENVELOPE_HYDROGEN_FRACTION: f64 = 0.75;

/// Hard regression window for GRAND-349 closure in CI.
pub const LI7_STELLAR_CLOSURE_DELTA_ABS_MAX: f64 = 3.0e-4;

#[derive(Debug, Clone, Copy)]
pub struct PopIiCaseInput {
    pub label: &'static str,
    pub mass_solar: f64,
    pub metallicity_z: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct PopIiCaseResult {
    pub input: PopIiCaseInput,
    pub core_temperature_k: f64,
    pub convective_base_temperature_k: f64,
    pub t9: f64,
    pub li7_burn_rate_per_year: f64,
    pub survival_factor: f64,
    pub depletion_percent: f64,
    pub closure_delta: f64,
}

#[derive(Debug, Clone)]
pub struct Lithium7StellarDepletionReport {
    pub eta10: f64,
    pub required_survival_factor: f64,
    pub required_depletion_percent: f64,
    pub convective_exposure_factor: f64,
    pub cases: Vec<PopIiCaseResult>,
    pub best_match: PopIiCaseResult,
    pub agreement_with_required: f64,
}

pub fn popii_core_temperature_k(mass_solar: f64, metallicity_z: f64) -> f64 {
    // Structural mass scaling from shared 3/5 exponent.
    let mass_scale = mass_solar.max(0.1).powf(HELIUM3_ETA_EXP);
    // Lower metallicity gives a modest compression/temperature boost.
    let metallicity_boost = (SOLAR_METALLICITY_Z / metallicity_z.max(1.0e-6))
        .powf(1.0 / CLIFFORD_STATE_COUNT_STRUCTURAL);
    SOLAR_CORE_TEMPERATURE_K * mass_scale * metallicity_boost
}

pub fn popii_convective_base_temperature_k(mass_solar: f64, metallicity_z: f64) -> f64 {
    // Convective-envelope base proxy from finite dark-sector partitioning.
    popii_core_temperature_k(mass_solar, metallicity_z) / DARK_STATE_COUNT_STRUCTURAL
}

pub fn convective_exposure_factor() -> f64 {
    // 12/6 = 2 from structural amplification over bivector channels.
    DARK_GEOMETRIC_AMPLIFICATION / BIVECTOR_TOTAL_COUNT
}

pub fn simulate_lithium_survival(
    burn_rate_per_year: f64,
    years: f64,
    step_years: f64,
    hydrogen_fraction: f64,
    exposure_factor: f64,
) -> f64 {
    if burn_rate_per_year <= 0.0 || years <= 0.0 {
        return 1.0;
    }
    let n_steps = ((years / step_years.max(1.0)).ceil() as usize).max(1);
    let dt = years / (n_steps as f64);
    let effective_rate = burn_rate_per_year * hydrogen_fraction.max(0.0) * exposure_factor.max(0.0);
    let mut li = 1.0;
    for _ in 0..n_steps {
        li *= (-effective_rate * dt).exp();
    }
    li.clamp(0.0, 1.0)
}

pub fn default_popii_cases() -> [PopIiCaseInput; 4] {
    [
        PopIiCaseInput {
            label: "popii_low_z",
            mass_solar: 0.75,
            metallicity_z: 2.0e-4,
        },
        PopIiCaseInput {
            label: "spite_anchor",
            mass_solar: 0.80,
            metallicity_z: 5.0e-4,
        },
        PopIiCaseInput {
            label: "popii_mid_z",
            mass_solar: 0.85,
            metallicity_z: 1.0e-3,
        },
        PopIiCaseInput {
            label: "metal_rich_control",
            mass_solar: 0.80,
            metallicity_z: 1.5e-3,
        },
    ]
}

pub fn evaluate_popii_case(
    input: PopIiCaseInput,
    required_survival: f64,
    rates: &RateEngine,
) -> PopIiCaseResult {
    let t_base = popii_convective_base_temperature_k(input.mass_solar, input.metallicity_z);
    let t9 = t_base * 1.0e-9;
    let burn_rate = rates.rate_for("li7_burn", t9).unwrap_or(0.0);
    let survival = simulate_lithium_survival(
        burn_rate,
        PRE_MAIN_SEQUENCE_YEARS,
        INTEGRATION_STEP_YEARS,
        ENVELOPE_HYDROGEN_FRACTION,
        convective_exposure_factor(),
    );
    let depletion_percent = 100.0 * (1.0 - survival);

    PopIiCaseResult {
        input,
        core_temperature_k: popii_core_temperature_k(input.mass_solar, input.metallicity_z),
        convective_base_temperature_k: t_base,
        t9,
        li7_burn_rate_per_year: burn_rate,
        survival_factor: survival,
        depletion_percent,
        closure_delta: survival - required_survival,
    }
}

pub fn evaluate_lithium7_stellar_depletion(
    cases: &[PopIiCaseInput],
) -> Lithium7StellarDepletionReport {
    let eta10 = eta10_from_baryogenesis();
    let required_survival = lithium7_residual_post_bbn_depletion_factor(eta10);
    let required_depletion_percent = 100.0 * (1.0 - required_survival);
    let rates = RateEngine::baseline_p1();

    let results: Vec<PopIiCaseResult> = cases
        .iter()
        .copied()
        .map(|c| evaluate_popii_case(c, required_survival, &rates))
        .collect();

    let best_match = results
        .iter()
        .min_by(|a, b| {
            a.closure_delta
                .abs()
                .partial_cmp(&b.closure_delta.abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .copied()
        .expect("at least one case");

    let agreement_with_required =
        1.0 - (best_match.closure_delta.abs() / required_survival.max(1.0e-12));

    Lithium7StellarDepletionReport {
        eta10,
        required_survival_factor: required_survival,
        required_depletion_percent,
        convective_exposure_factor: convective_exposure_factor(),
        cases: results,
        best_match,
        agreement_with_required,
    }
}

pub fn evaluate_lithium7_stellar_depletion_default() -> Lithium7StellarDepletionReport {
    evaluate_lithium7_stellar_depletion(&default_popii_cases())
}

pub fn lithium7_stellar_closure_pass(report: &Lithium7StellarDepletionReport) -> bool {
    report.best_match.closure_delta.abs() <= LI7_STELLAR_CLOSURE_DELTA_ABS_MAX
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn li7_depletion_default_has_close_best_match() {
        let report = evaluate_lithium7_stellar_depletion_default();
        assert!(report.best_match.closure_delta.abs() < 1.0e-3);
    }
}
