/*!
 * GUTOE Physics - Abiogenesis Threshold Lane (Kauffman Closure)
 * Copyright (C) 2026  Riff Labs
 *
 * This lane encodes a theorem-style closure check:
 *   Step 1: derive prebiotic network counts + kinetic lower bounds,
 *   Step 2: evaluate Kauffman closure control `N * p`,
 *   Step 3: apply uncertainty margin and emit binary inevitability gate.
 */

use crate::constants::{
    lambda_micro_finite_mode_rescale, ALPHA, ALPHA_LEADING_ORDER, CLIFFORD_STATE_COUNT_STRUCTURAL,
    DARK_FRACTION_GEOMETRIC_STRUCTURAL,
};
use crate::homochirality::amino_acid_enantiomer_split_ev;

/// Boltzmann constant in eV/K.
pub const KB_EV_PER_K: f64 = 8.617_333_262_145e-5;

/// Kauffman critical closure control threshold (mean catalytic branching).
pub const KAUFFMAN_CLOSURE_THRESHOLD: f64 = 1.0;

/// Structural monomer pool: `16 + 4 = 20`.
pub const STRUCTURAL_MONOMER_COUNT: f64 = 20.0;

/// Structural polymer-length scale used in this lane: `16 + 1 = 17`.
pub const STRUCTURAL_POLYMER_LENGTH: f64 = 17.0;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AbiogenesisWindows {
    pub closure_threshold: f64,
    pub catalytic_probability_min: f64,
    pub robust_margin_min: f64,
}

impl Default for AbiogenesisWindows {
    fn default() -> Self {
        Self {
            closure_threshold: KAUFFMAN_CLOSURE_THRESHOLD,
            catalytic_probability_min: 0.05,
            robust_margin_min: 0.25,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PrebioticNetworkScorecard {
    pub feedstock_species: usize,
    pub amino_acid_pool_left: usize,
    pub nucleotide_pool: usize,
    pub peptide_channels: usize,
    pub nucleotide_synthesis_channels: usize,
    pub phosphodiester_channels: usize,
    pub k_peptide: f64,
    pub k_nucleotide: f64,
    pub k_phosphodiester: f64,
    pub catalytic_probability_lower_bound: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct KauffmanClosureScorecard {
    pub monomer_count: f64,
    pub catalytic_probability: f64,
    pub n_times_p: f64,
    pub threshold: f64,
    pub closure_excess: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InevitabilityScorecard {
    pub pved_delta_e_ev: f64,
    pub thermal_chirality_bias: f64,
    pub alpha_rel_uncertainty: f64,
    pub micro_rel_uncertainty: f64,
    pub network_rel_uncertainty: f64,
    pub total_rel_uncertainty: f64,
    pub n_times_p_lower_3sigma: f64,
    pub robust_margin: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AbiogenesisScorecard {
    pub prebiotic: PrebioticNetworkScorecard,
    pub closure: KauffmanClosureScorecard,
    pub inevitability: InevitabilityScorecard,
    pub prebiotic_ok: bool,
    pub closure_ok: bool,
    pub inevitability_ok: bool,
}

impl AbiogenesisScorecard {
    pub const fn passes_all(self) -> bool {
        self.prebiotic_ok && self.closure_ok && self.inevitability_ok
    }
}

/// Thermal chirality bias from PVED: `tanh(ΔE/(2k_B T))`.
pub fn thermal_chirality_bias(temperature_k: f64) -> f64 {
    let de = amino_acid_enantiomer_split_ev();
    if !(de.is_finite() && temperature_k.is_finite()) || temperature_k <= 0.0 {
        return f64::NAN;
    }
    (de / (2.0 * KB_EV_PER_K * temperature_k)).tanh()
}

/// Step 1: derive prebiotic reaction-network counts and kinetic lower bound.
pub fn evaluate_prebiotic_network() -> PrebioticNetworkScorecard {
    // Canonical prebiotic feedstocks: CH4, NH3, H2O, H2.
    let feedstock_species = 4usize;

    // Chiral lane provides left-biased amino chemistry; use the canonical 20 pool.
    let amino_acid_pool_left = STRUCTURAL_MONOMER_COUNT as usize;
    let nucleotide_pool = 4usize;

    let peptide_channels = amino_acid_pool_left * amino_acid_pool_left;
    let nucleotide_synthesis_channels = nucleotide_pool * nucleotide_pool;
    let phosphodiester_channels = nucleotide_pool * nucleotide_pool;

    // Derived kinetic scales from alpha and geometric-contact factor (60/71).
    // We keep a conservative lower-bound lane by taking the minimum channel scale.
    let contact = DARK_FRACTION_GEOMETRIC_STRUCTURAL;
    let k_peptide = ALPHA_LEADING_ORDER * STRUCTURAL_POLYMER_LENGTH * contact;
    let k_nucleotide = ALPHA_LEADING_ORDER * (STRUCTURAL_POLYMER_LENGTH - 4.0) * contact;
    let k_phosphodiester = ALPHA_LEADING_ORDER * (STRUCTURAL_POLYMER_LENGTH - 6.0) * contact;

    let catalytic_probability_lower_bound = k_peptide.min(k_nucleotide).min(k_phosphodiester);

    PrebioticNetworkScorecard {
        feedstock_species,
        amino_acid_pool_left,
        nucleotide_pool,
        peptide_channels,
        nucleotide_synthesis_channels,
        phosphodiester_channels,
        k_peptide,
        k_nucleotide,
        k_phosphodiester,
        catalytic_probability_lower_bound,
    }
}

/// Step 2: Kauffman closure control `N * p`.
pub fn evaluate_kauffman_closure(
    monomer_count: f64,
    catalytic_probability: f64,
    threshold: f64,
) -> KauffmanClosureScorecard {
    let n_times_p = monomer_count * catalytic_probability;
    let closure_excess = n_times_p - threshold;
    KauffmanClosureScorecard {
        monomer_count,
        catalytic_probability,
        n_times_p,
        threshold,
        closure_excess,
    }
}

/// Step 3: uncertainty-aware inevitability margin.
pub fn evaluate_inevitability(
    closure: KauffmanClosureScorecard,
    temperature_k: f64,
) -> InevitabilityScorecard {
    let pved_delta_e_ev = amino_acid_enantiomer_split_ev();
    let thermal_bias = thermal_chirality_bias(temperature_k).abs();

    let alpha_rel_uncertainty = ((ALPHA - ALPHA_LEADING_ORDER) / ALPHA_LEADING_ORDER).abs();
    let micro_rel_uncertainty = (lambda_micro_finite_mode_rescale() - 1.0).abs();
    let network_rel_uncertainty =
        1.0 / (CLIFFORD_STATE_COUNT_STRUCTURAL * STRUCTURAL_MONOMER_COUNT);

    let total_rel_uncertainty = (alpha_rel_uncertainty.powi(2)
        + micro_rel_uncertainty.powi(2)
        + network_rel_uncertainty.powi(2))
    .sqrt();

    let n_times_p_lower_3sigma = closure.n_times_p * (1.0 - 3.0 * total_rel_uncertainty);
    let robust_margin = n_times_p_lower_3sigma - closure.threshold;

    InevitabilityScorecard {
        pved_delta_e_ev,
        thermal_chirality_bias: thermal_bias,
        alpha_rel_uncertainty,
        micro_rel_uncertainty,
        network_rel_uncertainty,
        total_rel_uncertainty,
        n_times_p_lower_3sigma,
        robust_margin,
    }
}

pub fn evaluate_abiogenesis_gate(
    w: AbiogenesisWindows,
    temperature_k: f64,
) -> AbiogenesisScorecard {
    let prebiotic = evaluate_prebiotic_network();
    let closure = evaluate_kauffman_closure(
        STRUCTURAL_MONOMER_COUNT,
        prebiotic.catalytic_probability_lower_bound,
        w.closure_threshold,
    );
    let inevitability = evaluate_inevitability(closure, temperature_k);

    let prebiotic_ok = prebiotic.feedstock_species == 4
        && prebiotic.amino_acid_pool_left == 20
        && prebiotic.catalytic_probability_lower_bound >= w.catalytic_probability_min
        && prebiotic.peptide_channels > 0
        && prebiotic.nucleotide_synthesis_channels > 0
        && prebiotic.phosphodiester_channels > 0;

    let closure_ok = closure.n_times_p > closure.threshold;

    let inevitability_ok = inevitability.pved_delta_e_ev > 0.0
        && inevitability.n_times_p_lower_3sigma > closure.threshold
        && inevitability.robust_margin >= w.robust_margin_min;

    AbiogenesisScorecard {
        prebiotic,
        closure,
        inevitability,
        prebiotic_ok,
        closure_ok,
        inevitability_ok,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn prebiotic_network_is_enumerated_and_positive() {
        let p = evaluate_prebiotic_network();
        assert_eq!(p.feedstock_species, 4);
        assert_eq!(p.amino_acid_pool_left, 20);
        assert_eq!(p.nucleotide_pool, 4);
        assert!(p.k_peptide > 0.0 && p.k_nucleotide > 0.0 && p.k_phosphodiester > 0.0);
        assert!(p.catalytic_probability_lower_bound > 0.0);
    }

    #[test]
    fn kauffman_threshold_exceeded_in_structural_lane() {
        let p = evaluate_prebiotic_network();
        let c = evaluate_kauffman_closure(
            STRUCTURAL_MONOMER_COUNT,
            p.catalytic_probability_lower_bound,
            KAUFFMAN_CLOSURE_THRESHOLD,
        );
        assert!(
            c.n_times_p > c.threshold,
            "N*p should exceed closure threshold: {c:#?}"
        );
    }

    #[test]
    fn inevitability_gate_passes_default() {
        let s = evaluate_abiogenesis_gate(AbiogenesisWindows::default(), 298.15);
        assert!(s.passes_all(), "abiogenesis gate failed: {s:#?}");
    }
}
