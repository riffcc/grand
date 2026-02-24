//! Energy and neutrino-loss accounting for stellar reactions.

use crate::{RateEngine, ReactionNetwork, Species};

#[derive(Debug, Clone, PartialEq)]
pub struct ReactionPowerRow {
    pub reaction_id: String,
    pub channel: String,
    pub rate: f64,
    pub q_mev: f64,
    pub neutrino_loss_fraction: f64,
    pub thermal_power: f64,
    pub neutrino_power: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PowerBudget {
    pub rows: Vec<ReactionPowerRow>,
    pub total_thermal_power: f64,
    pub total_neutrino_power: f64,
    pub total_power: f64,
}

pub fn neutrino_loss_fraction(reaction_products: &[crate::stellar_reactions::Stoich]) -> f64 {
    if reaction_products
        .iter()
        .any(|s| s.species == Species::ElectronNeutrino && s.coeff > 0)
    {
        // Baseline placeholder until spectrum-integrated neutrino transport is wired.
        0.35
    } else {
        0.0
    }
}

pub fn compute_power_budget(net: &ReactionNetwork, rates: &RateEngine, t9: f64) -> PowerBudget {
    let mut rows = Vec::with_capacity(net.reactions.len());
    let mut total_thermal_power = 0.0_f64;
    let mut total_neutrino_power = 0.0_f64;

    for r in &net.reactions {
        let rate = rates.rate_for(r.id, t9).unwrap_or(0.0) * r.branching_weight.max(0.0);
        let q = r.q_mev.max(0.0);
        let nu_frac = neutrino_loss_fraction(&r.products).clamp(0.0, 1.0);
        let gross = rate * q;
        let neutrino_power = gross * nu_frac;
        let thermal_power = gross - neutrino_power;

        total_thermal_power += thermal_power;
        total_neutrino_power += neutrino_power;
        rows.push(ReactionPowerRow {
            reaction_id: r.id.to_string(),
            channel: r.channel.to_string(),
            rate,
            q_mev: q,
            neutrino_loss_fraction: nu_frac,
            thermal_power,
            neutrino_power,
        });
    }

    PowerBudget {
        rows,
        total_thermal_power,
        total_neutrino_power,
        total_power: total_thermal_power + total_neutrino_power,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn budget_is_non_negative_and_balanced() {
        let net = ReactionNetwork::baseline_p1();
        let rates = RateEngine::baseline_p1();
        let b = compute_power_budget(&net, &rates, 0.02);
        assert!(b.total_power >= 0.0);
        assert!(b.total_thermal_power >= 0.0);
        assert!(b.total_neutrino_power >= 0.0);
        let lhs = b.total_thermal_power + b.total_neutrino_power;
        assert!((lhs - b.total_power).abs() < 1e-12);
    }

    #[test]
    fn neutrino_channels_have_nonzero_loss() {
        let net = ReactionNetwork::baseline_p1();
        let rates = RateEngine::baseline_p1();
        let b = compute_power_budget(&net, &rates, 0.02);
        let any_nu = b.rows.iter().any(|r| {
            (r.reaction_id == "pp_1" || r.reaction_id == "cno_2" || r.reaction_id == "cno_5")
                && r.neutrino_loss_fraction > 0.0
        });
        assert!(any_nu);
    }
}
