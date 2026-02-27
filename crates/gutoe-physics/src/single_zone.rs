//! Single-zone stellar core burn toy model for baseline validation.

use std::collections::HashMap;

use crate::{RateEngine, ReactionNetwork, Species};

#[derive(Debug, Clone, PartialEq)]
pub struct ZoneState {
    pub abund: HashMap<Species, f64>,
    pub thermal_power: f64,
}

impl ZoneState {
    pub fn solar_like_seed() -> Self {
        let mut abund = HashMap::new();
        abund.insert(Species::P1, 0.70);
        abund.insert(Species::He4, 0.28);
        abund.insert(Species::C12, 0.01);
        abund.insert(Species::N14, 0.005);
        abund.insert(Species::O15, 0.005);
        Self {
            abund,
            thermal_power: 0.0,
        }
    }

    pub fn get(&self, s: Species) -> f64 {
        *self.abund.get(&s).unwrap_or(&0.0)
    }
}

#[derive(Debug, Clone)]
pub struct SingleZoneBurn {
    pub net: ReactionNetwork,
    pub rates: RateEngine,
}

impl SingleZoneBurn {
    pub fn baseline() -> Self {
        Self {
            net: ReactionNetwork::baseline_p1(),
            rates: RateEngine::baseline_p1(),
        }
    }

    pub fn step(&self, state: &mut ZoneState, t9: f64, dt: f64) {
        let dt = dt.max(0.0);
        let mut delta: HashMap<Species, f64> = HashMap::new();
        let mut power = 0.0_f64;

        for r in &self.net.reactions {
            let base_rate = self.rates.rate_for(r.id, t9).unwrap_or(0.0) * r.branching_weight;
            let abund_factor = r.reactants.iter().fold(1.0_f64, |acc, st| {
                acc * state.get(st.species).powi(st.coeff.max(0))
            });
            let flux = base_rate * abund_factor * dt;

            for st in &r.reactants {
                *delta.entry(st.species).or_insert(0.0) -= flux * st.coeff as f64;
            }
            for st in &r.products {
                *delta.entry(st.species).or_insert(0.0) += flux * st.coeff as f64;
            }

            let neutrino_frac = if r
                .products
                .iter()
                .any(|s| s.species == Species::ElectronNeutrino)
            {
                0.35
            } else {
                0.0
            };
            power += flux * r.q_mev * (1.0 - neutrino_frac);
        }

        for (s, d) in delta {
            let next = (state.get(s) + d).max(0.0);
            if next > 0.0 {
                state.abund.insert(s, next);
            } else {
                state.abund.remove(&s);
            }
        }

        // Renormalize abundance sum to 1 for this toy single-zone model.
        let total = state.abund.values().sum::<f64>();
        if total > 1e-12 {
            for v in state.abund.values_mut() {
                *v /= total;
            }
        }
        state.thermal_power = power.max(0.0);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn solar_like_core_has_positive_power() {
        let burn = SingleZoneBurn::baseline();
        let mut s = ZoneState::solar_like_seed();
        burn.step(&mut s, 0.02, 1.0e6);
        assert!(s.thermal_power >= 0.0);
    }

    #[test]
    fn hydrogen_depletes_over_many_steps() {
        let burn = SingleZoneBurn::baseline();
        let mut s = ZoneState::solar_like_seed();
        let h0 = s.get(Species::P1);
        for _ in 0..100 {
            burn.step(&mut s, 0.02, 1.0e6);
        }
        let h1 = s.get(Species::P1);
        assert!(h1 <= h0);
    }
}
