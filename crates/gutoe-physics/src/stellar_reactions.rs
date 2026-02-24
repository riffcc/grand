//! Baseline stellar reaction graph model for physically simulated stars.
//!
//! Scope:
//! - P1 network scaffold with pp-chain, CNO-cycle, and triple-alpha anchors.
//! - Explicit conservation checks (baryon number and charge) per reaction.
//! - Graph-like query helpers for downstream rate engines/integrators.

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Species {
    P1,
    H2,
    He3,
    He4,
    Be7,
    B8,
    C12,
    N13,
    C13,
    N14,
    O15,
    N15,
    Positron,
    ElectronNeutrino,
    Gamma,
}

impl Species {
    pub fn baryon_number(self) -> i32 {
        match self {
            Species::P1 => 1,
            Species::H2 => 2,
            Species::He3 => 3,
            Species::He4 => 4,
            Species::Be7 => 7,
            Species::B8 => 8,
            Species::C12 => 12,
            Species::N13 => 13,
            Species::C13 => 13,
            Species::N14 => 14,
            Species::O15 => 15,
            Species::N15 => 15,
            Species::Positron | Species::ElectronNeutrino | Species::Gamma => 0,
        }
    }

    pub fn charge(self) -> i32 {
        match self {
            Species::P1 => 1,
            Species::H2 => 1,
            Species::He3 => 2,
            Species::He4 => 2,
            Species::Be7 => 4,
            Species::B8 => 5,
            Species::C12 => 6,
            Species::N13 => 7,
            Species::C13 => 6,
            Species::N14 => 7,
            Species::O15 => 8,
            Species::N15 => 7,
            Species::Positron => 1,
            Species::ElectronNeutrino | Species::Gamma => 0,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Stoich {
    pub species: Species,
    pub coeff: i32,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Reaction {
    pub id: &'static str,
    pub channel: &'static str,
    pub reactants: Vec<Stoich>,
    pub products: Vec<Stoich>,
    pub branching_weight: f64,
    pub q_mev: f64,
}

impl Reaction {
    pub fn baryon_balance(&self) -> i32 {
        let lhs = self
            .reactants
            .iter()
            .map(|s| s.coeff * s.species.baryon_number())
            .sum::<i32>();
        let rhs = self
            .products
            .iter()
            .map(|s| s.coeff * s.species.baryon_number())
            .sum::<i32>();
        rhs - lhs
    }

    pub fn charge_balance(&self) -> i32 {
        let lhs = self
            .reactants
            .iter()
            .map(|s| s.coeff * s.species.charge())
            .sum::<i32>();
        let rhs = self
            .products
            .iter()
            .map(|s| s.coeff * s.species.charge())
            .sum::<i32>();
        rhs - lhs
    }

    pub fn is_conserved(&self) -> bool {
        self.baryon_balance() == 0 && self.charge_balance() == 0
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct ReactionNetwork {
    pub reactions: Vec<Reaction>,
}

impl ReactionNetwork {
    pub fn baseline_p1() -> Self {
        // Minimal physically-grounded anchors for the first network pass:
        // pp chain (dominant for solar-like stars), CNO-I catalyst loop,
        // and triple-alpha for helium burning.
        let reactions = vec![
            Reaction {
                id: "pp_1",
                channel: "pp",
                reactants: vec![sto(Species::P1, 2)],
                products: vec![
                    sto(Species::H2, 1),
                    sto(Species::Positron, 1),
                    sto(Species::ElectronNeutrino, 1),
                ],
                branching_weight: 1.0,
                q_mev: 1.442,
            },
            Reaction {
                id: "pp_2",
                channel: "pp",
                reactants: vec![sto(Species::H2, 1), sto(Species::P1, 1)],
                products: vec![sto(Species::He3, 1), sto(Species::Gamma, 1)],
                branching_weight: 1.0,
                q_mev: 5.494,
            },
            Reaction {
                id: "pp_3",
                channel: "pp",
                reactants: vec![sto(Species::He3, 2)],
                products: vec![sto(Species::He4, 1), sto(Species::P1, 2)],
                branching_weight: 1.0,
                q_mev: 12.860,
            },
            Reaction {
                id: "cno_1",
                channel: "cno",
                reactants: vec![sto(Species::C12, 1), sto(Species::P1, 1)],
                products: vec![sto(Species::N13, 1), sto(Species::Gamma, 1)],
                branching_weight: 1.0,
                q_mev: 1.944,
            },
            Reaction {
                id: "cno_2",
                channel: "cno",
                reactants: vec![sto(Species::N13, 1)],
                products: vec![
                    sto(Species::C13, 1),
                    sto(Species::Positron, 1),
                    sto(Species::ElectronNeutrino, 1),
                ],
                branching_weight: 1.0,
                q_mev: 2.221,
            },
            Reaction {
                id: "cno_3",
                channel: "cno",
                reactants: vec![sto(Species::C13, 1), sto(Species::P1, 1)],
                products: vec![sto(Species::N14, 1), sto(Species::Gamma, 1)],
                branching_weight: 1.0,
                q_mev: 7.551,
            },
            Reaction {
                id: "cno_4",
                channel: "cno",
                reactants: vec![sto(Species::N14, 1), sto(Species::P1, 1)],
                products: vec![sto(Species::O15, 1), sto(Species::Gamma, 1)],
                branching_weight: 1.0,
                q_mev: 7.297,
            },
            Reaction {
                id: "cno_5",
                channel: "cno",
                reactants: vec![sto(Species::O15, 1)],
                products: vec![
                    sto(Species::N15, 1),
                    sto(Species::Positron, 1),
                    sto(Species::ElectronNeutrino, 1),
                ],
                branching_weight: 1.0,
                q_mev: 2.754,
            },
            Reaction {
                id: "cno_6",
                channel: "cno",
                reactants: vec![sto(Species::N15, 1), sto(Species::P1, 1)],
                products: vec![sto(Species::C12, 1), sto(Species::He4, 1)],
                branching_weight: 1.0,
                q_mev: 4.966,
            },
            Reaction {
                id: "triple_alpha",
                channel: "triple_alpha",
                reactants: vec![sto(Species::He4, 3)],
                products: vec![sto(Species::C12, 1), sto(Species::Gamma, 1)],
                branching_weight: 1.0,
                q_mev: 7.275,
            },
        ];
        Self { reactions }
    }

    pub fn all_conserved(&self) -> bool {
        self.reactions.iter().all(Reaction::is_conserved)
    }

    pub fn channel_reactions<'a>(&'a self, channel: &'a str) -> impl Iterator<Item = &'a Reaction> {
        self.reactions.iter().filter(move |r| r.channel == channel)
    }

    pub fn tracked_species(&self) -> Vec<Species> {
        vec![
            Species::P1,
            Species::H2,
            Species::He3,
            Species::He4,
            Species::C12,
            Species::N13,
            Species::C13,
            Species::N14,
            Species::O15,
            Species::N15,
            Species::Positron,
            Species::ElectronNeutrino,
            Species::Gamma,
        ]
    }

    /// Returns a reaction-by-species stoichiometric matrix (products - reactants).
    pub fn stoichiometric_matrix(&self) -> Vec<Vec<i32>> {
        let sp = self.tracked_species();
        self.reactions
            .iter()
            .map(|r| {
                sp.iter()
                    .map(|s| {
                        let prod = r
                            .products
                            .iter()
                            .filter(|x| x.species == *s)
                            .map(|x| x.coeff)
                            .sum::<i32>();
                        let react = r
                            .reactants
                            .iter()
                            .filter(|x| x.species == *s)
                            .map(|x| x.coeff)
                            .sum::<i32>();
                        prod - react
                    })
                    .collect()
            })
            .collect()
    }
}

fn sto(species: Species, coeff: i32) -> Stoich {
    Stoich { species, coeff }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn baseline_network_contains_required_channels() {
        let net = ReactionNetwork::baseline_p1();
        assert!(net.channel_reactions("pp").count() >= 3);
        assert!(net.channel_reactions("cno").count() >= 6);
        assert!(net.channel_reactions("triple_alpha").count() >= 1);
    }

    #[test]
    fn baseline_network_conserves_baryon_and_charge() {
        let net = ReactionNetwork::baseline_p1();
        assert!(net.all_conserved());
    }

    #[test]
    fn reaction_ids_are_unique() {
        let net = ReactionNetwork::baseline_p1();
        let mut ids: Vec<&str> = net.reactions.iter().map(|r| r.id).collect();
        ids.sort_unstable();
        ids.dedup();
        assert_eq!(ids.len(), net.reactions.len());
    }

    #[test]
    fn stoichiometric_matrix_dimensions_match_network() {
        let net = ReactionNetwork::baseline_p1();
        let m = net.stoichiometric_matrix();
        assert_eq!(m.len(), net.reactions.len());
        assert_eq!(m[0].len(), net.tracked_species().len());
    }
}
