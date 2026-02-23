use std::sync::Arc;

use gutoe_em::LatticeConfig;

use crate::actor::Actor;
use crate::gate::Gate;
use crate::metrics::{ExperimentResult, RoundMetrics};
use crate::state::{CowLattice, LatticeState};
use crate::store::StateStore;

/// Orchestrator: spawn N actors, apply gates to subsets, dedup, measure compression.
///
/// All actors start from a single shared `Arc<LatticeState>` (one insert in the store).
/// Each round: apply gate to targeted actors, dedup ALL actors, record metrics.
pub struct Orchestrator {
    actors: Vec<Actor>,
    store: StateStore,
    cfg: LatticeConfig,
    round_metrics: Vec<RoundMetrics>,
}

impl Orchestrator {
    /// Create N actors all sharing one initial state.
    pub fn new(initial: LatticeState, n_actors: usize, cfg: LatticeConfig) -> Self {
        let mut store = StateStore::new();
        let arc = store.insert(initial);

        let actors = (0..n_actors)
            .map(|i| Actor::new(i as u64, CowLattice::from_arc(Arc::clone(&arc))))
            .collect();

        Self {
            actors,
            store,
            cfg,
            round_metrics: Vec::new(),
        }
    }

    pub fn actors(&self) -> &[Actor] {
        &self.actors
    }

    pub fn store(&self) -> &StateStore {
        &self.store
    }

    pub fn round_metrics(&self) -> &[RoundMetrics] {
        &self.round_metrics
    }

    pub fn cfg(&self) -> &LatticeConfig {
        &self.cfg
    }

    /// Apply a gate to the specified actor indices, then dedup ALL actors.
    ///
    /// Returns per-round metrics and prints a summary line.
    pub fn apply_round(&mut self, gate: &dyn Gate, actor_indices: &[usize]) -> RoundMetrics {
        let round = self.round_metrics.len();
        let unique_before = self.store.unique_count();

        // Apply gate to targeted actors
        for &idx in actor_indices {
            self.actors[idx].apply_gate(gate, &self.cfg);
        }

        // Dedup ALL actors against the store
        for actor in &mut self.actors {
            actor.dedup(&mut self.store);
        }

        let unique_after = self.store.unique_count();
        let new_states = unique_after - unique_before;
        let compression = self.actors.len() as f64 / unique_after as f64;

        let metrics = RoundMetrics {
            round,
            gate_name: gate.name().to_string(),
            actors_targeted: actor_indices.len(),
            unique_states_before: unique_before,
            unique_states_after: unique_after,
            new_states_created: new_states,
            compression_factor: compression,
        };

        println!("  {metrics}");
        self.round_metrics.push(metrics);

        // Return a copy of the metrics for the caller
        let last = self.round_metrics.last().unwrap();
        RoundMetrics {
            round: last.round,
            gate_name: last.gate_name.clone(),
            actors_targeted: last.actors_targeted,
            unique_states_before: last.unique_states_before,
            unique_states_after: last.unique_states_after,
            new_states_created: last.new_states_created,
            compression_factor: last.compression_factor,
        }
    }

    /// Run a complete experiment: apply the same gate to all actors for N rounds.
    pub fn run_uniform(
        &mut self,
        gate: &dyn Gate,
        n_rounds: usize,
    ) -> ExperimentResult {
        let all_indices: Vec<usize> = (0..self.actors.len()).collect();

        for _ in 0..n_rounds {
            self.apply_round(gate, &all_indices);
        }

        self.result()
    }

    /// Compile final results from accumulated round metrics.
    pub fn result(&self) -> ExperimentResult {
        let growth_curve: Vec<usize> = self
            .round_metrics
            .iter()
            .map(|m| m.unique_states_after)
            .collect();

        let peak_compression = self
            .round_metrics
            .iter()
            .map(|m| m.compression_factor)
            .fold(f64::NEG_INFINITY, f64::max);

        let exponent = ExperimentResult::fit_growth_exponent(&growth_curve);

        ExperimentResult {
            total_actors: self.actors.len(),
            total_rounds: self.round_metrics.len(),
            final_unique_states: self.store.unique_count(),
            peak_compression,
            growth_curve,
            growth_exponent: exponent,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gate::{Z3CycleGate, XorProductGate};
    use gutoe_em::{LatticeConfig, VOID};

    #[test]
    fn z3_compression_4000_actors() {
        let cfg = LatticeConfig::default();
        let n = cfg.n_sites();
        // Start with a non-trivial state (all QUARK_SEED)
        let sites = vec![3u8; n];
        let initial = LatticeState::new(sites);

        let mut orch = Orchestrator::new(initial, 4000, cfg);
        let gate = Z3CycleGate;

        println!("test z3_compression_4000_actors ...");
        let result = orch.run_uniform(&gate, 3);

        // Z3 is order-3: identity + Z3 + Z3² = 3 distinct states, then saturates.
        // Round 0: initial -> Z3(initial) = 2 unique
        // Round 1: + Z3²(initial) = 3 unique
        // Round 2: Z3³ = identity (already in store) = 3 unique, 0 new
        assert_eq!(
            result.final_unique_states, 3,
            "Z3 order-3 should produce exactly 3 unique states, got {}",
            result.final_unique_states
        );

        let last_metrics = orch.round_metrics().last().unwrap();
        assert_eq!(last_metrics.new_states_created, 0, "round 2 should create 0 new states (Z3 saturated)");

        // Compression: 4000 actors / 3 states ≈ 1333x
        assert!(
            result.peak_compression > 1000.0,
            "compression should be >1000x, got {:.1}",
            result.peak_compression
        );
    }

    #[test]
    fn z3_half_actors_20_rounds() {
        let cfg = LatticeConfig::default();
        let sites = vec![3u8; cfg.n_sites()];
        let initial = LatticeState::new(sites);

        let mut orch = Orchestrator::new(initial, 4000, cfg);
        let gate = Z3CycleGate;

        println!("test z3_half_actors_20_rounds ...");

        // Apply to first half each round
        let half: Vec<usize> = (0..2000).collect();
        for _ in 0..20 {
            orch.apply_round(&gate, &half);
        }

        let result = orch.result();

        // Even with half-targeting, Z3 has only 3 states in its orbit.
        // Some actors stay at initial, some cycle: at most 3 unique states.
        assert!(
            result.final_unique_states <= 3,
            "Z3 half-apply should still saturate at ≤3 states, got {}",
            result.final_unique_states
        );
        assert!(
            result.peak_compression > 200.0,
            "compression should be >200x, got {:.1}",
            result.peak_compression
        );
    }

    #[test]
    fn xor_creates_diversity() {
        let cfg = LatticeConfig::default();
        let n = cfg.n_sites();

        // Start with a mixed lattice: alternating states 3 and 5
        let mut sites = vec![VOID; n];
        for i in 0..n {
            if i % 2 == 0 {
                sites[i] = 3;
            } else {
                sites[i] = 5;
            }
        }
        let initial = LatticeState::new(sites.clone());

        // XOR pairs: every even site with its right neighbour
        let mut pairs = Vec::new();
        for i in (0..n).step_by(2) {
            if i + 1 < n {
                pairs.push((i, i + 1));
            }
        }

        let mut orch = Orchestrator::new(initial, 100, cfg);
        let gate = XorProductGate::new(pairs);

        println!("test xor_creates_diversity ...");
        let result = orch.run_uniform(&gate, 10);

        // XOR should create some diversity but not exponential
        assert!(
            result.final_unique_states > 1,
            "XOR should create some diversity"
        );
        assert!(
            result.final_unique_states < 100,
            "XOR should not create exponential diversity, got {}",
            result.final_unique_states
        );
    }

    #[test]
    fn all_void_stays_compressed() {
        let cfg = LatticeConfig::default();
        let initial = LatticeState::new(vec![VOID; cfg.n_sites()]);

        let mut orch = Orchestrator::new(initial, 1000, cfg);
        let gate = Z3CycleGate;

        println!("test all_void_stays_compressed ...");
        let result = orch.run_uniform(&gate, 10);

        // Z3 on VOID is identity → all 1000 actors always have the same state
        assert_eq!(
            result.final_unique_states, 1,
            "all-void Z3 should stay at 1 unique state"
        );
        assert_eq!(result.peak_compression, 1000.0);
    }
}
