use std::sync::Arc;

use gutoe_em::LatticeConfig;

use crate::gate::Gate;
use crate::state::{CowLattice, LatticeState};
use crate::store::StateStore;

/// Independent simulation branch holding a CoW lattice handle.
///
/// Each actor starts as an O(1) clone of the initial state (shared Arc).
/// Applying a gate produces a new state; dedup may collapse it back to a
/// canonical Arc already in the store.
pub struct Actor {
    id: u64,
    lattice: CowLattice,
    history: Vec<[u8; 32]>,
}

impl Actor {
    pub fn new(id: u64, lattice: CowLattice) -> Self {
        let hash = *lattice.hash();
        Self {
            id,
            lattice,
            history: vec![hash],
        }
    }

    #[inline]
    pub fn id(&self) -> u64 {
        self.id
    }

    #[inline]
    pub fn hash(&self) -> &[u8; 32] {
        self.lattice.hash()
    }

    #[inline]
    pub fn lattice(&self) -> &CowLattice {
        &self.lattice
    }

    pub fn history(&self) -> &[[u8; 32]] {
        &self.history
    }

    /// Apply a gate to produce a new lattice state.
    /// The old state is replaced; the new hash is appended to history.
    pub fn apply_gate(&mut self, gate: &dyn Gate, cfg: &LatticeConfig) {
        let new_sites = gate.apply(self.lattice.sites(), cfg);
        self.lattice = CowLattice::new(LatticeState::new(new_sites));
        self.history.push(*self.lattice.hash());
    }

    /// Attempt to deduplicate this actor's state against the store.
    /// If the hash already exists, the actor's inner Arc is swapped to
    /// the canonical one (saving memory). Returns true if dedup occurred.
    pub fn dedup(&mut self, store: &mut StateStore) -> bool {
        let canonical = store.insert_arc(Arc::clone(self.lattice.arc()));
        if Arc::ptr_eq(self.lattice.arc(), &canonical) {
            false // already canonical or first insertion
        } else {
            self.lattice.swap_arc(canonical);
            true
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gate::Z3CycleGate;
    use gutoe_em::{LatticeConfig, VOID};

    #[test]
    fn actor_tracks_history() {
        let cfg = LatticeConfig::default();
        let sites = vec![3u8; cfg.n_sites()]; // all QUARK_SEED
        let cow = CowLattice::new(LatticeState::new(sites));
        let mut actor = Actor::new(0, cow);

        assert_eq!(actor.history().len(), 1);

        let gate = Z3CycleGate;
        actor.apply_gate(&gate, &cfg);
        assert_eq!(actor.history().len(), 2);
        assert_ne!(actor.history()[0], actor.history()[1], "Z3 should change hash");

        actor.apply_gate(&gate, &cfg);
        assert_eq!(actor.history().len(), 3);

        // Z3³ = identity, so applying gate 3 times total should return to start
        actor.apply_gate(&gate, &cfg);
        assert_eq!(actor.history().len(), 4);
        assert_eq!(
            actor.history()[0], actor.history()[3],
            "Z3³ should return to initial state"
        );
    }

    #[test]
    fn actor_dedup_shares_arc() {
        let cfg = LatticeConfig::default();
        let sites = vec![VOID; cfg.n_sites()];
        let state = LatticeState::new(sites.clone());
        let cow1 = CowLattice::new(LatticeState::new(sites));

        let mut store = StateStore::new();
        let canonical = store.insert(state);

        let mut actor = Actor::new(0, cow1);
        let deduped = actor.dedup(&mut store);
        assert!(deduped, "should dedup against store entry");
        assert!(
            Arc::ptr_eq(actor.lattice().arc(), &canonical),
            "actor should now hold canonical Arc"
        );
    }
}
