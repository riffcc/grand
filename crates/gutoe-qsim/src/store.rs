use std::collections::HashMap;
use std::sync::Arc;

use crate::state::LatticeState;

/// Content-addressed state store with BLAKE3-keyed deduplication.
///
/// Inserting a state whose BLAKE3 hash already exists returns the existing
/// `Arc<LatticeState>` — zero-copy dedup. BLAKE3 collision probability at
/// 10^7 states is ~10^-54, so no collision handling is needed.
pub struct StateStore {
    map: HashMap<[u8; 32], Arc<LatticeState>>,
    total_inserts: u64,
    dedup_hits: u64,
}

impl StateStore {
    pub fn new() -> Self {
        Self {
            map: HashMap::new(),
            total_inserts: 0,
            dedup_hits: 0,
        }
    }

    /// Insert a `LatticeState`, returning the canonical `Arc`.
    /// If the hash already exists, the input is dropped and the existing Arc returned (dedup hit).
    pub fn insert(&mut self, state: LatticeState) -> Arc<LatticeState> {
        self.total_inserts += 1;
        let hash = *state.hash();
        if let Some(existing) = self.map.get(&hash) {
            self.dedup_hits += 1;
            Arc::clone(existing)
        } else {
            let arc = Arc::new(state);
            self.map.insert(hash, Arc::clone(&arc));
            arc
        }
    }

    /// Insert from an existing `Arc<LatticeState>`.
    pub fn insert_arc(&mut self, arc: Arc<LatticeState>) -> Arc<LatticeState> {
        self.total_inserts += 1;
        let hash = *arc.hash();
        if let Some(existing) = self.map.get(&hash) {
            self.dedup_hits += 1;
            Arc::clone(existing)
        } else {
            self.map.insert(hash, Arc::clone(&arc));
            arc
        }
    }

    #[inline]
    pub fn unique_count(&self) -> usize {
        self.map.len()
    }

    #[inline]
    pub fn total_inserts(&self) -> u64 {
        self.total_inserts
    }

    #[inline]
    pub fn dedup_hits(&self) -> u64 {
        self.dedup_hits
    }

    /// Fraction of inserts that were dedup hits (0.0 = no dedup, 1.0 = all dedup).
    pub fn dedup_ratio(&self) -> f64 {
        if self.total_inserts == 0 {
            0.0
        } else {
            self.dedup_hits as f64 / self.total_inserts as f64
        }
    }

    /// Check if a hash already exists in the store.
    pub fn contains(&self, hash: &[u8; 32]) -> bool {
        self.map.contains_key(hash)
    }

    /// Get canonical Arc for a hash, if present.
    pub fn get(&self, hash: &[u8; 32]) -> Option<&Arc<LatticeState>> {
        self.map.get(hash)
    }
}

impl Default for StateStore {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dedup_returns_same_arc() {
        let mut store = StateStore::new();
        let s1 = LatticeState::new(vec![1, 2, 3]);
        let s2 = LatticeState::new(vec![1, 2, 3]);

        let a1 = store.insert(s1);
        let a2 = store.insert(s2);

        assert!(Arc::ptr_eq(&a1, &a2), "dedup must return same Arc");
        assert_eq!(store.unique_count(), 1);
        assert_eq!(store.dedup_hits(), 1);
        assert_eq!(store.total_inserts(), 2);
    }

    #[test]
    fn different_states_stored_separately() {
        let mut store = StateStore::new();
        let a1 = store.insert(LatticeState::new(vec![0; 10]));
        let a2 = store.insert(LatticeState::new(vec![1; 10]));

        assert!(!Arc::ptr_eq(&a1, &a2));
        assert_eq!(store.unique_count(), 2);
        assert_eq!(store.dedup_hits(), 0);
    }

    #[test]
    fn insert_arc_dedup() {
        let mut store = StateStore::new();
        let arc1 = store.insert(LatticeState::new(vec![7; 50]));
        let arc2 = Arc::new(LatticeState::new(vec![7; 50]));

        let result = store.insert_arc(arc2);
        assert!(Arc::ptr_eq(&arc1, &result), "insert_arc must dedup against existing");
    }

    #[test]
    fn dedup_ratio_correct() {
        let mut store = StateStore::new();
        assert_eq!(store.dedup_ratio(), 0.0);

        store.insert(LatticeState::new(vec![1]));
        store.insert(LatticeState::new(vec![1]));
        store.insert(LatticeState::new(vec![1]));
        store.insert(LatticeState::new(vec![2]));

        // 4 inserts, 2 dedup hits → 0.5
        assert!((store.dedup_ratio() - 0.5).abs() < 1e-10);
    }
}
