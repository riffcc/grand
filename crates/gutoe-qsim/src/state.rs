use std::sync::Arc;

/// Immutable, content-addressed lattice snapshot.
///
/// BLAKE3 hash is computed eagerly at construction — the hash invariant
/// is enforced structurally because `sites` is private and there are no
/// mutable accessors.
pub struct LatticeState {
    sites: Vec<u8>,
    hash: [u8; 32],
}

impl LatticeState {
    pub fn new(sites: Vec<u8>) -> Self {
        let hash = *blake3::hash(&sites).as_bytes();
        Self { sites, hash }
    }

    #[inline]
    pub fn hash(&self) -> &[u8; 32] {
        &self.hash
    }

    #[inline]
    pub fn sites(&self) -> &[u8] {
        &self.sites
    }

    #[inline]
    pub fn site_at(&self, idx: usize) -> u8 {
        self.sites[idx]
    }
}

/// O(1)-clone lattice handle with copy-on-write mutation.
///
/// Cloning a `CowLattice` is `Arc::clone` — no site data is copied.
/// Mutation always produces a NEW `CowLattice`; the original is unchanged.
#[derive(Clone)]
pub struct CowLattice {
    inner: Arc<LatticeState>,
}

impl CowLattice {
    pub fn new(state: LatticeState) -> Self {
        Self {
            inner: Arc::new(state),
        }
    }

    pub fn from_arc(arc: Arc<LatticeState>) -> Self {
        Self { inner: arc }
    }

    #[inline]
    pub fn hash(&self) -> &[u8; 32] {
        self.inner.hash()
    }

    #[inline]
    pub fn sites(&self) -> &[u8] {
        self.inner.sites()
    }

    #[inline]
    pub fn arc(&self) -> &Arc<LatticeState> {
        &self.inner
    }

    /// Copy-on-write mutation: clones sites, applies `f`, returns a new
    /// `CowLattice` with a fresh BLAKE3 hash. Original is unchanged.
    pub fn mutate<F: FnOnce(&mut Vec<u8>)>(&self, f: F) -> CowLattice {
        let mut sites = self.inner.sites.clone();
        f(&mut sites);
        CowLattice::new(LatticeState::new(sites))
    }

    /// Replace inner Arc with a canonical one (for dedup).
    pub fn swap_arc(&mut self, arc: Arc<LatticeState>) {
        self.inner = arc;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hash_determinism() {
        let sites = vec![0u8, 1, 2, 3, 4, 5, 6, 7, 8, 9];
        let s1 = LatticeState::new(sites.clone());
        let s2 = LatticeState::new(sites);
        assert_eq!(s1.hash(), s2.hash(), "same sites must produce same hash");
    }

    #[test]
    fn hash_sensitivity() {
        let s1 = LatticeState::new(vec![0; 100]);
        let s2 = LatticeState::new(vec![1; 100]);
        assert_ne!(
            s1.hash(),
            s2.hash(),
            "different sites must produce different hash"
        );
    }

    #[test]
    fn cow_clone_shares_arc() {
        let cow = CowLattice::new(LatticeState::new(vec![42; 1728]));
        let cow2 = cow.clone();
        assert!(Arc::ptr_eq(cow.arc(), cow2.arc()), "clone must share Arc");
    }

    #[test]
    fn mutate_produces_new_state() {
        let original = CowLattice::new(LatticeState::new(vec![0; 100]));
        let mutated = original.mutate(|sites| sites[50] = 7);

        assert_eq!(original.sites()[50], 0, "original must be unchanged");
        assert_eq!(mutated.sites()[50], 7, "mutated must reflect change");
        assert_ne!(original.hash(), mutated.hash(), "hashes must differ");
        assert!(
            !Arc::ptr_eq(original.arc(), mutated.arc()),
            "Arcs must differ"
        );
    }

    #[test]
    fn mutate_noop_preserves_hash() {
        let original = CowLattice::new(LatticeState::new(vec![5; 100]));
        let mutated = original.mutate(|_sites| { /* no-op */ });
        assert_eq!(
            original.hash(),
            mutated.hash(),
            "no-op mutate must preserve hash"
        );
    }
}
