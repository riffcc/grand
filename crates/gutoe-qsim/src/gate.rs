use std::collections::HashSet;

use gutoe_em::{mesh_neighbours, site_coords, LatticeConfig, LEPTON_SEED, VOID};

/// Z₃ cycle table: bit rotation `b0b1b2b3 -> b0|b3<<1|b1<<2|b2<<3` in Cl(1,3).
/// Order 3: Z3(Z3(Z3(s))) = s. VOID (0) is a fixed point.
/// Identical to `gutoe_em::sim::Z3_TABLE` (which is `const` not `pub`).
const Z3_TABLE: [u8; 17] = {
    let mut t = [0u8; 17];
    let mut s = 1u8;
    while s <= 16 {
        let mi = s - 1;
        let b0 = (mi >> 0) & 1;
        let b1 = (mi >> 1) & 1;
        let b2 = (mi >> 2) & 1;
        let b3 = (mi >> 3) & 1;
        t[s as usize] = (b0 | (b3 << 1) | (b1 << 2) | (b2 << 3)) + 1;
        s += 1;
    }
    t
};

/// Deterministic lattice transformation.
///
/// Contract: `apply(sites, cfg)` is a pure function — same input always produces
/// same output. This is what makes content-addressed dedup work: identical starting
/// states + same deterministic gate = one new state, not many.
pub trait Gate: Send + Sync {
    fn name(&self) -> &str;
    fn apply(&self, sites: &[u8], cfg: &LatticeConfig) -> Vec<u8>;
}

/// Z₃ cycle on all non-void, non-lepton sites.
///
/// Order 3: applying 3 times returns to identity. Starting from one state,
/// Z3 produces at most 3 unique states (identity + Z3 + Z3²), then saturates.
pub struct Z3CycleGate;

impl Gate for Z3CycleGate {
    fn name(&self) -> &str {
        "Z3Cycle"
    }

    fn apply(&self, sites: &[u8], _cfg: &LatticeConfig) -> Vec<u8> {
        sites
            .iter()
            .map(|&s| {
                if s == VOID || s == LEPTON_SEED {
                    s
                } else {
                    Z3_TABLE[s as usize]
                }
            })
            .collect()
    }
}

/// Clifford XOR product between specified site pairs.
///
/// For each `(a, b)` pair: `new[a] = ((sites[a]-1) ^ (sites[b]-1)) + 1`.
/// Both sites must be non-void, non-lepton; otherwise pair is skipped.
/// Creates more diversity than Z3: quark closure under Z3+XOR is
/// {3,5,7,9,11,13,15} (7 states).
pub struct XorProductGate {
    pairs: Vec<(usize, usize)>,
}

impl XorProductGate {
    pub fn new(pairs: Vec<(usize, usize)>) -> Self {
        Self { pairs }
    }
}

impl Gate for XorProductGate {
    fn name(&self) -> &str {
        "XorProduct"
    }

    fn apply(&self, sites: &[u8], _cfg: &LatticeConfig) -> Vec<u8> {
        let mut out = sites.to_vec();
        for &(a, b) in &self.pairs {
            let sa = sites[a];
            let sb = sites[b];
            if sa != VOID && sa != LEPTON_SEED && sb != VOID && sb != LEPTON_SEED {
                out[a] = ((sa - 1) ^ (sb - 1)) + 1;
            }
        }
        out
    }
}

/// Deterministic EM hop: each lepton moves toward the max-φ neighbour.
///
/// Ties broken by lowest site index (deterministic). Proton sites and
/// other lepton sites are excluded as hop targets (same rules as sim::step Pass 2).
/// Uses a frozen `phi` snapshot — no gauge evolution within the gate.
pub struct EmHopGate {
    phi: Vec<f64>,
    proton_sites: HashSet<usize>,
}

impl EmHopGate {
    pub fn new(phi: Vec<f64>, proton_sites: HashSet<usize>) -> Self {
        Self { phi, proton_sites }
    }
}

impl Gate for EmHopGate {
    fn name(&self) -> &str {
        "EmHop"
    }

    fn apply(&self, sites: &[u8], cfg: &LatticeConfig) -> Vec<u8> {
        let mut out = sites.to_vec();

        // Process leptons in ascending site order for determinism
        for site in 0..sites.len() {
            if sites[site] != LEPTON_SEED {
                continue;
            }
            let (r, c, z) = site_coords(site, cfg);
            let nbrs = mesh_neighbours(r, c, z, cfg);

            // Find max-φ non-lepton, non-proton neighbour (ties: lowest index)
            let mut best_phi = f64::NEG_INFINITY;
            let mut best_nb = None;
            for nb in nbrs {
                if out[nb] == LEPTON_SEED || self.proton_sites.contains(&nb) {
                    continue;
                }
                let p = self.phi[nb];
                if p > best_phi || (p == best_phi && best_nb.map_or(true, |prev| nb < prev)) {
                    best_phi = p;
                    best_nb = Some(nb);
                }
            }

            if let Some(target) = best_nb {
                let displaced = out[target];
                out[site] = displaced;
                out[target] = LEPTON_SEED;
            }
        }
        out
    }
}

/// Chains multiple gates in sequence: out = g_n(g_{n-1}(...g_1(sites)...)).
pub struct CompositeGate {
    gates: Vec<Box<dyn Gate>>,
    name: String,
}

impl CompositeGate {
    pub fn new(gates: Vec<Box<dyn Gate>>) -> Self {
        let name = gates
            .iter()
            .map(|g| g.name())
            .collect::<Vec<_>>()
            .join("+");
        Self { gates, name }
    }
}

impl Gate for CompositeGate {
    fn name(&self) -> &str {
        &self.name
    }

    fn apply(&self, sites: &[u8], cfg: &LatticeConfig) -> Vec<u8> {
        let mut current = sites.to_vec();
        for gate in &self.gates {
            current = gate.apply(&current, cfg);
        }
        current
    }
}

/// One complete deterministic sim step.
///
/// Equivalent to `sim::step` but fully deterministic (no RNG):
/// 1. Z₃ cycle on all quarks
/// 2. XOR product on pre-determined pairs
/// 3. Alignment: majority vote on all quarks (>k votes)
/// 4. EM hop using frozen φ
///
/// Site lists are pre-determined at construction time.
pub struct FullStepGate {
    xor_pairs: Vec<(usize, usize)>,
    phi: Vec<f64>,
    proton_sites: HashSet<usize>,
}

impl FullStepGate {
    /// Build a full step gate from lattice state + gauge fields.
    ///
    /// `xor_pairs`: deterministic pairs for Clifford XOR (e.g. every quark
    /// paired with its first active hex neighbour).
    pub fn new(
        xor_pairs: Vec<(usize, usize)>,
        phi: Vec<f64>,
        proton_sites: HashSet<usize>,
    ) -> Self {
        Self {
            xor_pairs,
            phi,
            proton_sites,
        }
    }

    /// Build xor_pairs deterministically from a lattice: each quark pairs
    /// with its lowest-index active (non-void, non-lepton) hex neighbour.
    pub fn deterministic_xor_pairs(sites: &[u8], cfg: &LatticeConfig) -> Vec<(usize, usize)> {
        let mut pairs = Vec::new();
        for (i, &s) in sites.iter().enumerate() {
            if s == VOID || s == LEPTON_SEED {
                continue;
            }
            let (r, c, z) = site_coords(i, cfg);
            let nbrs = mesh_neighbours(r, c, z, cfg);
            for nb in nbrs {
                let ns = sites[nb];
                if ns != VOID && ns != LEPTON_SEED {
                    pairs.push((i, nb));
                    break; // first active neighbour only
                }
            }
        }
        pairs
    }
}

impl Gate for FullStepGate {
    fn name(&self) -> &str {
        "FullStep"
    }

    fn apply(&self, sites: &[u8], cfg: &LatticeConfig) -> Vec<u8> {
        // Pass 1a: Z₃ cycle on all quarks
        let z3 = Z3CycleGate;
        let after_z3 = z3.apply(sites, cfg);

        // Pass 1b: XOR product
        let xor = XorProductGate::new(self.xor_pairs.clone());
        let after_xor = xor.apply(&after_z3, cfg);

        // Pass 1c: Alignment (majority vote with k=void_votes threshold)
        let after_align = alignment_step(&after_xor, cfg);

        // Pass 2: EM hop
        let em = EmHopGate::new(self.phi.clone(), self.proton_sites.clone());
        em.apply(&after_align, cfg)
    }
}

/// Deterministic alignment step: majority vote among active hex neighbours.
/// If the most popular neighbour state has > cfg.void_votes votes, adopt it.
fn alignment_step(sites: &[u8], cfg: &LatticeConfig) -> Vec<u8> {
    let mut out = sites.to_vec();
    for i in 0..sites.len() {
        let s = sites[i];
        if s == VOID || s == LEPTON_SEED {
            continue;
        }
        let (r, c, z) = site_coords(i, cfg);
        let nbrs = mesh_neighbours(r, c, z, cfg);

        // Count votes for each non-void, non-lepton neighbour state
        let mut counts = [0u8; 17]; // index by state value 1..=16
        for nb in nbrs {
            let ns = sites[nb];
            if ns != VOID && ns != LEPTON_SEED {
                counts[ns as usize] += 1;
            }
        }

        // Find winner (ties: lowest state value — deterministic)
        let mut best_state = 0u8;
        let mut best_count = 0u8;
        for state_val in 1u8..=16 {
            if counts[state_val as usize] > best_count {
                best_count = counts[state_val as usize];
                best_state = state_val;
            }
        }

        if best_count as usize > cfg.void_votes && best_state != 0 {
            out[i] = best_state;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use gutoe_em::config::QUARK_SEED;

    #[test]
    fn z3_table_matches_gutoe_em() {
        // Verify our local Z3_TABLE is identical to gutoe_em's by testing key properties
        assert_eq!(Z3_TABLE[0], VOID, "VOID is fixed point");
        assert_eq!(Z3_TABLE[QUARK_SEED as usize], 5, "QUARK_SEED=3 → 5");

        // Z3 is order-3 on all non-void states (same assertion as gutoe_em::sim::tests)
        for s in 1u8..=16 {
            let s2 = Z3_TABLE[s as usize];
            let s3 = Z3_TABLE[s2 as usize];
            let s4 = Z3_TABLE[s3 as usize];
            assert_eq!(s4, s, "Z3 not order-3 on state {s}");
        }
    }

    #[test]
    fn z3_gate_deterministic() {
        let cfg = LatticeConfig::default();
        let gate = Z3CycleGate;
        let sites = vec![QUARK_SEED; cfg.n_sites()];

        let out1 = gate.apply(&sites, &cfg);
        let out2 = gate.apply(&sites, &cfg);
        assert_eq!(out1, out2, "Z3 gate must be deterministic");
    }

    #[test]
    fn z3_gate_order_3() {
        let cfg = LatticeConfig::default();
        let gate = Z3CycleGate;
        let sites = vec![QUARK_SEED; cfg.n_sites()];

        let s1 = gate.apply(&sites, &cfg);
        let s2 = gate.apply(&s1, &cfg);
        let s3 = gate.apply(&s2, &cfg);
        assert_eq!(s3, sites, "Z3³ must be identity");
    }

    #[test]
    fn z3_preserves_void_and_lepton() {
        let cfg = LatticeConfig::default();
        let gate = Z3CycleGate;
        let mut sites = vec![VOID; cfg.n_sites()];
        sites[0] = LEPTON_SEED;
        sites[100] = QUARK_SEED;

        let out = gate.apply(&sites, &cfg);
        assert_eq!(out[0], LEPTON_SEED, "lepton must be preserved");
        for i in 1..cfg.n_sites() {
            if i == 100 {
                assert_ne!(out[i], VOID, "quark must be transformed");
            } else {
                assert_eq!(out[i], VOID, "void must be preserved");
            }
        }
    }

    #[test]
    fn xor_product_deterministic() {
        let cfg = LatticeConfig::default();
        let mut sites = vec![VOID; cfg.n_sites()];
        sites[0] = 3; // QUARK_SEED
        sites[1] = 5;
        let gate = XorProductGate::new(vec![(0, 1)]);

        let out1 = gate.apply(&sites, &cfg);
        let out2 = gate.apply(&sites, &cfg);
        assert_eq!(out1, out2);
        // (3-1)^(5-1) + 1 = 2^4 + 1 = 6 + 1 = 7
        assert_eq!(out1[0], ((3 - 1) ^ (5 - 1)) + 1);
        assert_eq!(out1[1], 5, "site_b unchanged");
    }

    #[test]
    fn xor_product_skips_void_lepton() {
        let cfg = LatticeConfig::default();
        let mut sites = vec![VOID; cfg.n_sites()];
        sites[0] = VOID;
        sites[1] = 5;
        sites[2] = LEPTON_SEED;
        sites[3] = 7;

        let gate = XorProductGate::new(vec![(0, 1), (2, 3)]);
        let out = gate.apply(&sites, &cfg);
        assert_eq!(out[0], VOID, "void site_a must be skipped");
        assert_eq!(out[2], LEPTON_SEED, "lepton site_a must be skipped");
    }

    #[test]
    fn composite_chains_gates() {
        let cfg = LatticeConfig::default();
        let sites = vec![QUARK_SEED; cfg.n_sites()];

        // Z3 twice = Z3²
        let gate = CompositeGate::new(vec![Box::new(Z3CycleGate), Box::new(Z3CycleGate)]);
        let out = gate.apply(&sites, &cfg);

        // Verify it equals Z3 applied twice
        let z3 = Z3CycleGate;
        let s1 = z3.apply(&sites, &cfg);
        let s2 = z3.apply(&s1, &cfg);
        assert_eq!(out, s2);
        assert_eq!(gate.name(), "Z3Cycle+Z3Cycle");
    }

    #[test]
    fn em_hop_deterministic() {
        let cfg = LatticeConfig::default();
        let n = cfg.n_sites();
        let mut sites = vec![QUARK_SEED; n];
        sites[100] = LEPTON_SEED;

        // φ gradient pointing toward site 101
        let mut phi = vec![0.0; n];
        let (r, c, z) = site_coords(100, &cfg);
        let nbrs = mesh_neighbours(r, c, z, &cfg);
        for (i, &nb) in nbrs.iter().enumerate() {
            phi[nb] = (i + 1) as f64;
        }

        let gate = EmHopGate::new(phi.clone(), HashSet::new());
        let out1 = gate.apply(&sites, &cfg);
        let out2 = gate.apply(&sites, &cfg);
        assert_eq!(out1, out2, "EM hop must be deterministic");

        // Lepton should have moved to max-phi neighbour
        assert_ne!(out1[100], LEPTON_SEED, "lepton should have hopped away");
        let max_nb = nbrs.last().copied().unwrap();
        assert_eq!(out1[max_nb], LEPTON_SEED, "lepton should be at max-phi neighbour");
    }

    #[test]
    fn em_hop_excludes_proton_sites() {
        let cfg = LatticeConfig::default();
        let n = cfg.n_sites();
        let mut sites = vec![QUARK_SEED; n];
        sites[100] = LEPTON_SEED;

        let (r, c, z) = site_coords(100, &cfg);
        let nbrs = mesh_neighbours(r, c, z, &cfg);

        // Make all neighbours have positive phi, but mark the max-phi one as proton
        let mut phi = vec![0.0; n];
        for (i, &nb) in nbrs.iter().enumerate() {
            phi[nb] = (i + 1) as f64;
        }

        let max_nb = *nbrs.last().unwrap();
        let proton_sites: HashSet<usize> = [max_nb].into_iter().collect();

        let gate = EmHopGate::new(phi, proton_sites);
        let out = gate.apply(&sites, &cfg);
        assert_ne!(out[max_nb], LEPTON_SEED, "lepton must not hop to proton site");
    }
}
