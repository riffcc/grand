//! GUTOE Cl(1,3) lattice simulation — Rust port of Python gutoe_2d5_toroid.py
//!
//! Two-pass step:
//!   Pass 1 — VOID differentiation, quark Z3/Clifford/alignment (leptons skipped)
//!   Pass 2 — lepton EM hops (grade-1 accessible, proton sites excluded)

use rand::Rng;

// ── Constants ──────────────────────────────────────────────────────────────────

pub const VOID: u8 = 0;
pub const LEPTON_SEED: u8 = 2; // γ⁰  (mi = 0b0001, grade 1)
pub const QUARK_SEED: u8 = 3; // γ¹  (mi = 0b0010, grade 1)

// ── Config ─────────────────────────────────────────────────────────────────────

#[derive(Clone, Debug)]
pub struct LatticeConfig {
    pub hex_rows: usize,
    pub hex_cols: usize,
    pub layers: usize,
    pub differentiation_prob: f64,
    pub cycle_prob: f64,
    pub clifford_prob: f64,
    pub alignment_strength: f64,
    pub quark_threshold: f64,
    pub void_votes: usize,
    pub em_prob: f64,
    pub photon_c: f64,
    pub photon_coupling: f64,
    pub poisson_iters: usize,
}

impl Default for LatticeConfig {
    fn default() -> Self {
        Self {
            hex_rows: 12,
            hex_cols: 12,
            layers: 12,
            differentiation_prob: 0.02,
            cycle_prob: 0.05,
            clifford_prob: 0.03,
            alignment_strength: 0.15,
            quark_threshold: 0.6,
            void_votes: 4,
            em_prob: 0.5,
            photon_c: 0.4,
            photon_coupling: 0.05,
            poisson_iters: 80,
        }
    }
}

// ── Pre-computed tables ────────────────────────────────────────────────────────

/// grade_of(s) for s in 0..=16; s=0 → -1 (VOID sentinel), s=1..=16 → popcount(s-1)
#[allow(dead_code)]
pub fn grade_of(s: u8) -> i8 {
    if s == 0 {
        -1
    } else {
        (s - 1).count_ones() as i8
    }
}

/// Z3 rotation table: b0,b1,b2,b3 → b3,b0,b1,b2 (cyclic on bits)
/// Index: state 0..=16, output: new state
pub fn make_z3_table() -> [u8; 17] {
    let mut t = [VOID; 17];
    for s in 1u8..=16 {
        let mi = s - 1;
        let b0 = (mi >> 0) & 1;
        let b1 = (mi >> 1) & 1;
        let b2 = (mi >> 2) & 1;
        let b3 = (mi >> 3) & 1;
        t[s as usize] = (b0 | (b3 << 1) | (b1 << 2) | (b2 << 3)) + 1;
    }
    t
}

/// Veracity table: similarity between two Clifford states
/// 1.0 if equal, √3/2 if Hamming-1, 0.5 if Hamming-2, 0.0 otherwise
/// VOID(0) with anything → 0.0
pub fn make_veracity_table() -> [[f32; 17]; 17] {
    let sqrt3_half = (3.0f32).sqrt() / 2.0;
    let mut t = [[0.0f32; 17]; 17];
    for s1 in 0usize..17 {
        for s2 in 0usize..17 {
            if s1 == 0 || s2 == 0 {
                t[s1][s2] = 0.0;
            } else if s1 == s2 {
                t[s1][s2] = 1.0;
            } else {
                let d = ((s1 - 1) ^ (s2 - 1)).count_ones();
                t[s1][s2] = match d {
                    1 => sqrt3_half,
                    2 => 0.5,
                    _ => 0.0,
                };
            }
        }
    }
    t
}

pub const Z3_ORBITS: [[u8; 3]; 4] = [[3, 5, 9], [4, 6, 10], [7, 11, 13], [8, 12, 14]];

// ── Geometry ───────────────────────────────────────────────────────────────────

pub fn site_coords(site: usize, cfg: &LatticeConfig) -> (usize, usize, usize) {
    let z = site / (cfg.hex_rows * cfg.hex_cols);
    let rem = site % (cfg.hex_rows * cfg.hex_cols);
    let r = rem / cfg.hex_cols;
    let c = rem % cfg.hex_cols;
    (r, c, z)
}

fn wrap(v: isize, n: usize) -> usize {
    ((v % n as isize + n as isize) as usize) % n
}

fn flat_idx(r: usize, c: usize, z: usize, cfg: &LatticeConfig) -> usize {
    (z * cfg.hex_rows + r) * cfg.hex_cols + c
}

pub fn hex_neighbours_rc(r: usize, c: usize, cfg: &LatticeConfig) -> [(usize, usize); 6] {
    let offsets: [(isize, isize); 6] = if r % 2 == 0 {
        [(-1, 0), (-1, 1), (0, -1), (0, 1), (1, 0), (1, 1)]
    } else {
        [(-1, -1), (-1, 0), (0, -1), (0, 1), (1, -1), (1, 0)]
    };
    let mut result = [(0usize, 0usize); 6];
    for (i, (dr, dc)) in offsets.iter().enumerate() {
        result[i] = (
            wrap(r as isize + dr, cfg.hex_rows),
            wrap(c as isize + dc, cfg.hex_cols),
        );
    }
    result
}

/// Intra-layer hex neighbours (6 sites, same z)
pub fn mesh_neighbours(r: usize, c: usize, z: usize, cfg: &LatticeConfig) -> [usize; 6] {
    let rc_nbrs = hex_neighbours_rc(r, c, cfg);
    let mut out = [0usize; 6];
    for (i, (nr, nc)) in rc_nbrs.iter().enumerate() {
        out[i] = flat_idx(*nr, *nc, z, cfg);
    }
    out
}

// ── State ──────────────────────────────────────────────────────────────────────

pub struct GutoeState {
    pub lattice: Vec<u8>,
    pub step_count: u64,
    pub phase: u8, // 1 = quarks only, 2 = EM active
    // cached tables
    z3_table: [u8; 17],
    veracity: [[f32; 17]; 17],
}

impl GutoeState {
    pub fn new(cfg: &LatticeConfig) -> Self {
        let n = cfg.hex_rows * cfg.hex_cols * cfg.layers;
        Self {
            lattice: vec![VOID; n],
            step_count: 0,
            phase: 1,
            z3_table: make_z3_table(),
            veracity: make_veracity_table(),
        }
    }

    pub fn n(&self, cfg: &LatticeConfig) -> usize {
        cfg.hex_rows * cfg.hex_cols * cfg.layers
    }

    fn local_fields(&self, site: usize, nbrs: &[usize; 6]) -> (f32, f32, f32) {
        let state = self.lattice[site] as usize;
        if state == 0 {
            return (0.0, 0.0, 0.0);
        }
        let mut total_v = 0.0f32;
        let mut grad = 0.0f32;
        let mut nbr_set = [false; 17];
        for &ni in nbrs {
            let ns = self.lattice[ni] as usize;
            let v = self.veracity[state][ns];
            total_v += v;
            grad += 1.0 - v;
            if ns != 0 {
                nbr_set[ns] = true;
            }
        }
        let n = 6.0f32;
        // z3 curvature: max over orbits of (orbit members seen - 1) / 2
        let z3_curv = Z3_ORBITS
            .iter()
            .map(|orbit| {
                let count = orbit.iter().filter(|&&s| nbr_set[s as usize]).count();
                if count > 0 {
                    (count - 1) as f32 / 2.0
                } else {
                    0.0
                }
            })
            .fold(0.0f32, f32::max);
        (total_v / n, z3_curv, grad / n)
    }

    /// Two-pass step. proton_sites must be pre-computed by caller.
    pub fn step<R: Rng>(
        &mut self,
        rng: &mut R,
        cfg: &LatticeConfig,
        phi: &[f64],
        proton_sites: &std::collections::HashSet<usize>,
    ) {
        let n = self.n(cfg);
        let mut new = self.lattice.clone();

        // ── Pass 1: VOID and quark dynamics (skip leptons) ──────────────────────
        for site in 0..n {
            let (r, c, z) = site_coords(site, cfg);
            let nbrs = mesh_neighbours(r, c, z, cfg);
            let state = self.lattice[site];

            if state == VOID {
                // Differentiation
                if rng.gen::<f64>() < cfg.differentiation_prob {
                    new[site] = QUARK_SEED;
                    continue;
                }
                let active: usize = nbrs.iter().filter(|&&ni| self.lattice[ni] != VOID).count();
                let total = nbrs.len();
                if active >= 2.max(total / 4)
                    && rng.gen::<f64>() < active as f64 / total as f64 * 0.4
                {
                    new[site] = QUARK_SEED;
                }
            } else if state == LEPTON_SEED {
                // Skip leptons in pass 1
                continue;
            } else {
                // Quark dynamics
                let r_val: f64 = rng.gen();
                if r_val < cfg.cycle_prob {
                    new[site] = self.z3_table[state as usize];
                } else if r_val < cfg.cycle_prob + cfg.clifford_prob {
                    let active_nbrs: Vec<u8> = nbrs
                        .iter()
                        .filter_map(|&ni| {
                            let ns = self.lattice[ni];
                            if ns != VOID {
                                Some(ns)
                            } else {
                                None
                            }
                        })
                        .collect();
                    if !active_nbrs.is_empty() {
                        let partner = active_nbrs[rng.gen_range(0..active_nbrs.len())];
                        new[site] = ((state - 1) ^ (partner - 1)) + 1;
                    }
                } else if r_val < cfg.cycle_prob + cfg.clifford_prob + cfg.alignment_strength {
                    // Alignment with void votes
                    let nbr_states: Vec<u8> = nbrs
                        .iter()
                        .filter_map(|&ni| {
                            let ns = self.lattice[ni];
                            if ns != VOID {
                                Some(ns)
                            } else {
                                None
                            }
                        })
                        .collect();
                    if !nbr_states.is_empty() {
                        // Count votes
                        let mut counts = [0u8; 17];
                        for &s in &nbr_states {
                            counts[s as usize] += 1;
                        }
                        // Find winner
                        let (winner, winner_count) = counts
                            .iter()
                            .enumerate()
                            .max_by_key(|&(_, &c)| c)
                            .map(|(i, &c)| (i as u8, c as usize))
                            .unwrap();
                        if winner_count > cfg.void_votes {
                            new[site] = winner;
                        }
                    }
                }
            }
        }

        // ── Pass 2: lepton EM hops ─────────────────────────────────────────────
        if self.phase == 2 {
            // Collect lepton sites from new (pass-1 result)
            let lepton_sites: Vec<usize> = (0..n).filter(|&s| new[s] == LEPTON_SEED).collect();

            for site in lepton_sites {
                if new[site] != LEPTON_SEED {
                    continue; // already moved
                }
                if rng.gen::<f64>() >= cfg.em_prob {
                    continue;
                }
                let (r, c, z) = site_coords(site, cfg);
                let nbrs = mesh_neighbours(r, c, z, cfg);

                // Candidates: not another lepton, not a proton site, grade-1 or VOID accessible
                let mut best_phi = f64::NEG_INFINITY;
                let mut best_nb = None;
                for &nb in &nbrs {
                    if new[nb] == LEPTON_SEED {
                        continue;
                    }
                    if proton_sites.contains(&nb) {
                        continue;
                    }
                    let p = phi[nb];
                    if p > best_phi {
                        best_phi = p;
                        best_nb = Some(nb);
                    }
                }
                if let Some(target) = best_nb {
                    let target_state = new[target];
                    new[site] = target_state;
                    new[target] = LEPTON_SEED;
                }
            }
        }

        self.lattice = new;
        self.step_count += 1;
    }
}

// ── Quark detection ────────────────────────────────────────────────────────────

#[derive(Clone, Debug)]
pub struct Quark {
    pub site: usize,
    pub r: usize,
    pub c: usize,
    pub z: usize,
    pub quark_type: QuarkType,
    #[allow(dead_code)]
    pub binding_coherence: f32,
    #[allow(dead_code)]
    pub veracity: f32,
    #[allow(dead_code)]
    pub curvature: f32,
}

#[derive(Clone, Debug, PartialEq)]
pub enum QuarkType {
    Up,
    Down,
}

pub fn detect_quarks(state: &GutoeState, cfg: &LatticeConfig) -> Vec<Quark> {
    let n = state.n(cfg);
    let mut quarks = Vec::new();
    for site in 0..n {
        let s = state.lattice[site];
        if s == VOID {
            continue;
        }
        let (r, c, z) = site_coords(site, cfg);
        let nbrs = mesh_neighbours(r, c, z, cfg);
        let (v, curv, grad) = state.local_fields(site, &nbrs);
        let bc = v / (1.0 + grad);
        if bc >= cfg.quark_threshold as f32 {
            quarks.push(Quark {
                site,
                r,
                c,
                z,
                quark_type: if v > curv {
                    QuarkType::Up
                } else {
                    QuarkType::Down
                },
                binding_coherence: bc,
                veracity: v,
                curvature: curv,
            });
        }
    }
    quarks
}

/// Returns list of (DOWN_site, UP1_site, UP2_site) proton triplets.
pub fn find_proton_triplets(quarks: &[Quark], cfg: &LatticeConfig) -> Vec<(usize, usize, usize)> {
    let quark_map: std::collections::HashMap<usize, &Quark> =
        quarks.iter().map(|q| (q.site, q)).collect();

    let nbr_cache: std::collections::HashMap<usize, std::collections::HashSet<usize>> = quarks
        .iter()
        .map(|q| {
            let nbrs = mesh_neighbours(q.r, q.c, q.z, cfg);
            (q.site, nbrs.iter().cloned().collect())
        })
        .collect();

    let mut used = std::collections::HashSet::new();
    let mut triplets = Vec::new();

    for q in quarks {
        if q.quark_type != QuarkType::Down || used.contains(&q.site) {
            continue;
        }
        let up_nbrs: Vec<&Quark> = nbr_cache[&q.site]
            .iter()
            .filter_map(|ni| quark_map.get(ni))
            .filter(|nq| nq.quark_type == QuarkType::Up && !used.contains(&nq.site))
            .cloned()
            .collect();
        if up_nbrs.len() < 2 {
            continue;
        }
        'outer: for i in 0..up_nbrs.len() {
            for j in (i + 1)..up_nbrs.len() {
                let p1 = up_nbrs[i].site;
                let p2 = up_nbrs[j].site;
                if nbr_cache.get(&p1).map_or(false, |s| s.contains(&p2)) {
                    triplets.push((q.site, p1, p2));
                    used.insert(q.site);
                    used.insert(p1);
                    used.insert(p2);
                    break 'outer;
                }
            }
        }
    }
    triplets
}

/// Inject n_inject leptons into sites within proton-containing layers (avoiding proton sites).
pub fn inject_leptons<R: Rng>(
    state: &mut GutoeState,
    triplets: &[(usize, usize, usize)],
    n_inject: usize,
    rng: &mut R,
    cfg: &LatticeConfig,
) {
    let layer_stride = cfg.hex_rows * cfg.hex_cols;
    let proton_sites: std::collections::HashSet<usize> = triplets
        .iter()
        .flat_map(|&(d, u1, u2)| [d, u1, u2])
        .collect();
    let proton_layers: std::collections::HashSet<usize> =
        triplets.iter().map(|&(d, _, _)| d / layer_stride).collect();

    let n = state.n(cfg);
    let mut candidates: Vec<usize> = (0..n)
        .filter(|s| !proton_sites.contains(s) && proton_layers.contains(&(s / layer_stride)))
        .filter(|&s| state.lattice[s] != LEPTON_SEED)
        .collect();

    if candidates.is_empty() {
        candidates = (0..n)
            .filter(|s| !proton_sites.contains(s) && state.lattice[*s] != LEPTON_SEED)
            .collect();
    }

    let count = n_inject.min(candidates.len());
    // Fisher-Yates partial shuffle to pick `count` random candidates
    for i in 0..count {
        let j = rng.gen_range(i..candidates.len());
        candidates.swap(i, j);
    }
    for &s in &candidates[..count] {
        state.lattice[s] = LEPTON_SEED;
    }
}

/// Proton shell: sites adjacent to any proton quark site (excluding proton sites themselves)
pub fn proton_shell(
    triplets: &[(usize, usize, usize)],
    cfg: &LatticeConfig,
) -> std::collections::HashSet<usize> {
    let proton_sites: std::collections::HashSet<usize> = triplets
        .iter()
        .flat_map(|&(d, u1, u2)| [d, u1, u2])
        .collect();
    let mut shell = std::collections::HashSet::new();
    for &ps in &proton_sites {
        let (r, c, z) = site_coords(ps, cfg);
        for nb in mesh_neighbours(r, c, z, cfg) {
            if !proton_sites.contains(&nb) {
                shell.insert(nb);
            }
        }
    }
    shell
}
