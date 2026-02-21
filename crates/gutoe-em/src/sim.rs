// GUTOE EM — Cl(1,3) lattice dynamics (two-pass step)
// Copyright (C) 2026 Riff Labs, AGPL-3.0-or-later
//
// Ports gutoe_em_hydrogen.py exactly:
//   Pass 1 — VOID differentiation + quark Z3/Clifford/alignment (skip leptons)
//   Pass 2 — lepton EM hops (grade-1 accessible, proton sites excluded)

use std::collections::{HashMap, HashSet};

use rand::Rng;

use crate::config::{LatticeConfig, LEPTON_SEED, QUARK_SEED, VOID};
use crate::gauge::GaugeFields;
use crate::geometry::{mesh_neighbours, site_coords};

// ── Cl(1,3) algebra tables ────────────────────────────────────────────────────

const SQRT3_HALF: f64 = 0.866_025_403_784_438_6; // √3/2

/// Veracity v(s1, s2) — alignment measure in Cl(1,3).
/// s ∈ 0..=16 where 0=VOID, states 1..=16 map to multivectors mi=(s−1)∈0..15.
///   VOID pair  → 0
///   equal      → 1
///   hamming-1  → √3/2
///   hamming-2  → 1/2
///   otherwise  → 0
#[inline]
pub fn veracity(s1: u8, s2: u8) -> f64 {
    if s1 == VOID || s2 == VOID {
        return 0.0;
    }
    if s1 == s2 {
        return 1.0;
    }
    let d = ((s1 - 1) ^ (s2 - 1)).count_ones();
    match d {
        1 => SQRT3_HALF,
        2 => 0.5,
        _ => 0.0,
    }
}

/// Z₃ cycle on Cl(1,3) states: bit rotation b₀b₁b₂b₃ → b₃b₀b₁b₂ on mi=(s−1).
const Z3_TABLE: [u8; 17] = {
    let mut t = [0u8; 17]; // t[0] = VOID stays VOID
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

/// Grade of a Clifford state (number of set bits in multivector index).
/// Returns -1 for VOID, 0..=4 for states 1..=16.
#[allow(dead_code)]
#[inline]
pub fn grade_of(s: u8) -> i32 {
    if s == VOID {
        -1
    } else {
        (s - 1).count_ones() as i32
    }
}

/// Z₃ quark orbits in Cl(1,3) — three-element groups under Z₃ rotation.
pub const Z3_ORBITS: [[u8; 3]; 4] = [[3, 5, 9], [4, 6, 10], [7, 11, 13], [8, 12, 14]];

// ── One-loop running coupling ─────────────────────────────────────────────────

/// Running Z₃ color coupling via the one-loop beta function.
///
/// α_s(t) = α_UV / (1 − (b₀/2π) × α_UV × ln(t+1))
///
/// Physics:
/// - UV (t=0): α_s = α_UV (quarks nearly free, cycle easily)
/// - IR (t → t_*): α_s → ∞ (Landau pole = confinement transition)
/// - t_* = exp(2π / (b₀ × α_UV)) − 1  (with default params: t_* ≈ 149)
///
/// b₀ comes from the Clifford grade structure:
///   b₀ = (11/3) × N_grade2 − (2/3) × N_grade1 = 58/3 ≈ 19.33
pub fn running_alpha_s(t: usize, cfg: &LatticeConfig) -> f64 {
    let a = cfg.coupling_uv;
    let b0_2pi = cfg.beta_coeff / (2.0 * std::f64::consts::PI);
    let denom = 1.0 - b0_2pi * a * ((t + 1) as f64).ln();
    if denom <= 0.0 {
        f64::INFINITY
    } else {
        a / denom
    }
}

/// Effective Z₃ cycle probability at timestep t.
///
/// As α_s grows (IR), the cycle rate DECREASES: quarks freeze into
/// color-singlet configurations (confinement). cycle_prob → 0 at the
/// Landau pole.
pub fn cycle_prob_rg(t: usize, cfg: &LatticeConfig) -> f64 {
    let alpha_s = running_alpha_s(t, cfg);
    if alpha_s.is_infinite() {
        0.0 // Fully confined: quarks cannot change color
    } else {
        // cycle_prob scales INVERSELY with α_s: more coupling = less cycling
        (cfg.cycle_prob * cfg.coupling_uv / alpha_s).min(1.0)
    }
}

/// Effective confinement energy scale (alignment strength) at timestep t.
///
/// The color binding energy INCREASES with α_s: stronger coupling →
/// tighter confinement → higher mass. This is the mechanism for the
/// proton-to-lepton mass ratio growing toward 1836.
pub fn alignment_rg(t: usize, cfg: &LatticeConfig) -> f64 {
    let alpha_s = running_alpha_s(t, cfg);
    if alpha_s.is_infinite() {
        cfg.alignment_strength * 1e4 // Past Landau pole: maximum confinement
    } else {
        cfg.alignment_strength * alpha_s / cfg.coupling_uv
    }
}

/// The Landau pole timestep: confinement transition.
/// t_* = exp(2π / (beta_coeff × coupling_uv)) − 1
pub fn landau_pole(cfg: &LatticeConfig) -> f64 {
    let b0_2pi = cfg.beta_coeff / (2.0 * std::f64::consts::PI);
    (1.0 / (b0_2pi * cfg.coupling_uv)).exp() - 1.0
}

// ── Lattice initialisation ────────────────────────────────────────────────────

pub fn init_lattice(cfg: &LatticeConfig) -> Vec<u8> {
    vec![VOID; cfg.n_sites()]
}

// ── Sample without replacement (Fisher-Yates) ─────────────────────────────────

pub fn sample_without_replacement<R: Rng>(rng: &mut R, pool: &[usize], n: usize) -> Vec<usize> {
    let mut pool = pool.to_vec();
    let take = n.min(pool.len());
    for i in 0..take {
        let j = rng.gen_range(i..pool.len());
        pool.swap(i, j);
    }
    pool[..take].to_vec()
}

// ── Single simulation step ────────────────────────────────────────────────────

/// One simulation step with two-pass design and RG-running Z₃ coupling.
///
/// `t` — current timestep, used to compute the running color coupling.
/// `gauge = None` → Phase 1 (quarks only).
/// `proton_sites` → set of proton quark sites excluded from lepton hop candidates.
///
/// RG dynamics:
///   cycle_prob_rg(t) → 0 as t → t_*   (quarks freeze: confinement)
///   alignment_rg(t)  → ∞ as t → t_*   (binding energy grows: mass)
///   t_* = exp(2π / (b₀ × α_UV)) ≈ 149 (end of Phase 1)
pub fn step<R: Rng>(
    lattice: &[u8],
    rng: &mut R,
    cfg: &LatticeConfig,
    gauge: Option<&GaugeFields>,
    proton_sites: &HashSet<usize>,
    t: usize,
) -> Vec<u8> {
    let n = cfg.n_sites();
    let mut new = lattice.to_vec();

    // Running coupling: cycle_prob decreases, alignment increases toward t_*
    let cp = cycle_prob_rg(t, cfg);
    let al = alignment_rg(t, cfg).min(1.0 - cp - cfg.clifford_prob);
    // al is capped so probabilities stay valid: cp + clifford + al ≤ 1

    // ── Pass 1: VOID differentiation + quark dynamics ─────────────────────────
    for site in 0..n {
        let (r, c, z) = site_coords(site, cfg);
        let state = lattice[site];

        if state == VOID {
            // Spontaneous differentiation
            if rng.gen::<f64>() < cfg.differentiation_prob {
                new[site] = QUARK_SEED;
                continue;
            }
            // Neighbourhood-driven activation (k=4 void votes)
            let nbrs = mesh_neighbours(r, c, z, cfg);
            let active = nbrs.iter().filter(|&&ni| lattice[ni] != VOID).count();
            let total = nbrs.len();
            if active >= 2.max(total / 4)
                && rng.gen::<f64>() < (active as f64 / total as f64) * 0.4
            {
                new[site] = QUARK_SEED;
            }
        } else if state == LEPTON_SEED && gauge.is_some() {
            // Skip leptons in Pass 1 — they are handled in Pass 2
        } else {
            // Quarks and all other non-void, non-lepton states
            let r_val: f64 = rng.gen();

            if r_val < cp {
                // Z₃ cycle: bit rotation in Cl(1,3)
                // Rate DECREASES with time: quarks freeze at confinement
                new[site] = Z3_TABLE[state as usize];
            } else if r_val < cp + cfg.clifford_prob {
                // Clifford XOR with a random active (non-lepton) neighbour
                let nbrs = mesh_neighbours(r, c, z, cfg);
                let partners: Vec<u8> = nbrs
                    .iter()
                    .filter_map(|&ni| {
                        let ns = lattice[ni];
                        if ns != VOID && ns != LEPTON_SEED {
                            Some(ns)
                        } else {
                            None
                        }
                    })
                    .collect();
                if !partners.is_empty() {
                    let partner = partners[rng.gen_range(0..partners.len())];
                    new[site] = ((state - 1) ^ (partner - 1)) + 1;
                }
            } else if r_val < cp + cfg.clifford_prob + al {
                // Alignment: majority vote among non-void, non-lepton neighbours
                // Rate INCREASES with time: stronger confinement → more alignment
                let nbrs = mesh_neighbours(r, c, z, cfg);
                let nbr_states: Vec<u8> = nbrs
                    .iter()
                    .filter_map(|&ni| {
                        let ns = lattice[ni];
                        if ns != VOID && ns != LEPTON_SEED {
                            Some(ns)
                        } else {
                            None
                        }
                    })
                    .collect();
                if !nbr_states.is_empty() {
                    // Count votes for each state
                    let mut votes = HashMap::new();
                    for &ns in &nbr_states {
                        *votes.entry(ns).or_insert(0usize) += 1;
                    }
                    let (&winner, &cnt) = votes.iter().max_by_key(|(_, &v)| v).unwrap();
                    if cnt > cfg.void_votes {
                        new[site] = winner;
                    }
                }
            }
        }
    }

    // ── Pass 2: lepton EM hops ────────────────────────────────────────────────
    // Reads lattice[] but writes to new[] to prevent race conditions.
    // Accessible = any non-lepton site that is not a proton quark.
    // Grade-1 accessible (not restricted to VOID/grade-2): the lattice fully
    // saturates to grade-1 by t≈100 so VOID restriction leaves leptons frozen.
    // Proton sites excluded: φ peaks at proton quarks; without exclusion the
    // lepton hops INTO the proton (destroying it) instead of orbiting the shell.
    if let Some(g) = gauge {
        let phi = &g.phi;
        for site in 0..n {
            if lattice[site] != LEPTON_SEED {
                continue;
            }
            if rng.gen::<f64>() < cfg.em_prob {
                let (r, c, z) = site_coords(site, cfg);
                let nbrs = mesh_neighbours(r, c, z, cfg);
                let candidates: Vec<(f64, usize)> = nbrs
                    .iter()
                    .filter_map(|&nb| {
                        if new[nb] != LEPTON_SEED && !proton_sites.contains(&nb) {
                            Some((phi[nb], nb))
                        } else {
                            None
                        }
                    })
                    .collect();
                if let Some(&(_, target)) =
                    candidates.iter().max_by(|a, b| a.0.partial_cmp(&b.0).unwrap())
                {
                    let displaced = new[target];
                    new[site] = displaced;
                    new[target] = LEPTON_SEED;
                }
            }
        }
    }

    new
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::{LatticeConfig, VOID, LEPTON_SEED};
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    fn small_cfg() -> LatticeConfig {
        LatticeConfig {
            hex_rows: 8,
            hex_cols: 8,
            layers: 1,
            ..Default::default()
        }
    }

    #[test]
    fn z3_table_is_order_3_on_quark_seed() {
        // QUARK_SEED=3 → orbit {3,5,9}
        let s = QUARK_SEED;
        let s2 = Z3_TABLE[s as usize];
        let s3 = Z3_TABLE[s2 as usize];
        let s4 = Z3_TABLE[s3 as usize];
        assert_eq!(s4, s, "Z3 is not order-3 on QUARK_SEED");
    }

    #[test]
    fn z3_table_void_fixed_point() {
        assert_eq!(Z3_TABLE[VOID as usize], VOID);
    }

    #[test]
    fn veracity_self_is_one() {
        for s in 1u8..=16 {
            assert_eq!(veracity(s, s), 1.0, "veracity({s},{s}) ≠ 1");
        }
    }

    #[test]
    fn veracity_void_is_zero() {
        for s in 0u8..=16 {
            assert_eq!(veracity(VOID, s), 0.0);
            assert_eq!(veracity(s, VOID), 0.0);
        }
    }

    #[test]
    fn grade_of_lepton_seed() {
        // LEPTON_SEED=2, mi=1=0b0001, grade=1
        assert_eq!(grade_of(LEPTON_SEED), 1);
    }

    #[test]
    fn step_void_lattice_grows() {
        let cfg = small_cfg();
        let lattice = init_lattice(&cfg);
        let mut rng = StdRng::seed_from_u64(42);
        let next = step(&lattice, &mut rng, &cfg, None, &HashSet::new(), 0);
        let active = next.iter().filter(|&&s| s != VOID).count();
        // With differentiation_prob=0.02 and N=64, expect some active sites
        assert!(active > 0, "Lattice should grow from VOID after one step");
    }

    #[test]
    fn lepton_conservation_within_grade1_orbit() {
        // γ⁰ (LEPTON_SEED=2, mi=0b0001) is NOT producible by Clifford XOR
        // within the grade-1 quark orbit {3,5,9} (γ¹,γ²,γ³; mi∈{2,4,8}).
        //
        // Proof: bits 1,2,3 (mi values 2,4,8) XOR to 6, 10, or 12 —
        // none of which have bit 0 set (mi=1 = LEPTON_SEED=2).
        //
        // Physical meaning: once the lattice saturates to the grade-1 orbit
        // by t≈100 (k=4 fully activates all sites), no Clifford quark–quark
        // interaction within that orbit can spontaneously create a lepton.
        // Lepton number is conserved by the Clifford algebra structure of {3,5,9}.
        let grade1_orbit = [3u8, 5, 9]; // the primary stable quark Z3 orbit
        for &s1 in &grade1_orbit {
            for &s2 in &grade1_orbit {
                if s1 == s2 {
                    continue; // identity XOR → 0, not a valid state
                }
                let xor_result = ((s1 - 1) ^ (s2 - 1)) + 1;
                assert_ne!(
                    xor_result, LEPTON_SEED,
                    "Clifford XOR within grade-1 orbit: {s1}^{s2} = {xor_result} = LEPTON_SEED — \
                     lepton number violated in primary quark orbit!"
                );
            }
        }
    }
}
