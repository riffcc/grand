// GUTOE EM — Particle detection, proton triplets, enrichment analysis
// Copyright (C) 2026 Riff Labs, AGPL-3.0-or-later

use std::collections::{HashMap, HashSet};

use crate::config::{LatticeConfig, QuarkType, LEPTON_SEED, VOID};
use crate::gauge::GaugeFields;
use crate::geometry::{mesh_neighbours, site_coords};
use crate::sim::veracity;

// ── Quark detection ───────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct Quark {
    pub site: usize,
    pub r: usize,
    pub c: usize,
    pub z: usize,
    pub quark_type: QuarkType,
}

/// Detect quark sites: sites where binding coherence bc = v/(1+grad) ≥ threshold.
/// Leptons and VOID are excluded. Grade-2+ states participate (stable under k=4).
pub fn detect_quarks(lattice: &[u8], cfg: &LatticeConfig) -> Vec<Quark> {
    let n = cfg.n_sites();
    let mut quarks = Vec::new();

    for site in 0..n {
        let state = lattice[site];
        if state == VOID || state == LEPTON_SEED {
            continue;
        }

        let (r, c, z) = site_coords(site, cfg);
        let nbrs = mesh_neighbours(r, c, z, cfg);
        let mut total_v = 0.0f64;
        let mut grad = 0.0f64;
        let mut nbr_set = HashSet::new();

        for &ni in &nbrs {
            let ns = lattice[ni];
            let v = veracity(state, ns);
            total_v += v;
            grad += 1.0 - v;
            if ns != VOID {
                nbr_set.insert(ns);
            }
        }

        let n_nbrs = nbrs.len() as f64;
        let v_mean = total_v / n_nbrs;

        // Z₃ curvature: max over 4 orbits
        const Z3_ORBITS: [[u8; 3]; 4] = [[3, 5, 9], [4, 6, 10], [7, 11, 13], [8, 12, 14]];
        let z3_curv = Z3_ORBITS
            .iter()
            .map(|orbit| {
                let cnt = orbit.iter().filter(|&&s| nbr_set.contains(&s)).count();
                (cnt as f64 - 1.0) / 2.0
            })
            .fold(f64::NEG_INFINITY, f64::max);

        let bc = v_mean / (1.0 + grad / n_nbrs);

        if bc >= cfg.quark_threshold {
            let qtype = if v_mean > z3_curv {
                QuarkType::Up
            } else {
                QuarkType::Down
            };
            quarks.push(Quark { site, r, c, z, quark_type: qtype });
        }
    }

    quarks
}

/// Find proton triplets: (DOWN, UP, UP) triangles in the hex lattice.
/// Returns Vec of [down_site, up1_site, up2_site].
pub fn find_proton_triplets(quarks: &[Quark], cfg: &LatticeConfig) -> Vec<[usize; 3]> {
    let quark_map: HashMap<usize, &Quark> = quarks.iter().map(|q| (q.site, q)).collect();
    let mut triplets = Vec::new();
    let mut used = HashSet::new();

    // Pre-cache neighbour sets for each quark site
    let nbr_cache: HashMap<usize, HashSet<usize>> = quarks
        .iter()
        .map(|q| {
            let nbrs: HashSet<usize> =
                mesh_neighbours(q.r, q.c, q.z, cfg).into_iter().collect();
            (q.site, nbrs)
        })
        .collect();

    for q in quarks {
        if q.quark_type != QuarkType::Down || used.contains(&q.site) {
            continue;
        }

        // Find unused UP neighbours of this DOWN quark
        let up_nbrs: Vec<&Quark> = nbr_cache
            .get(&q.site)
            .map(|s| {
                s.iter()
                    .filter_map(|&ni| {
                        quark_map
                            .get(&ni)
                            .filter(|qn| qn.quark_type == QuarkType::Up && !used.contains(&qn.site))
                            .copied()
                    })
                    .collect()
            })
            .unwrap_or_default();

        if up_nbrs.len() < 2 {
            continue;
        }

        // Find two UP neighbours that are also neighbours of each other (triangle)
        'outer: for i in 0..up_nbrs.len() {
            for j in (i + 1)..up_nbrs.len() {
                let p1 = up_nbrs[i].site;
                let p2 = up_nbrs[j].site;
                if nbr_cache.get(&p2).map_or(false, |s| s.contains(&p1)) {
                    triplets.push([q.site, p1, p2]);
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

// ── Enrichment analysis ───────────────────────────────────────────────────────

/// Result of one analysis snapshot.
pub struct AnalysisResult {
    pub protons: usize,
    pub leptons: usize,
    pub hydrogen: usize,
    /// Layer-restricted γ⁰ enrichment (lepton density in shell / background density).
    /// Capped at 20× to keep averages finite (rb=0 means all leptons in shell).
    pub enrich: f64,
    /// Δφ: mean φ at lepton sites minus mean φ at non-lepton sites.
    pub phi_ratio: f64,
}

/// Analyse one lattice snapshot for proton count, hydrogen count, and enrichment.
///
/// Layer-restricted metric: EM is intra-layer only (mesh_neighbours never crosses
/// layers), so comparing shell vs. background across all 12 layers dilutes the
/// signal with 7-8 empty layers.  We restrict to layers that contain ≥1 proton.
pub fn analyze(lattice: &[u8], gauge: Option<&GaugeFields>, cfg: &LatticeConfig) -> AnalysisResult {
    let n = cfg.n_sites();
    let layer_stride = cfg.layer_stride();

    let quarks = detect_quarks(lattice, cfg);
    let trips = find_proton_triplets(&quarks, cfg);

    let p_sites: HashSet<usize> = trips.iter().flat_map(|&[d, u1, u2]| [d, u1, u2]).collect();

    // Shell: non-proton sites adjacent to any proton quark
    let mut p_shell = HashSet::new();
    for &s in &p_sites {
        let (r, c, z) = site_coords(s, cfg);
        for nb in mesh_neighbours(r, c, z, cfg) {
            if !p_sites.contains(&nb) {
                p_shell.insert(nb);
            }
        }
    }

    let n_lep = lattice.iter().filter(|&&s| s == LEPTON_SEED).count();

    // Hydrogen: at least one adjacent γ⁰ in the proton shell
    let n_h = trips
        .iter()
        .filter(|&[d, u1, u2]| {
            let mut shell = HashSet::new();
            for &s in &[*d, *u1, *u2] {
                let (r, c, z) = site_coords(s, cfg);
                for nb in mesh_neighbours(r, c, z, cfg) {
                    if !p_sites.contains(&nb) {
                        shell.insert(nb);
                    }
                }
            }
            shell.iter().any(|&nb| lattice[nb] == LEPTON_SEED)
        })
        .count();

    // Layer-restricted enrichment
    let proton_layers: HashSet<usize> = trips
        .iter()
        .map(|&[d, _, _]| d / layer_stride)
        .collect();

    let lep_shell = p_shell.iter().filter(|&&s| lattice[s] == LEPTON_SEED).count();
    let shell_sz = p_shell.len().max(1);

    let bg_sites: Vec<usize> = (0..n)
        .filter(|&s| {
            proton_layers.contains(&(s / layer_stride))
                && !p_sites.contains(&s)
                && !p_shell.contains(&s)
        })
        .collect();
    let lep_bg = bg_sites.iter().filter(|&&s| lattice[s] == LEPTON_SEED).count();
    let bg_sz = bg_sites.len().max(1);

    let rs = lep_shell as f64 / shell_sz as f64;
    let rb = lep_bg as f64 / bg_sz as f64;
    let enrich = if rb > 1e-9 {
        (rs / rb).min(20.0)
    } else if rs > 0.0 {
        20.0 // all accessible leptons are in the shell
    } else {
        0.0
    };

    // φ-tracking: do leptons sit in higher-φ regions than background?
    let phi_ratio = if let (Some(g), true) = (gauge, n_lep > 0) {
        let phi = &g.phi;
        let lep_sum: f64 = (0..n)
            .filter(|&s| lattice[s] == LEPTON_SEED)
            .map(|s| phi[s])
            .sum();
        let bg_sum: f64 = (0..n)
            .filter(|&s| lattice[s] != LEPTON_SEED)
            .map(|s| phi[s])
            .sum();
        let n_bg = n - n_lep;
        let phi_lep = lep_sum / n_lep as f64;
        let phi_bg = if n_bg > 0 { bg_sum / n_bg as f64 } else { phi_lep };
        phi_lep - phi_bg
    } else {
        0.0
    };

    AnalysisResult { protons: trips.len(), leptons: n_lep, hydrogen: n_h, enrich, phi_ratio }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::{LatticeConfig, LEPTON_SEED, VOID};

    fn small_cfg() -> LatticeConfig {
        LatticeConfig {
            hex_rows: 8,
            hex_cols: 8,
            layers: 1,
            ..Default::default()
        }
    }

    #[test]
    fn detect_quarks_all_void() {
        let cfg = small_cfg();
        let lattice = vec![VOID; cfg.n_sites()];
        let quarks = detect_quarks(&lattice, &cfg);
        assert!(quarks.is_empty(), "All-VOID lattice should have no quarks");
    }

    #[test]
    fn detect_quarks_excludes_leptons() {
        let cfg = small_cfg();
        let mut lattice = vec![VOID; cfg.n_sites()];
        lattice[0] = LEPTON_SEED;
        let quarks = detect_quarks(&lattice, &cfg);
        assert!(
            quarks.iter().all(|q| q.site != 0),
            "Lepton at site 0 should not be detected as quark"
        );
    }

    #[test]
    fn find_proton_triplets_empty_quarks() {
        let cfg = small_cfg();
        let trips = find_proton_triplets(&[], &cfg);
        assert!(trips.is_empty());
    }

    #[test]
    fn analyze_all_void_gives_zeros() {
        let cfg = small_cfg();
        let lattice = vec![VOID; cfg.n_sites()];
        let result = analyze(&lattice, None, &cfg);
        assert_eq!(result.protons, 0);
        assert_eq!(result.leptons, 0);
        assert_eq!(result.hydrogen, 0);
        assert_eq!(result.enrich, 0.0);
    }
}
