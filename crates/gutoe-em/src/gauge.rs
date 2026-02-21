// GUTOE EM — U(1) gauge fields: Poisson solver, Maxwell wave, EM force
// Copyright (C) 2026 Riff Labs, AGPL-3.0-or-later
//
// Ports gutoe_gauge.py exactly:
//   phi[N]    — Coulomb potential  (∇²φ = −ρ, Jacobi)
//   A[N]      — photon wave field  (∂²A/∂t² = c²∇²A + J, leapfrog)
//   A_prev[N] — previous A for leapfrog

use std::collections::HashMap;

use crate::config::{LatticeConfig, QuarkType, LEPTON_CHARGE, LEPTON_SEED, UP_CHARGE, DOWN_CHARGE};
use crate::geometry::{site_coords, mesh_neighbours};

// ── Gauge field container ─────────────────────────────────────────────────────

pub struct GaugeFields {
    pub phi: Vec<f64>,
    pub a: Vec<f64>,
    pub a_prev: Vec<f64>,
}

impl GaugeFields {
    pub fn new(n: usize) -> Self {
        Self {
            phi: vec![0.0; n],
            a: vec![0.0; n],
            a_prev: vec![0.0; n],
        }
    }
}

// ── Charge density ────────────────────────────────────────────────────────────

pub fn compute_charge_density(
    lattice: &[u8],
    quark_types: &HashMap<usize, QuarkType>,
    cfg: &LatticeConfig,
) -> Vec<f64> {
    let mut rho = vec![0.0f64; cfg.n_sites()];
    for (&site, qtype) in quark_types {
        rho[site] = match qtype {
            QuarkType::Up => UP_CHARGE,
            QuarkType::Down => DOWN_CHARGE,
        };
    }
    for (site, &state) in lattice.iter().enumerate() {
        if state == LEPTON_SEED {
            rho[site] = LEPTON_CHARGE;
        }
    }
    rho
}

// ── Poisson solver (Jacobi) ───────────────────────────────────────────────────

/// Cached neighbour lists — build once, use for all Jacobi iterations.
fn build_nbr_cache(cfg: &LatticeConfig) -> Vec<Vec<usize>> {
    let n = cfg.n_sites();
    (0..n)
        .map(|site| {
            let (r, c, z) = site_coords(site, cfg);
            mesh_neighbours(r, c, z, cfg)
        })
        .collect()
}

/// Solve ∇²φ = −ρ on the hex lattice via Jacobi iteration.
///
/// Update rule: φ_new[i] = (Σⱼ∈nbrs φ[j]  +  n · ρ[i]) / n
///
/// Positive source (proton quark, +2/3 or combined +1) → positive φ.
/// Lepton drifts toward max φ → EM attraction.
pub fn jacobi_poisson(rho: &[f64], cfg: &LatticeConfig, n_iter: usize) -> Vec<f64> {
    let n = cfg.n_sites();
    let nbr_cache = build_nbr_cache(cfg);
    let mut phi = vec![0.0f64; n];
    let mut phi_new = vec![0.0f64; n];

    for _ in 0..n_iter {
        for site in 0..n {
            let nbrs = &nbr_cache[site];
            let k = nbrs.len() as f64;
            let sum_nbr: f64 = nbrs.iter().map(|&nb| phi[nb]).sum();
            phi_new[site] = (sum_nbr + k * rho[site]) / k;
        }
        std::mem::swap(&mut phi, &mut phi_new);
    }
    phi
}

// ── Discrete Laplacian ────────────────────────────────────────────────────────

fn hex_laplacian(field: &[f64], cfg: &LatticeConfig) -> Vec<f64> {
    let n = cfg.n_sites();
    let mut lap = vec![0.0f64; n];
    for site in 0..n {
        let (r, c, z) = site_coords(site, cfg);
        let nbrs = mesh_neighbours(r, c, z, cfg);
        let k = nbrs.len() as f64;
        let sum: f64 = nbrs.iter().map(|&nb| field[nb] - field[site]).sum();
        lap[site] = sum / k;
    }
    lap
}

// ── Maxwell scalar wave equation (leapfrog) ───────────────────────────────────

/// Leapfrog step: A_new = 2·A − A_prev + c²·∇²A + coupling·ρ
pub fn maxwell_wave_step(gauge: &mut GaugeFields, rho: &[f64], cfg: &LatticeConfig) {
    let c2 = cfg.photon_c * cfg.photon_c;
    let lap = hex_laplacian(&gauge.a, cfg);
    let n = cfg.n_sites();
    let mut a_new = vec![0.0f64; n];
    for i in 0..n {
        a_new[i] = 2.0 * gauge.a[i] - gauge.a_prev[i] + c2 * lap[i] + rho[i] * cfg.photon_coupling;
    }
    gauge.a_prev.clone_from(&gauge.a);
    gauge.a = a_new;
}

// ── EM force on lepton ────────────────────────────────────────────────────────

/// Lepton (charge −1) feels F = −q∇φ = +∇φ.
/// Returns the flat index of the neighbour with maximum φ + 0.3·A.
pub fn em_force_on_lepton(phi: &[f64], a: &[f64], site: usize, cfg: &LatticeConfig) -> usize {
    let (r, c, z) = site_coords(site, cfg);
    let nbrs = mesh_neighbours(r, c, z, cfg);
    nbrs.into_iter()
        .max_by(|&ai, &bi| {
            let va = phi[ai] + 0.3 * a[ai];
            let vb = phi[bi] + 0.3 * a[bi];
            va.partial_cmp(&vb).unwrap()
        })
        .unwrap()
}

// ── Full gauge update ─────────────────────────────────────────────────────────

/// Complete gauge field update per timestep.
/// 1. Compute full ρ (quarks + leptons)
/// 2. Poisson solve with quark-only ρ (lepton excluded → sees proton's +1 field)
/// 3. Maxwell leapfrog with full ρ (correct radiation physics)
pub fn update_gauge(
    gauge: &mut GaugeFields,
    lattice: &[u8],
    quark_types: &HashMap<usize, QuarkType>,
    cfg: &LatticeConfig,
) {
    let rho_full = compute_charge_density(lattice, quark_types, cfg);
    let mut rho_quarks = rho_full.clone();
    for (site, &state) in lattice.iter().enumerate() {
        if state == LEPTON_SEED {
            rho_quarks[site] = 0.0;
        }
    }
    gauge.phi = jacobi_poisson(&rho_quarks, cfg, cfg.poisson_iters);
    maxwell_wave_step(gauge, &rho_full, cfg);
}

// ── Tests: port of test_gutoe_gauge.py (all 13 tests) ────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::LatticeConfig;
    use crate::geometry::{site_coords, mesh_neighbours};

    fn small_cfg() -> LatticeConfig {
        LatticeConfig {
            hex_rows: 8,
            hex_cols: 8,
            layers: 1,
            ..Default::default()
        }
    }

    fn center_site(cfg: &LatticeConfig) -> usize {
        cfg.n_sites() / 2
    }

    // ── Poisson / Coulomb ─────────────────────────────────────────────────────

    /// No charges → φ = 0 everywhere.
    #[test]
    fn poisson_zero_charge_gives_zero() {
        let cfg = small_cfg();
        let rho = vec![0.0f64; cfg.n_sites()];
        let phi = jacobi_poisson(&rho, &cfg, 20);
        for (i, &v) in phi.iter().enumerate() {
            assert!(v.abs() < 1e-10, "phi[{i}] = {v}, expected 0");
        }
    }

    /// Single +1 charge → φ > 0 at charge site and neighbours.
    #[test]
    fn poisson_single_positive_charge_positive_phi() {
        let cfg = small_cfg();
        let mut rho = vec![0.0f64; cfg.n_sites()];
        let c = center_site(&cfg);
        rho[c] = 1.0;
        let phi = jacobi_poisson(&rho, &cfg, 100);
        assert!(phi[c] > 0.0, "φ at charge site should be positive, got {}", phi[c]);
        let (r, col, z) = site_coords(c, &cfg);
        let nbrs = mesh_neighbours(r, col, z, &cfg);
        for &nb in &nbrs {
            assert!(
                phi[c] >= phi[nb],
                "φ should peak at charge site: phi[center]={} phi[nbr]={}",
                phi[c],
                phi[nb]
            );
        }
    }

    /// Single −1 charge → φ < 0 at charge site.
    #[test]
    fn poisson_negative_charge_negative_phi() {
        let cfg = small_cfg();
        let mut rho = vec![0.0f64; cfg.n_sites()];
        rho[center_site(&cfg)] = -1.0;
        let phi = jacobi_poisson(&rho, &cfg, 100);
        assert!(phi[center_site(&cfg)] < 0.0);
    }

    /// φ decays with distance: phi[center] > phi[1-hop] > phi[2-hop].
    #[test]
    fn poisson_phi_decays_with_distance() {
        let cfg = small_cfg();
        let mut rho = vec![0.0f64; cfg.n_sites()];
        let c = center_site(&cfg);
        rho[c] = 1.0;
        let phi = jacobi_poisson(&rho, &cfg, 200);
        let (r, col, z) = site_coords(c, &cfg);
        let hop1 = mesh_neighbours(r, col, z, &cfg)[0];
        let (r1, c1, z1) = site_coords(hop1, &cfg);
        let hop2 = mesh_neighbours(r1, c1, z1, &cfg)
            .into_iter()
            .find(|&nb| nb != c)
            .unwrap();
        assert!(
            phi[c] > phi[hop1] && phi[hop1] > phi[hop2],
            "φ should decay: {:.4} > {:.4} > {:.4}",
            phi[c],
            phi[hop1],
            phi[hop2]
        );
    }

    /// Poisson is linear: φ(2ρ) = 2φ(ρ).
    #[test]
    fn poisson_linearity() {
        let cfg = small_cfg();
        let mut rho = vec![0.0f64; cfg.n_sites()];
        rho[center_site(&cfg)] = 1.0;
        let phi1 = jacobi_poisson(&rho, &cfg, 100);
        let rho2: Vec<f64> = rho.iter().map(|&x| 2.0 * x).collect();
        let phi2 = jacobi_poisson(&rho2, &cfg, 100);
        for i in 0..cfg.n_sites() {
            let ratio = if phi1[i].abs() > 1e-12 {
                phi2[i] / phi1[i]
            } else {
                // Both should be ~0
                if phi2[i].abs() < 1e-10 { 2.0 } else { phi2[i] }
            };
            assert!(
                (ratio - 2.0).abs() < 0.02,
                "φ(2ρ)/φ(ρ) = {ratio:.4} at site {i}, expected 2.0 (Poisson must be linear)"
            );
        }
    }

    // ── Maxwell wave equation ─────────────────────────────────────────────────

    /// No source, no initial field → A stays zero.
    #[test]
    fn maxwell_zero_source_no_field() {
        let cfg = small_cfg();
        let mut gauge = GaugeFields::new(cfg.n_sites());
        let rho = vec![0.0f64; cfg.n_sites()];
        for _ in 0..10 {
            maxwell_wave_step(&mut gauge, &rho, &cfg);
        }
        for (i, &v) in gauge.a.iter().enumerate() {
            assert!(v.abs() < 1e-10, "A[{i}] = {v}, expected 0 with no source");
        }
    }

    /// A source at center drives A away from zero.
    #[test]
    fn maxwell_source_drives_field() {
        let cfg = small_cfg();
        let mut gauge = GaugeFields::new(cfg.n_sites());
        let mut rho = vec![0.0f64; cfg.n_sites()];
        rho[center_site(&cfg)] = 1.0;
        for _ in 0..5 {
            maxwell_wave_step(&mut gauge, &rho, &cfg);
        }
        assert_ne!(
            gauge.a[center_site(&cfg)],
            0.0,
            "Source should drive A field at center"
        );
    }

    /// Field driven at center should reach neighbours within ~20 steps.
    #[test]
    fn maxwell_wave_propagates() {
        let cfg = small_cfg();
        let mut gauge = GaugeFields::new(cfg.n_sites());
        let mut rho = vec![0.0f64; cfg.n_sites()];
        let c = center_site(&cfg);
        rho[c] = 1.0;
        let (r, col, z) = site_coords(c, &cfg);
        let hop1 = mesh_neighbours(r, col, z, &cfg)[0];
        for _ in 0..20 {
            maxwell_wave_step(&mut gauge, &rho, &cfg);
        }
        assert!(
            gauge.a[hop1].abs() > 0.0,
            "Wave should propagate to neighbours. A[hop1]={:.6}",
            gauge.a[hop1]
        );
    }

    // ── EM force on lepton ────────────────────────────────────────────────────

    /// γ⁰ at site X feels force toward the site with highest φ.
    #[test]
    fn lepton_attracted_to_positive_phi() {
        let cfg = small_cfg();
        let mut phi = vec![0.0f64; cfg.n_sites()];
        let a = vec![0.0f64; cfg.n_sites()];
        let c = center_site(&cfg);
        let (r, col, z) = site_coords(c, &cfg);
        let nbrs = mesh_neighbours(r, col, z, &cfg);
        let target = nbrs[2];
        phi[target] = 10.0;
        let force_target = em_force_on_lepton(&phi, &a, c, &cfg);
        assert_eq!(
            force_target, target,
            "Lepton should move toward max φ={} at {target}, got {force_target}",
            phi[target]
        );
    }

    /// γ⁰ flees very negative φ (another negative charge repelling it).
    #[test]
    fn lepton_moves_away_from_negative_phi() {
        let cfg = small_cfg();
        let mut phi = vec![0.0f64; cfg.n_sites()];
        let a = vec![0.0f64; cfg.n_sites()];
        let c = center_site(&cfg);
        let (r, col, z) = site_coords(c, &cfg);
        let nbrs = mesh_neighbours(r, col, z, &cfg);
        for &nb in &nbrs {
            phi[nb] = 1.0;
        }
        let repulsion_site = nbrs[0];
        phi[repulsion_site] = -10.0;
        let force_target = em_force_on_lepton(&phi, &a, c, &cfg);
        assert_ne!(
            force_target, repulsion_site,
            "Lepton should flee negative φ at {repulsion_site}"
        );
    }

    // ── Charge density ────────────────────────────────────────────────────────

    /// A proton (2 UP + 1 DOWN) has net charge exactly +1.
    #[test]
    fn charge_density_proton_net_positive() {
        let net = 2.0 * UP_CHARGE + 1.0 * DOWN_CHARGE;
        assert!(
            (net - 1.0).abs() < 1e-9,
            "Proton charge = {net}, expected 1.0"
        );
    }

    /// γ⁰ lepton has charge −1.
    #[test]
    fn charge_density_lepton_negative() {
        assert_eq!(LEPTON_CHARGE, -1.0);
    }

    /// compute_charge_density returns array of correct shape and values.
    #[test]
    fn compute_charge_density_shape_and_values() {
        let cfg = small_cfg();
        let lattice = vec![0u8; cfg.n_sites()];
        let mut quark_types = HashMap::new();
        quark_types.insert(0usize, QuarkType::Up);
        quark_types.insert(1usize, QuarkType::Up);
        quark_types.insert(2usize, QuarkType::Down);
        let rho = compute_charge_density(&lattice, &quark_types, &cfg);
        assert_eq!(rho.len(), cfg.n_sites());
        assert!(rho[0] > 0.0, "UP quark site should have positive charge");
        assert!(rho[2] < 0.0, "DOWN quark site should have negative charge");
    }
}
