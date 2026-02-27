//! GUTOE U(1) gauge fields — Rust port of gutoe_gauge.py
//!
//! phi[N]    — scalar Coulomb potential (∇²φ = −ρ)
//! A[N]      — scalar photon field       (∂²A/∂t² = c²∇²A + J)
//! A_prev[N] — previous step A for leapfrog

use crate::sim::{mesh_neighbours, site_coords, LatticeConfig, QuarkType};

// ── Charge constants ───────────────────────────────────────────────────────────

pub const UP_CHARGE: f64 = 2.0 / 3.0;
pub const DOWN_CHARGE: f64 = -1.0 / 3.0;
pub const LEPTON_CHARGE: f64 = -1.0;

const LEPTON_SEED: u8 = 2;

// ── Gauge fields ───────────────────────────────────────────────────────────────

pub struct GaugeFields {
    pub phi: Vec<f64>,
    pub a: Vec<f64>,
    pub a_prev: Vec<f64>,
}

impl GaugeFields {
    pub fn new(cfg: &LatticeConfig) -> Self {
        let n = cfg.hex_rows * cfg.hex_cols * cfg.layers;
        Self {
            phi: vec![0.0; n],
            a: vec![0.0; n],
            a_prev: vec![0.0; n],
        }
    }
}

// ── Neighbour cache ────────────────────────────────────────────────────────────

pub struct NbrCache(pub Vec<[usize; 6]>);

impl NbrCache {
    pub fn build(cfg: &LatticeConfig) -> Self {
        let n = cfg.hex_rows * cfg.hex_cols * cfg.layers;
        let mut cache = Vec::with_capacity(n);
        for site in 0..n {
            let (r, c, z) = site_coords(site, cfg);
            cache.push(mesh_neighbours(r, c, z, cfg));
        }
        NbrCache(cache)
    }
}

// ── Charge density ─────────────────────────────────────────────────────────────

pub fn compute_charge_density(
    lattice: &[u8],
    quark_map: &std::collections::HashMap<usize, QuarkType>,
    cfg: &LatticeConfig,
) -> Vec<f64> {
    let n = cfg.hex_rows * cfg.hex_cols * cfg.layers;
    let mut rho = vec![0.0f64; n];
    for (&site, qtype) in quark_map {
        rho[site] = match qtype {
            QuarkType::Up => UP_CHARGE,
            QuarkType::Down => DOWN_CHARGE,
        };
    }
    for site in 0..n {
        if lattice[site] == LEPTON_SEED {
            rho[site] = LEPTON_CHARGE;
        }
    }
    rho
}

// ── Discrete Laplacian ─────────────────────────────────────────────────────────

fn hex_laplacian(field: &[f64], nbr_cache: &NbrCache) -> Vec<f64> {
    let n = field.len();
    let mut lap = vec![0.0f64; n];
    for site in 0..n {
        let nbrs = &nbr_cache.0[site];
        let sum: f64 = nbrs.iter().map(|&nb| field[nb] - field[site]).sum();
        lap[site] = sum / nbrs.len() as f64;
    }
    lap
}

// ── Poisson solver (Jacobi) ───────────────────────────────────────────────────

/// Solve ∇²φ = −ρ on hex lattice via Jacobi iteration.
/// φ[i] = (Σⱼ φ[j] + n·ρ[i]) / n
pub fn jacobi_poisson(
    rho: &[f64],
    _cfg: &LatticeConfig,
    nbr_cache: &NbrCache,
    n_iter: usize,
) -> Vec<f64> {
    let n = rho.len();
    let mut phi = vec![0.0f64; n];
    let mut phi_new = vec![0.0f64; n];
    for _ in 0..n_iter {
        for site in 0..n {
            let nbrs = &nbr_cache.0[site];
            let k = nbrs.len() as f64;
            let sum: f64 = nbrs.iter().map(|&nb| phi[nb]).sum();
            phi_new[site] = (sum + k * rho[site]) / k;
        }
        std::mem::swap(&mut phi, &mut phi_new);
    }
    phi
}

// ── Maxwell wave equation (leapfrog) ─────────────────────────────────────────

/// A_new = 2·A − A_prev + c²·∇²A + coupling·ρ
pub fn maxwell_wave_step(
    gauge: &mut GaugeFields,
    rho: &[f64],
    cfg: &LatticeConfig,
    nbr_cache: &NbrCache,
) {
    let c2 = cfg.photon_c * cfg.photon_c;
    let lap = hex_laplacian(&gauge.a, nbr_cache);
    let n = gauge.a.len();
    let mut a_new = vec![0.0f64; n];
    for site in 0..n {
        a_new[site] = 2.0 * gauge.a[site] - gauge.a_prev[site]
            + c2 * lap[site]
            + rho[site] * cfg.photon_coupling;
    }
    std::mem::swap(&mut gauge.a_prev, &mut gauge.a);
    gauge.a = a_new;
}

// ── Full gauge update ─────────────────────────────────────────────────────────

/// Full gauge update per timestep.
/// phi is set from proton-only quark ρ (lepton self-energy excluded).
/// Maxwell wave uses full ρ.
pub fn update_gauge(
    gauge: &mut GaugeFields,
    lattice: &[u8],
    quark_map: &std::collections::HashMap<usize, QuarkType>,
    proton_sites: &std::collections::HashSet<usize>,
    cfg: &LatticeConfig,
    nbr_cache: &NbrCache,
) {
    let n = cfg.hex_rows * cfg.hex_cols * cfg.layers;

    // Full ρ for Maxwell
    let rho_full = compute_charge_density(lattice, quark_map, cfg);

    // Proton-only ρ for Coulomb φ
    let mut rho_proton = vec![0.0f64; n];
    for &site in proton_sites {
        if let Some(qtype) = quark_map.get(&site) {
            rho_proton[site] = match qtype {
                QuarkType::Up => UP_CHARGE,
                QuarkType::Down => DOWN_CHARGE,
            };
        }
    }

    gauge.phi = jacobi_poisson(&rho_proton, cfg, nbr_cache, cfg.poisson_iters);
    maxwell_wave_step(gauge, &rho_full, cfg, nbr_cache);
}
