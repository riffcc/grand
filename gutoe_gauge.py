#!/usr/bin/env python3
"""
GUTOE: U(1) Gauge Field Physics Core

Two-layer architecture:
  Matter layer : discrete Clifford states on lattice sites
  Field layer  : continuous U(1) gauge fields

Gauge fields:
  phi[N]    — scalar Coulomb potential  (Poisson: ∇²φ = −ρ)
  A[N]      — scalar photon field       (wave eq: ∂²A/∂t² = c²∇²A + J)
  A_prev[N] — previous step A for leapfrog

Charge assignments (units of elementary charge e):
  UP quark   : +2/3
  DOWN quark : −1/3   → net proton charge = 2(+2/3) + (−1/3) = +1  ✓
  γ⁰ lepton  : −1
  everything else: 0

Copyright (C) 2026 Riff Labs, AGPL-3.0-or-later
"""

import numpy as np
from dataclasses import dataclass

# ── Constants ─────────────────────────────────────────────────────────────────

UP_CHARGE     = +2/3
DOWN_CHARGE   = -1/3
LEPTON_CHARGE = -1.0

VOID        = 0
LEPTON_SEED = 2   # γ⁰

_GRADE_TABLE = [-1] + [bin(mi).count('1') for mi in range(16)]

# ── Config ────────────────────────────────────────────────────────────────────

@dataclass
class LatticeConfig:
    hex_rows: int   = 12
    hex_cols: int   = 12
    layers:   int   = 12
    differentiation_prob: float = 0.02
    cycle_prob:           float = 0.05
    clifford_prob:        float = 0.03
    alignment_strength:   float = 0.15
    quark_threshold:      float = 0.6
    void_votes:           int   = 4
    em_prob:              float = 0.5    # prob lepton steps via EM each tick
    photon_c:             float = 0.4    # wave speed (stable: c < 1/√6 ≈ 0.408)
    photon_coupling:      float = 0.05   # source coupling J = coupling × ρ
    poisson_iters:        int   = 80     # Jacobi iterations per gauge update

# ── Geometry ──────────────────────────────────────────────────────────────────

def site_coords(site, cfg):
    z   = site // (cfg.hex_rows * cfg.hex_cols)
    rem = site %  (cfg.hex_rows * cfg.hex_cols)
    return rem // cfg.hex_cols, rem % cfg.hex_cols, z

def _idx(r, c, z, cfg):
    return ((z % cfg.layers) * cfg.hex_rows + (r % cfg.hex_rows)) * cfg.hex_cols \
           + (c % cfg.hex_cols)

def hex_neighbours(r, c, cfg):
    if r % 2 == 0:
        offs = [(-1,0),(-1,1),(0,-1),(0,1),(1,0),(1,1)]
    else:
        offs = [(-1,-1),(-1,0),(0,-1),(0,1),(1,-1),(1,0)]
    return [((r+dr) % cfg.hex_rows, (c+dc) % cfg.hex_cols) for dr,dc in offs]

def mesh_neighbours(r, c, z, cfg):
    return [_idx(nr, nc, z, cfg) for nr,nc in hex_neighbours(r, c, cfg)]

# ── Gauge field factory ───────────────────────────────────────────────────────

def make_gauge_fields(cfg):
    N = cfg.hex_rows * cfg.hex_cols * cfg.layers
    return {
        'phi':    np.zeros(N, dtype=np.float64),
        'A':      np.zeros(N, dtype=np.float64),
        'A_prev': np.zeros(N, dtype=np.float64),
    }

# ── Charge density ────────────────────────────────────────────────────────────

def compute_charge_density(lattice, quark_types, cfg):
    """
    quark_types : dict {site_idx: 'UP' | 'DOWN'}
    Returns rho[N] — charge per site.
    """
    N   = cfg.hex_rows * cfg.hex_cols * cfg.layers
    rho = np.zeros(N, dtype=np.float64)
    for site, qtype in quark_types.items():
        rho[site] = UP_CHARGE if qtype == 'UP' else DOWN_CHARGE
    for site in range(N):
        if int(lattice[site]) == LEPTON_SEED:
            rho[site] = LEPTON_CHARGE
    return rho

# ── Discrete Laplacian ────────────────────────────────────────────────────────

def hex_laplacian(field, cfg):
    """
    ∇²f[i] ≈ (1/n) Σⱼ∈nbrs (f[j] − f[i])
    n = number of neighbours (6 for hex-6).
    """
    N   = cfg.hex_rows * cfg.hex_cols * cfg.layers
    lap = np.zeros(N, dtype=np.float64)
    for site in range(N):
        r, c, z = site_coords(site, cfg)
        nbrs    = mesh_neighbours(r, c, z, cfg)
        n       = len(nbrs)
        lap[site] = sum(field[nb] - field[site] for nb in nbrs) / n
    return lap

# ── Poisson solver (Jacobi) ───────────────────────────────────────────────────

def jacobi_poisson(rho, cfg, n_iter=50):
    """
    Solve ∇²φ = −ρ on hex lattice via Jacobi iteration.

    Discrete Laplacian at site i with n neighbours:
        Σⱼ (φ[j] − φ[i]) / n = −ρ[i]
        Σⱼ φ[j] − n·φ[i]     = −n·ρ[i]
        φ[i] = (Σⱼ φ[j] + n·ρ[i]) / n    ← Jacobi update

    Positive source (proton quark) → positive φ.
    Lepton drifts toward max φ → EM attraction.
    """
    N   = cfg.hex_rows * cfg.hex_cols * cfg.layers
    phi = np.zeros(N, dtype=np.float64)

    # pre-cache neighbours for speed
    nbr_cache = []
    for site in range(N):
        r, c, z = site_coords(site, cfg)
        nbr_cache.append(mesh_neighbours(r, c, z, cfg))

    for _ in range(n_iter):
        phi_new = np.zeros(N, dtype=np.float64)
        for site in range(N):
            nbrs = nbr_cache[site]
            n    = len(nbrs)
            phi_new[site] = (sum(phi[nb] for nb in nbrs) + n * rho[site]) / n
        phi = phi_new

    return phi

# ── Maxwell scalar wave equation (leapfrog) ───────────────────────────────────

def maxwell_wave_step(gauge, rho, cfg, c=None):
    """
    Leapfrog step for scalar photon field A:

        A_new = 2·A − A_prev + c²·∇²A + coupling·ρ

    Stability: c² × max_eigenvalue(Laplacian) ≤ 1.
    For hex-6: max eigenvalue ≈ 2 (normalised), so c ≤ 1/√2 ≈ 0.707.
    Default cfg.photon_c = 0.4 gives comfortable margin.

    Source J = coupling × ρ drives the photon field from charge sites.
    """
    if c is None:
        c = cfg.photon_c

    A      = gauge['A']
    A_prev = gauge['A_prev']
    lap_A  = hex_laplacian(A, cfg)
    J      = rho * cfg.photon_coupling

    A_new  = 2.0 * A - A_prev + c**2 * lap_A + J

    gauge['A_prev'] = A.copy()
    gauge['A']      = A_new

# ── Full gauge update ─────────────────────────────────────────────────────────

def update_gauge(gauge, lattice, quark_types, cfg):
    """
    Complete gauge field update per timestep:
      1. Charge density ρ from current matter state
      2. Poisson solve: ∇²φ_q = −ρ_quarks  (quark-only Coulomb for lepton dynamics)
         Lepton self-energy is excluded so leptons always see the proton's +1 field.
         In real QED the electron responds to the external (proton) potential, not
         the self-consistent electron+proton combined neutral potential.
      3. φ_full = Poisson of full ρ (quarks + leptons) for completeness
      4. Maxwell leapfrog: ∂²A/∂t² = c²∇²A + J  (photon wave, full ρ)
    """
    rho_full     = compute_charge_density(lattice, quark_types, cfg)
    # Quark-only potential for lepton EM dynamics
    rho_quarks   = rho_full.copy()
    for site in range(cfg.hex_rows * cfg.hex_cols * cfg.layers):
        if int(lattice[site]) == LEPTON_SEED:
            rho_quarks[site] = 0.0
    gauge['phi'] = jacobi_poisson(rho_quarks, cfg, n_iter=cfg.poisson_iters)
    maxwell_wave_step(gauge, rho_full, cfg)

# ── EM force on lepton ────────────────────────────────────────────────────────

def em_force_on_lepton(phi, A, site, cfg):
    """
    Lepton (charge −1) feels Coulomb force F = −q∇φ = +∇φ.
    Moves toward neighbour with maximum combined EM potential φ + 0.3·A.

    Returns flat site index of the target neighbour.
    """
    r, c, z = site_coords(site, cfg)
    nbrs    = mesh_neighbours(r, c, z, cfg)
    potentials = [phi[nb] + 0.3 * A[nb] for nb in nbrs]
    return nbrs[int(np.argmax(potentials))]
