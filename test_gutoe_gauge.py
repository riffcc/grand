#!/usr/bin/env python3
"""TDD tests for gutoe_gauge.py — U(1) bond gauge field physics."""

import numpy as np
import pytest
import sys
sys.path.insert(0, '/mnt/riffcastle/castle/garage/grand-2026')

# These will fail until gutoe_gauge.py exists — that's expected (red phase)
try:
    from gutoe_gauge import (
        make_gauge_fields, compute_charge_density, jacobi_poisson,
        maxwell_wave_step, update_gauge, em_force_on_lepton,
        LatticeConfig, site_coords, mesh_neighbours,
    )
    IMPORT_OK = True
except ImportError:
    IMPORT_OK = False

pytestmark = pytest.mark.skipif(not IMPORT_OK, reason="gutoe_gauge not yet implemented")

@pytest.fixture
def small_cfg():
    return LatticeConfig(hex_rows=8, hex_cols=8, layers=1)

@pytest.fixture
def cfg():
    return LatticeConfig(hex_rows=12, hex_cols=12, layers=12)

def center_site(cfg):
    return cfg.hex_rows * cfg.hex_cols * cfg.layers // 2

# ── Poisson / Coulomb ─────────────────────────────────────────────────────────

def test_poisson_zero_charge_gives_zero(small_cfg):
    """No charges → φ = 0 everywhere."""
    N = small_cfg.hex_rows * small_cfg.hex_cols * small_cfg.layers
    rho = np.zeros(N)
    phi = jacobi_poisson(rho, small_cfg, n_iter=20)
    assert np.allclose(phi, 0, atol=1e-10)

def test_poisson_single_positive_charge_positive_phi(small_cfg):
    """Single +1 charge → φ > 0 at charge site and neighbors."""
    N = small_cfg.hex_rows * small_cfg.hex_cols * small_cfg.layers
    rho = np.zeros(N)
    c = center_site(small_cfg)
    rho[c] = 1.0
    phi = jacobi_poisson(rho, small_cfg, n_iter=100)
    assert phi[c] > 0, "φ at charge site should be positive"
    r, col, z = site_coords(c, small_cfg)
    nbrs = mesh_neighbours(r, col, z, small_cfg)
    assert all(phi[c] >= phi[nb] for nb in nbrs), "φ peaks at charge site"

def test_poisson_negative_charge_negative_phi(small_cfg):
    """Single -1 charge → φ < 0 at charge site."""
    N = small_cfg.hex_rows * small_cfg.hex_cols * small_cfg.layers
    rho = np.zeros(N)
    c = center_site(small_cfg)
    rho[c] = -1.0
    phi = jacobi_poisson(rho, small_cfg, n_iter=100)
    assert phi[c] < 0

def test_poisson_phi_decays_with_distance(small_cfg):
    """φ at charge site > φ at 2-hop distance > φ at 4-hop distance."""
    N = small_cfg.hex_rows * small_cfg.hex_cols * small_cfg.layers
    rho = np.zeros(N)
    c = center_site(small_cfg)
    rho[c] = 1.0
    phi = jacobi_poisson(rho, small_cfg, n_iter=200)
    r, col, z = site_coords(c, small_cfg)
    hop1 = mesh_neighbours(r, col, z, small_cfg)[0]
    r1, c1, z1 = site_coords(hop1, small_cfg)
    hop2 = [nb for nb in mesh_neighbours(r1, c1, z1, small_cfg) if nb != c][0]
    assert phi[c] > phi[hop1] > phi[hop2], \
        f"φ should decay: {phi[c]:.4f} > {phi[hop1]:.4f} > {phi[hop2]:.4f}"

def test_poisson_linearity(small_cfg):
    """Poisson is linear: φ(2ρ) = 2φ(ρ)."""
    N = small_cfg.hex_rows * small_cfg.hex_cols * small_cfg.layers
    rho = np.zeros(N)
    rho[center_site(small_cfg)] = 1.0
    phi1 = jacobi_poisson(rho, small_cfg, n_iter=100)
    phi2 = jacobi_poisson(2*rho, small_cfg, n_iter=100)
    assert np.allclose(phi2, 2*phi1, rtol=0.01), "Poisson must be linear"

# ── Maxwell wave equation ─────────────────────────────────────────────────────

def test_maxwell_zero_source_no_field(small_cfg):
    """No source, no initial field → A stays zero."""
    N = small_cfg.hex_rows * small_cfg.hex_cols * small_cfg.layers
    gauge = make_gauge_fields(small_cfg)
    rho = np.zeros(N)
    for _ in range(10):
        maxwell_wave_step(gauge, rho, small_cfg)
    assert np.allclose(gauge['A'], 0, atol=1e-10)

def test_maxwell_source_drives_field(small_cfg):
    """A source at center drives A away from zero."""
    N = small_cfg.hex_rows * small_cfg.hex_cols * small_cfg.layers
    gauge = make_gauge_fields(small_cfg)
    rho = np.zeros(N)
    rho[center_site(small_cfg)] = 1.0
    for _ in range(5):
        maxwell_wave_step(gauge, rho, small_cfg)
    assert gauge['A'][center_site(small_cfg)] != 0, "Source should drive A field"

def test_maxwell_wave_propagates(small_cfg):
    """Field driven at center for 1 step should reach neighbors within ~5 steps."""
    N = small_cfg.hex_rows * small_cfg.hex_cols * small_cfg.layers
    gauge = make_gauge_fields(small_cfg)
    rho = np.zeros(N)
    c = center_site(small_cfg)
    rho[c] = 1.0
    r, col, z = site_coords(c, small_cfg)
    hop1 = mesh_neighbours(r, col, z, small_cfg)[0]
    # Run long enough for signal to propagate
    for _ in range(20):
        maxwell_wave_step(gauge, rho, small_cfg)
    # Neighbor should have non-zero A
    assert abs(gauge['A'][hop1]) > 0, \
        f"Wave should propagate to neighbors. A[hop1]={gauge['A'][hop1]}"

# ── EM force on lepton ────────────────────────────────────────────────────────

def test_lepton_attracted_to_positive_phi(small_cfg):
    """γ⁰ at site X feels force toward site with higher φ."""
    N = small_cfg.hex_rows * small_cfg.hex_cols * small_cfg.layers
    phi = np.zeros(N)
    A = np.zeros(N)
    c = center_site(small_cfg)
    r, col, z = site_coords(c, small_cfg)
    nbrs = mesh_neighbours(r, col, z, small_cfg)
    # Set φ high at one specific neighbor
    target = nbrs[2]
    phi[target] = 10.0
    # Lepton at c should be pushed toward target
    force_target = em_force_on_lepton(phi, A, c, small_cfg)
    assert force_target == target, \
        f"Lepton should move toward max φ={phi[target]} at {target}, got {force_target}"

def test_lepton_moves_away_from_negative_phi(small_cfg):
    """γ⁰ flees negative φ (another negative charge)."""
    N = small_cfg.hex_rows * small_cfg.hex_cols * small_cfg.layers
    phi = np.zeros(N)
    A = np.zeros(N)
    c = center_site(small_cfg)
    r, col, z = site_coords(c, small_cfg)
    nbrs = mesh_neighbours(r, col, z, small_cfg)
    # Set φ very negative at one neighbor (another lepton), high elsewhere
    for nb in nbrs:
        phi[nb] = 1.0
    phi[nbrs[0]] = -10.0  # repulsion from one direction
    force_target = em_force_on_lepton(phi, A, c, small_cfg)
    assert force_target != nbrs[0], "Lepton should flee negative φ"

# ── Charge density ────────────────────────────────────────────────────────────

def test_charge_density_proton_net_positive(small_cfg):
    """A proton (2 UP + 1 DOWN quarks) has net charge +1."""
    from gutoe_gauge import UP_CHARGE, DOWN_CHARGE
    net = 2 * UP_CHARGE + 1 * DOWN_CHARGE
    assert abs(net - 1.0) < 1e-9, f"Proton charge should be +1, got {net}"

def test_charge_density_lepton_negative(small_cfg):
    """γ⁰ lepton has charge -1."""
    from gutoe_gauge import LEPTON_CHARGE
    assert LEPTON_CHARGE == -1.0

def test_compute_charge_density_shape(small_cfg):
    """charge_density returns array of shape N."""
    N = small_cfg.hex_rows * small_cfg.hex_cols * small_cfg.layers
    lattice = np.zeros(N, dtype=np.int8)
    quark_types = {0: 'UP', 1: 'UP', 2: 'DOWN'}
    rho = compute_charge_density(lattice, quark_types, small_cfg)
    assert rho.shape == (N,)
    assert rho[0] > 0   # UP quark
    assert rho[2] < 0   # DOWN quark

if __name__ == '__main__':
    pytest.main([__file__, '-v'])
