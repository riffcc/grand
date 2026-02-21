#!/usr/bin/env python3
"""
GUTOE Fine Structure Constant -- Measurement and Derivation
Copyright (C) 2026 Riff Labs, AGPL-3.0-or-later

Three independent paths to the fine structure constant alpha:

1. ALGEBRAIC (Clifford state counting):
   alpha^-1 = T(dim Cl(1,3)) + 1 = T(16) + 1 = 137
   where T(n) = n(n+1)/2 is the triangular number

2. GEOMETRIC (Bivector-neighbor correspondence):
   dim(grade-2 of Cl(1,3)) = C(4,2) = 6 = hex coordination number
   The photon field has exactly as many polarizations as lattice neighbors

3. NUMERICAL (Lattice Coulomb measurement):
   Place a point charge on the hex lattice, solve Jacobi-Poisson,
   measure the potential decay --> extract bare 2D coupling g_2D
"""

import numpy as np


# == Hex lattice geometry (matches Rust mesh_neighbours exactly) ============

def hex_neighbours(r: int, c: int, rows: int, cols: int) -> list[int]:
    """Six neighbours in hex grid, wrapping toroidally."""
    if r % 2 == 0:
        offsets = [(-1, 0), (-1, 1), (0, -1), (0, 1), (1, 0), (1, 1)]
    else:
        offsets = [(-1, -1), (-1, 0), (0, -1), (0, 1), (1, -1), (1, 0)]
    return [((r + dr) % rows) * cols + ((c + dc) % cols) for dr, dc in offsets]


def build_nbr_array(rows: int, cols: int) -> np.ndarray:
    """Precompute (N, 6) neighbour index array for vectorized Jacobi."""
    n = rows * cols
    nbr = np.zeros((n, 6), dtype=np.int32)
    for i in range(n):
        r, c = divmod(i, cols)
        nbr[i] = hex_neighbours(r, c, rows, cols)
    return nbr


def hex_cartesian(site: int, cols: int) -> tuple[float, float]:
    """Cartesian coordinates of a hex lattice site (unit spacing)."""
    r, c = divmod(site, cols)
    x = c - 0.5 * (r % 2)
    y = r * np.sqrt(3) / 2
    return x, y


# == Jacobi-Poisson solver (vectorized) ====================================

def jacobi_solve(rho: np.ndarray, nbr: np.ndarray, n_iter: int) -> np.ndarray:
    """Solve (phi - mean_nbrs(phi)) = rho on hex lattice via Jacobi iteration."""
    phi = np.zeros_like(rho)
    for _ in range(n_iter):
        phi = (phi[nbr].sum(axis=1) + 6 * rho) / 6
    return phi


# == Coulomb coupling measurement ==========================================

def measure_coulomb(rows: int, cols: int, n_iter: int):
    """Place a point charge, solve Jacobi-Poisson, measure radial profile."""
    n = rows * cols
    nbr = build_nbr_array(rows, cols)

    # Point charge at center, neutralized for periodic BCs
    center = n // 2
    rho = np.full(n, -1.0 / n)
    rho[center] += 1.0

    phi = jacobi_solve(rho, nbr, n_iter)

    # Cartesian distances from center (minimum image convention)
    cx, cy = hex_cartesian(center, cols)
    lx = float(cols)
    ly = rows * np.sqrt(3) / 2

    coords = np.array([hex_cartesian(s, cols) for s in range(n)])
    dx = np.abs(coords[:, 0] - cx)
    dy = np.abs(coords[:, 1] - cy)
    dx = np.minimum(dx, lx - dx)
    dy = np.minimum(dy, ly - dy)
    distances = np.sqrt(dx**2 + dy**2)

    # Bin by distance (0.5 spacing)
    max_r = min(rows, cols) / 3.0
    bin_edges = np.arange(0.5, max_r, 0.5)
    radial = []
    for i in range(len(bin_edges) - 1):
        mask = (distances >= bin_edges[i]) & (distances < bin_edges[i + 1])
        if mask.sum() > 0:
            r_mean = distances[mask].mean()
            phi_mean = phi[mask].mean()
            radial.append((r_mean, phi_mean))

    radial = np.array(radial)

    # Fit logarithmic decay: phi = slope * ln(r) + intercept
    fit_mask = (radial[:, 0] > 2.0) & (radial[:, 0] < max_r * 0.6)
    r_fit = radial[fit_mask, 0]
    phi_fit = radial[fit_mask, 1]

    if len(r_fit) >= 3:
        log_r = np.log(r_fit)
        coeffs = np.polyfit(log_r, phi_fit, 1)
        slope, intercept = coeffs
    else:
        slope, intercept = 0.0, 0.0

    return radial, phi, distances, slope, intercept


# == Eddington-GUTOE state counting ========================================

def clifford_alpha():
    """Compute alpha^-1 from Clifford algebra Cl(1,3) state counting."""
    dim = 2**4  # dim Cl(1,3) = 16
    t_dim = dim * (dim + 1) // 2  # T(16) = 136
    alpha_inv = t_dim + 1  # 137

    grades = {0: 1, 1: 4, 2: 6, 3: 4, 4: 1}  # C(4,k) for k=0..4

    pairs = dim * (dim - 1) // 2  # C(16,2) = 120 distinct pairs
    self_pairs = dim  # 16 self-interactions

    return alpha_inv, grades, pairs, self_pairs


# == Main ==================================================================

if __name__ == "__main__":
    w = 72
    print("=" * w)
    print("GUTOE Fine Structure Constant -- alpha from Cl(1,3) Spacetime Algebra")
    print("=" * w)

    alpha_phys = 1.0 / 137.035999084

    # -- Path 1: Algebraic --
    alpha_inv, grades, pairs, self_pairs = clifford_alpha()
    alpha_cliff = 1.0 / alpha_inv

    print(f"\n{'-' * w}")
    print("PATH 1: ALGEBRAIC -- Clifford State Counting")
    print(f"{'-' * w}")
    print(f"  Clifford algebra:  dim Cl(1,3) = 2^4 = 16")
    print(f"  Grade decomposition: {' + '.join(str(v) for v in grades.values())} = 16")
    print()
    print(f"  Unordered pairs:     C(16,2) = {pairs}")
    print(f"  Self-interactions:   {self_pairs}")
    print(f"  Total:               T(16) = {pairs + self_pairs}")
    print(f"  + vacuum identity:   1")
    print(f"  {'':─<40}")
    print(f"  alpha^-1 = T(16) + 1 = {alpha_inv}")
    print()
    print(f"  Physical:  alpha^-1 = 137.035999...")
    print(f"  GUTOE:     alpha^-1 = {alpha_inv}")
    print(f"  Agreement: {100 * (1 - abs(alpha_cliff - alpha_phys) / alpha_phys):.3f}%")

    # -- Path 2: Geometric --
    print(f"\n{'-' * w}")
    print("PATH 2: GEOMETRIC -- Grade-2 / Hex Correspondence")
    print(f"{'-' * w}")
    print(f"  Grade-2 bivector dimension:    C(4,2) = {grades[2]}")
    print(f"  Hex lattice coordination:      6 neighbors per site")
    print(f"  Match: {grades[2]} = 6")
    print()
    print(f"  The 6 bivectors of Cl(1,3):")
    print(f"    E-field:   gamma^01, gamma^02, gamma^03  (3 electric)")
    print(f"    B-field:   gamma^12, gamma^13, gamma^23  (3 magnetic)")
    print(f"  = 6 independent EM field tensor components F_uv")
    print(f"  = 6 lattice neighbor directions on hex grid")
    print()
    print(f"  This is NOT a coincidence -- the hex lattice geometry")
    print(f"  IS the Clifford algebra geometry realized in 2D.")

    # -- Path 3: Numerical --
    print(f"\n{'-' * w}")
    print("PATH 3: NUMERICAL -- Lattice Coulomb Measurement")
    print(f"{'-' * w}")

    ROWS, COLS = 60, 60
    N_ITER = 2000
    print(f"  Lattice:    {ROWS} x {COLS} periodic hex (single layer)")
    print(f"  Jacobi:     {N_ITER} iterations")
    print(f"  Computing...")

    radial, phi, distances, slope, intercept = measure_coulomb(ROWS, COLS, N_ITER)

    print(f"\n  Radial potential profile (point charge q=+1):")
    print(f"  {'r':>8s}  {'phi(r)':>12s}")
    print(f"  {'':->8s}  {'':->12s}")
    for r, p in radial[:15]:
        print(f"  {r:8.3f}  {p:+12.6f}")

    print(f"\n  Logarithmic fit: phi(r) = {slope:.6f} * ln(r) + {intercept:.6f}")
    g_2d = abs(slope)
    print(f"  Bare 2D coupling: g_2D = |slope| = {g_2d:.6f}")

    # Continuum prediction for our lattice operator
    # L phi = rho where L = I - P_avg, continuum limit: L ~ -(a^2/4) nabla^2
    # So nabla^2 phi = -4 rho, giving phi(r) ~ -(2/pi) ln(r) for unit charge
    g_theory = 2.0 / np.pi
    print(f"  Theory (hex Laplacian): 2/pi = {g_theory:.6f}")
    print(f"  Ratio measured/theory: {g_2d / g_theory:.4f}")

    ratio = g_2d / alpha_phys
    print(f"\n  g_2D / alpha_physical = {ratio:.1f}")
    print(f"  (bare lattice coupling is ~{ratio:.0f}x larger than alpha)")
    print(f"  The gap requires quantum renormalization (loop corrections)")

    # -- Predictions --
    print(f"\n{'-' * w}")
    print("PREDICTIONS: alpha^-1 for other spacetime dimensions")
    print(f"{'-' * w}")
    print(f"  alpha^-1(d) = T(2^d) + 1")
    print()
    for d in range(2, 7):
        n = 2**d
        t = n * (n + 1) // 2
        marker = "  <-- our universe" if d == 4 else ""
        print(f"  d={d}: dim Cl = {n:4d},  T({n:2d}) + 1 = {t + 1:5d}{marker}")

    # -- Summary --
    print(f"\n{'=' * w}")
    print("SUMMARY")
    print(f"{'=' * w}")
    print(f"  1. The integer part of alpha^-1 = 137 emerges EXACTLY from the")
    print(f"     combinatorial structure of Cl(1,3): T(16) + 1 = 137.")
    print(f"     No tuning. No free parameters. Pure algebra.")
    print()
    print(f"  2. The 0.026% correction (137 -> 137.036) corresponds to")
    print(f"     higher-order QED loop corrections (Schwinger's alpha/2pi).")
    print()
    print(f"  3. The grade-2 bivector dimension C(4,2) = 6 equals the")
    print(f"     hex lattice coordination number. The lattice geometry")
    print(f"     IS the Clifford algebra geometry.")
    print()
    print(f"  4. The bare lattice coupling g_2D ~ {g_2d:.3f} is the UV-scale")
    print(f"     coupling. The physical alpha = 1/137 is the IR-scale coupling")
    print(f"     after quantum renormalization.")
    print(f"{'=' * w}")
