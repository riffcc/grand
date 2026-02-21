#!/usr/bin/env python3
"""
GUTOE: Hexadecimal lattice — Cl(1,3) state space on hex mesh.

HEXADECIMAL = 16 states per cell, each a Cl(1,3) basis multivector.
HEX lattice = triangular (6-fold) spatial topology.

The 16 states are the 16 basis elements of Cl(1,3):
  Grade 0 (1):  {1}                  — scalar
  Grade 1 (4):  {γ⁰,γ¹,γ²,γ³}       — vectors
  Grade 2 (6):  {γ⁰¹,...,γ²³}        — bivectors
  Grade 3 (4):  {γ⁰¹²,...,γ¹²³}      — trivectors
  Grade 4 (1):  {γ⁰¹²³}             — pseudoscalar

Encoded as 4-bit multi-indices (subsets of {0,1,2,3}):
  bit 0 = index 0 (timelike), bits 1-3 = indices 1-3 (spacelike)

Z₃ acts on spacelike indices: 1→2→3→1 (bit 0 fixed).
This partitions the 16 states into:
  4 fixed points: {1, γ⁰, γ¹²³, γ⁰¹²³}
  4 orbits of 3:  {γ¹,γ²,γ³}, {γ⁰¹,γ⁰²,γ⁰³}, {γ¹²,γ²³,γ¹³}, {γ⁰¹²,γ⁰²³,γ⁰¹³}

Veracity = f(Hamming distance of multi-indices):
  d=0: 1.0,  d=1: √3/2,  d=2: 1/2,  d≥3: 0

Curvature = max Z₃ triple-junction diversity across all orbits:
  1 orbit member present → 0    (domain interior)
  2 orbit members        → 0.5  (domain boundary)
  3 orbit members        → 1.0  (triple junction)
  Same {0, 0.5, 1.0} scale as original Z₃ system.

Dynamics:
  - Differentiation: VOID → γ¹ (spacelike vector = original SINE)
  - Z₃ cycle: permute spacelike indices within orbits
  - Clifford interaction: XOR of multi-indices (creates higher grades)
  - Alignment: local majority wins

TIME IS NOT A LATTICE AXIS. Time is the simulation evolution.

Copyright (C) 2026 Riff Labs, AGPL-3.0-or-later
"""

import numpy as np
from dataclasses import dataclass
from collections import Counter
import sys

# ── Cl(1,3) basis: 16 states as 4-bit multi-indices ────────────────────────

VOID = 0  # not a Cl(1,3) element

# State s ∈ {1,...,16} encodes multi-index mi = s - 1 ∈ {0,...,15}
#   bit 0 → index 0 (timelike)
#   bit 1 → index 1 (spacelike)
#   bit 2 → index 2 (spacelike)
#   bit 3 → index 3 (spacelike)

CL_NAMES = [
    "1",      # mi=0   grade 0  scalar
    "γ⁰",     # mi=1   grade 1  timelike vector
    "γ¹",     # mi=2   grade 1  spacelike vector
    "γ⁰¹",    # mi=3   grade 2  boost x
    "γ²",     # mi=4   grade 1  spacelike vector
    "γ⁰²",    # mi=5   grade 2  boost y
    "γ¹²",    # mi=6   grade 2  rotation xy
    "γ⁰¹²",   # mi=7   grade 3  timelike trivector
    "γ³",     # mi=8   grade 1  spacelike vector
    "γ⁰³",    # mi=9   grade 2  boost z
    "γ¹³",    # mi=10  grade 2  rotation xz
    "γ⁰¹³",   # mi=11  grade 3  timelike trivector
    "γ²³",    # mi=12  grade 2  rotation yz
    "γ⁰²³",   # mi=13  grade 3  timelike trivector
    "γ¹²³",   # mi=14  grade 3  spatial pseudovector
    "γ⁰¹²³",  # mi=15  grade 4  pseudoscalar
]

# ── Precomputed tables ──────────────────────────────────────────────────────

_GRADE_TABLE = [-1] + [bin(mi).count('1') for mi in range(16)]

def _make_z3_table():
    table = [VOID]
    for s in range(1, 17):
        mi = s - 1
        b0 = (mi >> 0) & 1  # timelike — stays
        b1 = (mi >> 1) & 1  # index 1 → 2
        b2 = (mi >> 2) & 1  # index 2 → 3
        b3 = (mi >> 3) & 1  # index 3 → 1
        new_mi = b0 | (b3 << 1) | (b1 << 2) | (b2 << 3)
        table.append(new_mi + 1)
    return table

_Z3_TABLE = _make_z3_table()

_SQRT3_HALF = np.sqrt(3) / 2

def _make_veracity_table():
    table = {}
    for s1 in range(17):
        for s2 in range(17):
            if s1 == VOID or s2 == VOID:
                table[(s1, s2)] = 0.0
            elif s1 == s2:
                table[(s1, s2)] = 1.0
            else:
                d = bin((s1 - 1) ^ (s2 - 1)).count('1')
                if d == 1:
                    table[(s1, s2)] = _SQRT3_HALF
                elif d == 2:
                    table[(s1, s2)] = 0.5
                else:
                    table[(s1, s2)] = 0.0
    return table

_VERACITY_TABLE = _make_veracity_table()

# Z₃ orbits (sets of state numbers that cycle into each other)
Z3_ORBITS = [
    frozenset({3, 5, 9}),    # {γ¹, γ², γ³}         spacelike vectors
    frozenset({4, 6, 10}),   # {γ⁰¹, γ⁰², γ⁰³}      boosts
    frozenset({7, 13, 11}),  # {γ¹², γ²³, γ¹³}       rotations
    frozenset({8, 14, 12}),  # {γ⁰¹², γ⁰²³, γ⁰¹³}    timelike trivectors
]

# State → which orbit index (0-3), or -1 for fixed point
_ORBIT_TABLE = {}
for oi, orbit in enumerate(Z3_ORBITS):
    for s in orbit:
        _ORBIT_TABLE[s] = oi
# Fixed points: 1 (scalar), 2 (γ⁰), 15 (γ¹²³), 16 (γ⁰¹²³)
for s in [1, 2, 15, 16]:
    _ORBIT_TABLE[s] = -1

# The seed state: γ¹ (spacelike vector, state 3) = the original SINE
SEED_STATE = 3


# ── Convenience functions ───────────────────────────────────────────────────

def cl_name(s: int) -> str:
    if s == VOID: return "VOID"
    return CL_NAMES[s - 1]

def cl_grade(s: int) -> int:
    return _GRADE_TABLE[s]

def cl_z3_cycle(s: int) -> int:
    return _Z3_TABLE[s]

def cl_product(s1: int, s2: int) -> int:
    """Clifford product on multi-indices: XOR (symmetric difference)."""
    if s1 == VOID or s2 == VOID: return VOID
    return ((s1 - 1) ^ (s2 - 1)) + 1


# ── Hex lattice (6-fold triangular) ─────────────────────────────────────────

@dataclass
class LatticeConfig:
    hex_rows: int = 12
    hex_cols: int = 12
    layers: int = 12
    differentiation_prob: float = 0.02
    cycle_prob: float = 0.05
    clifford_prob: float = 0.03
    alignment_strength: float = 0.15
    quark_threshold: float = 0.6

def idx(r, c, z, cfg):
    return ((z % cfg.layers) * cfg.hex_rows + (r % cfg.hex_rows)) * cfg.hex_cols + (c % cfg.hex_cols)

def hex_planar_neighbours(r, c, cfg):
    if r % 2 == 0:
        offsets = [(-1, 0), (-1, 1), (0, -1), (0, 1), (1, 0), (1, 1)]
    else:
        offsets = [(-1, -1), (-1, 0), (0, -1), (0, 1), (1, -1), (1, 0)]
    return [((r + dr) % cfg.hex_rows, (c + dc) % cfg.hex_cols) for dr, dc in offsets]

def mesh_neighbours(r, c, z, cfg):
    """Hex-6: in-plane hex neighbours. Triangular topology."""
    return [idx(nr, nc, z, cfg) for nr, nc in hex_planar_neighbours(r, c, cfg)]

def init_lattice(cfg):
    return np.zeros(cfg.hex_rows * cfg.hex_cols * cfg.layers, dtype=np.int8)


# ── Field computation ───────────────────────────────────────────────────────

def compute_local_fields(lattice, site, r, c, z, cfg):
    """Compute (veracity, curvature, gradient) for a cell.

    Veracity: Hamming-based coherence with neighbours.
    Curvature: max Z₃ triple-junction diversity across all 4 orbits.
      - 1 orbit member present → 0
      - 2 orbit members        → 0.5
      - 3 orbit members        → 1.0  (triple junction)
    Gradient: average veracity deficit.
    """
    state = int(lattice[site])
    if state == VOID:
        return 0.0, 0.0, 0.0

    nbrs = mesh_neighbours(r, c, z, cfg)
    total_veracity = 0.0
    grad_sum = 0.0
    nbr_state_set = set()

    for ni in nbrs:
        ns = int(lattice[ni])
        v = _VERACITY_TABLE[(state, ns)]
        total_veracity += v
        grad_sum += 1.0 - v
        if ns != VOID:
            nbr_state_set.add(ns)

    n_total = len(nbrs)
    avg_veracity = total_veracity / n_total
    field_grad = grad_sum / n_total

    # Z₃ orbit curvature: max triple-junction diversity
    # Same {0, 0.5, 1.0} scale as original Z₃ system
    z3_curvature = 0.0
    for orbit in Z3_ORBITS:
        n_present = len(orbit & nbr_state_set)
        orbit_curv = (n_present - 1) / 2
        if orbit_curv > z3_curvature:
            z3_curvature = orbit_curv

    return avg_veracity, z3_curvature, field_grad


# ── Dynamics ────────────────────────────────────────────────────────────────

def step(lattice, rng, cfg):
    new = lattice.copy()
    N = cfg.hex_rows * cfg.hex_cols * cfg.layers

    for site in range(N):
        z = site // (cfg.hex_rows * cfg.hex_cols)
        rem = site % (cfg.hex_rows * cfg.hex_cols)
        r = rem // cfg.hex_cols
        c = rem % cfg.hex_cols
        state = int(lattice[site])

        if state == VOID:
            # Differentiation: VOID → γ¹ (spacelike vector = original SINE)
            if rng.random() < cfg.differentiation_prob:
                new[site] = SEED_STATE
                continue
            # Spread from active neighbours
            nbrs = mesh_neighbours(r, c, z, cfg)
            active = sum(1 for ni in nbrs if lattice[ni] != VOID)
            total = len(nbrs)
            if active >= max(2, total // 4):
                if rng.random() < active / total * 0.4:
                    new[site] = SEED_STATE
        else:
            r_val = rng.random()
            if r_val < cfg.cycle_prob:
                # Z₃ cycle: rotate spacelike indices
                # Fixed points (scalar, γ⁰, γ¹²³, γ⁰¹²³) don't change
                new[site] = _Z3_TABLE[state]

            elif r_val < cfg.cycle_prob + cfg.clifford_prob:
                # Clifford interaction: XOR with random active neighbour
                # grade-1 × grade-1 → grade-0 or grade-2
                # grade-2 × grade-1 → grade-1 or grade-3
                # This is how higher grades emerge from vector seeds
                nbrs = mesh_neighbours(r, c, z, cfg)
                active_nbrs = [int(lattice[ni]) for ni in nbrs
                               if lattice[ni] != VOID]
                if active_nbrs:
                    partner = active_nbrs[rng.integers(len(active_nbrs))]
                    new[site] = cl_product(state, partner)

            elif r_val < cfg.cycle_prob + cfg.clifford_prob + cfg.alignment_strength:
                # Alignment: adopt majority neighbour state
                nbrs = mesh_neighbours(r, c, z, cfg)
                nbr_states = [int(lattice[ni]) for ni in nbrs
                              if lattice[ni] != VOID]
                if nbr_states:
                    counts = Counter(nbr_states)
                    new[site] = counts.most_common(1)[0][0]

    return new


# ── Quark detection ─────────────────────────────────────────────────────────

@dataclass
class Quark:
    site: int
    r: int
    c: int
    z: int
    quark_type: str
    binding_coherence: float
    veracity: float
    curvature: float
    grade: int
    orbit: int  # Z₃ orbit index (-1 for fixed points)

def detect_quarks(lattice, cfg):
    quarks = []
    N = cfg.hex_rows * cfg.hex_cols * cfg.layers

    for site in range(N):
        state = int(lattice[site])
        if state == VOID:
            continue

        z = site // (cfg.hex_rows * cfg.hex_cols)
        rem = site % (cfg.hex_rows * cfg.hex_cols)
        r = rem // cfg.hex_cols
        c = rem % cfg.hex_cols

        v, curv, grad = compute_local_fields(lattice, site, r, c, z, cfg)
        bc = v / (1 + grad)

        if bc >= cfg.quark_threshold:
            qtype = "UP" if v > curv else "DOWN"
            quarks.append(Quark(site, r, c, z, qtype, bc, v, curv,
                                cl_grade(state), _ORBIT_TABLE.get(state, -1)))

    return quarks


def find_protons(quarks, cfg):
    """Count proton-like (uud) clusters: closed triangles on hex mesh."""
    quark_set = {q.site: q for q in quarks}
    strict_count = 0
    loose_count = 0
    used_strict = set()
    used_loose = set()

    nbr_cache = {}
    for q in quarks:
        nbr_cache[q.site] = set(mesh_neighbours(q.r, q.c, q.z, cfg))

    for q in quarks:
        if q.quark_type != "DOWN" or q.site in used_loose:
            continue

        up_nbrs = [quark_set[ni] for ni in nbr_cache[q.site]
                    if ni in quark_set and quark_set[ni].quark_type == "UP"
                    and ni not in used_strict]

        if len(up_nbrs) < 2:
            continue

        found_strict = False
        for i in range(len(up_nbrs)):
            if found_strict:
                break
            for j in range(i + 1, len(up_nbrs)):
                p1, p2 = up_nbrs[i].site, up_nbrs[j].site
                if p1 in nbr_cache.get(p2, set()):
                    strict_count += 1
                    used_strict.update([q.site, p1, p2])
                    found_strict = True
                    break

        if q.site not in used_loose and len(up_nbrs) >= 2:
            loose_count += 1
            used_loose.add(q.site)

    return strict_count, loose_count


# ── Grade & orbit distribution ──────────────────────────────────────────────

def grade_distribution(lattice):
    counts = {-1: 0, 0: 0, 1: 0, 2: 0, 3: 0, 4: 0}
    for s in lattice:
        counts[_GRADE_TABLE[int(s)]] += 1
    return counts

def orbit_distribution(lattice):
    """Count cells by Z₃ orbit. -1=fixed point, 0-3=orbit index."""
    counts = {-2: 0, -1: 0, 0: 0, 1: 0, 2: 0, 3: 0}  # -2=VOID
    for s in lattice:
        s = int(s)
        if s == VOID:
            counts[-2] += 1
        else:
            counts[_ORBIT_TABLE[s]] += 1
    return counts


# ── Main ────────────────────────────────────────────────────────────────────

def run_simulation(seed: int = 42, steps: int = 200):
    cfg = LatticeConfig()
    rng = np.random.default_rng(seed)
    lattice = init_lattice(cfg)
    N = cfg.hex_rows * cfg.hex_cols * cfg.layers

    print(f"GUTOE: Hexadecimal lattice — Cl(1,3) states on hex mesh")
    print(f"{'='*80}")
    print(f"Lattice: {cfg.hex_rows}×{cfg.hex_cols}×{cfg.layers} = {N} cells")
    print(f"States: 16 Cl(1,3) basis multivectors (grades 0-4)")
    print(f"Topology: hex-6 (triangular, in-plane)")
    print(f"Seed state: {cl_name(SEED_STATE)} (spacelike vector, like original SINE)")
    print(f"Z₃ cycle: spacelike index rotation 1→2→3→1")
    print(f"Clifford interaction: {cfg.clifford_prob:.0%} (XOR → higher grades)")
    print(f"Curvature: Z₃ orbit triple-junction (0 / 0.5 / 1.0)")
    print(f"Veracity: Hamming distance (d=0:1, d=1:√3/2, d=2:½, d≥3:0)")
    print(f"RNG seed: {seed}")
    print(f"{'='*80}")
    print()

    milestones = list(range(1, 31)) + [40, 50, 75, 100, 150, 200]

    for t in range(1, steps + 1):
        lattice = step(lattice, rng, cfg)

        if t in milestones or t == steps:
            quarks = detect_quarks(lattice, cfg)
            up = sum(1 for q in quarks if q.quark_type == "UP")
            dn = sum(1 for q in quarks if q.quark_type == "DOWN")
            strict, loose = find_protons(quarks, cfg)
            active = sum(1 for s in lattice if s != 0) / N
            gd = grade_distribution(lattice)

            ratio_str = '∞' if dn == 0 else f'{up/dn:.2f}'
            g_str = (f"g0={gd[0]:3d} g1={gd[1]:4d} "
                     f"g2={gd[2]:4d} g3={gd[3]:3d} g4={gd[4]:2d}")
            print(f"t={t:4d} | act={active:5.1%} | "
                  f"UP={up:4d} DN={dn:4d} r={ratio_str:>6s} | "
                  f"p={strict:3d} L={loose:3d} | {g_str}")

    # ── Final ───────────────────────────────────────────────────────────
    print()
    print(f"{'='*80}")
    quarks = detect_quarks(lattice, cfg)
    up = sum(1 for q in quarks if q.quark_type == "UP")
    dn = sum(1 for q in quarks if q.quark_type == "DOWN")
    strict, loose = find_protons(quarks, cfg)
    gd = grade_distribution(lattice)
    od = orbit_distribution(lattice)

    print(f"FINAL (t={steps}): UP={up} DOWN={dn} "
          f"ratio={up/dn if dn > 0 else float('inf'):.3f}")
    print(f"  Strict protons: {strict}  Loose clusters: {loose}")
    print()
    print(f"  Grade distribution:")
    for g in range(5):
        bar = '█' * (gd[g] // 20)
        print(f"    Grade {g}: {gd[g]:5d}  {bar}")
    print()
    print(f"  Orbit distribution:")
    print(f"    Fixed points: {od[-1]:5d}  "
          f"(scalar={sum(1 for s in lattice if int(s)==1)}, "
          f"γ⁰={sum(1 for s in lattice if int(s)==2)}, "
          f"γ¹²³={sum(1 for s in lattice if int(s)==15)}, "
          f"γ⁰¹²³={sum(1 for s in lattice if int(s)==16)})")
    for oi, orbit in enumerate(Z3_ORBITS):
        members = sorted(orbit)
        names = ', '.join(cl_name(s) for s in members)
        count = od[oi]
        print(f"    Orbit {oi} ({names}): {count:5d}")

    # Quark orbit breakdown
    up_orbits = Counter(q.orbit for q in quarks if q.quark_type == "UP")
    dn_orbits = Counter(q.orbit for q in quarks if q.quark_type == "DOWN")
    print()
    print(f"  Quark distribution by orbit:")
    print(f"    UP:   {dict(sorted(up_orbits.items()))}")
    print(f"    DOWN: {dict(sorted(dn_orbits.items()))}")

    if strict > 0:
        print()
        print("  ████████████████████████████████████████████████████████████")
        print("  ██  PROTONS FROM Cl(1,3) HEXADECIMAL LATTICE DYNAMICS   ██")
        print("  ████████████████████████████████████████████████████████████")

    # ── Multi-seed sweep ────────────────────────────────────────────────
    print()
    print(f"{'='*80}")
    print("MULTI-SEED: 20 seeds, peak protons during percolation")
    print(f"{'='*80}")

    peak_list = []
    peak_ratios = []
    for s in range(20):
        r = np.random.default_rng(s * 137 + 7)
        lat = init_lattice(cfg)
        best_protons = 0
        best_t = 0
        best_up = 0
        best_dn = 0
        for t in range(1, steps + 1):
            lat = step(lat, r, cfg)
            qs = detect_quarks(lat, cfg)
            u = sum(1 for q in qs if q.quark_type == "UP")
            d = sum(1 for q in qs if q.quark_type == "DOWN")
            ps, _ = find_protons(qs, cfg)
            if ps > best_protons:
                best_protons = ps
                best_t = t
                best_up = u
                best_dn = d
        peak_list.append(best_protons)
        ratio = best_up / best_dn if best_dn > 0 else float('inf')
        peak_ratios.append(ratio)
        ratio_str = f'{ratio:.1f}' if ratio < 1e6 else '∞'
        print(f"  seed={s:2d}: peak p={best_protons:3d} at t={best_t:3d} "
              f"(UP={best_up} DN={best_dn} r={ratio_str})")

    finite_ratios = [r for r in peak_ratios if r < 1e6]
    mean_p = np.mean(peak_list)
    std_p = np.std(peak_list)
    mean_r = np.mean(finite_ratios) if finite_ratios else float('inf')
    n_with = sum(1 for p in peak_list if p > 0)

    print(f"\n  Protons: {mean_p:.1f} ± {std_p:.1f}")
    print(f"  UP/DOWN ratio: {mean_r:.2f}")
    print(f"  Seeds with protons: {n_with}/20")


if __name__ == "__main__":
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 42
    run_simulation(seed=seed)
