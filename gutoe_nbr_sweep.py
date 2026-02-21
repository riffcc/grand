#!/usr/bin/env python3
"""
GUTOE: Neighbour count sweep.

Tests whether the UP/DOWN ratio and proton count depend on
the number of neighbours, or if it's just Z₃ combinatorics.

Runs 6, 12, 16, and 20 neighbours with the same seeds.
Also runs a REAL square lattice control (no hex, no triangles).

Copyright (C) 2026 Riff Labs, AGPL-3.0-or-later
"""

import numpy as np
from enum import IntEnum
from dataclasses import dataclass
from collections import Counter
import sys

# ── TriState system ──────────────────────────────────────────────────────────

class TriState(IntEnum):
    VOID = 0
    COSINE = 1
    SINE = 2
    TANGENT = 3

def cycle(s: TriState) -> TriState:
    if s == TriState.SINE: return TriState.COSINE
    if s == TriState.COSINE: return TriState.TANGENT
    if s == TriState.TANGENT: return TriState.SINE
    return TriState.VOID

def veracity_val(a: TriState, b: TriState) -> float:
    if a == TriState.VOID or b == TriState.VOID:
        return 0.0
    if a == b:
        return 1.0
    pair = frozenset({a, b})
    if pair == frozenset({TriState.SINE, TriState.COSINE}):
        return np.sqrt(3) / 2
    return 0.5

# ── Lattice config ───────────────────────────────────────────────────────────

@dataclass
class LatticeConfig:
    hex_rows: int = 12
    hex_cols: int = 12
    layers: int = 12
    differentiation_prob: float = 0.02
    cycle_prob: float = 0.05
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

# ── Mesh topologies (parameterized) ──────────────────────────────────────────

def mesh_hex6(r, c, z, cfg):
    """6 neighbours: in-plane hex only. Has triangles."""
    return [idx(nr, nc, z, cfg) for nr, nc in hex_planar_neighbours(r, c, cfg)]

def mesh_hex12(r, c, z, cfg):
    """12 neighbours: HCP (6 in-plane + 3 above + 3 below). Has triangles."""
    planar = hex_planar_neighbours(r, c, cfg)
    nbrs = [idx(nr, nc, z, cfg) for nr, nc in planar]
    for i in (0, 2, 4):
        nr, nc = planar[i]
        nbrs.append(idx(nr, nc, (z + 1) % cfg.layers, cfg))
    for i in (1, 3, 5):
        nr, nc = planar[i]
        nbrs.append(idx(nr, nc, (z - 1) % cfg.layers, cfg))
    return nbrs

def mesh_hex16(r, c, z, cfg):
    """16 neighbours: HCP + axial + extended. Has triangles."""
    planar = hex_planar_neighbours(r, c, cfg)
    nbrs = [idx(nr, nc, z, cfg) for nr, nc in planar]
    for i in (0, 2, 4):
        nr, nc = planar[i]
        nbrs.append(idx(nr, nc, (z + 1) % cfg.layers, cfg))
    for i in (1, 3, 5):
        nr, nc = planar[i]
        nbrs.append(idx(nr, nc, (z - 1) % cfg.layers, cfg))
    nbrs.append(idx(r, c, (z + 1) % cfg.layers, cfg))
    nbrs.append(idx(r, c, (z - 1) % cfg.layers, cfg))
    nbrs.append(idx(r, c, (z + 2) % cfg.layers, cfg))
    nbrs.append(idx(r, c, (z - 2) % cfg.layers, cfg))
    return nbrs

def mesh_hex20(r, c, z, cfg):
    """20 neighbours: full hex prism (6+6+6+2). Has triangles."""
    planar = hex_planar_neighbours(r, c, cfg)
    nbrs = [idx(nr, nc, z, cfg) for nr, nc in planar]
    for nr, nc in planar:
        nbrs.append(idx(nr, nc, (z + 1) % cfg.layers, cfg))
    for nr, nc in planar:
        nbrs.append(idx(nr, nc, (z - 1) % cfg.layers, cfg))
    nbrs.append(idx(r, c, (z + 1) % cfg.layers, cfg))
    nbrs.append(idx(r, c, (z - 1) % cfg.layers, cfg))
    return nbrs

def mesh_square6(r, c, z, cfg):
    """6 neighbours: square lattice (4 in-plane + 2 axial). NO triangles."""
    nbrs = []
    for dr, dc in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
        nbrs.append(idx((r + dr) % cfg.hex_rows, (c + dc) % cfg.hex_cols, z, cfg))
    nbrs.append(idx(r, c, (z + 1) % cfg.layers, cfg))
    nbrs.append(idx(r, c, (z - 1) % cfg.layers, cfg))
    return nbrs

# ── Simulation engine (mesh function as parameter) ───────────────────────────

def compute_local_fields(lattice, site, r, c, z, cfg, mesh_fn):
    state = TriState(lattice[site])
    if state == TriState.VOID:
        return 0.0, 0.0, 0.0

    nbrs = mesh_fn(r, c, z, cfg)
    total_veracity = 0.0
    grad_sum = 0.0
    nbr_states = []

    for nbr_idx in nbrs:
        nbr_state = TriState(lattice[nbr_idx])
        v = veracity_val(state, nbr_state)
        total_veracity += v
        grad_sum += 1.0 - v
        if nbr_state != TriState.VOID:
            nbr_states.append(nbr_state)

    n_total = len(nbrs)
    avg_veracity = total_veracity / n_total
    field_grad = grad_sum / n_total
    n_distinct = len(set(nbr_states))
    curvature = (n_distinct - 1) / 2

    return avg_veracity, curvature, field_grad

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

def detect_quarks(lattice, cfg, mesh_fn):
    quarks = []
    N = cfg.hex_rows * cfg.hex_cols * cfg.layers
    for site in range(N):
        state = TriState(lattice[site])
        if state == TriState.VOID:
            continue
        z = site // (cfg.hex_rows * cfg.hex_cols)
        rem = site % (cfg.hex_rows * cfg.hex_cols)
        r = rem // cfg.hex_cols
        c = rem % cfg.hex_cols
        v, curv, grad = compute_local_fields(lattice, site, r, c, z, cfg, mesh_fn)
        bc = v / (1 + grad)
        if bc >= cfg.quark_threshold:
            qtype = "UP" if v > curv else "DOWN"
            quarks.append(Quark(site, r, c, z, qtype, bc, v, curv))
    return quarks

def find_protons(quarks, cfg, mesh_fn):
    quark_set = {q.site: q for q in quarks}
    strict_count = 0
    used_strict = set()
    nbr_cache = {}
    for q in quarks:
        nbr_cache[q.site] = set(mesh_fn(q.r, q.c, q.z, cfg))
    for q in quarks:
        if q.quark_type != "DOWN" or q.site in used_strict:
            continue
        up_nbrs = [quark_set[ni] for ni in nbr_cache[q.site]
                    if ni in quark_set and quark_set[ni].quark_type == "UP"
                    and ni not in used_strict]
        if len(up_nbrs) < 2:
            continue
        found = False
        for i in range(len(up_nbrs)):
            if found:
                break
            for j in range(i + 1, len(up_nbrs)):
                p1, p2 = up_nbrs[i].site, up_nbrs[j].site
                if p1 in nbr_cache.get(p2, set()):
                    strict_count += 1
                    used_strict.update([q.site, p1, p2])
                    found = True
                    break
    return strict_count

def step(lattice, rng, cfg, mesh_fn):
    new = lattice.copy()
    N = cfg.hex_rows * cfg.hex_cols * cfg.layers
    for site in range(N):
        z = site // (cfg.hex_rows * cfg.hex_cols)
        rem = site % (cfg.hex_rows * cfg.hex_cols)
        r = rem // cfg.hex_cols
        c = rem % cfg.hex_cols
        state = TriState(lattice[site])
        if state == TriState.VOID:
            if rng.random() < cfg.differentiation_prob:
                new[site] = int(TriState.SINE)
                continue
            nbrs = mesh_fn(r, c, z, cfg)
            active_nbrs = sum(1 for ni in nbrs if lattice[ni] != int(TriState.VOID))
            total_nbrs = len(nbrs)
            if active_nbrs >= max(2, total_nbrs // 4):
                spread_prob = active_nbrs / total_nbrs * 0.4
                if rng.random() < spread_prob:
                    new[site] = int(TriState.SINE)
        else:
            if rng.random() < cfg.cycle_prob:
                new[site] = int(cycle(state))
            elif rng.random() < cfg.alignment_strength:
                nbrs = mesh_fn(r, c, z, cfg)
                nbr_states = [TriState(lattice[ni]) for ni in nbrs
                              if lattice[ni] != int(TriState.VOID)]
                if nbr_states:
                    counts = Counter(nbr_states)
                    majority = counts.most_common(1)[0][0]
                    new[site] = int(majority)
        new[site] = new[site]
    return new

def init_lattice(cfg):
    return np.zeros(cfg.hex_rows * cfg.hex_cols * cfg.layers, dtype=np.int8)

# ── Run sweep ────────────────────────────────────────────────────────────────

def run_config(name, mesh_fn, cfg, n_seeds=20, steps=200):
    """Run multi-seed experiment for one mesh configuration."""
    peaks = []
    peak_ratios = []

    for s in range(n_seeds):
        rng = np.random.default_rng(s * 137 + 7)
        lat = init_lattice(cfg)
        best_protons = 0
        best_t = 0
        best_ratio = 0.0

        for t in range(1, steps + 1):
            lat = step(lat, rng, cfg, mesh_fn)
            qs = detect_quarks(lat, cfg, mesh_fn)
            u = sum(1 for q in qs if q.quark_type == "UP")
            d = sum(1 for q in qs if q.quark_type == "DOWN")
            ps = find_protons(qs, cfg, mesh_fn)
            if ps > best_protons:
                best_protons = ps
                best_t = t
                best_ratio = u / d if d > 0 else float('inf')

        peaks.append(best_protons)
        peak_ratios.append(best_ratio)

    mean_p = np.mean(peaks)
    std_p = np.std(peaks)
    mean_r = np.mean([r for r in peak_ratios if r < 1e6])
    n_with = sum(1 for p in peaks if p > 0)

    print(f"  {name:30s} | protons={mean_p:6.1f}±{std_p:4.1f} | "
          f"UP/DN={mean_r:5.2f} | {n_with}/20 seeds")
    return mean_p, mean_r


if __name__ == "__main__":
    cfg = LatticeConfig()
    N = cfg.hex_rows * cfg.hex_cols * cfg.layers

    print(f"GUTOE Neighbour Count Sweep")
    print(f"{'='*75}")
    print(f"Lattice: {cfg.hex_rows}×{cfg.hex_cols}×{cfg.layers} = {N} cells")
    print(f"Curvature: Z₃ state diversity (n_distinct states among nbrs)")
    print(f"20 seeds each, 200 steps, peak protons tracked")
    print(f"{'='*75}")
    print()

    configs = [
        ("hex-6 (planar only)",     mesh_hex6),
        ("hex-12 (HCP)",            mesh_hex12),
        ("hex-16 (hexadecimal)",    mesh_hex16),
        ("hex-20 (full prism)",     mesh_hex20),
        ("square-6 (NO triangles)", mesh_square6),
    ]

    results = {}
    for name, fn in configs:
        p, r = run_config(name, fn, cfg)
        results[name] = (p, r)

    print()
    print(f"{'='*75}")
    print("ANALYSIS")
    print(f"{'='*75}")

    hex_results = [(n, p, r) for n, (p, r) in results.items() if 'hex' in n]
    sq_result = results.get("square-6 (NO triangles)", (0, 0))

    print()
    print("If proton count is ~constant across hex-6/12/16/20:")
    print("  → Z₃ diversity is combinatorial, not tied to 16")
    print()
    print("If proton count differs significantly:")
    print("  → Neighbour count matters (look for 16 as special)")
    print()
    print("If square-6 matches hex-6:")
    print("  → Triangles don't matter (topology is decorative)")
    print()
    print("If square-6 << hex-6:")
    print("  → Triangles are essential (topology is load-bearing)")
