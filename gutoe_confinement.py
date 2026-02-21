#!/usr/bin/env python3
"""
GUTOE: Confinement experiment — curvature as confining potential.

Hypothesis: alignment kills protons by homogenizing triple junctions.
Fix: alignment strength WEAKENS at high curvature.

  alignment_eff = alignment_strength × (1 - curvature)

At triple junction (curvature=1.0): alignment=0  → junction is STABLE
At domain boundary (curvature=0.5): alignment=0.075 → weak
At domain interior (curvature=0.0): alignment=0.15  → full

This makes curvature the confining potential: high-curvature regions
(where quarks live) resist the homogenization that destroys them.

Runs four configurations on identical seeds:
  A) BASELINE:    flat alignment (original dynamics)
  B) CURVATURE:   alignment × (1 - curvature)
  C) CURV+REDIR:  curvature confinement + grade-0→γ¹²³ redirect
  D) CURV+STABLE: curvature confinement + redirect + grade-3 stability

Copyright (C) 2026 Riff Labs, AGPL-3.0-or-later
"""

import numpy as np
from dataclasses import dataclass
from collections import Counter

# ── Cl(1,3) machinery ──────────────────────────────────────────────────────

VOID = 0
_GRADE_TABLE = [-1] + [bin(mi).count('1') for mi in range(16)]

def _make_z3_table():
    table = [VOID]
    for s in range(1, 17):
        mi = s - 1
        b0 = (mi >> 0) & 1
        b1 = (mi >> 1) & 1
        b2 = (mi >> 2) & 1
        b3 = (mi >> 3) & 1
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

Z3_ORBITS = [
    frozenset({3, 5, 9}),
    frozenset({4, 6, 10}),
    frozenset({7, 13, 11}),
    frozenset({8, 14, 12}),
]

_ORBIT_TABLE = {}
for oi, orbit in enumerate(Z3_ORBITS):
    for s in orbit:
        _ORBIT_TABLE[s] = oi
for s in [1, 2, 15, 16]:
    _ORBIT_TABLE[s] = -1

SEED_STATE = 3
BARYON_STATE = 15


# ── Lattice ─────────────────────────────────────────────────────────────────

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
    return [idx(nr, nc, z, cfg) for nr, nc in hex_planar_neighbours(r, c, cfg)]

def init_lattice(cfg):
    return np.zeros(cfg.hex_rows * cfg.hex_cols * cfg.layers, dtype=np.int8)


# ── Field computation ───────────────────────────────────────────────────────

def compute_local_fields(lattice, site, r, c, z, cfg):
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

    z3_curvature = 0.0
    for orbit in Z3_ORBITS:
        n_present = len(orbit & nbr_state_set)
        orbit_curv = (n_present - 1) / 2
        if orbit_curv > z3_curvature:
            z3_curvature = orbit_curv

    return avg_veracity, z3_curvature, field_grad


# ── Dynamics ────────────────────────────────────────────────────────────────

def step(lattice, rng, cfg, *,
         curvature_confinement=False,
         redirect_to_baryon=False,
         confine_grade3=False):
    new = lattice.copy()
    N = cfg.hex_rows * cfg.hex_cols * cfg.layers

    for site in range(N):
        z = site // (cfg.hex_rows * cfg.hex_cols)
        rem = site % (cfg.hex_rows * cfg.hex_cols)
        r = rem // cfg.hex_cols
        c = rem % cfg.hex_cols
        state = int(lattice[site])

        if state == VOID:
            if rng.random() < cfg.differentiation_prob:
                new[site] = SEED_STATE
                continue
            nbrs = mesh_neighbours(r, c, z, cfg)
            active = sum(1 for ni in nbrs if lattice[ni] != VOID)
            total = len(nbrs)
            if active >= max(2, total // 4):
                if rng.random() < active / total * 0.4:
                    new[site] = SEED_STATE
        else:
            grade = _GRADE_TABLE[state]

            if confine_grade3 and grade == 3:
                if rng.random() < cfg.cycle_prob:
                    new[site] = _Z3_TABLE[state]
                continue

            # Compute curvature for this cell (needed for curvature confinement)
            local_curv = 0.0
            if curvature_confinement:
                _, local_curv, _ = compute_local_fields(lattice, site, r, c, z, cfg)

            r_val = rng.random()
            if r_val < cfg.cycle_prob:
                new[site] = _Z3_TABLE[state]

            elif r_val < cfg.cycle_prob + cfg.clifford_prob:
                nbrs = mesh_neighbours(r, c, z, cfg)
                active_nbrs = [int(lattice[ni]) for ni in nbrs
                               if lattice[ni] != VOID]
                if active_nbrs:
                    partner = active_nbrs[rng.integers(len(active_nbrs))]
                    result_mi = (state - 1) ^ (partner - 1)
                    result = result_mi + 1
                    if redirect_to_baryon and _GRADE_TABLE[result] == 0:
                        result = BARYON_STATE
                    new[site] = result

            elif r_val < cfg.cycle_prob + cfg.clifford_prob + cfg.alignment_strength:
                # THE KEY CHANGE: curvature-dependent alignment
                # alignment_eff = alignment_strength × (1 - curvature)
                # High curvature → no alignment → triple junctions are STABLE
                if curvature_confinement:
                    alignment_eff = cfg.alignment_strength * (1 - local_curv)
                    # Re-check: does the random value fall within the reduced window?
                    threshold = cfg.cycle_prob + cfg.clifford_prob + alignment_eff
                    if r_val >= threshold:
                        continue  # curvature blocked alignment

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
            quarks.append(Quark(site, r, c, z, qtype, bc, v, curv))
    return quarks

def find_protons(quarks, cfg):
    quark_set = {q.site: q for q in quarks}
    strict_count = 0
    used = set()
    nbr_cache = {}
    for q in quarks:
        nbr_cache[q.site] = set(mesh_neighbours(q.r, q.c, q.z, cfg))
    for q in quarks:
        if q.quark_type != "DOWN" or q.site in used:
            continue
        up_nbrs = [quark_set[ni] for ni in nbr_cache[q.site]
                    if ni in quark_set and quark_set[ni].quark_type == "UP"
                    and ni not in used]
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
                    used.update([q.site, p1, p2])
                    found = True
                    break
    return strict_count


# ── Experiment ──────────────────────────────────────────────────────────────

def run_multiseed(name, cfg, n_seeds, steps, **kwargs):
    all_protons = np.zeros(steps)
    all_g0 = np.zeros(steps)
    all_g1 = np.zeros(steps)
    all_g3 = np.zeros(steps)
    all_up = np.zeros(steps)
    all_dn = np.zeros(steps)

    for s in range(n_seeds):
        seed = s * 137 + 7
        rng = np.random.default_rng(seed)
        lattice = init_lattice(cfg)

        for t in range(steps):
            lattice = step(lattice, rng, cfg, **kwargs)
            quarks = detect_quarks(lattice, cfg)
            u = sum(1 for q in quarks if q.quark_type == "UP")
            d = sum(1 for q in quarks if q.quark_type == "DOWN")
            p = find_protons(quarks, cfg)
            all_protons[t] += p
            all_g0[t] += sum(1 for x in lattice if _GRADE_TABLE[int(x)] == 0)
            all_g1[t] += sum(1 for x in lattice if _GRADE_TABLE[int(x)] == 1)
            all_g3[t] += sum(1 for x in lattice if _GRADE_TABLE[int(x)] == 3)
            all_up[t] += u
            all_dn[t] += d

    return {k: v / n_seeds for k, v in {
        'protons': all_protons, 'grade0': all_g0, 'grade1': all_g1,
        'grade3': all_g3, 'up': all_up, 'down': all_dn,
    }.items()}


if __name__ == "__main__":
    cfg = LatticeConfig()
    N = cfg.hex_rows * cfg.hex_cols * cfg.layers
    n_seeds = 10
    steps = 300

    print(f"GUTOE Confinement Experiment: Curvature as Confining Potential")
    print(f"{'='*80}")
    print(f"Lattice: {cfg.hex_rows}×{cfg.hex_cols}×{cfg.layers} = {N} cells")
    print(f"Key change: alignment_eff = alignment × (1 - curvature)")
    print(f"  Triple junction (curv=1.0): alignment = 0   → STABLE")
    print(f"  Domain boundary (curv=0.5): alignment = 0.075")
    print(f"  Domain interior (curv=0.0): alignment = 0.15 → full")
    print(f"Seeds: {n_seeds},  Steps: {steps}")
    print(f"{'='*80}")
    print()

    configs = [
        ("A) BASELINE",
         dict()),
        ("B) CURV-CONFINE",
         dict(curvature_confinement=True)),
        ("C) CURV+REDIR",
         dict(curvature_confinement=True, redirect_to_baryon=True)),
        ("D) CURV+REDIR+G3",
         dict(curvature_confinement=True, redirect_to_baryon=True, confine_grade3=True)),
    ]

    results = {}
    for name, kwargs in configs:
        print(f"Running {name}...")
        results[name] = run_multiseed(name, cfg, n_seeds, steps, **kwargs)
        print(f"  done.")

    # ── Time series ─────────────────────────────────────────────────────
    print()
    print(f"{'='*80}")
    print("PROTON TIME SERIES (averaged over seeds)")
    print(f"{'='*80}")
    print()

    header = f"{'t':>4s}"
    for name, _ in configs:
        short = name.split(')')[1].strip()[:10]
        header += f" | {short:>10s}  g0   g1   g3  U/D"
    print(header)
    print('-' * len(header))

    report_times = ([10, 15, 20, 25, 30, 35, 40, 50, 60, 75, 100,
                     125, 150, 200, 250, 300])
    for t in report_times:
        if t > steps:
            break
        line = f"{t:4d}"
        for name, _ in configs:
            r = results[name]
            i = t - 1
            p = r['protons'][i]
            g0 = r['grade0'][i]
            g1 = r['grade1'][i]
            g3 = r['grade3'][i]
            u = r['up'][i]
            d = r['down'][i]
            ratio = f"{u/d:.1f}" if d > 0.5 else "∞"
            line += f" | {p:10.0f} {g0:4.0f} {g1:4.0f} {g3:4.0f} {ratio:>5s}"
        print(line)

    # ── Summary ─────────────────────────────────────────────────────────
    print()
    print(f"{'='*80}")
    print("SUMMARY")
    print(f"{'='*80}")
    print()

    for name, _ in configs:
        r = results[name]
        peak_p = max(r['protons'])
        peak_t = np.argmax(r['protons']) + 1
        p100 = r['protons'][min(99, steps-1)]
        p200 = r['protons'][min(199, steps-1)] if steps >= 200 else 0
        p300 = r['protons'][min(299, steps-1)] if steps >= 300 else 0

        half = peak_p / 2
        half_t = None
        for i in range(int(peak_t), steps):
            if r['protons'][i] < half:
                half_t = i + 1
                break

        print(f"  {name}:")
        print(f"    Peak: {peak_p:.0f} protons at t={peak_t}")
        print(f"    t=100: {p100:.0f}  t=200: {p200:.0f}  t=300: {p300:.0f}")
        print(f"    Half-life: {f't={half_t}' if half_t else f'>t={steps}'}")
        print()

    # ── Verdict ─────────────────────────────────────────────────────────
    a = results[configs[0][0]]
    b = results[configs[1][0]]

    a_half = max(a['protons']) / 2
    b_half = max(b['protons']) / 2
    a_hl = steps
    b_hl = steps
    for i in range(int(np.argmax(a['protons'])), steps):
        if a['protons'][i] < a_half:
            a_hl = i + 1
            break
    for i in range(int(np.argmax(b['protons'])), steps):
        if b['protons'][i] < b_half:
            b_hl = i + 1
            break

    print(f"{'='*80}")
    print("VERDICT")
    print(f"{'='*80}")
    b_p300 = b['protons'][min(299, steps-1)]
    a_p300 = a['protons'][min(299, steps-1)]

    if b_hl > a_hl * 2:
        print(f"  CURVATURE CONFINEMENT WORKS")
        print(f"  Baseline half-life:  t={a_hl}")
        print(f"  Confined half-life:  t={b_hl}")
        print(f"  Lifetime extension:  {b_hl/a_hl:.1f}×")
        if b_p300 > 10:
            print(f"  Protons at t=300:    {b_p300:.0f} (PERSISTENT)")
        print()
        print(f"  Curvature IS the confining potential.")
        print(f"  High curvature → strong coupling → resists alignment.")
        print(f"  Triple junctions are stable attractors, not transient features.")
    elif b_hl > a_hl * 1.3:
        print(f"  PARTIAL CONFINEMENT: lifetime extended but not stable")
        print(f"  Baseline half-life: t={a_hl}")
        print(f"  Confined half-life: t={b_hl}")
    else:
        print(f"  NO CONFINEMENT: curvature-dependent alignment doesn't help")
        print(f"  Baseline half-life: t={a_hl}")
        print(f"  Confined half-life: t={b_hl}")
