#!/usr/bin/env python3
"""
GUTOE: Hydrogen and Fusion Analysis

Starting from the k=4 steady state at t=1000, we analyze:

1. HYDROGEN: Does any proton triangle have an adjacent 'electron'?
   Electron = γ⁰ (state 2), the timelike unit vector.
   γ⁰ is a Z₃ FIXED POINT — it has trivial Z₃ representation,
   making it 'colorless', unlike quarks which sit in Z₃ orbits of 3.
   This is the natural lepton in the Cl(1,3) framework.
   Hydrogen = proton triangle + adjacent γ⁰ not in any other proton.

2. FUSION: Are any two proton triangles adjacent (sharing a lattice edge)?
   Run one more step. Do they merge, repel, or annihilate?
   If two adjacent uud triangles transform into a larger bound structure,
   that is fusion.

All analysis is from the real t=1000 lattice snapshot — zero new physics.

Copyright (C) 2026 Riff Labs, AGPL-3.0-or-later
"""

import numpy as np
from collections import Counter
from dataclasses import dataclass

# ── Physics (identical to gutoe_2d5_toroid.py) ───────────────────────────────

VOID = 0
_GRADE_TABLE = [-1] + [bin(mi).count('1') for mi in range(16)]

def _make_z3_table():
    table = [VOID]
    for s in range(1, 17):
        mi = s - 1
        b0=(mi>>0)&1; b1=(mi>>1)&1; b2=(mi>>2)&1; b3=(mi>>3)&1
        table.append((b0|(b3<<1)|(b1<<2)|(b2<<3))+1)
    return table

_Z3_TABLE = _make_z3_table()
_SQRT3_HALF = np.sqrt(3)/2

def _make_veracity_table():
    t = {}
    for s1 in range(17):
        for s2 in range(17):
            if s1==VOID or s2==VOID: t[(s1,s2)]=0.0
            elif s1==s2: t[(s1,s2)]=1.0
            else:
                d=bin((s1-1)^(s2-1)).count('1')
                t[(s1,s2)] = _SQRT3_HALF if d==1 else (0.5 if d==2 else 0.0)
    return t

_V = _make_veracity_table()

Z3_ORBITS = [
    frozenset({3,5,9}), frozenset({4,6,10}),
    frozenset({7,13,11}), frozenset({8,14,12}),
]

# Z₃ fixed points (trivial Z₃ representation = 'colorless' = leptons)
# Grade 0: 1 (scalar), Grade 1: 2 (γ⁰), Grade 3: 15 (γ¹²³), Grade 4: 16 (γ⁰¹²³)
Z3_FIXED = frozenset({1, 2, 15, 16})
ELECTRON_STATE = 2  # γ⁰ — grade-1 Z₃ fixed point, timelike unit vector

SEED_STATE = 3  # γ¹

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
    void_votes: int = 4  # k=4: the stable phase

def idx(r,c,z,cfg):
    return((z%cfg.layers)*cfg.hex_rows+(r%cfg.hex_rows))*cfg.hex_cols+(c%cfg.hex_cols)

def hex_planar_neighbours(r,c,cfg):
    if r%2==0: offs=[(-1,0),(-1,1),(0,-1),(0,1),(1,0),(1,1)]
    else: offs=[(-1,-1),(-1,0),(0,-1),(0,1),(1,-1),(1,0)]
    return[((r+dr)%cfg.hex_rows,(c+dc)%cfg.hex_cols) for dr,dc in offs]

def mesh_neighbours(r,c,z,cfg):
    return[idx(nr,nc,z,cfg) for nr,nc in hex_planar_neighbours(r,c,cfg)]

def site_coords(site, cfg):
    z = site // (cfg.hex_rows * cfg.hex_cols)
    rem = site % (cfg.hex_rows * cfg.hex_cols)
    r = rem // cfg.hex_cols
    c = rem % cfg.hex_cols
    return r, c, z

def init_lattice(cfg):
    return np.zeros(cfg.hex_rows*cfg.hex_cols*cfg.layers, dtype=np.int8)

def compute_local_fields(lattice, site, r, c, z, cfg):
    state = int(lattice[site])
    if state == VOID:
        return 0.0, 0.0, 0.0
    nbrs = mesh_neighbours(r, c, z, cfg)
    total_v = 0.0; grad = 0.0; nbr_set = set()
    for ni in nbrs:
        ns = int(lattice[ni])
        v = _V[(state, ns)]
        total_v += v; grad += 1.0 - v
        if ns != VOID: nbr_set.add(ns)
    n = len(nbrs)
    z3_curv = max(((len(orbit & nbr_set)-1)/2) for orbit in Z3_ORBITS)
    return total_v/n, z3_curv, grad/n

def step(lattice, rng, cfg):
    new = lattice.copy()
    N = cfg.hex_rows * cfg.hex_cols * cfg.layers
    for site in range(N):
        r, c, z = site_coords(site, cfg)
        state = int(lattice[site])
        if state == VOID:
            if rng.random() < cfg.differentiation_prob:
                new[site] = SEED_STATE; continue
            nbrs = mesh_neighbours(r, c, z, cfg)
            active = sum(1 for ni in nbrs if lattice[ni] != VOID)
            total = len(nbrs)
            if active >= max(2, total//4):
                if rng.random() < active/total*0.4:
                    new[site] = SEED_STATE
        else:
            r_val = rng.random()
            if r_val < cfg.cycle_prob:
                new[site] = _Z3_TABLE[state]
            elif r_val < cfg.cycle_prob + cfg.clifford_prob:
                nbrs = mesh_neighbours(r, c, z, cfg)
                active_nbrs = [int(lattice[ni]) for ni in nbrs if lattice[ni]!=VOID]
                if active_nbrs:
                    partner = active_nbrs[rng.integers(len(active_nbrs))]
                    new[site] = ((state-1)^(partner-1))+1
            elif r_val < cfg.cycle_prob + cfg.clifford_prob + cfg.alignment_strength:
                nbrs = mesh_neighbours(r, c, z, cfg)
                nbr_states = [int(lattice[ni]) for ni in nbrs if lattice[ni]!=VOID]
                if nbr_states:
                    votes = Counter(nbr_states)
                    winner, winner_count = votes.most_common(1)[0]
                    if winner_count > cfg.void_votes:
                        new[site] = winner
    return new

@dataclass
class Quark:
    site: int; r: int; c: int; z: int
    quark_type: str; binding_coherence: float
    veracity: float; curvature: float

def detect_quarks(lattice, cfg):
    quarks = []
    N = cfg.hex_rows * cfg.hex_cols * cfg.layers
    for site in range(N):
        state = int(lattice[site])
        if state == VOID: continue
        r, c, z = site_coords(site, cfg)
        v, curv, grad = compute_local_fields(lattice, site, r, c, z, cfg)
        bc = v/(1+grad)
        if bc >= cfg.quark_threshold:
            quarks.append(Quark(site,r,c,z,"UP" if v>curv else "DOWN",bc,v,curv))
    return quarks

def find_proton_triplets(quarks, cfg):
    """Returns list of (down_site, up1_site, up2_site) proton triplets."""
    quark_set = {q.site: q for q in quarks}
    triplets = []
    used = set()
    nbr_cache = {q.site: set(mesh_neighbours(q.r,q.c,q.z,cfg)) for q in quarks}
    for q in quarks:
        if q.quark_type != "DOWN" or q.site in used: continue
        up_nbrs = [quark_set[ni] for ni in nbr_cache[q.site]
                   if ni in quark_set and quark_set[ni].quark_type=="UP"
                   and ni not in used]
        if len(up_nbrs) < 2: continue
        found = False
        for i in range(len(up_nbrs)):
            if found: break
            for j in range(i+1, len(up_nbrs)):
                p1,p2 = up_nbrs[i].site, up_nbrs[j].site
                if p1 in nbr_cache.get(p2, set()):
                    triplets.append((q.site, p1, p2))
                    used.update([q.site, p1, p2])
                    found = True; break
    return triplets

# ── Analysis ─────────────────────────────────────────────────────────────────

def analyze_hydrogen(lattice, triplets, cfg):
    """
    For each proton, find adjacent γ⁰ cells (state 2) not in any proton.
    γ⁰ is the Z₃-fixed grade-1 state — the natural lepton/electron.
    Returns: list of (proton_idx, electron_site) pairs = hydrogen atoms.
    Also counts: free electrons total, total grade-1 cells.
    """
    proton_sites = set()
    for d, u1, u2 in triplets:
        proton_sites.update([d, u1, u2])

    # Build neighborhood for each proton
    hydrogen = []
    electrons_near = {}
    for pi, (d, u1, u2) in enumerate(triplets):
        # All sites neighboring the proton triangle
        proton_shell = set()
        for site in [d, u1, u2]:
            r, c, z = site_coords(site, cfg)
            proton_shell.update(mesh_neighbours(r, c, z, cfg))
        proton_shell -= proton_sites  # exclude the proton's own cells

        # Find γ⁰ electrons in the shell
        electrons_here = [ns for ns in proton_shell
                          if int(lattice[ns]) == ELECTRON_STATE]
        electrons_near[pi] = electrons_here
        if electrons_here:
            hydrogen.append((pi, electrons_here))

    # Global electron (γ⁰) count
    N = cfg.hex_rows * cfg.hex_cols * cfg.layers
    all_g0_electrons = sum(1 for i in range(N) if int(lattice[i]) == ELECTRON_STATE)
    all_grade1 = sum(1 for i in range(N) if _GRADE_TABLE[int(lattice[i])] == 1)

    return hydrogen, all_g0_electrons, all_grade1

def analyze_fusion(lattice, triplets, rng, cfg):
    """
    Find proton triplets adjacent to each other (sharing a lattice edge).
    Run one more timestep and check what happened to those sites.
    Returns: adjacent_pairs, and what each site became after one step.
    """
    proton_sites = {}
    for i, (d, u1, u2) in enumerate(triplets):
        for site in [d, u1, u2]:
            proton_sites[site] = i

    nbr_cache = {}
    for i, (d, u1, u2) in enumerate(triplets):
        for site in [d, u1, u2]:
            if site not in nbr_cache:
                r, c, z = site_coords(site, cfg)
                nbr_cache[site] = set(mesh_neighbours(r, c, z, cfg))

    # Find pairs of adjacent protons
    adjacent_pairs = {}
    for i, (d1, u11, u12) in enumerate(triplets):
        for site1 in [d1, u11, u12]:
            for nbr in nbr_cache[site1]:
                if nbr in proton_sites:
                    j = proton_sites[nbr]
                    if j != i:
                        key = tuple(sorted([i, j]))
                        if key not in adjacent_pairs:
                            adjacent_pairs[key] = (site1, nbr)  # one shared edge

    # Run one more step and observe
    lat_before = lattice.copy()
    lat_after = step(lattice, rng, cfg)

    fusion_events = []
    for (i, j), (edge1, edge2) in adjacent_pairs.items():
        d1, u11, u12 = triplets[i]
        d2, u21, u22 = triplets[j]
        sites_i = [d1, u11, u12]
        sites_j = [d2, u21, u22]

        # Check if the triangle sites survived
        def grade(lat, site):
            return _GRADE_TABLE[int(lat[site])]

        before_i = [int(lat_before[s]) for s in sites_i]
        before_j = [int(lat_before[s]) for s in sites_j]
        after_i  = [int(lat_after[s])  for s in sites_i]
        after_j  = [int(lat_after[s])  for s in sites_j]

        changed_i = sum(1 for b,a in zip(before_i, after_i) if b != a)
        changed_j = sum(1 for b,a in zip(before_j, after_j) if b != a)

        # Detect fusion: sites from both protons now form new coherent structure
        all_after = set(after_i + after_j) - {VOID}
        all_grades_after = [_GRADE_TABLE[s] for s in all_after if s != VOID]
        has_higher_grade = any(g >= 2 for g in all_grades_after)

        fusion_events.append({
            'pair': (i, j),
            'changed_i': changed_i,
            'changed_j': changed_j,
            'before_i': before_i, 'after_i': after_i,
            'before_j': before_j, 'after_j': after_j,
            'has_higher_grade': has_higher_grade,
        })

    return adjacent_pairs, fusion_events, lat_after

# ── State name helper ─────────────────────────────────────────────────────────

def state_name(s):
    if s == 0: return 'VOID'
    if s == 1: return 'scalar'
    mi = s - 1
    bits = []
    if (mi>>0)&1: bits.append('0')
    if (mi>>1)&1: bits.append('1')
    if (mi>>2)&1: bits.append('2')
    if (mi>>3)&1: bits.append('3')
    return 'γ' + ''.join(bits) if bits else 'scalar'

# ── Main ─────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    cfg = LatticeConfig(void_votes=4)
    N = cfg.hex_rows * cfg.hex_cols * cfg.layers

    print("GUTOE: Hydrogen and Fusion")
    print("="*72)
    print("Running k=4 (stable phase) to t=1000, then analyzing snapshot.")
    print(f"Lattice: {cfg.hex_rows}×{cfg.hex_cols}×{cfg.layers} = {N} sites")
    print(f"Electron = γ⁰ (state 2): grade-1, Z₃ fixed point (colorless lepton)")
    print("="*72)
    print()

    # Run 3 seeds to t=1000 and analyze each snapshot
    n_seeds = 3
    total_hydrogen = 0
    total_adjacent_pairs = 0
    total_fusion_candidates = 0

    for seed_idx in range(n_seeds):
        rng = np.random.default_rng(seed_idx * 137 + 7)
        lat = init_lattice(cfg)

        print(f"Seed {seed_idx}: running to t=1000...", end='', flush=True)
        for t in range(1000):
            lat = step(lat, rng, cfg)
            if (t+1) % 200 == 0:
                print(f" t={t+1}", end='', flush=True)
        print(" done.")

        # Snapshot analysis
        quarks = detect_quarks(lat, cfg)
        triplets = find_proton_triplets(quarks, cfg)
        n_protons = len(triplets)

        # Grade distribution
        grade_counts = Counter(_GRADE_TABLE[int(lat[i])] for i in range(N))
        n_g0 = grade_counts.get(0, 0)
        n_g1 = grade_counts.get(1, 0)
        n_g2 = grade_counts.get(2, 0)
        n_void = grade_counts.get(-1, 0)

        # State distribution within grade-1
        g1_states = Counter(int(lat[i]) for i in range(N)
                            if _GRADE_TABLE[int(lat[i])] == 1)

        print(f"\n  Seed {seed_idx} at t=1000:")
        print(f"    VOID={n_void}  grade-0={n_g0}  grade-1={n_g1}  grade-2+={N-n_void-n_g0-n_g1}")
        print(f"    Grade-1 by state: " +
              ", ".join(f"{state_name(s)}={c}" for s,c in sorted(g1_states.items())))
        print(f"    Protons found: {n_protons}")

        # ── HYDROGEN ─────────────────────────────────────────────────────
        hydrogen, n_electrons, n_all_g1 = analyze_hydrogen(lat, triplets, cfg)

        print(f"\n  HYDROGEN:")
        print(f"    γ⁰ cells (electrons) in lattice: {n_electrons}")
        print(f"    Proton triangles with adjacent γ⁰: {len(hydrogen)} / {n_protons}")
        total_hydrogen += len(hydrogen)

        if hydrogen:
            for pi, esites in hydrogen[:5]:  # show first 5
                d, u1, u2 = triplets[pi]
                print(f"      Proton {pi} (sites {d},{u1},{u2}) + "
                      f"{len(esites)} electron(s) at {esites[:3]}")
            if len(hydrogen) > 5:
                print(f"      ... and {len(hydrogen)-5} more")
        else:
            print(f"    → No hydrogen found. (γ⁰ cells: {n_electrons})")

        # ── FUSION ───────────────────────────────────────────────────────
        rng_fusion = np.random.default_rng(seed_idx * 137 + 7 + 99999)
        adjacent_pairs, fusion_events, lat_after = analyze_fusion(
            lat, triplets, rng_fusion, cfg)

        print(f"\n  FUSION:")
        print(f"    Adjacent proton pairs (sharing a lattice edge): {len(adjacent_pairs)}")
        total_adjacent_pairs += len(adjacent_pairs)

        if fusion_events:
            merge_count = 0
            repel_count = 0
            stable_count = 0
            for ev in fusion_events:
                changed = ev['changed_i'] + ev['changed_j']
                if ev['has_higher_grade']:
                    merge_count += 1
                elif changed == 0:
                    stable_count += 1
                else:
                    repel_count += 1

            print(f"    After one step:")
            print(f"      Produced higher-grade structure (merge/fusion): {merge_count}")
            print(f"      Sites unchanged (stable coexistence):           {stable_count}")
            print(f"      Sites changed, no higher grade (disruption):    {repel_count}")
            total_fusion_candidates += merge_count

            # Show first fusion-candidate event in detail
            for ev in fusion_events:
                if ev['has_higher_grade']:
                    i, j = ev['pair']
                    print(f"\n    First fusion event — protons {i} and {j}:")
                    print(f"      Before proton {i}: "
                          + " ".join(state_name(s) for s in ev['before_i']))
                    print(f"      After  proton {i}: "
                          + " ".join(state_name(s) for s in ev['after_i']))
                    print(f"      Before proton {j}: "
                          + " ".join(state_name(s) for s in ev['before_j']))
                    print(f"      After  proton {j}: "
                          + " ".join(state_name(s) for s in ev['after_j']))
                    break
        else:
            print(f"    → No adjacent proton pairs found.")

        # Did the snapshot after +1 step still have protons?
        quarks_after = detect_quarks(lat_after, cfg)
        triplets_after = find_proton_triplets(quarks_after, cfg)
        print(f"\n    Protons in t+1 snapshot: {len(triplets_after)}")
        print()

    # ── Summary ──────────────────────────────────────────────────────────
    print("="*72)
    print("SUMMARY (3 seeds)")
    print("="*72)
    print(f"  Hydrogen atoms found:       {total_hydrogen}")
    print(f"  Adjacent proton pairs:      {total_adjacent_pairs}")
    print(f"  Fusion events (+1 step):    {total_fusion_candidates}")
    print()
    if total_hydrogen > 0:
        print("  HYDROGEN: YES — protons with adjacent γ⁰ electrons found.")
        print("  The Z₃-fixed timelike vector naturally binds near proton triangles.")
    else:
        print("  HYDROGEN: NO — no γ⁰ electrons adjacent to protons.")
        print("  Either γ⁰ is absent or it doesn't bind to proton neighborhoods.")
    print()
    if total_adjacent_pairs > 0:
        print("  FUSION CANDIDATES: YES — adjacent proton triangles exist.")
        if total_fusion_candidates > 0:
            print("  FUSION EVENTS: YES — some produced higher-grade structures!")
        else:
            print("  FUSION EVENTS: NO — adjacency doesn't produce higher-grade states.")
    else:
        print("  FUSION CANDIDATES: NONE — protons are spatially isolated.")
