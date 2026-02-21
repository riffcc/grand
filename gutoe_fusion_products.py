#!/usr/bin/env python3
"""
GUTOE: Fusion Product Identification

At t=25 (dense phase), adjacent proton triangles fuse at ~15% rate.
This script captures the exact product states of every merge event.

Grade-2 states in Cl(1,3):
  Spacelike bivectors (Z₃ orbit {7,11,13}):
    γ¹²  = state  7  (mi=0110) — SO(3) rotation generator
    γ¹³  = state 11  (mi=1010) — SO(3) rotation generator
    γ²³  = state 13  (mi=1100) — SO(3) rotation generator
  Mixed bivectors (Z₃ orbit {4,6,10}):
    γ⁰¹  = state  4  (mi=0011) — Lorentz boost generator
    γ⁰²  = state  6  (mi=0101) — Lorentz boost generator
    γ⁰³  = state 10  (mi=1001) — Lorentz boost generator

Grade-3 states (trivectors):
    γ⁰¹² = state  8  (mi=0111)
    γ⁰¹³ = state 12  (mi=1011)
    γ⁰²³ = state 14  (mi=1101)
    γ¹²³ = state 15  (mi=1110) — Z₃ fixed point, baryon number

If the product is consistently one specific bivector: it's a named particle.
p + p → (grade-2 product) is the fusion reaction.

Copyright (C) 2026 Riff Labs, AGPL-3.0-or-later
"""

import numpy as np
from collections import Counter
from dataclasses import dataclass

VOID = 0
QUARK_SEED = 3

_GRADE_TABLE = [-1] + [bin(mi).count('1') for mi in range(16)]

STATE_NAMES = {
    0:  'VOID',
    1:  'scalar',
    2:  'γ⁰',
    3:  'γ¹',
    4:  'γ⁰¹',
    5:  'γ²',
    6:  'γ⁰²',
    7:  'γ¹²',
    8:  'γ⁰¹²',
    9:  'γ³',
    10: 'γ⁰³',
    11: 'γ¹³',
    12: 'γ⁰¹³',
    13: 'γ²³',
    14: 'γ⁰²³',
    15: 'γ¹²³',
    16: 'γ⁰¹²³',
}

GRADE2_NAMES = {
    4:  'γ⁰¹ (boost)',
    6:  'γ⁰² (boost)',
    7:  'γ¹²  (rotation)',
    10: 'γ⁰³ (boost)',
    11: 'γ¹³  (rotation)',
    13: 'γ²³  (rotation)',
}

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
    void_votes: int = 4

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
    return rem // cfg.hex_cols, rem % cfg.hex_cols, z

def init_lattice(cfg):
    return np.zeros(cfg.hex_rows*cfg.hex_cols*cfg.layers, dtype=np.int8)

def compute_local_fields(lattice, site, r, c, z, cfg):
    state = int(lattice[site])
    if state == VOID: return 0.0, 0.0, 0.0
    nbrs = mesh_neighbours(r, c, z, cfg)
    total_v=0.0; grad=0.0; nbr_set=set()
    for ni in nbrs:
        ns = int(lattice[ni])
        v = _V[(state, ns)]
        total_v += v; grad += 1.0-v
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
                new[site] = QUARK_SEED; continue
            nbrs = mesh_neighbours(r, c, z, cfg)
            active = sum(1 for ni in nbrs if lattice[ni] != VOID)
            total = len(nbrs)
            if active >= max(2, total//4):
                if rng.random() < active/total*0.4:
                    new[site] = QUARK_SEED
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
    site: int; r: int; c: int; z: int; quark_type: str

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
            quarks.append(Quark(site,r,c,z,"UP" if v>curv else "DOWN"))
    return quarks

def find_proton_triplets(quarks, cfg):
    quark_set = {q.site: q for q in quarks}
    triplets = []; used = set()
    nbr_cache = {q.site: set(mesh_neighbours(q.r,q.c,q.z,cfg)) for q in quarks}
    for q in quarks:
        if q.quark_type!="DOWN" or q.site in used: continue
        up_nbrs = [quark_set[ni] for ni in nbr_cache[q.site]
                   if ni in quark_set and quark_set[ni].quark_type=="UP"
                   and ni not in used]
        if len(up_nbrs)<2: continue
        found=False
        for i in range(len(up_nbrs)):
            if found: break
            for j in range(i+1,len(up_nbrs)):
                p1,p2=up_nbrs[i].site,up_nbrs[j].site
                if p1 in nbr_cache.get(p2,set()):
                    triplets.append((q.site,p1,p2))
                    used.update([q.site,p1,p2])
                    found=True; break
    return triplets

# ── Main ─────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    cfg = LatticeConfig()
    N = cfg.hex_rows * cfg.hex_cols * cfg.layers
    n_seeds = 20  # more seeds = more fusion events = better statistics

    print("GUTOE: Fusion Product Identification")
    print("="*72)
    print(f"Lattice: {cfg.hex_rows}×{cfg.hex_cols}×{cfg.layers}, k=4, {n_seeds} seeds at t=25")
    print()

    # Tally all changed states across all merge events
    # A merge event = adjacent proton pair where at least one site
    # becomes grade-2+ after one step.
    product_counter    = Counter()  # what states appear in merge events
    reactant_counter   = Counter()  # what states went IN (should be grade-1 quarks)
    event_signatures   = []         # (frozenset of before-states, frozenset of after-states)
    total_pairs        = 0
    total_merges       = 0
    total_protons_seen = 0

    for seed_idx in range(n_seeds):
        rng = np.random.default_rng(seed_idx * 137 + 7)
        lat = init_lattice(cfg)
        for _ in range(25):
            lat = step(lat, rng, cfg)

        quarks = detect_quarks(lat, cfg)
        triplets = find_proton_triplets(quarks, cfg)
        total_protons_seen += len(triplets)

        proton_map = {}
        for i, (d,u1,u2) in enumerate(triplets):
            for s in [d,u1,u2]: proton_map[s] = i

        nbr_cache = {}
        for i, (d,u1,u2) in enumerate(triplets):
            for s in [d,u1,u2]:
                if s not in nbr_cache:
                    r,c,z = site_coords(s, cfg)
                    nbr_cache[s] = set(mesh_neighbours(r,c,z,cfg))

        pairs = {}
        for i, (d,u1,u2) in enumerate(triplets):
            for s in [d,u1,u2]:
                for nb in nbr_cache[s]:
                    if nb in proton_map and proton_map[nb] != i:
                        key = tuple(sorted([i, proton_map[nb]]))
                        pairs[key] = True
        total_pairs += len(pairs)

        # Step forward with a fixed rng (deterministic from snapshot)
        rng2 = np.random.default_rng(seed_idx * 999983 + 42)
        lat_after = step(lat, rng2, cfg)

        for (i, j) in pairs:
            sites_i = list(triplets[i])
            sites_j = list(triplets[j])
            all_sites = sites_i + sites_j

            before = [int(lat[s]) for s in all_sites]
            after  = [int(lat_after[s]) for s in all_sites]

            grades_after = [_GRADE_TABLE[a] for a in after if a != VOID]
            is_merge = any(g >= 2 for g in grades_after)

            if is_merge:
                total_merges += 1
                for b in before: reactant_counter[b] += 1
                for a in after:
                    if a != VOID:
                        product_counter[a] += 1
                # Record event signature: before-multiset → after-multiset
                sig = (tuple(sorted(before)), tuple(sorted(after)))
                event_signatures.append(sig)

    # ── Results ──────────────────────────────────────────────────────────
    print(f"Statistics ({n_seeds} seeds at t=25):")
    print(f"  Total protons seen:    {total_protons_seen}")
    print(f"  Adjacent pairs:        {total_pairs}")
    print(f"  Merge events:          {total_merges}")
    rate = total_merges / total_pairs if total_pairs > 0 else 0
    print(f"  Fusion rate:           {rate:.1%} per adjacent pair per step")
    print()

    # ── What states go IN ────────────────────────────────────────────────
    print("Reactant states (what enters the fusion zone):")
    for state, count in sorted(reactant_counter.items(), key=lambda x: -x[1]):
        g = _GRADE_TABLE[state]
        print(f"  {STATE_NAMES[state]:8s} (grade-{g}): {count:4d}  "
              f"{'█'*min(count//2, 40)}")
    print()

    # ── What states come OUT ─────────────────────────────────────────────
    print("Product states (what the fusion zone becomes):")
    grade2_total = sum(c for s,c in product_counter.items()
                       if _GRADE_TABLE[s] == 2)
    grade3_total = sum(c for s,c in product_counter.items()
                       if _GRADE_TABLE[s] == 3)
    for state, count in sorted(product_counter.items(), key=lambda x: -x[1]):
        g = _GRADE_TABLE[state]
        name = GRADE2_NAMES.get(state, STATE_NAMES[state])
        bar = '█' * min(count//2, 40)
        print(f"  {name:20s} (grade-{g}): {count:4d}  {bar}")
    print()
    print(f"  Grade-1 survivors:  {sum(c for s,c in product_counter.items() if _GRADE_TABLE[s]==1)}")
    print(f"  Grade-2 products:   {grade2_total}")
    print(f"  Grade-3 products:   {grade3_total}")
    print()

    # ── Most common full event signatures ────────────────────────────────
    sig_counter = Counter(event_signatures)
    print(f"Top fusion event signatures (before → after):")
    for (before, after), count in sig_counter.most_common(10):
        before_names = ' '.join(STATE_NAMES[s] for s in before)
        after_names  = ' '.join(STATE_NAMES[s] for s in after)
        print(f"  [{count}×]  {before_names}")
        print(f"         → {after_names}")
        print()

    # ── Verdict ──────────────────────────────────────────────────────────
    print("="*72)
    print("VERDICT")
    print("="*72)
    if product_counter:
        dominant = product_counter.most_common(1)[0]
        dom_state, dom_count = dominant
        dom_frac = dom_count / sum(product_counter.values())
        dom_grade = _GRADE_TABLE[dom_state]
        print(f"\n  Dominant fusion product: {STATE_NAMES[dom_state]} "
              f"({dom_frac:.0%} of all product sites, grade-{dom_grade})")
        if dom_grade == 2:
            if dom_state in (7, 11, 13):
                print(f"  → Spacelike bivector: SO(3) rotation generator.")
                print(f"  → Interpretation: two protons merge into a SPIN-1 bound state.")
                print(f"  → Candidate: deuteron analogue (spin-aligned pp).")
            elif dom_state in (4, 6, 10):
                print(f"  → Mixed bivector: Lorentz boost generator.")
                print(f"  → Interpretation: two protons merge into a BOOST state.")
        elif dom_grade == 1:
            print(f"  → Grade-1 dominant: protons mostly survive, partial rearrangement.")
        elif dom_grade == 3:
            print(f"  → Trivector dominant: three-index structure, possible baryon merger.")

        # Check if grade-2 products cluster in one Z₃ orbit
        spacelike = sum(product_counter[s] for s in [7,11,13])
        mixed     = sum(product_counter[s] for s in [4,6,10])
        if grade2_total > 0:
            print(f"\n  Grade-2 breakdown:")
            print(f"    Spacelike bivectors {{γ¹², γ¹³, γ²³}} (rotation): {spacelike}  "
                  f"({spacelike/grade2_total:.0%})")
            print(f"    Mixed bivectors {{γ⁰¹, γ⁰², γ⁰³}} (boost):      {mixed}  "
                  f"({mixed/grade2_total:.0%})")
            if spacelike > 2 * mixed:
                print(f"  → Rotations dominate: product has angular momentum structure.")
            elif mixed > 2 * spacelike:
                print(f"  → Boosts dominate: product has relativistic boost structure.")
            else:
                print(f"  → Mixed: no clear geometric preference.")
