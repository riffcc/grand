#!/usr/bin/env python3
"""
GUTOE: Lepton Seeding — Hydrogen from First Principles

Add γ⁰ (state 2) as a second seed species alongside γ¹ (state 3).

γ⁰ is the timelike unit vector — a Z₃ FIXED POINT. It doesn't cycle.
Unlike quarks {γ¹,γ²,γ³} which form Z₃ orbits of size 3 (colored),
γ⁰ has trivial Z₃ representation (colorless) → the natural lepton.

Key questions:
  1. Does γ⁰ survive? (Z₃ cycle doesn't kill it, but alignment might)
  2. Is γ⁰ preferentially found near proton triangles? (binding)
  3. Does γ⁰ fraction near protons exceed background rate? (selective binding)

Why γ⁰ might bind to protons:
  - Protons live at triple junctions (3-way quark domain boundaries)
  - At triple junctions, no single quark state wins >4 votes → alignment blocked
  - γ⁰ is ALSO protected from alignment at mixed boundaries
  - → γ⁰ and protons co-inhabit the same protected niches

Clifford interaction: γ⁰ · γ¹ = γ⁰¹ (grade-2 bivector) — the EM vertex.

Copyright (C) 2026 Riff Labs, AGPL-3.0-or-later
"""

import numpy as np
from collections import Counter
from dataclasses import dataclass

# ── Physics (k=4 void votes, as in stable phase) ─────────────────────────────

VOID         = 0
QUARK_SEED   = 3   # γ¹ — spacelike, Z₃ orbit member, colored
LEPTON_SEED  = 2   # γ⁰ — timelike, Z₃ fixed point, colorless

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

@dataclass
class LatticeConfig:
    hex_rows: int = 12
    hex_cols: int = 12
    layers: int = 12
    differentiation_prob: float = 0.02
    lepton_fraction: float = 0.1    # fraction of seeds that become γ⁰ (vs γ¹)
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
                # Seed quark or lepton
                seed = LEPTON_SEED if rng.random() < cfg.lepton_fraction else QUARK_SEED
                new[site] = seed; continue
            nbrs = mesh_neighbours(r, c, z, cfg)
            active = sum(1 for ni in nbrs if lattice[ni] != VOID)
            total = len(nbrs)
            if active >= max(2, total//4):
                if rng.random() < active/total*0.4:
                    seed = LEPTON_SEED if rng.random() < cfg.lepton_fraction else QUARK_SEED
                    new[site] = seed
        else:
            r_val = rng.random()
            if r_val < cfg.cycle_prob:
                # Z₃ cycle: γ⁰ maps to itself (fixed point)
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
            quarks.append(Quark(site,r,c,z,"UP" if v>curv else "DOWN",bc))
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
                    triplets.append((q.site,p1,p2)); used.update([q.site,p1,p2])
                    found=True; break
    return triplets

# ── Hydrogen analysis ─────────────────────────────────────────────────────────

def analyze_snapshot(lattice, cfg):
    quarks = detect_quarks(lattice, cfg)
    triplets = find_proton_triplets(quarks, cfg)
    proton_sites = set()
    for d,u1,u2 in triplets:
        proton_sites.update([d,u1,u2])

    N = cfg.hex_rows * cfg.hex_cols * cfg.layers

    # Count γ⁰ (lepton) cells
    n_lepton = sum(1 for i in range(N) if int(lattice[i]) == LEPTON_SEED)
    n_g1_quark = sum(1 for i in range(N)
                     if _GRADE_TABLE[int(lattice[i])]==1 and int(lattice[i])!=LEPTON_SEED)

    if n_lepton == 0:
        return len(triplets), 0, 0, 0.0, 0.0

    # Build proton shell (all sites neighboring any proton quark)
    proton_shell = set()
    for site in proton_sites:
        r, c, z = site_coords(site, cfg)
        for ni in mesh_neighbours(r, c, z, cfg):
            if ni not in proton_sites:
                proton_shell.add(ni)

    # Count leptons in shell vs. background
    lepton_in_shell  = sum(1 for i in range(N)
                           if int(lattice[i])==LEPTON_SEED and i in proton_shell)
    lepton_elsewhere = n_lepton - lepton_in_shell

    shell_size = len(proton_shell)
    non_shell_size = N - len(proton_sites) - shell_size

    # Rates: fraction of shell/non-shell that is γ⁰
    rate_shell    = lepton_in_shell  / shell_size     if shell_size     > 0 else 0
    rate_bg       = lepton_elsewhere / non_shell_size if non_shell_size > 0 else 0

    # Hydrogen: protons with at least one adjacent γ⁰
    hydrogen_count = 0
    for d,u1,u2 in triplets:
        nbr_sites = set()
        for site in [d,u1,u2]:
            r,c,z = site_coords(site, cfg)
            nbr_sites.update(mesh_neighbours(r,c,z,cfg))
        nbr_sites -= proton_sites
        if any(int(lattice[ns])==LEPTON_SEED for ns in nbr_sites):
            hydrogen_count += 1

    return len(triplets), n_lepton, hydrogen_count, rate_shell, rate_bg

# ── Main ─────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    cfg = LatticeConfig()
    N = cfg.hex_rows * cfg.hex_cols * cfg.layers

    lepton_fractions = [0.05, 0.10, 0.20]
    steps = 1000
    n_seeds = 5

    print("GUTOE: Lepton Seeding — Hydrogen Test")
    print("="*72)
    print(f"Lattice: {cfg.hex_rows}×{cfg.hex_cols}×{cfg.layers} = {N} sites")
    print(f"Quark seed:  γ¹ (state 3) — Z₃ orbit, colored")
    print(f"Lepton seed: γ⁰ (state 2) — Z₃ fixed point, colorless")
    print(f"k=4 void votes (stable proton phase)")
    print(f"Lepton fraction: fraction of seeds that become γ⁰")
    print(f"{n_seeds} seeds × {steps} steps")
    print()
    print(f"  Veracity γ⁰↔γ¹ = 0.5  (Hamming distance 2 → not a quark)")
    print(f"  Z₃ cycle: γ⁰ → γ⁰  (fixed, doesn't rotate away)")
    print(f"  Clifford: γ⁰·γ¹ = γ⁰¹  (EM vertex)")
    print("="*72)
    print()

    for lf in lepton_fractions:
        cfg_lf = LatticeConfig(lepton_fraction=lf)
        print(f"lepton_fraction={lf:.0%}  ({lf:.0%} of seeds → γ⁰)")
        print(f"  Running {n_seeds} seeds to t={steps}...")

        total_protons = 0
        total_leptons = 0
        total_hydrogen = 0
        rate_shells = []
        rate_bgs = []

        for seed_idx in range(n_seeds):
            rng = np.random.default_rng(seed_idx * 137 + 7)
            lat = init_lattice(cfg_lf)
            for t in range(steps):
                lat = step(lat, rng, cfg_lf)

            n_p, n_l, n_h, rs, rb = analyze_snapshot(lat, cfg_lf)
            total_protons  += n_p
            total_leptons  += n_l
            total_hydrogen += n_h
            if n_l > 0:
                rate_shells.append(rs)
                rate_bgs.append(rb)

        avg_p = total_protons  / n_seeds
        avg_l = total_leptons  / n_seeds
        avg_h = total_hydrogen / n_seeds
        avg_rs = np.mean(rate_shells) if rate_shells else 0
        avg_rb = np.mean(rate_bgs)   if rate_bgs   else 0
        enrich = avg_rs / avg_rb if avg_rb > 1e-9 else float('inf')

        print(f"  Protons (avg):          {avg_p:.1f}")
        print(f"  γ⁰ leptons (avg):       {avg_l:.1f}")
        print(f"  Hydrogen atoms (avg):   {avg_h:.1f}  "
              f"({avg_h/avg_p*100:.0f}% of protons bound)" if avg_p>0 else "")
        print(f"  γ⁰ rate near protons:   {avg_rs:.4f}")
        print(f"  γ⁰ rate background:     {avg_rb:.4f}")
        print(f"  Enrichment ratio:       {enrich:.2f}×")
        if enrich > 1.5:
            print(f"  → γ⁰ preferentially binds to proton neighborhoods!")
        elif enrich > 1.0:
            print(f"  → Slight preference, but weak.")
        else:
            print(f"  → No preferential binding.")
        print()
