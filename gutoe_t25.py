#!/usr/bin/env python3
"""
GUTOE: t=25 Dense Phase Analysis + Lepton Injection

At t=25 (percolation peak), ~94 protons are packed in a 12×12×12 lattice.
This is the densest proton moment in the simulation.

Phase 1 — Snapshot analysis at t=25:
  - Adjacent proton pairs (sharing a lattice edge) → fusion candidates
  - What happens at t=26? Do adjacent protons merge, repel, or ignore?

Phase 2 — Lepton injection at t=25:
  - Inject N γ⁰ cells at random non-proton sites
  - Run 200 more steps
  - Track: γ⁰ survival, γ⁰ proximity to protons, hydrogen formation
  - γ⁰ is Z₃ fixed — it won't cycle away, must survive alignment pressure

Copyright (C) 2026 Riff Labs, AGPL-3.0-or-later
"""

import numpy as np
from collections import Counter
from dataclasses import dataclass

VOID        = 0
QUARK_SEED  = 3   # γ¹
LEPTON_SEED = 2   # γ⁰

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

def proton_shell(triplets, proton_sites, cfg):
    shell = set()
    for site in proton_sites:
        r, c, z = site_coords(site, cfg)
        for ni in mesh_neighbours(r, c, z, cfg):
            if ni not in proton_sites:
                shell.add(ni)
    return shell

# ── Main ─────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    cfg = LatticeConfig()
    N = cfg.hex_rows * cfg.hex_cols * cfg.layers
    n_seeds = 5

    print("GUTOE: t=25 Dense Phase — Fusion + Lepton Injection")
    print("="*72)
    print(f"Lattice: {cfg.hex_rows}×{cfg.hex_cols}×{cfg.layers} = {N} sites, k=4")
    print()

    # ── PHASE 1: Fusion analysis ──────────────────────────────────────────
    print("PHASE 1: FUSION — adjacent proton pairs at t=25")
    print("-"*72)

    total_protons = 0
    total_adjacent = 0
    total_merge = 0
    total_repel = 0
    total_stable = 0

    for seed_idx in range(n_seeds):
        rng = np.random.default_rng(seed_idx * 137 + 7)
        lat = init_lattice(cfg)
        for _ in range(25):
            lat = step(lat, rng, cfg)

        quarks = detect_quarks(lat, cfg)
        triplets = find_proton_triplets(quarks, cfg)
        n_p = len(triplets)
        total_protons += n_p

        # Build proton site map
        proton_map = {}  # site → triplet index
        for i, (d,u1,u2) in enumerate(triplets):
            for s in [d,u1,u2]: proton_map[s] = i

        nbr_cache = {}
        for i, (d,u1,u2) in enumerate(triplets):
            for s in [d,u1,u2]:
                if s not in nbr_cache:
                    r,c,z = site_coords(s, cfg)
                    nbr_cache[s] = set(mesh_neighbours(r,c,z,cfg))

        # Find adjacent proton pairs
        pairs = {}
        for i, (d,u1,u2) in enumerate(triplets):
            for s in [d,u1,u2]:
                for nb in nbr_cache[s]:
                    if nb in proton_map and proton_map[nb] != i:
                        key = tuple(sorted([i, proton_map[nb]]))
                        pairs[key] = True
        total_adjacent += len(pairs)

        # Run one step and check outcomes
        rng2 = np.random.default_rng(seed_idx * 137 + 7 + 50000)
        lat_after = step(lat, rng2, cfg)

        for (i, j) in pairs:
            sites_i = list(triplets[i])
            sites_j = list(triplets[j])
            before_i = [int(lat[s]) for s in sites_i]
            before_j = [int(lat[s]) for s in sites_j]
            after_i  = [int(lat_after[s]) for s in sites_i]
            after_j  = [int(lat_after[s]) for s in sites_j]
            grades_after = [_GRADE_TABLE[s] for s in after_i+after_j if s!=VOID]
            changed = sum(b!=a for b,a in zip(before_i+before_j, after_i+after_j))
            if any(g >= 2 for g in grades_after):
                total_merge += 1
            elif changed == 0:
                total_stable += 1
            else:
                total_repel += 1

        print(f"  seed {seed_idx}: {n_p} protons, {len(pairs)} adjacent pairs", end='')
        if pairs:
            print(f"  (merge={sum(1 for (i,j) in pairs if any(_GRADE_TABLE[int(lat_after[s])]>=2 for s in list(triplets[i])+list(triplets[j])))}, "
                  f"repel/other={len(pairs)-sum(1 for (i,j) in pairs if any(_GRADE_TABLE[int(lat_after[s])]>=2 for s in list(triplets[i])+list(triplets[j])))})")
        else:
            print()

    print(f"\n  Total: {total_protons/n_seeds:.0f} avg protons, "
          f"{total_adjacent/n_seeds:.1f} avg adjacent pairs")
    print(f"  After t=26: merge={total_merge}  repel/disrupted={total_repel}  stable={total_stable}")
    print()

    # ── PHASE 2: Lepton injection at t=25 ────────────────────────────────
    print("PHASE 2: LEPTON INJECTION — drop γ⁰ into dense proton phase")
    print("-"*72)
    print("Inject γ⁰ cells at t=25, run 200 more steps, track binding.")
    print()

    n_inject = 20   # γ⁰ cells to inject
    run_more = 200

    for seed_idx in range(n_seeds):
        # Run to t=25 (quark-only)
        rng = np.random.default_rng(seed_idx * 137 + 7)
        lat = init_lattice(cfg)
        for _ in range(25):
            lat = step(lat, rng, cfg)

        quarks_before = detect_quarks(lat, cfg)
        triplets_before = find_proton_triplets(quarks_before, cfg)
        proton_sites_before = set(s for d,u1,u2 in triplets_before for s in [d,u1,u2])

        # Inject γ⁰ at random non-proton sites
        non_proton = [i for i in range(N) if i not in proton_sites_before]
        inject_sites = rng.choice(non_proton, size=min(n_inject, len(non_proton)),
                                  replace=False)
        lat_with_leptons = lat.copy()
        for s in inject_sites:
            lat_with_leptons[s] = LEPTON_SEED

        # Track over time
        lat_run = lat_with_leptons.copy()
        lepton_counts = []
        proton_counts = []
        hydrogen_counts = []
        shell_rates = []
        bg_rates = []

        for t in range(run_more):
            lat_run = step(lat_run, rng, cfg)
            n_l = sum(1 for i in range(N) if int(lat_run[i]) == LEPTON_SEED)
            lepton_counts.append(n_l)

            if (t+1) % 25 == 0 or t == run_more-1:
                qs = detect_quarks(lat_run, cfg)
                trips = find_proton_triplets(qs, cfg)
                n_p = len(trips)
                proton_counts.append((t+1, n_p))

                p_sites = set(s for d,u1,u2 in trips for s in [d,u1,u2])
                shell = proton_shell(trips, p_sites, cfg)

                n_lep_shell = sum(1 for i in range(N)
                                  if int(lat_run[i])==LEPTON_SEED and i in shell)
                n_lep_bg    = n_l - n_lep_shell
                shell_sz    = len(shell)
                bg_sz       = N - len(p_sites) - shell_sz

                rs = n_lep_shell / shell_sz if shell_sz > 0 else 0
                rb = n_lep_bg   / bg_sz     if bg_sz    > 0 else 0
                shell_rates.append(rs)
                bg_rates.append(rb)

                # Hydrogen: protons with adjacent γ⁰
                n_h = 0
                for d,u1,u2 in trips:
                    nbrs = set()
                    for s in [d,u1,u2]:
                        r,c,z = site_coords(s, cfg)
                        nbrs.update(mesh_neighbours(r,c,z,cfg))
                    nbrs -= p_sites
                    if any(int(lat_run[ns])==LEPTON_SEED for ns in nbrs):
                        n_h += 1
                hydrogen_counts.append((t+1, n_h))

        final_l = lepton_counts[-1]
        survival = final_l / n_inject * 100
        avg_rs = np.mean(shell_rates) if shell_rates else 0
        avg_rb = np.mean(bg_rates) if bg_rates else 0
        enrich = avg_rs / avg_rb if avg_rb > 1e-9 else float('inf')

        print(f"  seed {seed_idx}: injected {n_inject} γ⁰")
        print(f"    Survival at t=225: {final_l}/{n_inject} ({survival:.0f}%)")
        proton_str = "  protons: " + "  ".join(f"t={t}:{p}" for t,p in proton_counts)
        print(f"    {proton_str}")
        hydro_str = "hydrogen: " + "  ".join(f"t={t}:{h}" for t,h in hydrogen_counts)
        print(f"    {hydro_str}")
        print(f"    γ⁰ rate shell={avg_rs:.4f}  bg={avg_rb:.4f}  enrichment={enrich:.2f}×")
        if enrich > 1.5:
            print(f"    → BINDING: γ⁰ preferentially found near protons!")
        elif enrich < 0.7:
            print(f"    → EXCLUSION: protons repel γ⁰ (like-charge?)")
        else:
            print(f"    → No spatial preference.")
        print()
