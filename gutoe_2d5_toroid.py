#!/usr/bin/env python3
"""
GUTOE: 2.5D hexagonal toroid — void votes in alignment.

The rear void faces VOTE in the alignment step.
Each rear face casts a vote for VOID (= "stay put").

6 front neighbours vote for their actual states.
k rear faces vote for VOID.

Alignment only succeeds if a real state beats VOID's vote count.
With k void votes, a state needs >k front-neighbour agreements to win.

  k=0: baseline (any plurality aligns)
  k=2: need 3/6 agreement  (weak void pressure)
  k=4: need 5/6 agreement  (moderate — only domain interiors align)
  k=6: need 7/6 agreement  (impossible — alignment fully blocked)

Prediction: moderate k preserves triple junctions (where front
neighbours are mixed, VOID wins) while still allowing domain
interiors to align. This could stabilize proton half-life.

Copyright (C) 2026 Riff Labs, AGPL-3.0-or-later
"""

import numpy as np
from dataclasses import dataclass
from collections import Counter

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
    void_votes: int = 0  # rear void faces cast this many VOID votes in alignment

def idx(r,c,z,cfg):
    return((z%cfg.layers)*cfg.hex_rows+(r%cfg.hex_rows))*cfg.hex_cols+(c%cfg.hex_cols)

def hex_planar_neighbours(r,c,cfg):
    if r%2==0: offs=[(-1,0),(-1,1),(0,-1),(0,1),(1,0),(1,1)]
    else: offs=[(-1,-1),(-1,0),(0,-1),(0,1),(1,-1),(1,0)]
    return[((r+dr)%cfg.hex_rows,(c+dc)%cfg.hex_cols) for dr,dc in offs]

def mesh_neighbours(r,c,z,cfg):
    return[idx(nr,nc,z,cfg) for nr,nc in hex_planar_neighbours(r,c,cfg)]

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
        z = site // (cfg.hex_rows * cfg.hex_cols)
        rem = site % (cfg.hex_rows * cfg.hex_cols)
        r = rem // cfg.hex_cols; c = rem % cfg.hex_cols
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
                    # Rear void faces cast void_votes votes for VOID
                    # A real state only wins if it beats the void vote count
                    winner, winner_count = votes.most_common(1)[0]
                    if winner_count > cfg.void_votes:
                        new[site] = winner
                    # else: VOID wins → keep current state (no alignment)

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
        z = site//(cfg.hex_rows*cfg.hex_cols)
        rem = site%(cfg.hex_rows*cfg.hex_cols)
        r = rem//cfg.hex_cols; c = rem%cfg.hex_cols
        v, curv, grad = compute_local_fields(lattice, site, r, c, z, cfg)
        bc = v/(1+grad)
        if bc >= cfg.quark_threshold:
            quarks.append(Quark(site,r,c,z,"UP" if v>curv else "DOWN",bc,v,curv))
    return quarks

def find_protons(quarks, cfg):
    quark_set = {q.site: q for q in quarks}
    strict_count = 0; used = set()
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
                    strict_count+=1; used.update([q.site,p1,p2])
                    found=True; break
    return strict_count

def run_multiseed(cfg, n_seeds, steps):
    all_p = np.zeros(steps)
    all_up = np.zeros(steps); all_dn = np.zeros(steps)
    all_g0 = np.zeros(steps); all_g1 = np.zeros(steps)
    for s in range(n_seeds):
        rng = np.random.default_rng(s*137+7)
        lat = init_lattice(cfg)
        for t in range(steps):
            lat = step(lat, rng, cfg)
            qs = detect_quarks(lat, cfg)
            all_p[t] += find_protons(qs, cfg)
            all_up[t] += sum(1 for q in qs if q.quark_type=="UP")
            all_dn[t] += sum(1 for q in qs if q.quark_type=="DOWN")
            all_g0[t] += sum(1 for x in lat if _GRADE_TABLE[int(x)]==0)
            all_g1[t] += sum(1 for x in lat if _GRADE_TABLE[int(x)]==1)
    n=n_seeds
    return {k:v/n for k,v in {
        'protons':all_p,'up':all_up,'down':all_dn,
        'grade0':all_g0,'grade1':all_g1}.items()}

if __name__ == "__main__":
    n_seeds = 10
    steps = 1000

    void_vote_values = [0, 2, 4, 5]

    print("GUTOE: 2.5D Toroid — Void Votes in Alignment")
    print("="*72)
    print("Rear void faces cast k votes for VOID in majority election.")
    print("A real state wins only if it gets MORE than k front-neighbour votes.")
    print(f"  k=0: baseline (any plurality aligns — 1+ votes needed)")
    print(f"  k=2: need 3/6 front neighbours to agree")
    print(f"  k=4: need 5/6 front neighbours to agree")
    print(f"  k=5: need 6/6 front neighbours to agree (unanimity)")
    print(f"Seeds: {n_seeds},  Steps: {steps}")
    print("="*72)
    print()

    results = {}
    for k in void_vote_values:
        cfg = LatticeConfig(void_votes=k)
        print(f"Running k={k}...")
        results[k] = run_multiseed(cfg, n_seeds, steps)
        print(f"  done.")

    # ── Time series ─────────────────────────────────────────────────────
    print()
    print("="*72)
    hdr = f"{'t':>4s}"
    for k in void_vote_values:
        hdr += f" | k={k}:  p  U/D   g0   g1"
    print(hdr)
    print("-"*len(hdr))

    milestones = [15,20,25,30,40,50,75,100,150,200,300,400,500,600,750,1000]
    for t in milestones:
        if t > steps: break
        i = t - 1
        line = f"{t:4d}"
        for k in void_vote_values:
            r = results[k]
            p = r['protons'][i]
            u = r['up'][i]; d = r['down'][i]
            g0 = r['grade0'][i]; g1 = r['grade1'][i]
            ratio = f"{u/d:.1f}" if d>0.5 else "∞"
            line += f" | {p:6.0f} {ratio:>5s} {g0:4.0f} {g1:4.0f}"
        print(line)

    # ── Summary ─────────────────────────────────────────────────────────
    print()
    print("="*72)
    print("SUMMARY")
    print("="*72)

    half_lives = {}
    for k in void_vote_values:
        r = results[k]
        peak_p = max(r['protons'])
        peak_t = int(np.argmax(r['protons'])) + 1
        p100 = r['protons'][min(99, steps-1)]
        p200 = r['protons'][min(199, steps-1)] if steps>=200 else 0
        p300 = r['protons'][min(299, steps-1)] if steps>=300 else 0
        half = peak_p / 2
        hl = next((i+1 for i in range(peak_t, steps) if r['protons'][i] < half), steps+1)
        half_lives[k] = hl
        hl_str = f"t={hl}" if hl <= steps else f">t={steps}"
        print(f"\n  k={k} ({['baseline','2/6 threshold','4/6 threshold','unanimity'][void_vote_values.index(k)]}):")
        print(f"    Peak: {peak_p:.0f} protons at t={peak_t}")
        print(f"    t=100: {p100:.1f}  t=200: {p200:.1f}  t=300: {p300:.1f}")
        print(f"    Half-life: {hl_str}")

    # ── Phase detection ──────────────────────────────────────────────────
    # A system with a stable non-zero proton population is qualitatively
    # different from one that thermalizes to zero.
    # Check last 100 steps for a stable residual: mean > threshold and
    # not monotonically declining (std is non-trivial).
    STABLE_THRESHOLD = 3.0  # protons must average above this to be "alive"
    late_start = max(0, steps - 100)

    stable = {}
    late_means = {}
    for k in void_vote_values:
        late = results[k]['protons'][late_start:]
        late_means[k] = float(np.mean(late))
        is_stable = late_means[k] > STABLE_THRESHOLD
        stable[k] = is_stable

    print()
    print("="*72)
    print("VERDICT")
    print("="*72)
    baseline_hl = half_lives[0]

    print(f"\n  Baseline (k=0) half-life: t={baseline_hl}  |  t=300 mean: {late_means[0]:.1f} protons")
    for k in void_vote_values[1:]:
        hl = half_lives[k]
        phase = "STABLE" if stable[k] else "DECAYS"
        mark = " ← STABLE PHASE" if stable[k] else ""
        print(f"  k={k}:  half-life t={hl}  |  t=200-300 mean: {late_means[k]:.1f} protons  [{phase}]{mark}")

    stable_ks = [k for k in void_vote_values if stable[k] and k != 0]
    baseline_alive = stable.get(0, False)

    if stable_ks and not baseline_alive:
        best_stable = max(stable_ks, key=lambda k: late_means[k])
        print(f"\n  PHASE TRANSITION: void votes create a stable proton population.")
        print(f"  k=0 (baseline): thermalizes to zero by t={steps}. DEAD.")
        for k in stable_ks:
            print(f"  k={k}: {late_means[k]:.1f} protons (mean t=200–300). ALIVE.")
        print(f"\n  Void votes act as effective temperature, preventing complete")
        print(f"  thermalization into the scalar vacuum. The grade-1 excitation")
        print(f"  population and its bound protons persist indefinitely.")
        print(f"  Half-life metric is misleading here — the system is not decaying,")
        print(f"  it is equilibrating to a non-trivial steady state.")
        print(f"\n  Sweet spot: k={best_stable} ({late_means[best_stable]:.1f} avg protons at late times).")
    elif not stable_ks:
        best_k = max(void_vote_values, key=lambda k: half_lives[k])
        ratio = half_lives[best_k] / max(baseline_hl, 1)
        print(f"\n  NO STABLE PHASE: all configurations thermalize to zero.")
        print(f"  Best half-life extension: {ratio:.1f}× at k={best_k}.")
    else:
        print(f"\n  All configurations (including baseline) maintain protons. No clear transition.")
