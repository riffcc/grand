#!/usr/bin/env python3
"""
GUTOE: Electroweak Extension

Adds electromagnetism and the weak force to the stable proton simulation.

ELECTROMAGNETISM (U(1)):
  The photon is a grade-2 bivector — we just proved this: pp fusion
  produces γ¹² (SO(3) rotation generator) = the spin-1 photon candidate.
  EM force: γ⁰ (lepton) surfs through the grade-2 photon field toward
  grade-1 quarks. γ⁰ swaps with an adjacent grade-2 site if that site
  is adjacent to a quark. This is photon-mediated attraction.
  γ⁰ is IMMUNE to the standard alignment step — it has its own EM dynamics.

WEAK FORCE (SU(2)):
  The W boson is γ⁰¹ (mixed bivector, mi=0011): γ⁰ · γ¹ = γ⁰¹.
  Weak vertex: γ⁰ (lepton) + γ¹ (quark) → γ⁰¹ (W boson) + VOID (neutrino)
  This converts: lepton + quark → W boson, lepton absorbed.
  W boson (γ⁰¹) can decay back: γ⁰¹ → γ⁰ + γ¹ (reverse weak vertex).
  Net effect on proton: a quark (γ¹) becomes W boson (γ⁰¹), changing
  the proton's structure — the weak force transmutes quarks.

NEUTRON:
  Proton  = triangle of 2 UP + 1 DOWN quarks.
  Neutron = triangle of 1 UP + 2 DOWN quarks.
  Detected the same way as protons but with flipped majority.

DEUTERON:
  Adjacent proton + neutron triangles sharing a lattice edge.
  This is the pp chain fusion product (via weak quark flip).

HYDROGEN:
  Proton triangle + adjacent γ⁰ cell = hydrogen atom.
  With EM attraction, γ⁰ should preferentially bind to proton shells.

TWO-PHASE PROTOCOL:
  Phase 1: Run k=4 to t=500 (stable proton phase, no leptons).
  Phase 2: Inject γ⁰ leptons, activate EM+weak, run 500 more steps.
  Compare: hydrogen count, γ⁰ enrichment near protons vs. baseline.

Copyright (C) 2026 Riff Labs, AGPL-3.0-or-later
"""

import numpy as np
from collections import Counter
from dataclasses import dataclass

VOID        = 0
QUARK_SEED  = 3   # γ¹ — spacelike grade-1, Z₃ orbit member
LEPTON_SEED = 2   # γ⁰ — timelike grade-1, Z₃ fixed point

_GRADE_TABLE = [-1] + [bin(mi).count('1') for mi in range(16)]

STATE_NAMES = {
    0:'VOID', 1:'scalar', 2:'γ⁰', 3:'γ¹', 4:'γ⁰¹', 5:'γ²', 6:'γ⁰²',
    7:'γ¹²', 8:'γ⁰¹²', 9:'γ³', 10:'γ⁰³', 11:'γ¹³', 12:'γ⁰¹³',
    13:'γ²³', 14:'γ⁰²³', 15:'γ¹²³', 16:'γ⁰¹²³',
}

W_BOSON_STATES = frozenset({4, 6, 10})   # γ⁰¹, γ⁰², γ⁰³ — mixed bivectors
PHOTON_STATES  = frozenset({7, 11, 13})  # γ¹², γ¹³, γ²³ — spacelike bivectors

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
    # Electroweak
    em_prob: float = 0.0     # probability γ⁰ drifts through photon field
    weak_prob: float = 0.0   # probability weak vertex fires
    w_decay_prob: float = 0.0  # probability W boson decays back to quark+lepton

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

        # ── VOID ──────────────────────────────────────────────────────
        if state == VOID:
            if rng.random() < cfg.differentiation_prob:
                new[site] = QUARK_SEED; continue
            nbrs = mesh_neighbours(r, c, z, cfg)
            active = sum(1 for ni in nbrs if lattice[ni] != VOID)
            total = len(nbrs)
            if active >= max(2, total//4):
                if rng.random() < active/total*0.4:
                    new[site] = QUARK_SEED

        # ── LEPTON (γ⁰) — EM + weak, IMMUNE to alignment ─────────────
        elif state == LEPTON_SEED:
            r_val = rng.random()

            if r_val < cfg.em_prob:
                # EM drift: γ⁰ surfs through grade-2 (photon) field toward grade-1 (quark)
                # Find grade-2 neighbors adjacent to grade-1 quarks
                nbrs = mesh_neighbours(r, c, z, cfg)
                photon_gateway = []
                for ni in nbrs:
                    ni_state = int(lattice[ni])
                    if _GRADE_TABLE[ni_state] == 2:  # photon site
                        r2,c2,z2 = site_coords(ni, cfg)
                        ni_nbrs = mesh_neighbours(r2, c2, z2, cfg)
                        if any(_GRADE_TABLE[int(lattice[nni])]==1 and nni!=site
                               for nni in ni_nbrs):
                            photon_gateway.append(ni)
                if photon_gateway:
                    target = photon_gateway[rng.integers(len(photon_gateway))]
                    target_state = int(lattice[target])
                    new[site] = target_state   # photon comes to lepton's old site
                    new[target] = LEPTON_SEED  # lepton moves to photon's site

            elif r_val < cfg.em_prob + cfg.weak_prob:
                # Weak vertex: γ⁰ (lepton) + γ¹ (quark) → γ⁰¹ (W boson) + VOID (ν)
                nbrs = mesh_neighbours(r, c, z, cfg)
                quark_nbrs = [ni for ni in nbrs
                              if _GRADE_TABLE[int(lattice[ni])]==1
                              and int(lattice[ni]) != LEPTON_SEED]
                if quark_nbrs:
                    q_site = quark_nbrs[rng.integers(len(quark_nbrs))]
                    q_state = int(lattice[q_site])
                    w_state = ((LEPTON_SEED-1) ^ (q_state-1)) + 1  # γ⁰ XOR quark = W
                    new[q_site] = w_state  # quark → W boson
                    new[site]   = VOID     # lepton → neutrino (escapes lattice)
            # else: lepton does nothing this step (no alignment applied)

        # ── W BOSON — can decay back ───────────────────────────────────
        elif state in W_BOSON_STATES and cfg.w_decay_prob > 0:
            if rng.random() < cfg.w_decay_prob:
                # W → quark + lepton (reverse weak vertex)
                # W (γ⁰¹, mi=0011) → find a VOID neighbor to emit lepton into
                nbrs = mesh_neighbours(r, c, z, cfg)
                void_nbrs = [ni for ni in nbrs if int(lattice[ni]) == VOID]
                if void_nbrs:
                    emit_site = void_nbrs[rng.integers(len(void_nbrs))]
                    # W decays: W → quark (recover original γ¹ via XOR with γ⁰)
                    q_state = ((state-1) ^ (LEPTON_SEED-1)) + 1
                    new[site]      = q_state       # W → quark
                    new[emit_site] = LEPTON_SEED   # emit lepton into void
            else:
                # W boson standard dynamics (Z₃ cycle + alignment)
                r_val = rng.random()
                if r_val < cfg.cycle_prob:
                    new[site] = _Z3_TABLE[state]
                elif r_val < cfg.cycle_prob + cfg.alignment_strength:
                    nbrs = mesh_neighbours(r, c, z, cfg)
                    nbr_states = [int(lattice[ni]) for ni in nbrs if lattice[ni]!=VOID]
                    if nbr_states:
                        votes = Counter(nbr_states)
                        winner, cnt = votes.most_common(1)[0]
                        if cnt > cfg.void_votes:
                            new[site] = winner

        # ── QUARKS + ALL OTHER STATES — standard dynamics ─────────────
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
                    winner, cnt = votes.most_common(1)[0]
                    if cnt > cfg.void_votes:
                        new[site] = winner
    return new

# ── Particle detection ────────────────────────────────────────────────────────

@dataclass
class Quark:
    site: int; r: int; c: int; z: int; quark_type: str

def detect_quarks(lattice, cfg):
    quarks = []
    N = cfg.hex_rows * cfg.hex_cols * cfg.layers
    for site in range(N):
        state = int(lattice[site])
        if state == VOID or state == LEPTON_SEED: continue
        if _GRADE_TABLE[state] != 1: continue  # quarks are grade-1 non-lepton
        r, c, z = site_coords(site, cfg)
        v, curv, grad = compute_local_fields(lattice, site, r, c, z, cfg)
        bc = v/(1+grad)
        if bc >= cfg.quark_threshold:
            quarks.append(Quark(site,r,c,z,"UP" if v>curv else "DOWN"))
    return quarks

def find_baryon_triplets(quarks, cfg, n_up_req, n_dn_req):
    """Find triangles with n_up_req UP and n_dn_req DOWN quarks.
       Proton: n_up=2, n_dn=1.  Neutron: n_up=1, n_dn=2."""
    quark_set = {q.site: q for q in quarks}
    triplets = []; used = set()
    nbr_cache = {q.site: set(mesh_neighbours(q.r,q.c,q.z,cfg)) for q in quarks}
    minority = "UP" if n_up_req < n_dn_req else "DOWN"
    majority = "DOWN" if minority=="UP" else "UP"
    for q in quarks:
        if q.quark_type != minority or q.site in used: continue
        maj_nbrs = [quark_set[ni] for ni in nbr_cache[q.site]
                    if ni in quark_set and quark_set[ni].quark_type==majority
                    and ni not in used]
        if len(maj_nbrs) < 2: continue
        found=False
        for i in range(len(maj_nbrs)):
            if found: break
            for j in range(i+1, len(maj_nbrs)):
                p1,p2 = maj_nbrs[i].site, maj_nbrs[j].site
                if p1 in nbr_cache.get(p2, set()):
                    triplets.append((q.site,p1,p2))
                    used.update([q.site,p1,p2])
                    found=True; break
    return triplets

def find_deuterons(proton_trips, neutron_trips, cfg):
    """Adjacent proton+neutron pairs sharing a lattice edge."""
    p_sites = {}
    for i,(d,u1,u2) in enumerate(proton_trips):
        for s in [d,u1,u2]: p_sites[s] = i
    deuterons = []; used = set()
    for j,(d,u1,u2) in enumerate(neutron_trips):
        if j in used: continue
        for s in [d,u1,u2]:
            r,c,z = site_coords(s, cfg)
            for nb in mesh_neighbours(r,c,z,cfg):
                if nb in p_sites:
                    pi = p_sites[nb]
                    key = (pi, j)
                    if key not in used:
                        deuterons.append(key)
                        used.add(key)
                        break
    return deuterons

def analyze(lattice, cfg):
    N = cfg.hex_rows * cfg.hex_cols * cfg.layers
    quarks = detect_quarks(lattice, cfg)
    protons  = find_baryon_triplets(quarks, cfg, n_up_req=2, n_dn_req=1)
    neutrons = find_baryon_triplets(quarks, cfg, n_up_req=1, n_dn_req=2)
    deuterons = find_deuterons(protons, neutrons, cfg)

    n_lepton = sum(1 for i in range(N) if int(lattice[i])==LEPTON_SEED)
    n_wboson = sum(1 for i in range(N) if int(lattice[i]) in W_BOSON_STATES)

    # Hydrogen: proton with adjacent γ⁰
    p_sites = set(s for d,u1,u2 in protons for s in [d,u1,u2])
    p_shell  = set()
    for s in p_sites:
        r,c,z = site_coords(s, cfg)
        for nb in mesh_neighbours(r,c,z,cfg):
            if nb not in p_sites: p_shell.add(nb)
    n_hydrogen = sum(1 for d,u1,u2 in protons
                     if any(int(lattice[nb])==LEPTON_SEED
                            for s in [d,u1,u2]
                            for nb in mesh_neighbours(*site_coords(s,cfg)[:], cfg)
                            if nb not in p_sites))
    # EM enrichment
    lep_shell = sum(1 for s in p_shell if int(lattice[s])==LEPTON_SEED)
    lep_bg    = n_lepton - lep_shell
    shell_sz  = len(p_shell)
    bg_sz     = N - len(p_sites) - shell_sz
    rs = lep_shell/shell_sz if shell_sz>0 else 0
    rb = lep_bg/bg_sz       if bg_sz>0    else 0
    enrich = rs/rb if rb>1e-9 else (float('inf') if rs>0 else 0.0)

    return {
        'protons': len(protons), 'neutrons': len(neutrons),
        'deuterons': len(deuterons), 'hydrogen': n_hydrogen,
        'leptons': n_lepton, 'wbosons': n_wboson,
        'enrich': enrich,
    }

# ── Main ─────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    N = 12*12*12

    print("GUTOE: Electroweak Extension")
    print("="*72)
    print("Photon = γ¹² (spacelike bivector, SO(3) rotation generator)")
    print("W boson = γ⁰¹ (mixed bivector, timelike×spacelike)")
    print("Electron = γ⁰ (timelike grade-1, Z₃ fixed point)")
    print()
    print("EM:   γ⁰ surfs through grade-2 photon field toward grade-1 quarks")
    print("Weak: γ⁰ + γ¹ → γ⁰¹ (W⁺) + VOID (ν);  W → γ⁰ + γ¹ (decay)")
    print()
    print("Protocol: Phase 1 = t=500 quarks only (stable protons).")
    print("          Phase 2 = inject 50 γ⁰, activate EM+weak, run 500 steps.")
    print("="*72)
    print()

    n_seeds = 5
    n_inject = 50
    phase1_steps = 500
    phase2_steps = 500
    report_every = 50

    # Sweep EM and weak coupling strengths
    configs = [
        ("no EM, no weak (control)",  0.00, 0.00, 0.00),
        ("EM only",                   0.10, 0.00, 0.00),
        ("weak only",                 0.00, 0.02, 0.10),
        ("EM + weak",                 0.10, 0.02, 0.10),
    ]

    for label, em_p, wk_p, wd_p in configs:
        cfg = LatticeConfig(em_prob=em_p, weak_prob=wk_p, w_decay_prob=wd_p)
        print(f"─── {label} ───")
        print(f"    em={em_p}  weak={wk_p}  w_decay={wd_p}")
        print()

        totals = {k: np.zeros(phase2_steps//report_every) for k in
                  ['protons','neutrons','deuterons','hydrogen','leptons','wbosons']}
        enrichments = np.zeros(phase2_steps//report_every)

        for seed_idx in range(n_seeds):
            # Phase 1: quark-only to t=500
            rng = np.random.default_rng(seed_idx * 137 + 7)
            lat = init_lattice(cfg)
            for _ in range(phase1_steps):
                lat = step(lat, rng, cfg)

            # Inject γ⁰ leptons at non-proton sites
            quarks = detect_quarks(lat, cfg)
            proton_sites = set(s for d,u1,u2 in
                               find_baryon_triplets(quarks, cfg, 2, 1)
                               for s in [d,u1,u2])
            candidates = [i for i in range(N) if i not in proton_sites
                          and int(lat[i]) != VOID]  # inject into existing matter
            inject_sites = rng.choice(candidates,
                                      size=min(n_inject, len(candidates)),
                                      replace=False)
            for s in inject_sites:
                lat[s] = LEPTON_SEED

            # Phase 2: with electroweak
            for ti in range(phase2_steps):
                lat = step(lat, rng, cfg)
                if (ti+1) % report_every == 0:
                    ri = (ti+1)//report_every - 1
                    a = analyze(lat, cfg)
                    for k in totals: totals[k][ri] += a[k]
                    enrichments[ri] += a['enrich']

        # Print time series
        col = f"{'t':>5s} | {'p':>5s} {'n':>5s} {'d':>5s} {'H':>5s} {'γ⁰':>5s} {'W':>5s} {'enrich':>7s}"
        print("  " + col)
        print("  " + "-"*len(col))
        for ri in range(phase2_steps//report_every):
            t = phase1_steps + (ri+1)*report_every
            p  = totals['protons'][ri]/n_seeds
            n  = totals['neutrons'][ri]/n_seeds
            d  = totals['deuterons'][ri]/n_seeds
            h  = totals['hydrogen'][ri]/n_seeds
            l  = totals['leptons'][ri]/n_seeds
            w  = totals['wbosons'][ri]/n_seeds
            e  = enrichments[ri]/n_seeds
            print(f"  {t:5d} | {p:5.1f} {n:5.1f} {d:5.1f} {h:5.1f} {l:5.1f} {w:5.1f} {e:7.2f}×")

        # Verdict
        final_h = totals['hydrogen'][-1]/n_seeds
        final_d = totals['deuterons'][-1]/n_seeds
        final_e = enrichments[-1]/n_seeds
        print()
        print(f"  Final: hydrogen={final_h:.1f}  deuterons={final_d:.1f}  "
              f"γ⁰-enrichment={final_e:.2f}×")
        if final_e > 1.5:
            print(f"  HYDROGEN BINDING: γ⁰ preferentially near protons (EM works).")
        elif final_e < 0.5:
            print(f"  LEPTON EXCLUSION: γ⁰ expelled from proton shell.")
        else:
            print(f"  NO PREFERENCE: EM coupling not strong enough or photons absent.")
        if final_d > 0.3:
            print(f"  DEUTERON FORMATION: weak force transmutes protons to neutrons.")
        print()
