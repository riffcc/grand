#!/usr/bin/env python3
"""
GUTOE: 20 quick experiments.
Running #1, #4, #5, #16, #17, #18 first (highest impact per minute).

Copyright (C) 2026 Riff Labs, AGPL-3.0-or-later
"""

import numpy as np
from enum import IntEnum
from dataclasses import dataclass, field
from collections import Counter
import sys

# ── Core (copied from gutoe_spontaneous_uud.py) ─────────────────────────────

class TriState(IntEnum):
    VOID = 0; COSINE = 1; SINE = 2; TANGENT = 3

def cycle(s): return [TriState.VOID, TriState.TANGENT, TriState.COSINE, TriState.SINE][s]
def veracity_val(a, b):
    if a == 0 or b == 0: return 0.0
    if a == b: return 1.0
    if frozenset({a,b}) == frozenset({2,1}): return np.sqrt(3)/2
    return 0.5

@dataclass
class Cfg:
    hex_rows: int = 12; hex_cols: int = 12; layers: int = 8
    reduced_mesh: bool = True
    diff_prob: float = 0.02; cycle_prob: float = 0.05
    align: float = 0.15; threshold: float = 0.6
    use_square: bool = False  # for experiment #3

def idx(r,c,z,cfg): return ((z%cfg.layers)*cfg.hex_rows+(r%cfg.hex_rows))*cfg.hex_cols+(c%cfg.hex_cols)

def hex_nbrs(r,c,cfg):
    if r%2==0: offs=[(-1,0),(-1,1),(0,-1),(0,1),(1,0),(1,1)]
    else: offs=[(-1,-1),(-1,0),(0,-1),(0,1),(1,-1),(1,0)]
    return [((r+dr)%cfg.hex_rows,(c+dc)%cfg.hex_cols) for dr,dc in offs]

def square_nbrs(r,c,cfg):
    return [((r-1)%cfg.hex_rows,c),(r,(c+1)%cfg.hex_cols),
            ((r+1)%cfg.hex_rows,c),(r,(c-1)%cfg.hex_cols)]

def planar_nbrs(r,c,cfg):
    return square_nbrs(r,c,cfg) if cfg.use_square else hex_nbrs(r,c,cfg)

def mesh_nbrs(r,c,z,cfg):
    nbrs = [(idx(nr,nc,z,cfg),'p') for nr,nc in planar_nbrs(r,c,cfg)]
    if not cfg.reduced_mesh:
        for nr,nc in planar_nbrs(r,c,cfg):
            nbrs.append((idx(nr,nc,(z+1)%cfg.layers,cfg),'ah'))
            nbrs.append((idx(nr,nc,(z-1)%cfg.layers,cfg),'bh'))
    nbrs.append((idx(r,c,(z+1)%cfg.layers,cfg),'v'))
    nbrs.append((idx(r,c,(z-1)%cfg.layers,cfg),'v'))
    return nbrs

def decompose(site,cfg):
    z=site//(cfg.hex_rows*cfg.hex_cols); rem=site%(cfg.hex_rows*cfg.hex_cols)
    return rem//cfg.hex_cols, rem%cfg.hex_cols, z

def fields(lat,site,cfg):
    state=TriState(lat[site])
    if state==0: return 0.,0.,0.
    r,c,z=decompose(site,cfg); nbrs=mesh_nbrs(r,c,z,cfg)
    tv=0.; sc=0.; tc=0.; gs=0.; ns=0; nt=0
    for ni,ct in nbrs:
        v=veracity_val(state,TriState(lat[ni])); tv+=v; d=1.-v
        if ct=='v': tc+=d; nt+=1
        else: sc+=d; ns+=1
        gs+=d
    n=len(nbrs); av=tv/n
    if ns>0: sc/=ns
    if nt>0: tc/=nt
    curv=max(abs(sc-tc),sc); grad=gs/n
    return av, curv, grad

def step(lat,rng,cfg):
    new=lat.copy(); N=len(lat)
    for s in range(N):
        st=TriState(lat[s])
        if st==0:
            if rng.random()<cfg.diff_prob: new[s]=2; continue
            r,c,z=decompose(s,cfg); nbrs=mesh_nbrs(r,c,z,cfg)
            act=sum(1 for ni,_ in nbrs if lat[ni]!=0); tn=len(nbrs)
            if act>=max(2,tn//4) and rng.random()<act/tn*0.4: new[s]=2
        else:
            if rng.random()<cfg.cycle_prob: new[s]=int(cycle(st))
            elif cfg.align>0 and rng.random()<cfg.align:
                r,c,z=decompose(s,cfg); nbrs=mesh_nbrs(r,c,z,cfg)
                ns=[TriState(lat[ni]) for ni,_ in nbrs if lat[ni]!=0]
                if ns: new[s]=int(Counter(ns).most_common(1)[0][0])
    return new

@dataclass
class Quark:
    site:int; qtype:str; bc:float; v:float; curv:float; grad:float

def detect(lat,cfg):
    qs=[]
    for s in range(len(lat)):
        if lat[s]==0: continue
        v,cu,gr=fields(lat,s,cfg); bc=v/(1+gr)
        if bc>=cfg.threshold:
            qs.append(Quark(s,"UP" if v>cu else "DOWN",bc,v,cu,gr))
    return qs

def find_baryons(qs,cfg):
    """Find protons (uud) and neutrons (udd)."""
    qset={q.site:q for q in qs}; ncache={}
    for q in qs:
        r,c,z=decompose(q.site,cfg)
        ncache[q.site]={ni for ni,_ in mesh_nbrs(r,c,z,cfg)}
    protons=0; neutrons=0; used=set()
    for q in qs:
        if q.site in used: continue
        adj=[qset[ni] for ni in ncache[q.site] if ni in qset and ni not in used]
        if len(adj)<2: continue
        for i in range(len(adj)):
            found=False
            for j in range(i+1,len(adj)):
                trio=[q,adj[i],adj[j]]
                p1,p2=adj[i].site,adj[j].site
                if p1 not in ncache.get(p2,set()): continue
                types=Counter(t.qtype for t in trio)
                if types.get("UP",0)==2 and types.get("DOWN",0)==1:
                    protons+=1; used.update(t.site for t in trio); found=True; break
                elif types.get("UP",0)==1 and types.get("DOWN",0)==2:
                    neutrons+=1; used.update(t.site for t in trio); found=True; break
            if found: break
    return protons, neutrons

def run_seeds(cfg, steps=200, n_seeds=20, label=""):
    peaks_p=[]; peaks_n=[]; peak_ts=[]
    for s in range(n_seeds):
        rng=np.random.default_rng(s*137+7)
        lat=np.zeros(cfg.hex_rows*cfg.hex_cols*cfg.layers,dtype=np.int8)
        pp=0; pn=0; pt=0
        for t in range(1,steps+1):
            lat=step(lat,rng,cfg)
            qs=detect(lat,cfg)
            p,n=find_baryons(qs,cfg)
            if p+n>pp+pn: pp=p; pn=n; pt=t
        peaks_p.append(pp); peaks_n.append(pn); peak_ts.append(pt)
    print(f"  {label}: protons={np.mean(peaks_p):.1f}±{np.std(peaks_p):.1f}  "
          f"neutrons={np.mean(peaks_n):.1f}±{np.std(peaks_n):.1f}  "
          f"seeds_with_baryons={sum(1 for p,n in zip(peaks_p,peaks_n) if p+n>0)}/{n_seeds}")
    return peaks_p, peaks_n, peak_ts

# ══════════════════════════════════════════════════════════════════════════════
print("="*70)
print("GUTOE: 20 EXPERIMENTS")
print("="*70)

# ── #1: Kill alignment ───────────────────────────────────────────────────────
print("\n#1: KILL ALIGNMENT (align=0 vs align=0.15)")
print("-"*50)
run_seeds(Cfg(align=0.15), label="align=0.15 (baseline)")
run_seeds(Cfg(align=0.0),  label="align=0.00 (no align)")

# ── #3: Square mesh vs hex ───────────────────────────────────────────────────
print("\n#3: SQUARE MESH vs HEX MESH")
print("-"*50)
run_seeds(Cfg(use_square=False), label="hex (6 planar)")
run_seeds(Cfg(use_square=True),  label="square (4 planar)")

# ── #5: Count neutrons ───────────────────────────────────────────────────────
print("\n#5: NEUTRONS (udd) — already tracked above!")
print("-"*50)
print("  (neutron counts shown in all experiments above)")

# ── #2: Sweep lattice layers ─────────────────────────────────────────────────
print("\n#2: SWEEP VERTICAL LAYERS")
print("-"*50)
for L in [2, 4, 8, 16]:
    run_seeds(Cfg(layers=L), n_seeds=10, label=f"layers={L:2d}")

# ── #4: Percolation threshold ────────────────────────────────────────────────
print("\n#4: PERCOLATION THRESHOLD (active fraction at first proton)")
print("-"*50)
first_fracs = []
cfg4 = Cfg()
for s in range(50):
    rng=np.random.default_rng(s*137+7)
    lat=np.zeros(cfg4.hex_rows*cfg4.hex_cols*cfg4.layers,dtype=np.int8)
    N=len(lat)
    for t in range(1,200):
        lat=step(lat,rng,cfg4)
        qs=detect(lat,cfg4)
        p,n=find_baryons(qs,cfg4)
        if p>0:
            frac=sum(1 for x in lat if x!=0)/N
            first_fracs.append(frac)
            break
if first_fracs:
    print(f"  First proton at active fraction: {np.mean(first_fracs):.3f} ± {np.std(first_fracs):.3f}")
    print(f"  Min: {min(first_fracs):.3f}  Max: {max(first_fracs):.3f}")
    print(f"  Known 2D hex percolation threshold: 0.6527...")
    print(f"  Seeds that produced protons: {len(first_fracs)}/50")

# ── #16: Mass ratio UP/DOWN ──────────────────────────────────────────────────
print("\n#16: MASS RATIO from field geometry")
print("-"*50)
cfg16 = Cfg()
all_up_mass = []; all_dn_mass = []
for s in range(20):
    rng=np.random.default_rng(s*137+7)
    lat=np.zeros(cfg16.hex_rows*cfg16.hex_cols*cfg16.layers,dtype=np.int8)
    for t in range(1,200):
        lat=step(lat,rng,cfg16)
        qs=detect(lat,cfg16)
        if any(q.qtype=="DOWN" for q in qs):
            for q in qs:
                # mass ∝ v * curv * grad / l * λ_QG²
                # l = 1 (lattice units), λ_QG = 1/12
                mass = q.v * max(q.curv, 0.01) * max(q.grad, 0.01) * (1/12)**2
                if q.qtype == "UP": all_up_mass.append(mass)
                else: all_dn_mass.append(mass)
if all_up_mass and all_dn_mass:
    ratio = np.mean(all_up_mass) / np.mean(all_dn_mass)
    print(f"  Mean UP mass:   {np.mean(all_up_mass):.6f} ± {np.std(all_up_mass):.6f}")
    print(f"  Mean DOWN mass: {np.mean(all_dn_mass):.6f} ± {np.std(all_dn_mass):.6f}")
    print(f"  Ratio m_UP/m_DOWN: {ratio:.3f}")
    print(f"  Observed m_u/m_d:  0.47")
else:
    print(f"  Not enough DOWN quarks for mass ratio (UP={len(all_up_mass)}, DN={len(all_dn_mass)})")

# ── #17: Charge conservation ─────────────────────────────────────────────────
print("\n#17: CHARGE CONSERVATION")
print("-"*50)
cfg17 = Cfg()
rng17 = np.random.default_rng(42)
lat17 = np.zeros(cfg17.hex_rows*cfg17.hex_cols*cfg17.layers, dtype=np.int8)
charges = []
for t in range(1, 101):
    lat17 = step(lat17, rng17, cfg17)
    qs = detect(lat17, cfg17)
    # UP = +2/3, DOWN = -1/3
    total_charge = sum(2/3 if q.qtype=="UP" else -1/3 for q in qs)
    charges.append(total_charge)
    if t in [1,5,10,20,50,100]:
        u=sum(1 for q in qs if q.qtype=="UP")
        d=sum(1 for q in qs if q.qtype=="DOWN")
        print(f"  t={t:3d}: UP={u:4d} DN={d:3d} total_charge={total_charge:+8.1f}")
print(f"  Charge range: [{min(charges):+.1f}, {max(charges):+.1f}]")
print(f"  Conserved? {'NO — charge grows with quark count' if max(charges)-min(charges)>10 else 'YES'}")

# ── #18: Energy conservation ─────────────────────────────────────────────────
print("\n#18: ENERGY (sum of veracity²) OVER TIME")
print("-"*50)
cfg18 = Cfg()
rng18 = np.random.default_rng(42)
lat18 = np.zeros(cfg18.hex_rows*cfg18.hex_cols*cfg18.layers, dtype=np.int8)
for t in range(1, 101):
    lat18 = step(lat18, rng18, cfg18)
    if t in [1,5,10,20,50,100]:
        energy = 0.0
        N18 = len(lat18)
        for s in range(N18):
            if lat18[s] == 0: continue
            v, _, _ = fields(lat18, s, cfg18)
            energy += v * v
        print(f"  t={t:3d}: E={energy:8.1f}  (active={sum(1 for x in lat18 if x!=0)})")

# ── #15: Lattice big bang ────────────────────────────────────────────────────
print("\n#15: LATTICE BIG BANG (single SINE, no thermal noise)")
print("-"*50)
cfg15 = Cfg(diff_prob=0.0, align=0.0)  # no thermal noise, no alignment
rng15 = np.random.default_rng(42)
lat15 = np.zeros(cfg15.hex_rows*cfg15.hex_cols*cfg15.layers, dtype=np.int8)
# Plant one SINE at the centre
centre = idx(cfg15.hex_rows//2, cfg15.hex_cols//2, cfg15.layers//2, cfg15)
lat15[centre] = int(TriState.SINE)
for t in range(1, 101):
    lat15 = step(lat15, rng15, cfg15)
    if t in [1,2,3,5,10,20,50,100]:
        active = sum(1 for x in lat15 if x != 0)
        print(f"  t={t:3d}: active={active:5d}/{len(lat15)}")

# ── #19: Reverse time ────────────────────────────────────────────────────────
print("\n#19: TIME REVERSAL (reverse update order)")
print("-"*50)
cfg19 = Cfg()
# Forward run
rng_f = np.random.default_rng(42)
lat_f = np.zeros(cfg19.hex_rows*cfg19.hex_cols*cfg19.layers, dtype=np.int8)
states_forward = []
for t in range(1, 51):
    lat_f = step(lat_f, rng_f, cfg19)
    active = sum(1 for x in lat_f if x != 0)
    states_forward.append(active)

# "Reverse": start from equilibrium, reverse cycle direction
def reverse_cycle(s):
    if s == TriState.COSINE: return TriState.SINE
    if s == TriState.TANGENT: return TriState.COSINE
    if s == TriState.SINE: return TriState.TANGENT
    return TriState.VOID

def step_reverse(lat, rng, cfg):
    new = lat.copy(); N = len(lat)
    for s in range(N-1, -1, -1):  # reverse order
        st = TriState(lat[s])
        if st == 0:
            if rng.random() < cfg.diff_prob: new[s] = 2; continue
            r,c,z = decompose(s, cfg); nbrs = mesh_nbrs(r,c,z,cfg)
            act = sum(1 for ni,_ in nbrs if lat[ni] != 0); tn = len(nbrs)
            if act >= max(2, tn//4) and rng.random() < act/tn*0.4: new[s] = 2
        else:
            if rng.random() < cfg.cycle_prob: new[s] = int(reverse_cycle(st))
            elif cfg.align > 0 and rng.random() < cfg.align:
                r,c,z = decompose(s, cfg); nbrs = mesh_nbrs(r,c,z,cfg)
                ns = [TriState(lat[ni]) for ni,_ in nbrs if lat[ni] != 0]
                if ns: new[s] = int(Counter(ns).most_common(1)[0][0])
    return new

rng_r = np.random.default_rng(99)
lat_r = lat_f.copy()  # start from final state
states_reverse = []
for t in range(1, 51):
    lat_r = step_reverse(lat_r, rng_r, cfg19)
    active = sum(1 for x in lat_r if x != 0)
    states_reverse.append(active)

print(f"  Forward:  active went {states_forward[0]:4d} → {states_forward[-1]:4d}")
print(f"  Reverse:  active went {states_reverse[0]:4d} → {states_reverse[-1]:4d}")
print(f"  Time-reversible? {'YES' if states_reverse[-1] < states_forward[-1]*0.5 else 'NO — arrow of time exists'}")

# ── #12: HexState↔fermion counting ──────────────────────────────────────────
print("\n#12: HEXSTATE ↔ STANDARD MODEL FERMION COUNTING")
print("-"*50)
print("  12 HexStates (from HexStates.lean):")
print("  Angles: 0°, 30°, 60°, 90°, 120°, 150°, 180°, 210°, 240°, 270°, 300°, 330°")
print()
print("  First generation SM fermions (12):")
print("  Quarks (6):  u_r u_g u_b d_r d_g d_b")
print("  Leptons (6): e⁻ νₑ (and their antiparticles, or 3 colors × 2)")
print()
print("  Z₆ symmetry of HexStates → color SU(3)?")
print("  If rotation by 120° = color rotation:")
print("    0°,120°,240° = three colors of one quark")
print("    60°,180°,300° = three colors of another quark")
print("  Z₂ (negation) = particle/antiparticle?")
print("    negation flips 0°↔180°, 60°↔240°, 120°↔300°")
print()
charge_map = {0: 0, 60: 1/3, 120: 2/3, 180: 1, 240: 2/3, 300: 1/3}
print("  If charge = angle/180° mod 1:")
for angle, q in sorted(charge_map.items()):
    print(f"    {angle:3d}° → charge {q:.2f}")
print(f"  Proton (3 quarks at 120°): 3 × 2/3 = {3*2/3:.1f}... too high!")
print(f"  Need: charge = angle/270° for u quarks? Or different mapping.")

# ── #20: Fine structure constant ─────────────────────────────────────────────
print("\n#20: FINE STRUCTURE CONSTANT HUNTING")
print("-"*50)
sqrt3_2 = np.sqrt(3)/2
candidates = [
    ("1/(4π·12·√3/2)", 1/(4*np.pi*12*sqrt3_2)),
    ("1/(4π·12)", 1/(4*np.pi*12)),
    ("1/(2π·12·√3)", 1/(2*np.pi*12*np.sqrt(3))),
    ("1/(12²·√3/2)", 1/(144*sqrt3_2)),
    ("1/(12·(12-1/12))", 1/(12*(12-1/12))),
    ("(1/12)²·6·√3/2", (1/12)**2*6*sqrt3_2),
    ("1/(4π·(12+1/12))", 1/(4*np.pi*(12+1/12))),
    ("√3/(4π·12²)", np.sqrt(3)/(4*np.pi*144)),
]
alpha_real = 1/137.036
print(f"  α = 1/137.036 = {alpha_real:.6f}")
print()
for name, val in sorted(candidates, key=lambda x: abs(x[1]-alpha_real)):
    pct = (val/alpha_real - 1)*100
    print(f"  {name:30s} = {val:.6f}  (1/{1/val:.1f})  [{pct:+.1f}%]")

print()
print("="*70)
print("ALL EXPERIMENTS COMPLETE")
print("="*70)
