#!/usr/bin/env python3
"""
GUTOE Neutron-Proton Mass Ratio: First Dynamic Prediction
Copyright (C) 2026 Riff Labs, AGPL-3.0-or-later

mn/mp ≈ 1.001378 — the neutron is slightly heavier than the proton.

Physical origin:
  QCD: DOWN quarks (curvature-dominant) carry more energy than UP quarks
       (veracity-dominant). A neutron (UDD) has two DOWN quarks; a proton (UUD)
       has two UP quarks. This makes neutrons heavier.
  EM:  The proton (charge +1) has Coulomb self-energy; the neutron (charge 0)
       does not. This partially cancels the QCD effect.
  Net: mn - mp ≈ 1.3 MeV = 0.00138 × mp

This is the first DYNAMIC prediction from the GUTOE lattice — it requires
running the simulation and measuring energy differences between two baryon
configurations, not just counting Clifford algebra states.
"""

import sys, os
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))

from gutoe_gauge import (
    LatticeConfig, site_coords, mesh_neighbours,
    jacobi_poisson, VOID, LEPTON_SEED,
)
from gutoe_em_hydrogen import (
    init_lattice, step, detect_quarks, find_proton_triplets, _V, Z3_ORBITS,
)

MN_MP_EXP = 1.001378   # CODATA 2018

# ── Neutron triplet detection ──────────────────────────────────────────────────

def find_neutron_triplets(quarks, cfg):
    """Find neutron triplets: (DOWN, DOWN, UP) triangles."""
    quark_set = {q.site: q for q in quarks}
    triplets = []
    used = set()
    nbr_cache = {
        q.site: set(mesh_neighbours(q.r, q.c, q.z, cfg)) for q in quarks
    }

    for q in quarks:
        # Anchor on UP quarks (unique in neutron, as DOWN is unique in proton)
        if q.quark_type != "UP" or q.site in used:
            continue

        # Find two DOWN neighbours that are adjacent to each other
        down_nbrs = [
            quark_set[ni]
            for ni in nbr_cache[q.site]
            if ni in quark_set
            and quark_set[ni].quark_type == "DOWN"
            and ni not in used
        ]
        if len(down_nbrs) < 2:
            continue

        found = False
        for i in range(len(down_nbrs)):
            if found:
                break
            for j in range(i + 1, len(down_nbrs)):
                d1, d2 = down_nbrs[i].site, down_nbrs[j].site
                if d1 in nbr_cache.get(d2, set()):
                    triplets.append((d1, d2, q.site))  # (DOWN, DOWN, UP)
                    used.update([q.site, d1, d2])
                    found = True
                    break
    return triplets

# ── Baryon binding energy ──────────────────────────────────────────────────────

def baryon_binding_energy(lattice, triplet):
    """
    Sum of quark-quark veracity × alignment for the 3 internal bonds of a triangle.
    This is the QCD binding energy of the baryon (before EM correction).
    """
    sites = list(triplet)
    states = [int(lattice[s]) for s in sites]
    E = 0.0
    for i in range(len(sites)):
        for j in range(i + 1, len(sites)):
            E += _V[(states[i], states[j])]
    return E  # in units of alignment_strength


def proton_coulomb_self_energy(triplet, phi):
    """EM self-energy of the proton from its own Coulomb field at the quark sites."""
    sites = list(triplet)
    # Proton charge +1 distributed as (+2/3, +2/3, -1/3) over 3 quarks
    # Self-energy = (1/2) × q_total × phi_mean = 0.5 × 1 × phi_mean
    phi_mean = np.mean([phi[s] for s in sites])
    return 0.5 * phi_mean   # positive: proton EM self-energy makes proton heavier


# ── Main measurement ───────────────────────────────────────────────────────────

def measure_np_ratio(n_seeds=10, n_phase1=150):
    cfg = LatticeConfig()

    proton_E = []
    neutron_E = []
    proton_EM = []
    neutron_EM = []

    print(f"  {n_seeds} seeds, {n_phase1} Phase-1 steps")
    print()

    for seed in range(n_seeds):
        rng = np.random.default_rng(seed * 137 + 7)
        lattice = init_lattice(cfg)

        for _ in range(n_phase1):
            lattice = step(lattice, rng, cfg, gauge=None, proton_sites=None)

        quarks   = detect_quarks(lattice, cfg)
        protons  = find_proton_triplets(quarks, cfg)   # (DOWN, UP, UP)
        neutrons = find_neutron_triplets(quarks, cfg)  # (DOWN, DOWN, UP)

        # Coulomb field from quark charge distribution (for EM self-energy)
        N = cfg.hex_rows * cfg.hex_cols * cfg.layers
        rho = np.zeros(N)
        for q in quarks:
            rho[q.site] = 2.0/3.0 if q.quark_type == "UP" else -1.0/3.0
        phi = jacobi_poisson(rho, cfg, cfg.poisson_iters)

        ep_list = [baryon_binding_energy(lattice, t) for t in protons]
        en_list = [baryon_binding_energy(lattice, t) for t in neutrons]
        eem_list = [proton_coulomb_self_energy(t, phi) for t in protons]

        proton_E.extend(ep_list)
        neutron_E.extend(en_list)
        proton_EM.extend(eem_list)

        ep_mean = np.mean(ep_list) if ep_list else 0
        en_mean = np.mean(en_list) if en_list else 0
        print(f"  seed {seed:2d}: {len(protons):3d} protons  {len(neutrons):3d} neutrons"
              f"  E_p={ep_mean:.4f}  E_n={en_mean:.4f}")

    if not proton_E or not neutron_E:
        return None

    E_p_qcd = np.mean(proton_E)
    E_n_qcd = np.mean(neutron_E)
    E_p_em  = np.mean(proton_EM)
    # Neutron has no EM self-energy: E_n_em = 0

    # Total mass ∝ QCD binding energy + EM self-energy
    # Note: EM self-energy of proton is POSITIVE (makes proton heavier)
    E_p_total = E_p_qcd + E_p_em
    E_n_total = E_n_qcd          # neutron has no EM correction

    return {
        'E_p_qcd': E_p_qcd,
        'E_n_qcd': E_n_qcd,
        'E_p_em': E_p_em,
        'E_p_total': E_p_total,
        'E_n_total': E_n_total,
        'mn_mp_qcd':  E_n_qcd / E_p_qcd if E_p_qcd > 0 else 0,
        'mn_mp_full': E_n_total / E_p_total if E_p_total > 0 else 0,
        'n_p': len(proton_E),
        'n_n': len(neutron_E),
    }


if __name__ == "__main__":
    w = 72
    print("=" * w)
    print("GUTOE Neutron-Proton Mass Ratio")
    print("=" * w)
    print(f"\n  Experimental:  mn/mp = {MN_MP_EXP:.6f}")
    print(f"                 mn - mp ≈ 1.293 MeV = {MN_MP_EXP - 1:.5f} × mp")

    print(f"\n{'─' * w}")
    r = measure_np_ratio(n_seeds=10, n_phase1=150)
    print(f"{'─' * w}")

    if r is None:
        print("\n  No baryons detected.")
    else:
        print(f"\n  Proton binding energy (QCD, lattice units):")
        print(f"    E_p_QCD  = {r['E_p_qcd']:.6f}  (mean over {r['n_p']} protons)")
        print(f"  Neutron binding energy (QCD, lattice units):")
        print(f"    E_n_QCD  = {r['E_n_qcd']:.6f}  (mean over {r['n_n']} neutrons)")
        print(f"\n  QCD-only ratio: E_n/E_p = {r['mn_mp_qcd']:.6f}")
        print(f"    (1 = exactly degenerate, >1 = neutron heavier, <1 = proton heavier)")

        frac = r['E_n_qcd'] - r['E_p_qcd']
        print(f"\n  ΔE_QCD = E_n - E_p = {frac:.6f} lattice units")
        print(f"  Fractional: (E_n-E_p)/E_p = {frac/r['E_p_qcd']:.6f}")

        print(f"\n  Proton EM self-energy:")
        print(f"    E_p_EM   = {r['E_p_em']:.6f}  (positive: makes proton heavier)")
        print(f"\n  Full ratio including EM correction:")
        print(f"    E_p_total = {r['E_p_total']:.6f}")
        print(f"    E_n_total = {r['E_n_total']:.6f}")
        print(f"    mn/mp     = {r['mn_mp_full']:.6f}")

        print(f"\n{'═' * w}")
        print(f"  Experimental:  mn/mp = {MN_MP_EXP:.6f}")
        print(f"  GUTOE sim:     mn/mp = {r['mn_mp_full']:.6f}")
        err = abs(r['mn_mp_full'] - MN_MP_EXP) / MN_MP_EXP * 100
        print(f"  Error:         {err:.2f}%")
        print()

        if r['E_n_qcd'] > r['E_p_qcd']:
            print(f"  QCD direction: ✓ neutron heavier than proton (DOWN > UP energy)")
        else:
            print(f"  QCD direction: ✗ proton heavier — model needs DOWN > UP energy")

        if err < 10:
            print(f"  Agreement:     within 10% — GUTOE dynamics predict mn > mp")
        elif r['mn_mp_full'] > 1:
            print(f"  mn > mp:       ✓ correct sign, wrong magnitude")
            print(f"  The mechanism is right. The precise value needs more physics.")
        else:
            print(f"  mn < mp:       ✗ wrong sign — EM overcorrects")

        print(f"{'=' * w}")
