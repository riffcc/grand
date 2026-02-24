# 005 — 2.73x Lepton Enrichment Around Protons

## Objective

Document and reproduce the reported lepton-shell enrichment signal in the proton-coupled EM simulation.

## Reproduction run

Command:

`python3 -u gutoe_em_hydrogen.py`

Configuration (from script defaults):

- Lattice: `12x12x12`
- Seeds: `10`
- Phase 1: `t=0..150` (quarks only)
- Phase 2: `500` steps with EM active
- Injection: `20` gamma0 leptons at `t=150`
- Report interval: `50` steps

## Summary output (averaged over 10 seeds)

At `t=350`:

- protons: `12.5`
- leptons: `20.0`
- hydrogen count: `0.90`
- enrichment: **`2.73x`**
- `Delta phi (lep-all)`: `+0.1124`

Script verdict:

- `HYDROGEN: YES — peak enrichment 2.73x at t=350 (EM binding confirmed).`
- `phi-TRACKING: YES — leptons at Delta phi=+0.113 above lattice mean.`

## Interpretation

This run reproduces the 2.73x peak enrichment as a transient feature in the phase-2 window, with positive phi-bias indicating preferential occupancy near proton-associated positive potential regions.
