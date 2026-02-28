# GRAND — Cl(1,3) Candidate Forms For Frozen Triangulation Constants

Date: 2026-02-28

## Context

Frozen constants were stamped in:

- `findings/154-triangulation-constants-freeze.json`
- `findings/154-triangulation-constants-freeze.md`

Targets:

- `p_ratio = 13.688110433760`
- `kappa_geo = 34.697396055505`
- `ew_coeff_required = 8.460487692308`

## Candidate Structural Forms

Evaluated in:

- `crates/gutoe-physics/src/bin/triangulation_clifford_candidates.rs`
- `crates/gutoe-physics/src/bin/triangulation_clifford_candidates_ci_gate.rs`

Candidate formulas:

1. `p = 137/10 - 1/(7*12)`
2. `kappa = (60/11) * (19/3 + 1/36 + 1/(7*13*136))`
3. `ew_coeff = 8 + 6/13 - 1/(7*136)`

Where counts map to Cl(1,3) lane integers:

- `7 = grade2 + 1`
- `12 = total gauge generators`
- `19 = 16 + 3`
- `13 = 16 - 3`
- `136 = T(16)`
- `60/11 = geometric dark/visible ratio lane factor`

## Numerical Match (relative error vs frozen targets)

- `p_rel = -1.110e-6`
- `kappa_rel = +4.144e-7`
- `ew_rel = +4.126e-8`

## Gate

`triangulation_clifford_candidates_ci_gate` added and wired into global gate run.

Windows:

- `|p_rel| <= 2e-6`
- `|kappa_rel| <= 1e-6`
- `|ew_rel| <= 1e-7`

Current status: pass.

## Honesty

These are high-fidelity candidate reconstructions, not final Lean-proven derivations yet.

They are now pinned as explicit hypotheses with CI protection while formal derivation work proceeds.
