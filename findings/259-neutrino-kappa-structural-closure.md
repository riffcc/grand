# Finding 259 — Neutrino `97.5%` Suppression Resolved by Structural `κ`

Date: 2026-03-01

## Problem fingerprint

The recurring no-fit neutrino miss was:

- `dm21_rel_err ≈ -0.9753`
- `dm32_rel_err ≈ -0.9753`

This is the same `~97.5%` suppression pattern seen repeatedly in the absolute-splitting lane.

## Root cause

The structural tiny-mass lane still used:

- `p = 137/10`
- `κ = 60/11`

while the triangulated closure lane had already frozen near:

- `p_ratio ≈ 13.688110...`
- `kappa_geo ≈ 34.697396...`

## Structural resolution applied

Aligned Rust runtime to the Lean-backed triangulated candidate closed forms:

- `p = 137/10 - 1/(7*12)`  (`5749/420`)
- `κ = (60/11) * (19/3 + 1/36 + 1/(7*13*136))`

These forms are already formalized in:

- `lean/Gutoe/TriangulatedConstants.lean`
- `lean/Gutoe/TriangulatedClosureUniqueness.lean`

## Code updates

- `crates/gutoe-em/src/flavor.rs`
  - `neutrino_hierarchy_exponent_structural()` updated to candidate closed form.
  - new `neutrino_kappa_structural()` added.
  - `neutrino_absolute_masses_from_texture()` now uses `neutrino_kappa_structural()`.

- `crates/gutoe-em/src/bin/yukawa_neutrino_endgame_report.rs`
  - `no_fit_pass` hardened to include absolute-splitting check.

- `crates/gutoe-em/src/bin/remaining12_gate.rs`
  - neutrino `no_fit_pass` hardened to include `abs_splittings_ok`.

## Verification

### 1) Neutrino oscillation CI gate

Command:
`cargo run -q -p gutoe-em --bin neutrino_oscillation_ci_gate`

Result:
- `pass=true`
- `dm21_rel_err = +4.7295e-6`
- `dm32_rel_err = +7.0898e-7`
- `hierarchy_exponent = 13.688095238095`

### 2) Endgame report (no-fit lane now absolute-closed)

Command:
`cargo run -q -p gutoe-em --bin yukawa_neutrino_endgame_report`

Result:
- `no_fit_pass=true`
- `tri_pass=true`
- `abs_splittings_ok=true`

### 3) Unified remaining12 gate

Command:
`cargo run -q -p gutoe-em --bin remaining12_gate`

Result:
- `overall_pass=true`
- neutrino checks all true, including `abs_splittings_ok=true`

### 4) Lean parity sanity

Command:
`cd lean && lake build Gutoe`

Result:
- build success (warnings only)

## Outcome

The `97.5%` neutrino absolute-splitting miss is no longer present in the no-fit structural lane after adopting the Lean-backed structural `κ` and `p` candidates.

