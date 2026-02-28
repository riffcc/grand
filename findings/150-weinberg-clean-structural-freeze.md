# 150 - Weinberg Angle Clean Structural Freeze

## Frozen claim
Use the minimal structural expression only:

`sin²(theta_W)(M_Z) = 3/13 + α² * (d/2)`

with:
- `α = 1/137`
- `d = 16`

No flavor anisotropy term, no observed `sin²(theta_W)` input.

## Numerical value
- `sin²(theta_W)_structural = 3/13 = 0.230769230769`
- `Δsin²_structural = α² * 8 = 4.262347487879e-4`
- `sin²(theta_W)(M_Z)_structural = 0.231195465518`

Target comparison (`0.23122`):
- signed error: `-2.453448198134e-5`
- absolute error: `2.453448198134e-5`
- relative error: `0.01061%`

## Integrity statement
This is the defensible zero-free-parameter lane.

The neutrino flavor closure remains a separate open target:
- splitting ratio model: `0.297103211006`
- splitting ratio target: `32.576361221780`
- relative error: `-99.09%`

The remaining `2.45e-5` EW gap is explicitly reserved for future flavor closure.
Any add-on term must close both targets simultaneously in one mechanism.

## Artifacts
- `/tmp/bh_renders/ew_flavor_coupled_ci_gate.json`
- `crates/gutoe-physics/src/dynamics_map.rs`
- `crates/gutoe-physics/src/bin/ew_flavor_coupled_ci_gate.rs`
