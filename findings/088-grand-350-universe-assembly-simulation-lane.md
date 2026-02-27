# Finding 088: GRAND-350 End-to-End Universe Assembly Simulation Lane

Date: 2026-02-27  
Status: GRAND-350 initial executable lane landed

## Goal

Assemble the already-closed physics lanes into one executable simulation that
runs a coherent universe history and enforces cross-era consistency gates.

## What landed

### New Rust module

- `crates/gutoe-physics/src/universe.rs`

This module composes existing lanes:

- inflation (`GRAND-347`)
- baryogenesis (`GRAND-348`)
- BBN (`GRAND-349`)
- dark matter falsification (`GRAND-346`, unified branch)
- FRW expansion from `Λ_full` and structural matter budget

### New binaries

- `crates/gutoe-physics/src/bin/universe_sim.rs`
- `crates/gutoe-physics/src/bin/universe_ci_gate.rs`

### Wiring

- `crates/gutoe-physics/src/lib.rs` now exports `universe`

## Model assembly used

Global budget:

- `Ω_b0 = 0.0493`
- `Ω_dm0 = Ω_b0 * (60/11)`
- `Ω_m0 = Ω_b0 + Ω_dm0`
- `Ω_r0 = 9e-5`
- `Ω_k0 = 0`
- `Ω_Λ0 = 1 - Ω_m0 - Ω_r0 - Ω_k0`

Expansion source:

- `Λ = lambda_cosmological_full_candidate()`
- `H0` solved from `Ω_Λ = Λ c² / (3 H0²)`

Time evolution:

- FRW integral with radiation+matter+curvature+Λ terms
- explicit epoch extraction (baryogenesis proxy, BBN, recombination,
  matter-Λ equality, today)
- dense history table (257 rows) from `z=0` to `z=1e9`

## Quantitative output

From `/tmp/bh_renders/universe_sim_report.txt`:

- `H0 = 68.0163 km/s/Mpc`
- `age = 13.6269 Gyr`
- `BBN age ≈ 177.61 s`
- `recombination age ≈ 368.68 kyr`
- full pipeline gate: pass

## Gate definition

`passes_all` requires all of:

- inflation gate pass
- baryogenesis gate pass
- BBN gate pass
- unified dark matter gate pass
- `H0` relative error within 3%
- age in `[13.0, 14.5]` Gyr
- recombination age in `[200, 500]` kyr
- BBN timing in `[10, 2000]` s

## Build/runtime sanity

- `cargo check -p gutoe-physics --bin universe_sim --bin universe_ci_gate` ✅
- `cargo test -p gutoe-physics universe -- --nocapture` ✅
- `cargo run -q -p gutoe-physics --bin universe_sim` ✅
- `cargo run -q -p gutoe-physics --bin universe_ci_gate` ✅

Artifacts written:

- `/tmp/bh_renders/universe_sim_report.txt`
- `/tmp/bh_renders/universe_sim_report.json`
- `/tmp/bh_renders/universe_ci_gate.json`

## Honest boundary

This is an assembled FRW+checkpoint simulation lane, not yet a full coupled
Boltzmann + perturbation solver.

Still open for hard-mode closure:

1. add perturbation growth transfer (`P(k,z)` lane) and BAO/CMB transfer checks.
2. couple BBN/recombination to explicit reaction/opacity network evolution.
3. add uncertainty propagation across all upstream gates (not only pass/fail
   windows).
