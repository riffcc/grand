# TODO — Parked State (Context Switch)

Timestamp: 2026-02-28 21:28:20 UTC

## Current External Benchmark (default path)
- phase_accuracy = 1.000000
- density_mae_g_cm3 = 2.502334
- melting_mae_k = 324.656297
- boiling_mae_k = 760.899932
- ionization_mae_ev = 0.307058
- elements_with_any_red = 75

## Work Just Completed
- Added ratio-locked holdout runner:
  - `crates/gutoe-physics/src/bin/refractory_ratio_holdout.rs`
  - Runs 4d-fit->5d-validate and 5d-fit->4d-validate under exact `g_f/g_v = 12/7`.
- Added env-overridable refractory gains in thermo lane:
  - `crates/gutoe-physics/src/chemical_thermo.rs`
  - `GUTOE_CHEM_REFRACTORY_FUSION_GAIN_Q`
  - `GUTOE_CHEM_REFRACTORY_VAPOR_GAIN_Q`
- Added Lean structural ratio lock (`12/7`) and odd-parity count theorem:
  - `lean/Gutoe/ThermalEntropyClosure.lean`
- Added findings doc:
  - `findings/188-refractory-ratio-holdout-lock.md`

## Key Holdout Result
- Under denominator <= 20 and exact ratio lock `12/7`, bidirectional winner was:
  - `g_f = 12/13`, `g_v = 7/13`
- Structural pair check preserved:
  - `(3/5)/(7/20) = 12/7` exactly.

## Uncommitted Files
- Modified:
  - `crates/gutoe-physics/src/chemical_thermo.rs`
  - `lean/Gutoe/ThermalEntropyClosure.lean`
- Untracked:
  - `crates/gutoe-physics/src/bin/refractory_ratio_holdout.rs`
  - `findings/186-refractory-crystal-cohesive-coupling-pass.md`
  - `findings/187-red-burndown-molecular-floor-and-radius-cap.md`
  - `findings/188-refractory-ratio-holdout-lock.md`

## Resume Plan (after context switch)
1. Decide whether to keep structural defaults (`3/5`, `7/20`) vs holdout-optimal (`12/13`, `7/13`).
2. If keeping structural defaults, keep holdout as validation-only evidence lane.
3. If switching defaults, rerun strict external benchmark + global gate and compare full regressions.
