# Finding 042 — Uniform MB Thermal-Average Witness for pp Rate

## Scope
Advance GRAND-282 from a fixed 3-point witness to a general `(n+1)`-sample uniform Maxwell-Boltzmann thermal average over an energy ladder.

## Lean
- Added `ppThermalAverageUniform` in `lean/Gutoe/StellarFusion.lean`.
- Proved `pp_thermal_average_uniform_positive` under:
  - `g ≠ 0`, `f0 ≠ 0`
  - `protonDensity > 0`, `mReduced > 0`
  - `E0 > 0`, `dE ≥ 0`
- Result: strict positivity of the finite uniform thermal average for any `n : ℕ`.

## Rust parity
- Added `pp_thermal_average_uniform(...)` in `crates/gutoe-physics/src/equations.rs` mirroring Lean semantics.
- Added tests:
  - `pp_thermal_average_uniform_is_strictly_positive_under_physical_inputs`
  - `pp_thermal_average_uniform_rejects_nonphysical_grid`

## Verification
- `lake build Gutoe.StellarFusion` ✅
- `lake build Gutoe` ✅
- `cargo test -p gutoe-physics pp_thermal_average_uniform -- --nocapture` ✅

## Note
This is still a finite-sample witness (not yet a full continuous integral limit). It replaces a hardcoded 3-point special case with a general finite quadrature family and keeps Lean/Rust parity.
