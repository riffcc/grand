# Finding 041 — Maxwell-Boltzmann Weighted pp Thermal Kernel

Date: 2026-02-25  
Scope: GRAND-282 (in progress)

## What Landed

We extended the strict-positive pp weak-rate result to include thermal weighting.

### Lean (`lean/Gutoe/StellarFusion.lean`)

Added:
- `maxwellBoltzmannWeight T E = exp(-E/T)` with positivity theorem.
- `ppThermalKernel = ppWeakRateFromSU2 * maxwellBoltzmannWeight`.
- `pp_thermal_kernel_positive` under the same finite/positive physical assumptions used for the weak+Gamow kernel.
- `ppThermalAverage3` and theorem `pp_thermal_average3_positive` as a positive 3-point quadrature witness for thermal averaging.

### Rust (`crates/gutoe-physics/src/equations.rs`)

Added matching runtime functions:
- `maxwell_boltzmann_weight`
- `pp_thermal_kernel`
- `pp_thermal_average3`

Added tests verifying positivity under physical inputs.

## Verification

- `lake build Gutoe.StellarFusion` ✅
- `lake build Gutoe` ✅
- `cargo test -p gutoe-physics maxwell_boltzmann_weight_is_positive_for_positive_temperature -- --nocapture` ✅
- `cargo test -p gutoe-physics pp_thermal_kernel_is_strictly_positive_under_physical_inputs -- --nocapture` ✅
- `cargo test -p gutoe-physics pp_thermal_average3_is_strictly_positive_under_physical_inputs -- --nocapture` ✅

## Status

This is a thermal-kernel and positive quadrature witness, not yet the full continuous Maxwell-Boltzmann integral theorem.

- GRAND-282 moved to `In Progress` for the full integral/bounds completion.
