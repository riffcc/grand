# 181 — Crystal Blend Weights via Lean Finite Elimination

Date: 2026-02-28

## Why
User directive: stop float fitting and force coefficients by formal elimination.

## What changed
### Lean (formal constraint closure)
Added new module:
- `lean/Gutoe/CrystalStructureWeights.lean`

Key formal objects/theorems:
- `familyCount := grade1_4d.card`
- `crystalCount := magneticTriplet.card`
- `familyWeightQ`, `crystalWeightQ`
- `blend_weights_exact` proves:
  - `familyWeightQ = 4/7`
  - `crystalWeightQ = 3/7`
- `goodNumeratorPairs` finite candidate set over `0..totalCount`
- `good_numerator_pairs_unique` proves unique survivor `(4,3)`
- `blend_numerators_unique` lifts uniqueness propositionally

No `sorry` added.

### Lean build integration
- Added `Gutoe.CrystalStructureWeights` to `lean/lakefile.lean` roots.
- Verified:
  - `cd lean && lake build Gutoe.CrystalStructureWeights` ✅
  - `cd lean && lake build Gutoe` ✅

### Rust wiring (no free decimal for blend)
In `crates/gutoe-physics/src/chemical_thermo.rs`:
- Crystal structure lane remains active (`CrystalPrototype` + structure-based packing prior).
- Blend updated from tuned decimal to Lean-forced rational:
  - `base = (4/7) * base_family + (3/7) * base_crystal`

This removes one major free/tuned scalar from this lane.

## Bench impact (with Lean-constrained weights)
External benchmark (`Z=1..94`):
- phase accuracy: `1.000000`
- density MAE: `2.829502 g/cm³`
- condensed density MAE: `3.209043 g/cm³`
- ionization MAE: `0.400367 eV`

## Notes
- Remaining coefficients in density/thermal channels are still heuristic and should be moved to the same finite-grammar + Lean-elimination workflow incrementally.
- This pass establishes the exact pattern for doing that safely.
