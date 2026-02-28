# 182 — Crystal Coefficients via Lean Finite Elimination

Date: 2026-02-28

## Objective
Replace additional crystal-channel floating-point heuristics with Lean-constrained exact rationals and prove finite uniqueness/elimination over structural denominator families.

## Lean additions
New module:
- `lean/Gutoe/CrystalCoefficientClosure.lean`

Registered in:
- `lean/lakefile.lean`

### Structural counts used
- `d = 2^4 = 16`
- `|grade1| = 4`
- `|grade2| = C(4,2) = 6`
- `|SU(2)| = |magneticTriplet| = 3`
- `Z3 grade-1 fixed count = 1`

### Coefficients closed in Lean
- corridor d-weight = `13/20`
- corridor v-weight = `7/20`
- transition pack gain = `1/2`
- post-transition pack gain = `5/12`
- lanthanide pack gain = `3/10`
- actinide pack gain = `3/25`
- lanthanide radius gain = `9/50`
- actinide radius gain = `2/25`

### Finite elimination theorems
For each coefficient family, numerator candidates over the structural denominator are finitely enumerated and filtered; uniqueness theorems prove only the selected numerator survives.

## Rust wiring
Updated:
- `crates/gutoe-physics/src/chemical_thermo.rs`

Replaced decimal literals with exact rational forms aligned to Lean closures for the coefficients above.

## Verification
- `cd lean && lake build Gutoe.CrystalCoefficientClosure` ✅
- `cd lean && lake build Gutoe` ✅
- `cargo check -q -p gutoe-physics --lib` ✅
- `cargo test -q -p gutoe-physics --lib chemical_thermo` ✅ (8 passed)
- `GUTOE_CHEM_PHASE_OVERRIDE=1 cargo run -q -p gutoe-physics --bin mass_periodic_report` ✅
- `cargo run -q -p gutoe-physics --bin element_unified_external_benchmark` ✅
- `cargo run -q -p gutoe-physics --bin chemical_thermo_calibrate` ✅

## Metrics
External benchmark (`Z=1..94`):
- phase accuracy: `1.000000`
- density MAE: `2.829566 g/cm³`
- condensed-only density MAE: `3.209115 g/cm³`
- ionization MAE: `0.400367 eV`

Holdout (baseline -> fitted MAE):
- period split: `3.227652 -> 3.553069`
- s/p -> d/f: `3.235277 -> 3.293001`
- d/f -> s/p: `2.345299 -> 2.218300`

Interpretation: this pass primarily improves epistemic integrity (formal closure, no decimal tuning for selected coefficients) while preserving the existing density/phase/IE frontier.
