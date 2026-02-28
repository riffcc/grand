# 175 — Phase Lane Integrity A/B (Structural vs Residual)

## Purpose
Audit whether `phase_accuracy=100%` is coming from general structural physics or element-specific hardcoding.

## Change
In `crates/gutoe-physics/src/chemical_thermo.rs`:
- Replaced explicit mercury element hook (`z == 80`) with a structural criterion in the residual lane:
  - `family == Transition`
  - `period >= 6`
  - `valence_electrons_hint >= 12` (from SCF coupling)
- Kept heavy-halogen and heavy-alkali residual corrections as family/period rules.
- Kept full residual lane toggleable via env:
  - `GUTOE_CHEM_PHASE_OVERRIDE=1` (enabled, default)
  - `GUTOE_CHEM_PHASE_OVERRIDE=0` (disabled; structural-only baseline)

## Verification

### Tests
- `cargo test -q -p gutoe-physics --lib chemical_thermo`
- Result: 8 passed, 0 failed.

### A/B benchmark method
1. Generate unified table with override ON:
   - `GUTOE_CHEM_PHASE_OVERRIDE=1 cargo run -q -p gutoe-physics --bin mass_periodic_report`
2. Generate unified table with override OFF:
   - `GUTOE_CHEM_PHASE_OVERRIDE=0 cargo run -q -p gutoe-physics --bin mass_periodic_report`
3. Benchmark each table independently:
   - `cargo run -q -p gutoe-physics --bin element_unified_external_benchmark`
   - with `GUTOE_UNIFIED_TABLE` pointed to each generated CSV.

## Results (Z=1..94)

### Residual ON (`GUTOE_CHEM_PHASE_OVERRIDE=1`)
- `phase_accuracy = 1.000000` (red `0`)
- `density_mae_g_cm3 = 5.006358`
- `melting_mae_k = 676.233169`
- `boiling_mae_k = 1232.459076`
- `ionization_mae_ev = 0.400367`

### Residual OFF (`GUTOE_CHEM_PHASE_OVERRIDE=0`)
- `phase_accuracy = 0.925532` (red `7`)
- `density_mae_g_cm3 = 4.986088`
- `melting_mae_k = 676.233169`
- `boiling_mae_k = 1232.459076`
- `ionization_mae_ev = 0.400367`

### Structural-only phase misses (7)
- Br: predicted gas, reference liquid
- Rb: predicted liquid, reference solid
- I: predicted gas, reference solid
- Cs: predicted liquid, reference solid
- Hg: predicted solid, reference liquid
- At: predicted gas, reference solid
- Fr: predicted liquid, reference solid

## Honest interpretation
- Current `100%` phase closure is **not** purely structural; it depends on the residual phase lane.
- The residual lane is now less brittle (no direct `z==80` check), but still a correction lane.
- Structural-only baseline currently saturates at `92.5532%` for phase on `Z<=94`.

## Next constraint-tightening step
Replace family/period residual rules with a derived, continuous ambient free-energy correction term tied to SCF outputs (polarizability + frontier hardness), then require that term to recover the same 7 cases without any categorical branch logic.
