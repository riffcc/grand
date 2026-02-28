# 183 — Thermal Entropy Lean Sweep (Bonding/Valence-Coupled)

## Scope
Targeted the two remaining red thermal lanes (`melting_k`, `boiling_k`) using:
- Lean-constrained rational coefficients only
- No element-specific targeting/hardcoded per-element overrides
- Bonding/valence structure from existing coupled chemistry diagnostics

## Lean closure added
New module:
- `lean/Gutoe/ThermalEntropyClosure.lean`

Closed coefficients:
- `dFusionGainQ = 5/6`
- `dVaporGainQ = 3/5`
- `covalentGainQ = 4/5`
- `metalloidPenaltyQ = 7/4`
- `molecularPenaltyQ = 5/4`

Finite numerator elimination theorems added for each denominator family.

## Rust wiring
File:
- `crates/gutoe-physics/src/chemical_thermo.rs`

Changes:
- Added `thermal_entropy_scales(...)` driven by:
  - family/period/crystal prototype
  - valence-electron topology
  - d-band occupancy proxy (hints + valence fallback)
- Applied entropy scales in thermal transduction:
  - `T_m ~ ΔH_fus / (S_fus * fusion_scale)`
  - `T_b ~ ΔH_vap / (S_vap * vapor_scale)`
- Enforced physical ordering gate:
  - `boiling_temperature_k >= melting_temperature_k + 2 K`

## Verification
- `cargo test -q -p gutoe-physics --lib chemical_thermo` ✅
- `cargo test -q -p gutoe-physics --lib ab_initio_qchem` ✅
- `lake build Gutoe.ThermalEntropyClosure` ✅
- `lake build Gutoe` ✅
- `mass_periodic_report` + `element_unified_external_benchmark` ✅

## External benchmark impact (Z=1..94)
Before this sweep:
- `melting_mae_k = 676.938609`
- `boiling_mae_k = 1229.031166`
- `elements_with_any_red = 94`

After pass 1 (thermal entropy layer):
- `melting_mae_k = 622.583256`
- `boiling_mae_k = 1128.266698`
- `elements_with_any_red = 93`

After pass 2 (transition valence gate + metalloid topology penalty):
- `melting_mae_k = 448.703068`
- `boiling_mae_k = 973.868170`
- `elements_with_any_red = 90`

Unchanged/held:
- `phase_accuracy = 1.000000`
- `ionization_mae_ev = 0.307058`
- `density_mae_g_cm3 = 2.829566`

## Notes
- The largest remaining thermal misses cluster in high-melting 4d/5d transition metals and some p-block allotropy-sensitive systems.
- This lane is still reduced-order; crystal/allotropy resolution (bcc/fcc/hcp + polymorph fractions) remains the dominant missing physics for thermal closure.
