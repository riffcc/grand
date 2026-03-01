# 186 — Refractory Crystal+Cohesive Coupling (Transition Rows, Period ≥ 5)

## Scope
Strengthened the refractory transition-metal lane by coupling crystal prototype + d-shell structure into both:
- thermal entropy scaling (`thermal_entropy_scales`)
- cohesive transduction (`raw_cohesive` / cohesive clamp window)

No element-specific branches were added.

## Structural changes
File: `crates/gutoe-physics/src/chemical_thermo.rs`

1) Added `transition_refractory_strength(...)`:
- Inputs: family, period, crystal prototype, orbital hints
- Signals used:
  - d-open-shell peak `4 d_fill (1-d_fill)`
  - valence corridor
  - valence edge gate
  - crystal prototype gate (bcc/hcp/fcc/other)
  - period gate (`period>=5`)
  - blended 4d corridor contribution

2) Thermal lane update:
- In `thermal_entropy_scales`, transition `period>=5` now applies moderated refractory suppression:
  - fusion scale gain: `1/2`
  - vapor scale gain: `3/10`

3) Cohesive lane update:
- Added `refractory_cohesion_gate` multiplying `raw_cohesive`
- Added refractory-dependent upper clamp scale for `coupled_cohesive`

## Verification
- `cargo test -q -p gutoe-physics --lib chemical_thermo` -> pass
- `cargo run -q -p gutoe-physics --bin mass_periodic_report` -> pass
- `cargo run -q -p gutoe-physics --bin element_unified_external_benchmark` -> pass
- `cd lean && lake build Gutoe` -> pass
- Red canary (strict) -> pass (`elements_with_any_red` did not increase)

## External benchmark delta (Z=1..94)
Baseline before this pass:
- melting_mae_k = `448.703068`
- boiling_mae_k = `973.868170`
- melting_red = `72`
- boiling_red = `67`
- elements_with_any_red = `90`

After this pass:
- melting_mae_k = `436.480955` (improved)
- boiling_mae_k = `972.843967` (improved)
- melting_red = `72` (unchanged)
- boiling_red = `65` (improved)
- elements_with_any_red = `90` (no regression)
- phase_accuracy = `1.000000` (unchanged)
- density_mae_g_cm3 = `2.829566` (unchanged)
- ionization_mae_ev = `0.307058` (unchanged)

## Honest assessment
This closes a stronger refractory coupling pass without canary regression.
Thermal MAE improved modestly and boiling red count improved by 2, but global red-element count remains 90.
The lane is now better coupled physically; next gain likely requires explicit crystal-prototype-dependent latent channels beyond scalar cohesion boosts.
