# 185 — Refractory Transition Thermal Entropy Lever (General, Non-Targeted)

## Scope
Implemented a general refractory correction in `chemical_thermo.rs` to raise transition-metal thermal locking via entropy suppression, with no element-specific branching.

## Change
Added a new branch in `thermal_entropy_scales(...)` for `ChemicalFamily::Transition && period >= 4`:

- Refractory strength uses only structural signals:
  - `d_frac`
  - half-filled shell peak `1 - 2|d_fill - 0.5|`
  - open-shell gate
  - crystal prototype gate (`bcc/hcp/fcc/...`)
  - period weighting
- Rational gains (documented as Lean-constrained style coefficients):
  - fusion suppression: `7/20`
  - vapor suppression: `1/4`

## Verification
- `cargo test -q -p gutoe-physics --lib chemical_thermo` -> pass (8/8)
- `cd lean && lake build Gutoe` -> pass (warnings only)

## External benchmark impact (Z=1..94)
After regenerating `mass_periodic_report` and rerunning `element_unified_external_benchmark`:

- phase_accuracy: `1.000000` (unchanged)
- density_mae_g_cm3: `2.829566` (unchanged)
- melting_mae_k: `448.703068 -> 448.363222` (improved by `0.339846` K)
- boiling_mae_k: `973.868170 -> 973.058846` (improved by `0.809324` K)
- ionization_mae_ev: `0.307058` (unchanged)
- elements_with_any_red: `90` (unchanged)
- red canary: pass (`strict=true`, no regression)

## Honest assessment
This lever is directionally correct but low-amplitude. It improves thermal MAE slightly without any regression, but it does not reduce `elements_with_any_red` and does not materially change the dominant red clusters.

Likely next step: stronger physics coupling (not coefficient inflation), specifically a crystal-prototype + cohesive split for refractory transition rows.
