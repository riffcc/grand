# Finding 122: GRAND-354 Parameter Degeneracy Audit + CI Gate

Status: GRAND-354 complete
Date: 2026-02-28

## Scope

Close `GRAND-354` by adding an explicit degree-of-freedom audit for the cosmology lane with:
- parameter provenance inventory,
- numeric sensitivity/Jacobian structure,
- rank/condition summary,
- hidden-ΛCDM re-encoding check,
- machine-readable CI artifact.

## Implementation

Added new binary:
- `crates/gutoe-physics/src/bin/parameter_degeneracy_ci_gate.rs`

Integrated into global aggregate gate:
- `crates/gutoe-physics/src/bin/global_gate_report.rs`

The gate now writes:
- `/tmp/bh_renders/parameter_degeneracy_ci_gate.txt`
- `/tmp/bh_renders/parameter_degeneracy_ci_gate.json`

and exits nonzero on audit failure.

## Parameter Inventory (Provenance)

The artifact now explicitly classifies parameters as:
- `derived_constant`
- `observational_anchor`
- `fixed_assumption`
- `tunable_runtime`

Current core cosmology inventory totals:
- total entries: `10`
- free/tunable entries: `2`
  - `pmns_theta23_alpha2_c`
  - `leptogenesis_pmns_gain`

## Sensitivity Matrix + Degeneracy Structure

Outputs audited:
- `H0`, `Omega_m`, `Omega_lambda`, `r_s`, `theta*`, `l1`, `l2`, `Yp`, `D/H`, `eta_B`

Analyzed knobs:
- `pmns_theta23_alpha2_c` (tunable)
- `leptogenesis_pmns_gain` (tunable)
- `omega_r0` (fixed assumption)
- `omega_k0` (fixed assumption)

From `/tmp/bh_renders/parameter_degeneracy_ci_gate.json`:
- all-knob rank: `3`
- all-knob condition number: `1.3513579932701172e4`
- tunable-only rank: `1`
- tunable→transfer max sensitivity: `0.0`
- tunable→baryo max sensitivity: `2.9143897999010107e-3`

Interpretation:
- tunable freedom is low-dimensional and effectively one-directional for baryogenesis (`rank=1` in tunable subspace),
- tunable knobs do not couple into transfer/CMB geometry outputs in the assembled lane,
- the lane does not have enough hidden freedom to emulate a generic ΛCDM refit through these knobs.

## Verdict

Gate verdict:
- `no_hidden_lcdm_reencoding_in_core_lane`
- `overall_pass = true`

Global CI aggregation now includes this gate and remains green (`global_gate overall_pass=true`).
