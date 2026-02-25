# 034 — SM × QCD Unification Gate

Status: Lean closure + runtime parity rows landed.

## Lean additions

- New module: `lean/Gutoe/SMQCDUnification.lean`
- Included in build roots: `lean/lakefile.lean`

Key theorems:

- `qcd_core_gate_holds`
- `sm_qcd_unified_structural_holds`
- `sm_qcd_general_case_bundle_holds`

## What is unified

This bundle now ties together, in one theorem chain:

1. one-generation SM closure (`CanonicalRepConstraints`, anomaly cancellation,
   Witten parity),
2. QCD structural core (`|quarkOrbit|=3`, gluons `=8`, `beta0Clifford>0`,
   `thetaQcdStructural=0`),
3. Strong-CP general split from GRAND-267:
   - emergent-image GUTOE route ⇒ `θ` unphysical,
   - nonzero topological-sector route ⇒ `θ` physical.

## Runtime parity

`crates/gutoe-physics/src/bin/theorem_parity.rs` now includes:

- `qcd_beta0_structural = 58/3`
- `qcd_su3_generators_structural = 8`
- `sm_total_gauge_generators_structural = 12`

Latest run (`/tmp/bh_renders/theorem_runtime_parity.csv`) reports all rows `ok`.
