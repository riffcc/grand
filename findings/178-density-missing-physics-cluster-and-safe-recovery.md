# 178 — Density Missing-Physics Cluster + Safe Recovery

## Objective
Push density MAE down by adding physics terms (not element IDs), while preserving:
- phase closure at 100%,
- ionization lane integrity,
- holdout diagnostics.

## Kept changes (physics-side, stable)
File: `crates/gutoe-physics/src/chemical_thermo.rs`

Added structural mechanisms:
1. **d half-fill compaction** in packing/radius (targets dense 4d/5d metallic rows).
2. **p-network porosity** in packing (targets directional covalent solids).
3. **lanthanide f-core contraction** + tighter lanthanide/actinide radius upper scales.

Resulting benchmark (Z=1..94):
- phase_accuracy: `1.000000`
- density_mae_g_cm3: `3.266714`
- density_mae_g_cm3_condensed_only: `3.704905`
- ionization_mae_ev: `0.400367`

Progression context:
- previous: `3.434563`
- now: `3.266714`

## Rejected change (rolled back)
File attempted: `crates/gutoe-physics/src/ab_initio_qchem.rs`

Attempted poorer d/f screening in SCF kernel (angular screening weights).
Observed impact was catastrophic to global chemistry:
- density_mae_g_cm3 worsened to `4.745467`
- ionization_mae_ev exploded to `33.541276`

Decision: **reverted** that SCF-kernel patch.
This preserves chemistry integrity and keeps the model on the best-known frontier.

## Holdout calibration diagnostic (still active)
Using `chemical_thermo_calibrate` lane:
- period holdout baseline improved vs earlier (`holdout ~3.83`),
- but coefficient fitting still overfits across block splits.

Interpretation: missing physics remains real; optimization alone is not sufficient.

## Current dominant residual families (density)
- Underpredicted: post-transition, lanthanide, transition, actinide
- Overpredicted: nonmetal (notably C, S)

This points to next missing-physics layers:
1. metallic crystal-structure channel (FCC/BCC/HCP/post-transition polymorphs),
2. nonmetal allotropy channel (graphitic/ring/chain topology),
3. d-band occupancy-to-volume mapping for 4d cluster (Nb–Pd).
