# 179 — Density Missing-Physics Pass (4d Corridor + Allotropy Network)

Date: 2026-02-28
Scope: `crates/gutoe-physics/src/chemical_thermo.rs`

## Goal
Apply the 3 requested structural channels without element hardcoding:
1. crystal-structure channel,
2. nonmetal allotropy channel,
3. d-band occupancy map (Nb–Pd corridor).

## Changes Implemented
- Added/extended crystal packing and radius multipliers with explicit period-5 transition corridor logic.
- Broadened 4d occupancy map from narrow half-fill peak to a full Nb–Pd corridor strength function:
  - `transition_4d_corridor_strength(...)`
- Strengthened open-network allotropy response for nonmetal/metalloid solids:
  - higher porosity reduction in packing
  - larger radius expansion for directional/ring/layer networks
- Relaxed period-5 transition radius floor dynamically by corridor strength (general rule; no element IDs).

## Verification
- `cargo check -q -p gutoe-physics --lib` ✅
- `cargo test -q -p gutoe-physics --lib chemical_thermo` ✅ (8 passed)
- `GUTOE_CHEM_PHASE_OVERRIDE=1 cargo run -q -p gutoe-physics --bin mass_periodic_report` ✅
- `cargo run -q -p gutoe-physics --bin element_unified_external_benchmark` ✅
- `cargo run -q -p gutoe-physics --bin chemical_thermo_calibrate` ✅

## External Benchmark Movement (Z=1..94)
Baseline before this pass: phase=1.000000, density MAE=3.254827, density condensed MAE=3.691423, ionization MAE=0.400367

After this pass:
- phase_accuracy: **1.000000** (unchanged, still closed)
- density_mae_g_cm3: **2.923209** (improved by 0.331618)
- density_mae_g_cm3_condensed_only: **3.315320** (improved by 0.376103)
- ionization_mae_ev: **0.400367** (unchanged)
- melting/boiling lanes: unchanged in this pass

## Holdout Integrity (No-Fit vs Fitted)
- period holdout (train p1–p4, hold p5–p7): 3.376855 -> 3.725264 (fit worsens holdout)
- block holdout (train s/p, hold d/f): 3.419792 -> 6.052889 (fit overfits strongly)
- block holdout (train d/f, hold s/p): 2.339349 -> 2.241959 (small transfer gain)
- phase holdout accuracy: 1.000000 for all splits

Interpretation: This pass improved baseline structural physics; coefficient-fitting still fails to generalize in two of three splits, so the remaining density error is still physics-limited (not calibration-limited).

## Largest Remaining Density Residuals (absolute)
Top cluster now includes:
- Tl, In, Sn, Pb (post-transition heavy p-block)
- S (open-network allotropy still too dense)
- Nd/Sm/Pm/Np (f-block/actinide contraction gap)
- Nb (early 4d corridor still underdense)

## Next Structural Targets
1. Post-transition heavy p-block crystal channel (In/Sn/Tl/Pb packing + relativistic valence effects).
2. f-block screening correction beyond current f-core scalar (lanthanide/actinide volume compression).
3. nonmetal phase-specific allotrope map extension (C/S/Se family) with structural toggles and holdout checks.
