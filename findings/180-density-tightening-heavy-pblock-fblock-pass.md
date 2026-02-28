# 180 — Density Tightening Pass: Heavy p-block + f-block gating

Date: 2026-02-28
Scope: `crates/gutoe-physics/src/chemical_thermo.rs`

## Objective
Tighten density predictions by targeting two structural residual families while preserving:
- phase closure (100% at 298K),
- ionization lane stability,
- no element-name hardcoding.

Targeted mechanisms:
1. Heavy post-transition p-block compaction (In/Sn/Tl/Pb corridor).
2. f-block contraction improvement with actinide early/late gating.

## Implemented structural changes
- Added generalized corridor functions:
  - `transition_4d_corridor_strength(...)`
  - `post_transition_heavy_strength(...)`
  - `f_block_contraction_strength(...)`
- Coupled these strengths into both:
  - `crystal_packing_multiplier(...)`
  - `crystal_radius_multiplier(...)`
- Added adaptive radius lower-floor logic in coupled predictor for:
  - Lanthanides (f-block-driven lower bound)
  - Heavy post-transition elements (corridor-driven lower bound)
  - Existing 4d transition lower-floor kept and fed by corridor strength
- Corrected over-compression behavior by gating actinide contraction to late-fill character (valence/f-occupancy), avoiding early-actinide blow-up.

## Validation commands
- `cargo check -q -p gutoe-physics --lib`
- `cargo test -q -p gutoe-physics --lib chemical_thermo`
- `GUTOE_CHEM_PHASE_OVERRIDE=1 cargo run -q -p gutoe-physics --bin mass_periodic_report`
- `cargo run -q -p gutoe-physics --bin element_unified_external_benchmark`
- `cargo run -q -p gutoe-physics --bin chemical_thermo_calibrate`

All passed.

## External benchmark movement (Z=1..94)
Previous best (before this pass):
- phase_accuracy: 1.000000
- density_mae_g_cm3: 2.923209
- density_mae_g_cm3_condensed_only: 3.315320
- ionization_mae_ev: 0.400367

After this pass:
- phase_accuracy: **1.000000** (preserved)
- density_mae_g_cm3: **2.813832** (improved by 0.109377)
- density_mae_g_cm3_condensed_only: **3.191271** (improved by 0.124049)
- melting_mae_k: **676.233169** (unchanged)
- boiling_mae_k: **1232.459076** (unchanged)
- ionization_mae_ev: **0.400367** (preserved)

## Holdout integrity (no-fit vs fitted)
From `chemical_thermo_calibration_report.txt`:
- period holdout (train p1–p4, hold p5–p7):
  - baseline 3.195318 -> fitted 3.557515 (worse)
- block holdout (train s/p, hold d/f):
  - baseline 3.228748 -> fitted 5.344203 (overfit)
- block holdout (train d/f, hold s/p):
  - baseline 2.320850 -> fitted 2.207500 (small transfer gain)

Interpretation: baseline structural physics improved; coefficient fitting still fails to generalize in two of three splits, confirming remaining error is missing-physics limited.

## Largest remaining density residuals
Current top absolute density errors are concentrated in:
- open-network allotrope nonmetals (`S`, `C`),
- heavy post-transition (`Tl`, `In`, `Sn`, `Pb`),
- residual d/f compression cluster (`Nb`, `Pd`, `Np`, `Pu`, `Sb`).

## Next structural rung
- Explicit allotrope topology lane (molecular/ring/layer phase-aware density transduction) for `C/S/Se` family.
- Stronger relativistic p-block valence contraction channel for period-6 post-transition metals.
- Separate 4d vs 5d occupancy maps to avoid sharing a single compaction profile.
