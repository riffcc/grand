# 177 — Chemical Thermo Empirical Calibration with Holdout Gates

## Goal
Add an auditable empirical calibration lane for density-tail coefficients while enforcing anti-overfit discipline:
- always report **train** vs **holdout** side-by-side,
- include block and period holdout regimes,
- no silent refit on holdout.

## Implemented

### 1) Calibration-ready thermo API
File: `crates/gutoe-physics/src/chemical_thermo.rs`

- Added `ChemicalThermoCalibration` (explicit coefficient vector + defaults).
- Added env-based calibration loading (`GUTOE_CHEM_CAL_*`) via `chemical_thermo_calibration_from_env()`.
- Added prefetch/cache path to avoid repeated SCF in optimizer:
  - `CoupledThermoPrefetch`
  - `prefetch_element_thermo_coupled(z, a)`
  - `predict_element_thermo_coupled_from_prefetch_calibrated(...)`
- Existing runtime path remains:
  - `predict_element_thermo_coupled_with_diagnostics(...)` uses default/env calibration.

### 2) New calibrator binary
File: `crates/gutoe-physics/src/bin/chemical_thermo_calibrate.rs`

- Reads external reference CSV.
- Builds per-element prefetch cache once.
- Runs constrained coordinate-descent calibration.
- Emits side-by-side train/holdout metrics for:
  1. `period_holdout_train_p1_p4_hold_p5_p7`
  2. `block_holdout_train_sp_hold_df`
  3. `block_holdout_train_df_hold_sp`

Outputs:
- `/tmp/nuclear_chart/chemical_thermo_calibration_report.txt`
- `/tmp/nuclear_chart/chemical_thermo_calibration_report.json`

## Verification
- `cargo check -q -p gutoe-physics --bin chemical_thermo_calibrate` ✅
- `cargo test -q -p gutoe-physics --lib chemical_thermo` ✅ (8 passed)

## Key holdout results (density MAE)

### Period holdout (train periods 1..4, hold 5..7)
- Baseline train: `2.526741`
- Baseline holdout: `4.061768`
- Fitted train: `2.127343`
- Fitted holdout: `4.441324`
- Read: train improves, holdout degrades (overfit across periods).

### Block holdout A (train s/p, hold d/f)
- Baseline train: `2.560025`
- Baseline holdout: `4.248047`
- Fitted train: `2.246664`
- Fitted holdout: `6.937170`
- Read: severe overfit to s/p, poor d/f transfer.

### Block holdout B (train d/f, hold s/p)
- Baseline train: `4.248047`
- Baseline holdout: `2.560025`
- Fitted train: `4.021180`
- Fitted holdout: `2.545188`
- Read: modest improvement with acceptable transfer in this direction.

## Honest interpretation
The calibration lane is now fully auditable and fast, but current objective still tends to overfit in two of three structural splits. This is exactly why train/holdout gating is necessary before accepting empirically derived coefficients as structural.
