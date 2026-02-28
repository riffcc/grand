# Cannabinoid Panel (CB1 Neuron Lane) — Multi-Compound Expansion

Date: 2026-02-28  
Status: Implemented (panel report + panel CI gate)

## Goal
Scale the single-THC lane to a comparative cannabinoid panel:
- multiple cannabinoids,
- shared CB1 thermodynamic/decomposition framework,
- shared neuron-response transduction (occupancy -> release/firing).

## What shipped

- Module:
  - `crates/gutoe-physics/src/cannabinoid_panel.rs`
- Binaries:
  - `crates/gutoe-physics/src/bin/cannabinoid_panel_report.rs`
  - `crates/gutoe-physics/src/bin/cannabinoid_panel_ci_gate.rs`
- Reuse:
  - `crates/gutoe-physics/src/neuro_thc.rs` generic CB1 wrappers

## Panel coverage

13 compounds:
- `11_oh_thc`
- `delta9_thc`
- `delta8_thc`
- `anandamide`
- `thcv`
- `cbn`
- `cbg`
- `2_ag`
- `thca`
- `cbdv`
- `cbc`
- `cbd`
- `cbda`

## Executed result

Command:
- `cargo run -q -p gutoe-physics --bin cannabinoid_panel_report`

Output:
- `count = 13`
- `mean_explained_fraction_of_abs_delta_g = 0.1868`
- `mean_abs_residual_closure_error_kj_mol = 5.496`

Top potency (lowest CB1 Ki):
1. `11_oh_thc` (20 nM)
2. `delta9_thc` (40 nM)
3. `delta8_thc` (44 nM)
4. `anandamide` (61 nM)
5. `thcv` (75 nM)

Top occupancy at 100 nM:
1. `11_oh_thc` (`0.833`)
2. `delta9_thc` (`0.714`)
3. `delta8_thc` (`0.694`)
4. `anandamide` (`0.621`)
5. `thcv` (`0.571`)

Artifacts:
- `/tmp/bh_renders/cannabinoid_panel/cannabinoid_panel_report.txt`
- `/tmp/bh_renders/cannabinoid_panel/cannabinoid_panel_report.csv`
- `/tmp/bh_renders/cannabinoid_panel/cannabinoid_panel_report.json`

## CI gate

Command:
- `cargo run -q -p gutoe-physics --bin cannabinoid_panel_ci_gate`

Default windows:
- panel count >= 10
- mean explained fraction >= 0.10
- mean abs residual closure error <= 8.0 kJ/mol
- potency sanity: `delta9_thc occupancy_100nM > cbd occupancy_100nM`

## Honesty statement

This panel is a first-pass comparative scaffold.
Ki priors are literature-scale seeds and can vary by assay/protocol.
Some compounds still show large residual closure errors and need second-pass
descriptor refinement for tighter quantitative closure.
