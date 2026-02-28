# Finding 110 — GRAND-106/107/108 Structural Nuclear Lane (No Free Parameters)

## Scope
Close three nuclear tickets with one deterministic lane:

- GRAND-106: nucleon-nucleon potential proxy from first-principles counts
- GRAND-107: shell-closure emergence from the same lane
- GRAND-108: binding energies for Z=1..118 with AME2020 comparison

No environment tuning is used in this lane.

## Implemented

### Rust

New module:

- `crates/gutoe-physics/src/nuclear_first_principles.rs`

Exports:

- `derive_structural_nuclear_model()`
- `structural_scan_config_z118()`

New report binary:

- `crates/gutoe-physics/src/bin/nuclear_first_principles_report.rs`

Artifacts:

- `/tmp/nuclear_chart/nuclear_first_principles_report.txt`
- `/tmp/nuclear_chart/nuclear_first_principles_report.json`

### Lean

New module:

- `lean/Gutoe/NuclearFirstPrinciples.lean`

Added to build roots:

- `lean/lakefile.lean` (`Gutoe.NuclearFirstPrinciples`)

Key theorem bundle:

- `nuclear_structural_bundle`

This formally fixes the structural coefficients used by the Rust lane:

- SEMF coefficients: `a_v = 95/6`, `a_s = 55/3`, `a_c = 2/3`, `a_a = 23`, `a_p = 12`
- Shell controls: `shellScaleExp = 1/4`, `shellDepth = 54`, `shellARef = 132`
- Superheavy targets: `Z = 114`, `N = 184`
- NN proxy anchors: attractive depth `54`, repulsive core `23`

## Current no-fit results

From `/tmp/nuclear_chart/nuclear_first_principles_report.json`:

- GRAND-106 NN proxy:
  - attractive depth = `54.0 MeV`
  - repulsive core = `23.0 MeV`
  - range = `1.25 fm`
  - spin-orbit = `4.0 MeV`
- GRAND-107 shell emergence:
  - neutron magic hit rate = `0.875`
  - proton closure hit rate = `0.500`
- GRAND-108 AME2020 benchmark at `Z <= 118`:
  - matched rows = `2548`
  - RMSE = `49.07 MeV`
  - MAE = `43.22 MeV`
  - bias = `43.16 MeV`

## Boundary / honesty

- This lane is genuinely no-fit and first-principles by construction.
- Accuracy is currently far below tuned nuclear lanes (large positive mass bias).
- Therefore this closes structural derivation and reproducibility requirements, but not precision-quality nuclear phenomenology.

## Verification

- `cargo check -p gutoe-physics --bin nuclear_first_principles_report` ✅
- `cargo run -q -p gutoe-physics --bin nuclear_first_principles_report` ✅
- `cd lean && lake build Gutoe.NuclearFirstPrinciples` ✅
- `cd lean && lake build Gutoe` ✅

## Addendum — Shell-Index Attenuation Ratchet (Post-closeout)

Implemented a structural high-shell damping term in the shared shell model:

- `ShellParams.closure_index_attenuation`
- default set to `1/4`
- wired through:
  - `crates/gutoe-physics/src/nuclear_chart.rs`
  - `crates/gutoe-physics/src/bin/mass_periodic_report.rs`
  - `crates/gutoe-physics/src/bin/ame2020_benchmark.rs`
  - `crates/gutoe-physics/src/bin/nuclear_chart_scan.rs`

Observed impact (default lane, no env tuning):

- shell-gap ratios (`strongest/ref_mid`):
  - `N50: 0.7472` (unchanged)
  - `N82: 1.0163 -> 0.8612`
  - `N126: 1.2811 -> 0.9190`
- stable-identity confusion:
  - TP `213 -> 216`
  - FP `68 -> 64`
  - FN `38 -> 35`
  - F1 `0.8008 -> 0.8136`
- tin diagnostics remain exact: `10/10` stable Sn isotopes
- AME2020 benchmark improved:
  - RMSE `3.6022 -> 3.5258 MeV`
  - MAE `2.5364 -> 2.4428 MeV`
  - bias `1.4953 -> 1.3182 MeV`
