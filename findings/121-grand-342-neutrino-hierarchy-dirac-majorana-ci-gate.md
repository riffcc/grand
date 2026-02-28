# Finding 121: GRAND-342 Neutrino Hierarchy + Dirac/Majorana CI Gate

Status: GRAND-342 complete
Date: 2026-02-28

## Scope

Close `GRAND-342` by making the neutrino hierarchy prediction explicit (`normal` vs `inverted`) and wiring a falsifiable Dirac/Majorana determination into CI, with hard regression gates and machine-readable artifacts.

## Implementation

### 1) Structural Majorana-symmetry residual

In `crates/gutoe-em/src/flavor.rs`:
- added `neutrino_majorana_symmetry_residual() -> f64`
  - computes `max |M_ij - M_ji|` on the texture matrix.
  - this is a direct symmetry residual for the Majorana condition (`M = M^T`).
- added `neutrino_dirac_majorana_prediction() -> &'static str`
  - returns `"majorana_like"` when residual <= `1e-12`, else `"dirac"`.

### 2) Public exports

In `crates/gutoe-em/src/lib.rs`:
- exported the two new functions above.

### 3) Dedicated CI gate binary

Added `crates/gutoe-em/src/bin/neutrino_ci_gate.rs`:
- evaluates hierarchy prediction,
- evaluates Dirac/Majorana prediction,
- enforces residual exclusion threshold (`majorana_symmetry_residual > 1e-12`),
- enforces absolute-mass lane bounds (`m3 < 0.8 eV`, `sum(mν) < 0.12 eV`),
- writes `/tmp/bh_renders/neutrino_ci_gate.json`,
- exits nonzero on any gate failure.

### 4) Global gate integration

In `crates/gutoe-physics/src/bin/global_gate_report.rs`:
- runs `gutoe-em` neutrino gate binary,
- ingests `/tmp/bh_renders/neutrino_ci_gate.json`,
- includes neutrino lane in both text and JSON global reports,
- folds neutrino pass into `overall_pass`.

## Tests and Run Evidence

Commands executed:
- `cargo test -q -p gutoe-em neutrino_ -- --nocapture`
- `cargo run -q -p gutoe-em --bin neutrino_ci_gate`
- `cargo run -q -p gutoe-physics --bin global_gate_report`

All passed.

## Physics Outputs (from gate artifacts)

From `/tmp/bh_renders/neutrino_ci_gate.json`:
- hierarchy prediction: `normal`
- mass-character prediction: `dirac`
- Majorana symmetry residual: `9.948906419733e-1`
- absolute masses:
  - `m1 = 8.497119214462e-4 eV`
  - `m2 = 6.952371162470e-3 eV`
  - `m3 = 7.904528763902e-3 eV`
  - `sum(mν) = 1.570661184782e-2 eV`
- checks: all true

From `/tmp/bh_renders/global_gate/global_gate_report.json`:
- neutrino lane `pass = true`
- `overall_pass = true`

## Falsifiability Surface

- Hierarchy prediction (`normal`) is directly testable by long-baseline/JUNO-style hierarchy analyses.
- Dirac/Majorana lane is encoded as an explicit structural discriminator (`Majorana symmetry residual`) with a hard CI threshold; future texture changes cannot silently flip interpretation.

## Result

`GRAND-342` now has explicit runtime predictions, hard CI enforcement, and report-level integration into the unified global gate.
