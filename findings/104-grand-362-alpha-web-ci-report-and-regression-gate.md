# GRAND-362 — Alpha Web CI Report + Regression Gate

## Goal

Make `alpha^{-1}=137` closure operationally reproducible in one command, with
cross-lane propagation into lepton masses and the `G` bridge.

## Implemented

1. Unified CI-style report binary

- `crates/gutoe-physics/src/bin/alpha_web_ci_report.rs`
- Outputs:
  - `/tmp/bh_renders/alpha_web_ci_report/alpha_web_ci_report.txt`
  - `/tmp/bh_renders/alpha_web_ci_report/alpha_web_ci_report.json`

Report includes:
- alpha identity block (`T(16)+1=137`, structural vs physical alpha)
- structural-vs-physical alpha lepton-lane deltas (`mu`, `tau`)
- `G` bridge values from the electron transduction factor
- CI gate booleans and overall `passes_all`

2. Hard regression test in `gutoe-em`

- File: `crates/gutoe-em/src/alpha.rs`
- Test: `structural_alpha_identity_and_lane_regression_gate`
- Gate asserts:
  - exact identity: `triangular(2^4) + 1 = 137`
  - structural lane sanity: `|mu_rel|<1%`, `|tau_rel|<1%`

## Verification

Commands run:

- `cargo test -p gutoe-em structural_alpha_identity_and_lane_regression_gate -- --nocapture`
- `cargo run -q -p gutoe-physics --bin alpha_web_ci_report`

Results:

- test passes
- report generated
- `ci_gate passes_all=true`

## Why this matters

This turns theorem closure (`alpha^{-1}=137`) into a one-button reproducible
artifact with downstream consistency checks, so review does not depend on chat
or manual reconstruction.
