# Finding 089 — GRAND-351 CMB/BAO Transfer Lane

Date: 2026-02-27
Status: COMPLETE

## Summary
Implemented a new transfer lane that predicts CMB/BAO observables from derived GUTOE background + inflation inputs, with hard CI gates.

New files:
- `crates/gutoe-physics/src/cosmo_transfer.rs`
- `crates/gutoe-physics/src/bin/cmb_transfer_report.rs`
- `crates/gutoe-physics/src/bin/cmb_transfer_ci_gate.rs`

Integrated into assembled universe gate:
- `crates/gutoe-physics/src/universe.rs`
- `crates/gutoe-physics/src/bin/universe_sim.rs`
- `crates/gutoe-physics/src/bin/universe_ci_gate.rs`

## What the lane computes
- Drag redshift `z_drag`
- Sound horizon `r_s`
- Comoving CMB acoustic angle `theta_*`
- Acoustic scale and peak proxies (`l1`, `l2`)
- Growth suppression and `P(k,z)` pivot proxy

## Verified outputs (default run)
From `cargo run -q -p gutoe-physics --bin cmb_transfer_report`:
- `r_s = 149.228 Mpc`
- `theta_* = 1.090493e-2`
- `l1 = 217.61`
- `l2 = 530.96`
- Gate pass: `true`

Artifacts:
- `/tmp/bh_renders/cmb_transfer_report.txt`
- `/tmp/bh_renders/cmb_transfer_report.json`
- `/tmp/bh_renders/cmb_transfer_ci_gate.json`

## Notes
- Corrected acoustic-angle computation to use comoving angular-diameter distance (`D_M`) for CMB peak geometry.
- Transfer lane now contributes directly to `passes_late_universe()` and full universe CI pass/fail.
