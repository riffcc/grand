# Finding 094 — GRAND-354 BBN Likelihood Contours

Date: 2026-02-27
Status: COMPLETE

## Summary
Added a dedicated BBN contour scanner to quantify sensitivity in `(Omega_b, eta10)` instead of pass/fail-only gates.

New file:
- `crates/gutoe-physics/src/bin/bbn_likelihood_scan.rs`

## Output artifacts
From:
- `cargo run -q -p gutoe-physics --bin bbn_likelihood_scan`

Artifacts:
- `/tmp/bh_renders/bbn_likelihood_scan.csv`
- `/tmp/bh_renders/bbn_likelihood_scan.json`

## What it computes
- Grid scan over:
  - `Omega_b` factor in `[0.90, 1.10]`
  - `eta10` factor in `[0.80, 1.20]`
- Joint `chi2` against BBN anchors:
  - `Y_p = 0.245 ± 0.003`
  - `D/H = 2.547e-5 ± 0.050e-5`
- 1σ/2σ/3σ contour occupancy fractions
- Local logarithmic sensitivities:
  - `dln(D/H)/dln(Omega_b)`
  - `dln(D/H)/dln(eta10)`

## Baseline run summary
- Best scan point:
  - `chi2 = 0.2928`
  - `Omega_b factor = 0.9000`
  - `eta10 factor = 0.8300`
  - `Y_p = 0.24649`
  - `D/H = 2.536e-5`
- Sensitivity (baseline neighborhood):
  - `dln(D/H)/dln(Omega_b) = 0.000`
  - `dln(D/H)/dln(eta10) = -1.462`

## Note
This scanner is an analysis tool to expose tension/sensitivity surfaces explicitly. It does not introduce fitted parameters into the core theory lane.
