# Finding 091 — GRAND-353 Upstream Uncertainty Propagation

Date: 2026-02-27
Status: COMPLETE

## Summary
Implemented quantitative uncertainty propagation across upstream cosmology lanes (not only pass/fail windows), with report + CI gate outputs.

New files:
- `crates/gutoe-physics/src/uncertainty.rs`
- `crates/gutoe-physics/src/bin/uncertainty_report.rs`
- `crates/gutoe-physics/src/bin/uncertainty_ci_gate.rs`

## What is propagated
Sampled perturbations feed through:
- Inflation observables (`n_s`, `A_s`)
- Baryogenesis (`eta_B` via `eta10`)
- BBN abundance formulas
- Dark-sector matter fraction lane
- CMB/BAO transfer lane
- Microphysics lane
- FRW background outputs (`H0`, age, timing anchors)

Outputs include distributions (`mean/std/p05/p50/p95/min/max`) for:
- `n_s`, `A_s`, `eta10`, dark-matter fraction
- `H0`, age
- `r_s`, `theta_*`, `l1`, `l2`
- microphysics outputs (`Y_p`, `D/H`, `z_visibility_peak`)
- component pass fractions + overall pass fraction

## Verified outputs (default run, 768 samples)
From `cargo run -q -p gutoe-physics --bin uncertainty_report`:
- `valid_samples = 768`
- `overall pass fraction = 0.819`
- `H0 p50 = 68.024`
- `theta* p50 = 1.090468e-2`
- `Yp p50 = 0.24925`

CI gate (`uncertainty_ci_gate`) checks:
- minimum overall pass fraction
- 95% relative span constraints for `H0` and `theta_*`
- 95% absolute span constraint for `Y_p`

Gate result: pass `true`.

Artifacts:
- `/tmp/bh_renders/uncertainty_report.txt`
- `/tmp/bh_renders/uncertainty_report.json`
- `/tmp/bh_renders/uncertainty_ci_gate.json`
