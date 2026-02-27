# Finding 093 — GRAND-344 Derived Silk Envelope (First Pass)

Date: 2026-02-27
Status: IN PROGRESS

## Summary
Implemented a parameter-free microphysics-derived damping module for CMB TT using in-framework recombination and diffusion integrals.

New files:
- `crates/gutoe-physics/src/cmb_damping.rs`
- `crates/gutoe-physics/src/bin/cmb_derived_damping_report.rs`

Updated:
- `crates/gutoe-physics/src/lib.rs`

## What is derived (no CMB fit knobs)
From `MicrophysicsAssumptions`:
- Visibility peak redshift `z_*`
- Visibility width `sigma_z`
- Diffusion scale `k_D` (and `ell_diff = k_D * D_M(z_*)`)
- Visibility-width scale `ell_vis`

Then applied envelope to TT spectrum:
- `D_ell -> D_ell * exp(-(ell/ell_diff)^2) * exp(-(ell/ell_vis)^2)`

## First-pass output
From:
- `GUTOE_CLASS_BIN=/tmp/class_public/class cargo run -q -p gutoe-physics --bin cmb_derived_damping_report`

Derived scales:
- `z_* = 1063.00`
- `ell_diff = 1598.6`
- `ell_vis = inf` (visibility-width term negligible in this pass)

Fit impact:
- Binned TT: `chi2 1984.5 -> 122247.9`
- Full TT: `chi2 4014.7 -> 91072.7`

## Interpretation
This confirms the key boundary cleanly:
- CLASS baseline already includes physical diffusion damping.
- Applying an additional absolute diffusion envelope post-CLASS double-counts damping and degrades fit.

So the right next step is not an absolute multiplier. It is a **derived differential correction** (microphysics-vs-lane mismatch operator), then cross-channel validation on TE/EE.

## Next step (unchanged objective)
- Build a differential envelope from explicit microphysics lane against transfer/Boltzmann lane assumptions.
- Re-test on TT.
- Apply unchanged operator to TE/EE for cross-channel validation.
