# Finding 095 — GRAND-343 Derived `tau_reio` from Structure Timing

Date: 2026-02-27
Status: COMPLETE

## Summary
Implemented a structural `tau_reio` derivation lane and wired it into the CLASS report/three-way CMB runs as the default when `GUTOE_TAU_REIO` is not explicitly set.

New files:
- `crates/gutoe-physics/src/cmb_reionization.rs`
- `crates/gutoe-physics/src/bin/cmb_tau_derived_report.rs`

Updated:
- `crates/gutoe-physics/src/lib.rs`
- `crates/gutoe-physics/src/bin/cmb_class_report.rs`
- `crates/gutoe-physics/src/bin/cmb_three_way_compare.rs`

## Structural chain
- Compute structural reionization redshift:
  - `z_reion ~ (12-4) * ((4+6)/(6+3)) * (eta10/6)^(1/3)`
- Integrate optical depth from `z=0..z_reion` with expansion history + baryon density:
  - `tau_reio = ∫ c sigma_T n_e(z)/((1+z)H(z)) dz`
- Use derived `tau_reio` directly in CLASS input when no override is provided.

## Baseline output
From:
- `GUTOE_CLASS_BIN=/tmp/class_public/class cargo run -q -p gutoe-physics --bin cmb_tau_derived_report`

Derived:
- `z_reion = 9.035`
- `tau_reio = 0.067531`
- `exp(-2 tau) = 0.873661`

Fit impact versus explicit assumption `tau=0.054`:
- Binned TT: `chi2 1984.5 -> 940.7` (`Δ=-1043.8`)
- Full TT: `chi2 4014.7 -> 3168.9` (`Δ=-845.8`)
- Full reduced `chi2`: `1.607 -> 1.268`

## Notes
- This derivation is upstream-coupled (depends on baryogenesis `eta10` and background cosmology), not a CMB residual fit.
- `GUTOE_TAU_REIO` remains available as an explicit override for ablation/testing.
