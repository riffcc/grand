# Finding 092 — GRAND-355 CLASS Full-Spectrum CMB Lane

Date: 2026-02-27
Status: COMPLETE

## Summary
Added a full Boltzmann-lane harness using CLASS/classy and compared derived GUTOE cosmology inputs against Planck TT datasets (binned and unbinned), with reproducible scans and diagnostics.

New files:
- `crates/gutoe-physics/src/cmb_class.rs`
- `crates/gutoe-physics/src/bin/cmb_class_report.rs`
- `crates/gutoe-physics/src/bin/cmb_likelihood_scan.rs`
- `crates/gutoe-physics/src/bin/cmb_three_way_compare.rs`
- `scripts/cmb_pull_denoise.py`
- `crates/gutoe-physics/data/COM_PowerSpect_CMB-TT-binned_R3.01.txt`
- `crates/gutoe-physics/data/COM_PowerSpect_CMB-TT-full_R3.01.txt`

Updated:
- `crates/gutoe-physics/src/lib.rs`
- `crates/gutoe-physics/data/README.md`

## Key results
With derived baseline inputs (`A_s=2.219e-9`, `n_s=0.9654`, `h=0.6802`, `Omega_b h^2=0.02281`, `Omega_cdm h^2=0.12440`) and explicit `tau_reio` assumption:

- Binned TT: `chi2=1984.47` over `n=83` (`red=24.20`)
- Full TT (unbinned): `chi2=4014.70` over `n=2499` (`red=1.607`)

Three-way consistency:
- Prediction vs binned: `chi2=1984.47`, `n=83`
- Prediction vs full: `chi2=4014.70`, `n=2499`
- Binned vs rebinned-from-full: `chi2=7.63`, `n=83`, `red=0.093`

Likelihood scans:
- Coarse full-spectrum `A_s/tau` scan best:
  - `A_s=2.136061e-9`, `tau=0.055`, `chi2=3057.968`, `red=1.224`
- Refined scan best:
  - `A_s=2.116643e-9`, `tau=0.051`, `chi2=3057.114`, `red=1.224`

## Pull-structure diagnostics
Denoise analysis found strong low-dimensional structure:
- Global affine pull map `p_t ≈ 0.724 p_b - 2.621`, `R^2=0.779`
- Local affine residual compression: RMSE `3.909 -> 0.267` (~`14.65x`)
- Full delta smooth-component extraction: RMSE `0.6219 -> 0.0085` (~`73.3x`)

Interpretation: mismatch is dominated by smooth multipole-envelope structure, not random per-bin noise.

## Known boundary
- `tau_reio` is still explicit assumption in this lane.
- This finding is diagnostic and likelihood-focused; no fitted correction was promoted into theory parameters.
