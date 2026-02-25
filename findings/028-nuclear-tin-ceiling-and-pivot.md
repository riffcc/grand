# Finding 028 — Tin Ceiling Confirmed, Pivot To EHT Pipeline

Date: 2026-02-25

## Summary

We ran focused nuclear calibration sweeps and confirmed the same outcome repeatedly:

- Within the current model class, `Sn` stabilizes at **7/10** max.
- Pushing tin count above that forces shell-anchor drift and/or bias regression.
- Additional broad sweeps are low-yield; the bottleneck is missing model structure, not search coverage.

This is the stopping signal for sweep-first iteration in this lane.

## Key Numbers

From reproducible runs (`/tmp/nuclear_find028/summary.json`):

### 1) Baseline tuned anchor-preserving lane

- `RMSE = 3.2065 MeV`
- `MAE = 2.1579 MeV`
- `Bias = +0.0045 MeV`
- Magic-N residual means:
  - `N=50: -1.3556`
  - `N=82: -0.2160`
  - `N=126: +1.9438`
- Tin:
  - Predicted stable-like: `6/10`
  - Missing: `[112, 115, 119, 124]`

### 2) Best global RMSE in sweep (non-anchor-preserving)

- `RMSE = 3.1453 MeV`
- `MAE = 2.0965 MeV`
- `Bias = -0.1925 MeV`
- Magic-N residual means:
  - `N=50: -2.6135`
  - `N=82: -1.7858`
  - `N=126: +0.1685`
- Tin:
  - Predicted stable-like: `7/10`
  - Missing: `[112, 115, 124]`

### 3) With new localized Z≈50 isovector valley term

- `RMSE = 3.1642 MeV`
- `MAE = 2.1474 MeV`
- `Bias = +0.2381 MeV`
- Magic-N residual means:
  - `N=50: -1.1601`
  - `N=82: +0.3212`
  - `N=126: +1.9438`
- Tin:
  - Predicted stable-like: `7/10`
  - Missing: `[112, 115, 124]`

Interpretation: the new term is useful as a local degree of freedom, but it does not break the `7/10` ceiling without paying elsewhere.

## Diagnostic Evidence

Tin-local margins (neighbor binding deltas) remain positive where closure is missing:

- `A=112` vs `Z-1` remains `> 0`
- `A=115` vs `Z-1` remains `> 0`
- `A=124` vs `Z+1` remains `> 0`

These are now directly emitted in:

- `mass_periodic_report.json` → `tin_diagnostics.neighbor_binding_deltas_mev`
- `tin_isotope_diagnostics.csv`

## Decision

Stop large parameter sweeps for this nuclear lane for now.  
The model has reached a credible scaffold baseline:

- ~`3.15–3.21 MeV` RMSE over `2,548` AME2020 nuclides
- clear residual structure and quantified failure modes

Next high-leverage move is **not more tuning**.

## Next Energy (Selected)

Proceed to **EHT pipeline / Point 27 lane**:

1. lock synthetic-observation path to measurable shadow-diameter output,
2. generate a direct GUTOE-vs-GR diameter comparison artifact,
3. keep nuclear lane parked at current baseline until structural upgrade (e.g., deeper shell/mean-field term or Cl(1,3)-derived spin-orbit closure physics).
