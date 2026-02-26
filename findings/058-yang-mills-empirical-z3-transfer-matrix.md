# Finding 058 — Empirical Z3-Driven Transfer Matrix and GEVP Gap Extraction

Date: 2026-02-26

Scope: GRAND-296/297/298

## What was done

The Yang-Mills mass-gap scaffold was upgraded from toy matrices to real lattice-driven estimators from `gutoe-em::step_counted` trajectories.

Two extraction lanes now run per volume:

1. **Surrogate transition transfer matrix**
- Build scalar observable from Z3 cycle density + orbit skew.
- Bin into 3 states.
- Build empirical transition kernel with Laplace smoothing.
- Construct reversible symmetric PSD transfer surrogate.

2. **Operator-basis correlator GEVP lane (new hardening)**
- Operator basis per sample: `[cycle_density, skew01, skew12]`.
- Build covariance/correlator matrices `C(t0)` and `C(t1)` from real time series.
- Build generalized-eigenvalue transfer proxy `M = L^{-1} C(t1) L^{-T}` using Cholesky of `C(t0)`.
- Extract `(lambda0, lambda1)` and `m_gap = -ln(lambda1/lambda0) / ((t1 - t0) a_t)`.

The selected reporting lane is GEVP when available, with surrogate as fallback/cross-check.

## Artifacts

- `/tmp/bh_renders/ym_mass_gap_report.txt`
- `/tmp/bh_renders/ym_mass_gap_report.json`

## Key run configuration

- `a_t = 1.0`
- `L = 6, 8, 10, 12`
- `burn_in = 24`
- `steps = 420`
- `seeds_per_volume = 6`
- `GEVP: t0=1, t1=2, reg=1e-8`

## Selected (GEVP) gap results

- `L=6`: `lambda0=0.98859`, `lambda1=0.97121`, `gap=0.01774`
- `L=8`: `lambda0=0.99102`, `lambda1=0.97363`, `gap=0.01770`
- `L=10`: `lambda0=0.99422`, `lambda1=0.97757`, `gap=0.01690`
- `L=12`: `lambda0=0.99455`, `lambda1=0.97934`, `gap=0.01542`

Trend checks:
- `monotone_nonincreasing_in_volume = true` (tolerance `0.01`)
- Continuum stability band (using seed/statistical envelope):
  - `[0.001317, 0.029524]`
- Continuum proxy fit vs `a^2`:
  - intercept `0.01556` (positive)

## Interpretation

This is now a real-data extraction from Z3 lattice dynamics with an operator-basis GEVP path in place. Under this estimator, the finite-volume gap remains positive and the continuum intercept proxy is positive.

## Honesty / limitations

- This is still a coarse 3-operator basis, not yet a full Wilson-loop operator tower.
- The continuum statement is still a numerical phenomenology proxy fit, not a theorem-level continuum proof.
- Next step: extend operator basis and solve multi-operator GEVP with plateau diagnostics for excited-state contamination control.
