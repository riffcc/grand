# Finding 253 — Gap Intercept Probe (`gap_lg` vs `alpha_s`, quadratic forms)

## Scope
Probed the requested next-order model for the ladder gap

- `gap_lg = 3 - s²(L_g)`

using the full-dynamics scan points (`m_t ... 1e19 GeV`):

1. `gap_lg = c1 * alpha_s + c2 * alpha_s²`
2. `gap_lg = alpha_s + c2 * alpha_s²` (fixed `c1 = 1`)
3. reference baselines:
   - origin linear: `gap_lg = c1 * alpha_s`
   - affine linear: `gap_lg = b0 + c1 * alpha_s`

## Point estimates

From direct least squares on current 7-point scan:

- origin linear:
  - `c1 = 0.997641667`
  - `RMSE = 0.017541056`

- affine linear:
  - `b0 = 0.028325941`
  - `c1 = 0.566651146`
  - `RMSE = 0.006481320`

- quadratic, no intercept:
  - `c1 = 1.916340718`
  - `c2 = -10.888164740`
  - `RMSE = 0.003042030`

- quadratic with fixed `c1=1`:
  - `c2 = -1.212140610`
  - `RMSE = 0.016546424`

## Stability checks

### Leave-one-scale-out (LOO), quadratic no-intercept

- `c1` range: `[1.874689, 2.072949]`
- `c2` range: `[-14.272545, -10.471485]`
- RMSE range: `[0.001236, 0.003279]`

### Bootstrap (49,999 usable resamples)

95% intervals:

- `c1`: `[1.729534, 2.125691]`
- `c2`: `[-15.066038, -8.843911]`
- RMSE: `[0.000434, 0.003804]`

For fixed-`c1=1` model:

- `c2` median near `-1.161`
- broad/high-tail instability in bootstrap (`97.5%` quantile large), indicating poor conditioning under fixed-`c1` with this short dataset.

## Interpretation

What survives:

- A pure linear origin model (`c1≈1`) captures first-order behavior.
- Adding an `alpha_s²` term materially improves residuals.

What does **not** yet survive:

- The specific constrained hypothesis `gap_lg = alpha_s + c2 * alpha_s²` is not strongly preferred by fit quality and is unstable under resampling compared to the unconstrained quadratic lane.

Current status:

- Intercept/leakage is real and can be absorbed by higher-order terms.
- Coefficients are not yet at a “clean structural constant” stage with current 7-point depth.
- More UV/IR anchors or a theory-prior-constrained fit is needed before claiming a closed-form group-factor identity for `c2`.
