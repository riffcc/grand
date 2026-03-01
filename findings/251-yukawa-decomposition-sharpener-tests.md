# Finding 251 — Yukawa Decomposition Sharpeners (QCD / Yukawa / Casimir)

## Scope
Executed the sharpened decomposition tests requested for the full-dynamics Yukawa lane:

1. **Equation 2 sharpener**: test `gap_lg = 3 - s²(L_g)` against `alpha_s(μ)`.
2. **Equation 3 sharpener**: test whether `d(Δ_ud)/dlnμ` correlates with `y_t²`.
3. **Equation 1 sharpener**: check `Δ_QCD = s²_down - 2` against `4/9 = C_F/N_c`.

Implementation:
- `crates/gutoe-em/src/bin/yukawa_mode_decomp_fit.rs`

Outputs:
- `/tmp/bh_renders/yukawa_mode_decomp_fit.txt`
- `/tmp/bh_renders/yukawa_mode_decomp_fit.csv`
- `/tmp/bh_renders/yukawa_mode_decomp_fit.json`

## Core fitted results

From run summary:

- `s2_down mean = 2.440686475`
- `s2_lg(1e19) = 2.968411312`
- `delta_ud(1e19) = 0.866437062`

### (A) Gap-to-3 vs strong coupling

- Free affine fit:
  - `gap_lg = intercept + slope * alpha_s`
  - `intercept = 0.028325941`
  - `slope = 0.566651152`
  - `R² = 0.879701`
  - `RMSE = 0.006481320`

- Origin-constrained fit:
  - `gap_lg = c1 * alpha_s`
  - `c1 = 0.997641667`
  - `RMSE = 0.017541055`

Interpretation:
- A first-order `gap_lg ~ alpha_s` law is supported.
- Strict through-origin form is plausible but leaves larger residuals than free affine.
- Nonzero intercept indicates subleading/non-`alpha_s` dressing remains at this scan depth.

### (B) Split derivative vs Yukawa strength

Segment-level derivative test:

- `d(delta_ud)/dlnμ = intercept + slope * y_t²(mid)`
- `intercept = -0.000236557`
- `slope = 0.003362222`
- `R² = 0.999786`
- `RMSE = 0.000006633`

Interpretation:
- Extremely strong correlation supports an anomalous-dimension style cumulative Yukawa mechanism.
- This is consistent with the observed behavior: `delta_ud` grows with scale while pointwise `y_t` decreases.

### (C) Casimir proximity check for plateau shift

- `Δ_QCD = s²_down_mean - 2 = 0.440686475`
- `4/9 = 0.444444444`
- Difference:
  - `Δ_QCD - 4/9 = -0.003757970`
  - Relative = `-0.846%`

Interpretation:
- The plateau shift is very close to the SU(3) Casimir ratio candidate `C_F/N_c = 4/9`.
- This is a high-quality structural hit, pending derivation of sub-percent correction terms.

## Additional check: naive `4 alpha_s / pi` scaling

If forced into `Δ_QCD = kappa * (4 alpha_s / pi)`:

- `kappa(mt) = 3.227381254`
- `kappa(1e19) = 18.456963086`

Interpretation:
- `kappa` is not stable; naive one-factor `4 alpha_s / pi` closure does not explain the plateau shift.

## Verdict

All three sharpeners produced discriminating signal:

- **Equation 2**: `gap_lg` tracks `alpha_s` strongly, with evidence for additional offset terms.
- **Equation 3**: derivative-level test strongly supports Yukawa-integral dynamics.
- **Equation 1**: `Δ_QCD` sits within ~0.85% of `4/9`, strongly favoring a Casimir-structured plateau.

This materially strengthens the “vacuum-dressed base instanton + RG dressing” interpretation and provides concrete coefficients for next-order derivation work.
