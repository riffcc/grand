# Finding 254 — Rational-Lock Candidate for `gap_lg(α_s)`

## Claim tested
For the ladder gap

- `gap_lg = 3 - s²(L_g)`

test whether the fitted quadratic coefficients are structurally:

- `c1 = 23/12`
- `c2 = -98/9`

in

- `gap_lg = c1 * alpha_s + c2 * alpha_s²`.

## Numerical result

From the full-dynamics 7-point scan (`m_t` to `1e19 GeV`):

- Free fit:
  - `c1 = 1.916340718`
  - `c2 = -10.888164740`
  - `RMSE = 0.003042030`

- Structural forced fit:
  - `c1 = 23/12 = 1.916666667`
  - `c2 = -98/9 = -10.888888889`
  - `RMSE = 0.003042064`

Difference in RMSE:
- `ΔRMSE ≈ 3.35e-8` (negligible at current scan precision).

Relative coefficient offsets:
- `|c1_free - 23/12| / (23/12) ≈ 0.017%`
- `|c2_free + 98/9| / (98/9) ≈ 0.0067%`

## Interpretation

At current resolution, the structural rational form is effectively indistinguishable from the free quadratic fit.

This is a strong candidate lock, but not yet a formal proof:
- data depth is short (7 anchors),
- higher-order and threshold effects remain visible in alternate parameterizations,
- coefficient-lock should be re-tested after scan densification.

## Structural mapping candidate

- `23 = 11*N_c - 2*n_f` at `N_c=3, n_f=5`
- `98 = 2*(7²)` with `7 = beta0(n_f=6)` one-loop coefficient
- denominators `12` and `9` suggest gauge/count and color-square normalizations.

This keeps the working hypothesis physically grounded in asymptotic-freedom coefficients while remaining falsifiable.
