# Finding 062 — Explicit Doeblin Mass-Gap Lower Bound (Lean)

Date: 2026-02-26

Scope: GRAND-297 / GRAND-298 bridge hardening

## Added proofs

File:
- `lean/Gutoe/YangMillsMassGap.lean`

New definitions/theorems:
- `doeblinGapLowerBound`
- `mass_gap_ge_doeblin_bound`
- `doeblin_bound_positive`
- `mass_gap_positive_of_doeblin_ratio`

## Result

Lean now proves the explicit inequality chain:

- If `0 < eps < 1` and `lambda1/lambda0 ≤ 1 - eps`, then
  `doeblinGapLowerBound a_t eps ≤ massGapFromEigenRatio a_t lambda0 lambda1`.
- Since `doeblinGapLowerBound a_t eps = -log(1-eps)/a_t` and `0 < eps < 1`,
  this lower bound is strictly positive.
- Therefore `massGapFromEigenRatio > 0`.

This is the exact analytic form:

`m_gap ≥ -ln(1-ε)/a_t > 0`.

## Why this matters

We now have an explicit, computable lower bound object in Lean tied to the minorization parameter.
This turns the final step from qualitative positivity into a quantitative inequality target.

## Build verification

- `cd lean && lake build Gutoe.YangMillsMassGap` ✅
- `cd lean && lake build Gutoe` ✅
