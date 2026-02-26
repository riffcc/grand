# Finding 063 — Doeblin Subdominant Eigenvalue Bound (Lean)

Date: 2026-02-26

Scope: GRAND-297 / final analytic bridge to explicit gap bound

Note: this finding captures the contraction inequality stage. The external
`|mu| ≤ 1` assumption used here is discharged structurally in
`findings/064-yang-mills-theorem-a-structural-gap-closure.md`.

## Added theorem

File:
- `lean/Gutoe/YangMillsStructuralGap.lean`
- `lean/Gutoe/YangMillsMassGap.lean`

New theorems:
- `uniformKernel_mulVec_eq_avg`
- `uniformKernel_mulVec_zero_of_sum_zero`
- `abs_eigenvalue_le_one_sub_eps_of_decomposition`
- `eigenvalue_le_one_sub_eps_of_decomposition`
- `mass_gap_positive_of_doeblin_mode`

## Statement proved

For a decomposition

`P = eps * U + (1-eps) * R`

with:
- `0 < eps < 1` (implemented as `eps < 1`; nonnegativity follows in the target application),
- `v` a zero-sum mode (`v 0 + v 1 + v 2 = 0`),
- `P * v = lam * v`,
- `R * v = mu * v`,
- `|mu| ≤ 1`,

Lean now proves:

`|lam| ≤ 1 - eps`.

This is the exact subdominant-mode contraction inequality needed for the Doeblin mass-gap route.

## Connection to existing gap theorem

`lean/Gutoe/YangMillsMassGap.lean` already had:
- `mass_gap_ge_doeblin_bound`
- `doeblin_bound_positive`
- `mass_gap_positive_of_doeblin_ratio`

Together, the chain is now explicit:

`|lam| ≤ 1-eps  ⇒  m_gap ≥ -log(1-eps)/a_t > 0`.

And in Lean this is directly packaged in:

- `mass_gap_positive_of_doeblin_mode`

which consumes:
- a Doeblin decomposition on `P`,
- a zero-sum mode eigenpair for `P`,
- a bounded residual-mode eigenvalue `|mu| ≤ 1`,
- `0 < eps < 1`, `0 < lam`, `0 < a_t`,

and returns:

`0 < massGapFromEigenRatio a_t 1 lam`.

## Build verification

- `cd lean && lake build Gutoe.YangMillsStructuralGap` ✅
- `cd lean && lake build Gutoe.YangMillsMassGap` ✅
- `cd lean && lake build Gutoe` ✅
