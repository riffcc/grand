# Finding 064 — Yang-Mills Theorem A Structural Gap Closure (Lean)

Date: 2026-02-26

Scope: GRAND-297 / GRAND-298 (Theorem A closure)

## Result

The Yang-Mills structural gap chain is now closed in Lean without any extra spectral hypothesis on the residual mode.

Specifically, the prior external assumption `|mu| ≤ 1` has been discharged from matrix structure (`R ∈ rowStochastic`) using Gershgorin.

## New/updated theorem path

Files:
- `lean/Gutoe/YangMillsStructuralGap.lean`
- `lean/Gutoe/YangMillsMassGap.lean`

Structural lemmas and closure:
- `abs_eigenvalue_le_one_of_rowStochastic`
- `abs_eigenvalue_le_one_sub_eps_of_decomposition_stochastic`
- `eigenvalue_le_one_sub_eps_of_decomposition_stochastic`
- `mass_gap_positive_of_doeblin_mode` (now uses stochastic specialization directly)

## Closed chain (Lean)

`Cl(1,3)` basis cardinality
→ Z₃ transfer basis (`dim = 3`)
→ Laplace-smoothed positive kernel
→ row-stochastic normalization
→ Doeblin decomposition `P = eps*U + (1-eps)*R`, `eps > 0`
→ `R ∈ rowStochastic`
→ `|mu| ≤ 1` (Gershgorin, proven)
→ `|lam| ≤ 1-eps`
→ `m_gap ≥ -log(1-eps)/a_t`
→ `m_gap > 0`

## What is now removed

- No standalone `|mu| ≤ 1` input assumption in the structural mass-gap bridge.
- No numerical eigenvalue fixtures are needed for this theorem path.
- No `sorry`.

## What remains outside Theorem A

This finding closes the structural positivity theorem (Theorem A path).

Remaining macro-bridges are unchanged:
- Theorem B: continuum-survival proof of the bound.
- Theorem C: formal Wilson-action equivalence bridge.

## Build verification

- `cd lean && lake build Gutoe.YangMillsStructuralGap` ✅
- `cd lean && lake build Gutoe.YangMillsMassGap` ✅
- `cd lean && lake build Gutoe` ✅

