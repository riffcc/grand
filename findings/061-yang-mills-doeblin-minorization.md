# Finding 061 — Yang-Mills Doeblin Minorization Layer (Lean)

Date: 2026-02-26

Scope: GRAND-297 / Theorem A structural closure track

## Added proofs (no `sorry`)

All additions are in:
- `lean/Gutoe/YangMillsStructuralGap.lean`

### 1. Row-stochastic normalization from counted transitions

- `smoothed_transition_row_sum_one`

Given:
- Laplace smoothing with `alpha > 0`
- `rowTotals i = ∑ j counts i j`

Proves each smoothed row sums to exactly `1`.

### 2. Uniform global positivity floor

- `maxRowTotal`
- `rowTotal_le_maxRowTotal`
- `laplaceGlobalFloor`
- `laplace_global_floor_pos`
- `smoothed_transition_entry_ge_global_floor`

Proves every smoothed entry has an explicit, construction-derived lower bound:

`P(i,j) ≥ alpha / (maxRowTotal + 3*alpha) > 0`.

### 3. Explicit minorization constant range

- `laplace_global_floor_le_one_third`
- `minorizationEps`
- `minorization_eps_range`

Defines `minorizationEps = 3 * laplaceGlobalFloor` and proves:

- `0 < minorizationEps`
- `minorizationEps ≤ 1`

This is a formal Doeblin-style one-step minorization witness for the 3-state transfer basis.

### 4. Full algebraic Doeblin decomposition witness

- `uniformKernel`
- `residualKernel`
- `doeblin_decomposition`

For `alpha > 0`, normalized rows, and strict `minorizationEps < 1`, Lean now proves existence of
`R` such that:

- `P = εU + (1-ε)R`
- `R` is entrywise nonnegative
- each row of `R` sums to `1`

This gives an explicit finite-state decomposition suitable for downstream contraction/spectral bounds.

## Why this matters

This moves the YM gap story from:
- “entrywise positive matrices are primitive”

to:
- “the smoothed kernel is row-stochastic and has an explicit uniform minorization constant.”

That is the right pre-spectral structure for contraction/mixing and downstream spectral-gap arguments.

## Build verification

- `cd lean && lake build Gutoe.YangMillsStructuralGap` ✅
- `cd lean && lake build Gutoe` ✅
