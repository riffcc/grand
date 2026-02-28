# 161 — Infinite Elimination Canonical Closure

Date: 2026-02-28

## Result

The infinite-to-finite elimination lane is now closed end-to-end in Lean for all three triangulated constants.

## New closure theorems

In `Gutoe.TriangulatedGrammarUniverse`:

- `p_index_window_xl_singleton`
- `p_xl_bounds_and_index_select_candidate`
- `p_infinite_elimination_xl_complete`

- `kappa_index_window_xl_singleton`
- `kappa_xl_bounds_and_index_select_candidate`
- `kappa_infinite_elimination_xl_complete`

- `ew_index_exact_xl_singleton`
- `ew_xl_bounds_and_index_select_candidate`
- `ew_infinite_elimination_xl_complete`

## Meaning

Given XL structural bounds plus near-target condition:

- `p` collapses to `(-1, 0)`
- `kappa` collapses to `(1, (1, 1))`
- `ew` collapses to `(1, (1, -1))`

The proofs now explicitly pass through index elimination windows, then canonical tuple selection.

## Verification

- `lake build Gutoe.TriangulatedGrammarUniverse` passes
- `lake build Gutoe` passes

No theorem-direction/sign changes were made.
