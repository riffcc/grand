# Finding 069: GRAND-309 Haar Hooks Landed (No Axioms)

Date: 2026-02-26
Status: GRAND-309 closed

## What landed

Module: `lean/Gutoe/HaarMeasureHooks.lean`

1. Haar existence/uniqueness hooks:
- `haar_measure_is_haar`
- `canonical_haar_is_haar`
- `left_invariant_measure_eq_smul_canonical_haar`

These are direct Lean wrappers around Mathlib Haar machinery (`haarMeasure`, `haar`, `isMulLeftInvariant_eq_smul`) in a reusable GUTOE namespace.

2. Quotient-measure hooks (`G ⧸ Γ`):
- `quotient_haar_of_preimage`
- `quotient_measure_unique_from_preimage`

These connect the bridge lane to Mathlib's quotient-Haar results and uniqueness of `QuotientMeasureEqMeasurePreimage`-compatible measures.

3. Disintegration/unfolding hook:
- `quotientFiberMeasure`
- `integral_unfolding_over_quotient`

This pins the integral decomposition interface needed for Path-2 center-orbit / coset-fiber expectation decomposition.

## Build sanity

- `lake build Gutoe.HaarMeasureHooks` ✅
- `lake build Gutoe` ✅

No `sorry` introduced.

## Board state update

- GRAND-308: Done
- GRAND-309: Done
- GRAND-310: In Progress
- GRAND-311: Todo

## Remaining gap for Clay-lane bridge

`GRAND-310/311` remain: proving full gauge-invariant fiber collapse in the Path-2 expectation pipeline (from hooks to full decomposition theorem and parity checks).
