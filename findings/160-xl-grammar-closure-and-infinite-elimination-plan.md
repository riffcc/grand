# 160 — XL Grammar Closure and Infinite-Space Elimination Plan

Date: 2026-02-28

## What was expanded

Lean module `Gutoe.TriangulatedGrammarUniverse` now includes an XL supergrammar tier:

- `p`: coefficient band expanded to `[-32, 32] × [-32, 32]`
- `kappa`: coefficient band expanded to `[-32, 32]^3`
- `ew`: coefficient band expanded to `m ∈ [0, 5]`, `b,c ∈ [-48, 48]`

New closure theorem:

- `all_subgrammars_unique_if_candidate_included_xl`

This proves unique near-solution selection for any subgrammar included in the XL universes, provided the candidate tuple is present.

Also added:

- `legacy_universes_embed_in_xl`

## Validation

- `lake build Gutoe.TriangulatedGrammarUniverse` passes
- `lake build Gutoe` passes
- Rust triangulation gates still pass:
  - `triangulation_ci_gate`
  - `triangulation_clifford_candidates_ci_gate`

## Honest remaining gaps (unchanged)

1. Supergrammar closure is still over declared finite universes, not all infinite integer grammars.
2. Neutrino absolute mass scale remains off (kappa identified but not yet fully propagated).
3. Texture phase drift remains open (`GRAND-291` partial).
4. Axiom-to-reality mapping (`why Cl(1,3) is nature`) is philosophical/empirical, not formally provable inside Lean.

## Infinite-space tightening plan (elimination + ranges)

To represent infinite candidate spaces safely in finite proof effort:

1. Eliminate to lattice-index forms (linear combinations of integer coefficients).
2. Convert near-target constraints into narrow integer index intervals.
3. Apply structural range constraints from Cl(1,3) combinatorics (grade counts, orbit bounds, generator-count limits).
4. Prove uniqueness in each constrained slice, then compose slices.

This avoids brute-force infinite enumeration while keeping the proof constructive and auditable.

## Delivered now (infinite -> finite certificates)

Implemented in `Gutoe.TriangulatedGrammarUniverse`:

- `pExprQ_eq_index` with lattice index `pIndex = 442*a + 3*b`
- `p_near_forces_index_window`: near-target implies `pIndex ∈ {-442, -441}`
- `kappaExprQ_eq_index` with lattice index `kappaIndex = 3527160*a + 15470*b + 45*c`
- `kappa_near_forces_index_window`: near-target implies `3542672 ≤ kappaIndex ≤ 3542675`
- `ewExprQ_eq_index` with lattice index `ewIndex = 99008*m + 5712*b + 13*c`
- `ew_near_forces_index_exact`: near-target implies `ewIndex = 104707`

This is the core “break infinity” move: all integer tuples collapse to tiny finite index windows before grammar-range uniqueness is applied.
