# GRAND — Grammar-Universe Closure

Date: 2026-02-28

## Claim Proven

For each constant lane (`p`, `kappa`, `ew`), uniqueness now holds for **every subgrammar**
of an explicitly declared finite supergrammar universe (provided the subgrammar includes the candidate tuple).

Formalized in:

- `lean/Gutoe/TriangulatedGrammarUniverse.lean`

Core theorem:

- `all_subgrammars_unique_if_candidate_included`

This theorem gives, for any subgrammar `G` under each universe:

1. existence of a near-solution (if candidate included),
2. uniqueness of that near-solution inside `G`.

## Supergrammar Universes

- `pUniverse`: coefficients in `[-8,8] × [-8,8]` over
  `137/10 + a*(1/(7*12)) + b*(1/(7*13*136))`
- `kappaUniverse`: coefficients in `[-8,8]^3` over
  `(60/11)*(a*(19/3) + b*(1/36) + c*(1/(7*13*136)))`
- `ewUniverse`: `(m,b,c)` with `m∈[0,2]`, `b,c∈[-8,8]` over
  `m*8 + b*(6/13) + c*(1/(7*136))`

Singleton near sets were proven in each universe.

## Build Verification

- `lake build Gutoe.TriangulatedGrammarUniverse` passes.
- `lake build Gutoe` was previously green with this module integrated.

## Boundary

This is “every grammar” quantified over all subgrammars of the declared supergrammar.
If the supergrammar is expanded, closure must be re-proven on the expanded universe.
