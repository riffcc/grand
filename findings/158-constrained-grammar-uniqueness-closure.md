# GRAND — Constrained Grammar Uniqueness Closure

Date: 2026-02-28

## What closed

Added formal finite-grammar uniqueness closure in:

- `lean/Gutoe/TriangulatedClosureUniqueness.lean`

Registered in:

- `lean/lakefile.lean`

## Theorem content

Defined finite coefficient grammars and proved uniqueness of the selected combination inside frozen windows:

- `p` grammar: `s ∈ {-1,0,1}` for `137/10 + s*(1/(7*12))`
- `kappa` grammar: `(a,b,c) ∈ {0,1}^3` for
  `(60/11) * (a*(19/3) + b*(1/36) + c*(1/(7*13*136)))`
- `ew` grammar: `(b,s) ∈ {0,1} × {-1,0,1}` for
  `8 + b*(6/13) + s*(1/(7*136))`

Uniqueness proofs:

- `p_good_signs_unique`
- `kappa_good_tuples_unique`
- `ew_good_tuples_unique`
- `constrained_grammar_uniqueness_closure`
- `constrained_grammar_selected_coefficients`

Recovered unique coefficients:

- `s_p = -1`
- `(a,b,c) = (1,1,1)`
- `(b_ew,s_ew) = (1,-1)`

## Verification

- `cd lean && lake build Gutoe`
- Success: `Gutoe.TriangulatedClosureUniqueness` built
- Full build green (`8131 jobs`)

## Honesty boundary

This is formal uniqueness **within the explicitly declared finite grammar**.
If grammar is expanded, uniqueness must be reproven under the enlarged search space.
