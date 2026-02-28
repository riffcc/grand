# 162 — Lean Native Coupled Triangulation Closure

Date: 2026-02-28

## What was added

In `Gutoe.TriangulatedGrammarUniverse`:

- `coupled_triangulation_xl_complete`
- `coupled_triangulation_values_xl`

These theorems couple `p`, `kappa`, and `ew` natively in Lean.

## Statement-level effect

Under XL structural bounds + near-anchor assumptions for all three lanes, Lean now proves in one coupled theorem that:

- `(ap, bp) = pCandidateTuple`
- `(ak, (bk, ck)) = kappaCandidateTuple`
- `(aew, (bew, cew)) = ewCandidateTuple`

Then the second theorem lifts tuple closure to value closure:

- `pExprQ = pCandidateQ`
- `kappaExprQ = kappaCandidateQ`
- `ewExprQ = ewCoeffCandidateQ`

## Verification

- `lake build Gutoe.TriangulatedGrammarUniverse` passed.
- `lake build Gutoe` passed.

No theorem direction/sign flips were introduced.
