# Finding 076: GRAND-303 Compact-Simple Gauge Scope Generalization

Date: 2026-02-26
Status: GRAND-303 complete

## Updated module

- `lean/Gutoe/YangMillsGaugeScope.lean`

## What landed

### 1) GRAND-312 normalization bridge lifted to finite-center compact-group scope

New theorem:
- `normalized_expectation_reduce_to_center_of_finite_center`

This lifts the normalized expectation reduction from the SU(3)-named lane into an
abstract group scope with finite center, by deriving the required countability
instance and reusing the GRAND-312 theorem.

### 2) Full Path-2 package theorem at compact-simple scope

New theorem:
- `compact_simple_scope_supports_full_path2`

Under `CompactSimpleGaugeScope G` plus explicit measurable/quotient hypotheses,
it jointly returns:

1. expectation decomposition over `G ⧸ Z(G)`
2. normalized expectation reduction (from quotient-normalization data)
3. non-vanishing continuum Doeblin gap lower bound

This gives a single theorem-level generalization artifact for Clay-scope review.

## Why this matters

- The argument is no longer framed only in SU(3)-specific naming.
- The load-bearing Haar + normalization + continuum-gap pieces are now stated at
  compact-simple finite-center group scope.
- Residual assumptions are explicit theorem arguments (no hidden axioms).

## Build sanity

- `lake build Gutoe.YangMillsGaugeScope` ✅
- `lake build Gutoe` ✅

No `sorry` introduced.
