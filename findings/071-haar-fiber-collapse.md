# Finding 071: GRAND-311 Fiber Collapse + Parity Closure

Date: 2026-02-26
Status: GRAND-311 complete

## New module

- `lean/Gutoe/HaarFiberCollapse.lean`

## What landed

### 1) Gauge invariance -> fiber constancy

- `GaugeInvariantUnder`
- `gauge_invariant_implies_fiber_constant_of_surjective_action`

This formalizes the Path-2 claim in reusable form:
if subgroup fiber moves are covered by gauge transformations (`ρ : H →* G`), then
any gauge-invariant observable is constant on subgroup/coset fibers.

Also included:
- `FactorsThroughQuotient`
- `fiber_constant_of_factors_through_quotient`
- `center_fiber_constant_of_factorization`

### 2) Fiber/coset normalization collapse

- `normalizedExpectation`
- `normalized_expectation_scale_cancel`
- `normalized_expectation_collapse_of_common_factor`
- `normalized_expectation_reduce_to_center`

This proves the exact cancellation rule used in Path-2:
if both observable integral and total mass carry the same fiber scalar factor `c`,
normalized expectations are identical between full and center sectors.

### 3) Parity closure with transfer lane

- `finite_parity_bridge`

This theorem is an explicit closure hook back to:
- `HaarExpectationDecomposition.finite_parity_with_transfer_lane`
- `YangMillsWilsonBridge.finite_fiber_expectation_collapse`

So the continuous collapse lane and finite transfer lane share one collapse structure.

## Explicit assumptions remaining (recorded)

The normalized collapse theorem requires explicit hypotheses:
- common-factor integral relation `hInt`
- common-factor total-mass relation `hMass`
- nondegeneracy `c ≠ 0` and quotient mass nonzero

These are now isolated assumptions, not hidden axioms.

## Build sanity

- `lake build Gutoe.HaarFiberCollapse` ✅
- `lake build Gutoe` ✅

No `sorry` introduced.
