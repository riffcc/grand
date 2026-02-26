# Finding 066 — Yang-Mills Uniform Epsilon Floor from SC Coordination (Lean)

Date: 2026-02-26

Scope: GRAND-305

## Result

Added a structural epsilon-floor lane that ties the Doeblin minorization constant
directly to Cl(1,3)-forced simple-cubic coordination (`6`).

Files:
- `lean/Gutoe/YangMillsStructuralGap.lean`
- `lean/Gutoe/YangMillsContinuumSurvival.lean`

## New Lean theorems

In `YangMillsStructuralGap`:
- `SCRegularRowTotals`
- `z3LocalMultiplicity`
- `z3_local_multiplicity_eq_two`
- `z3CanonicalLocalCounts`
- `rowTotalsFromCounts`
- `Z3SCLocalRegularCounts`
- `z3_canonical_local_counts_regular`
- `z3CanonicalRowTotals`
- `z3_canonical_row_totals_sc_regular`
- `sc_regular_row_totals_of_z3_local_regular_counts`
- `z3_canonical_row_total_eq_coordination`
- `maxRowTotal_le_coordination_of_sc_bound`
- `maxRowTotal_eq_coordination_of_sc_regular`
- `minorization_eps_ge_sc_coordination_floor`
- `minorization_eps_ge_bounded_max_row_total_floor`
- `minorization_eps_eq_sc_regular`

In `YangMillsContinuumSurvival`:
- `uniform_eps_floor_of_sc_regular_schedule`
- `uniform_eps_floor_of_bounded_max_row_total_schedule`
- `continuum_hypotheses_of_sc_regular_schedule`
- `continuum_hypotheses_of_bounded_max_row_total_schedule`
- `rowTotalsScheduleFromCountsSchedule`
- `sc_regular_schedule_of_z3_local_regular_counts`
- `uniform_eps_floor_of_z3_local_regular_schedule`
- `z3CanonicalRowTotalsSchedule`
- `z3_canonical_schedule_sc_regular`
- `uniform_eps_floor_of_z3_canonical_schedule`
- `continuum_hypotheses_of_z3_canonical_schedule`
- `continuum_survival_gap_nonvanishing_of_z3_canonical_schedule`
- `empiricalMaxRowTotalBound` (= `846`)
- `uniform_eps_floor_of_empirical_bound_846`
- `continuum_survival_gap_nonvanishing_of_empirical_bound_846`

## What is proven

1. If row totals are SC-bounded (`rowTotals i ≤ coordinationNumber`), then:

`(3*alpha)/((coordinationNumber:ℝ)+3*alpha) ≤ minorizationEps rowTotals alpha`.

2. If row totals are SC-regular (`rowTotals i = coordinationNumber`), then:

`minorizationEps rowTotals alpha = (3*alpha)/((coordinationNumber:ℝ)+3*alpha)`.

3. For a refinement schedule with fixed `alpha > 0` and SC-regular rows at each
step, there exists a **uniform positive epsilon floor** independent of refinement
index.

## Interpretation

This discharges the previous generic epsilon-floor gap at the theorem level for
the SC-regular transfer schedule class, and plugs directly into the
continuum-survival hypothesis package, including a canonical Z3/SC non-vanishing
mass-gap lower-bound corollary.

The remaining practical bridge is proving the empirical transfer construction
satisfies `Z3SCLocalRegularCounts` at each refinement step (or is bounded by a
schedule that does), so the local-regular bridge theorems can be applied
directly to measured lanes.

Layer-1 empirical lane now has an explicit theorem target: if measured
refinement rows satisfy `maxRowTotal ≤ 846` and remain positive, the
continuum-survival non-vanishing gap corollary follows with a computable floor.

## Build verification

- `cd lean && lake build Gutoe.YangMillsStructuralGap` ✅
- `cd lean && lake build Gutoe.YangMillsContinuumSurvival` ✅
- `cd lean && lake build Gutoe` ✅
