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
- `maxRowTotal_le_coordination_of_sc_bound`
- `maxRowTotal_eq_coordination_of_sc_regular`
- `minorization_eps_ge_sc_coordination_floor`
- `minorization_eps_eq_sc_regular`

In `YangMillsContinuumSurvival`:
- `uniform_eps_floor_of_sc_regular_schedule`
- `continuum_hypotheses_of_sc_regular_schedule`

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
continuum-survival hypothesis package.

The remaining practical bridge is proving the empirical transfer construction
satisfies `SCRegularRowTotals` at each refinement step.

## Build verification

- `cd lean && lake build Gutoe.YangMillsStructuralGap` ✅
- `cd lean && lake build Gutoe.YangMillsContinuumSurvival` ✅
- `cd lean && lake build Gutoe` ✅

