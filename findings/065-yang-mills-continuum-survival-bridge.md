# Finding 065 — Yang-Mills Continuum-Survival Bridge (Lean)

Date: 2026-02-26

Scope: GRAND-299

## Result

Added a formal continuum-survival bridge that turns the structural Doeblin
mass-gap inequality into a **uniform non-vanishing lower bound** across a
joint refinement schedule.

File:
- `lean/Gutoe/YangMillsContinuumSurvival.lean`

## New Lean artifacts

- `ContinuumSurvivalHypotheses`
- `doeblin_gap_lower_bound_of_floor_cap`
- `continuum_survival_gap_nonvanishing`

## What is proven

Given a refinement schedule `n ↦ (a_t n, eps n)`, if:
- `a_t n > 0` for all `n`,
- `0 < eps n < 1` for all `n`,
- there exists a uniform positive floor `epsFloor` with `epsFloor ≤ eps n`,
- there exists a uniform positive cap `aCap` with `a_t n ≤ aCap`,

then Lean proves:

`∃ c > 0, ∀ n, c ≤ doeblinGapLowerBound (a_t n) (eps n)`.

So the Doeblin gap lane cannot collapse to zero under these explicit structural
bounds.

## Interpretation

This closes GRAND-299 as an explicit theorem-level bridge from Theorem A
structural positivity to a non-vanishing continuum-survival statement, while
keeping all assumptions visible and auditable.

## Build verification

- `cd lean && lake build Gutoe.YangMillsContinuumSurvival` ✅
- `cd lean && lake build Gutoe` ✅

