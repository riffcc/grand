# Finding 067 — Yang-Mills Wilson→Z3 Structural Bridge (Lean)

Date: 2026-02-26

Scope: GRAND-300 (Theorem C bridge lane)

## Result

Added a structural Wilson-action bridge module that maps a Wilson schedule into
Z3 nearest-neighbor transfer targets and inherits the SC-regular continuum
survival mass-gap lane with no empirical row-total hypothesis.

Files:
- `lean/Gutoe/YangMillsWilsonBridge.lean`
- `lean/lakefile.lean`

## New Lean definitions/theorems

In `YangMillsWilsonBridge`:
- `WilsonZ3Action`
- `wilsonRowTotalsSchedule`
- `wilson_row_totals_sc_regular`
- `wilson_max_row_total_eq_coordination`
- `wilson_minorization_eps_closed_form`
- `wilson_action_bridge_nonvanishing_gap`

## What is proven

1. Any Wilson schedule represented as refinement-indexed Z3 nearest-neighbor
   targets induces SC-regular row totals at every step.

2. The induced max row total is exactly the SC coordination number (`6`) for
   each refinement index.

3. The induced Doeblin minorization constant has the closed form

   `(3*alpha)/((coordinationNumber:ℝ)+3*alpha)`.

4. Therefore the non-vanishing continuum mass-gap lower-bound theorem from the
   structural Yang-Mills chain applies directly to this Wilson bridge lane.

## Interpretation

This is a structural Theorem-C bridge slice: it removes dependence on empirical
`maxRowTotal` certificates for the represented Wilson schedule class and routes
through the already-proven Z3 nearest-neighbor SC-regular schedule theorems.

Remaining work is still the full semantic equivalence proof between this Wilson
representation lane and the conventional Wilson-action SU(3) lattice path
integral (operator/action-level correspondence), but the transfer-gap bridge is
now theoremized end-to-end in this lane.

## Build verification

- `cd lean && lake build Gutoe.YangMillsWilsonBridge` ✅
- `cd lean && lake build Gutoe` ✅
