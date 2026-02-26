# Finding 070: GRAND-310 Expectation Decomposition Landed

Date: 2026-02-26
Status: GRAND-310 complete

## New module

- `lean/Gutoe/HaarExpectationDecomposition.lean`

## What was proven

### Continuous decomposition (Path-2 core)

1. `expectation_decomposition_over_subgroup`
   - The formal decomposition theorem:
     - `E_G[f] = E_{G⧸Γ}[E_fiber[f]]`
   - Lean form uses:
     - `expectation`
     - `fiberExpectation = QuotientGroup.automorphize`
     - `quotientFiberMeasure`
   - Proof is an exact bridge to `integral_unfolding_over_quotient` from GRAND-309.

2. `expectation_decomposition_over_center`
   - Center-specialized corollary:
     - `E_G[f] = E_{G⧸Z(G)}[E_fiber[f]]`
   - This is the statement needed for the SU(3)/Z3 bridge lane.

### Finite analog parity check

3. `finite_parity_with_transfer_lane`
   - Defines finite LHS/RHS expectation forms and proves exact equality by reducing to:
     - `YangMillsWilsonBridge.finite_fiber_expectation_collapse`
   - Confirms continuous decomposition shape is coherent with the already-proven transfer-lane collapse theorem.

## Build sanity

- `lake build Gutoe.HaarExpectationDecomposition` ✅
- `lake build Gutoe` ✅

No `sorry` introduced.

## Board implication

- GRAND-310: done
- GRAND-311 remains: finalize explicit gauge-invariant fiber-constancy/coset-collapse theorem package and parity closure narrative.
