# Finding 274 — RH Convergence Contract End-to-End Instance (Lean)

Date: 2026-03-01  
Scope: Prove the new convergence-transfer lane is executable end-to-end by constructing a fully discharged instance.

## Updated

- `lean/Gutoe/RiemannFiniteXiModel.lean`

## New constructions/theorems

Added in `RiemannFiniteXiModel`:

- `XiFiniteConst`
- `specConst`
- `tolZero`
- `finiteBridgeFamily_XiFiniteConst`
- `zeroTol_tolZero`
- `approxZero_XiFiniteConst`
- `rigidity_XiFiniteConst`
- `XiFiniteConvergenceContract`
- `rh_XiFinite_via_convergence_contract`

## What this proves

The convergence-transfer path is no longer only a contract API:

- an explicit `RHConvergenceTransferContract` instance is constructed,
- all obligations are discharged with no `sorry`,
- RH-for-`XiFinite` follows through the convergence/rigidity route.

This validates that the new Step-4 formalization is operational in Lean.

## Build verification

Executed:

```bash
cd lean
lake build Gutoe
```

Result: **passes** (`8169` jobs, warnings only).

