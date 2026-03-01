# Finding 273 — RH Convergence-Transfer Contract (Lean)

Date: 2026-03-01  
Scope: Encode “finite spectrum converges to `Xi` zeros” as explicit theorem obligations that derive `ZeroForwardTransfer`.

## Added

New module:

- `lean/Gutoe/RiemannConvergenceTransfer.lean`

Wired into roots:

- `lean/lakefile.lean` now includes `Gutoe.RiemannConvergenceTransfer`.

Integrated into closure:

- `lean/Gutoe/RiemannRHClosure.lean` imports the module and adds wrapper theorem.

## New obligation layer

In `RiemannConvergenceTransfer`:

- `ApproxZeroConvergence Xi XiN tol`
  - each `Xi` zero is tolerance-close to some finite-level evaluation.
- `SpectralRigidity XiN tol`
  - tolerance-close at level `N` implies exact zero at level `N`.
- `zeroForward_of_convergence_and_rigidity`
  - derives `ZeroForwardTransfer Xi XiN`.

Contract:

- `RHConvergenceTransferContract Xi`
  - fields:
    - finite bridge family (`FiniteBridgeFamily`)
    - tolerance profile
    - approximate zero convergence
    - spectral rigidity

Closure:

- `rh_of_convergence_transfer_contract`
- and closure API in `RiemannRHClosure`:
  - `rh_from_convergence_transfer_contract`

## Why this matters

This turns the Step-4 frontier into a machine-checkable interface:

- no need to assume `ZeroForwardTransfer` directly;
- it is now a theorem consequence of two explicit analytic obligations
  (convergence + rigidity).

That is exactly the “prove finite operator spectrum converges to `ξ` zeros” lane, formalized.

## Build verification

Executed:

```bash
cd lean
lake build Gutoe
```

Result: **passes** (`8168` jobs, warnings only).

