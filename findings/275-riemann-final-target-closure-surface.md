# Finding 275 — RH Final-Target Closure Surface (Lean)

Date: 2026-03-01  
Scope: Bind the contract lane to Mathlib’s actual `RiemannHypothesis` statement with explicit remaining obligation.

## Added

New module:

- `lean/Gutoe/RiemannFinalTarget.lean`

Wired into roots:

- `lean/lakefile.lean` now includes `Gutoe.RiemannFinalTarget`.

## Final target objects

- `XiTarget : ℂ → ℂ := completedRiemannZeta₀`
- `rhXiTarget_of_contract`
  - RH-for-`XiTarget` from `RHConvergenceTransferContract XiTarget`.

## Explicit final obligation

- `NontrivialZeroTransferToXiTarget`
  - transfer from nontrivial `riemannZeta` zeros to `XiTarget` zeros.

## Final closure theorem surface

- `mathlibRH_of_contract_and_transfer`
  - assumptions:
    1. convergence-transfer contract for `XiTarget`,
    2. nontrivial-zero transfer.
  - conclusion:
    - Mathlib’s `RiemannHypothesis`.

This is the direct bridge from the contract stack to the canonical RH proposition in mathlib.

## Build verification

Executed:

```bash
cd lean
lake build Gutoe
```

Result: **passes** (`8170` jobs, warnings only).

