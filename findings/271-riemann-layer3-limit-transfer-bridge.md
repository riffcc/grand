# Finding 271 — RH Layer-3 Limit-Transfer Bridge (Lean)

Date: 2026-03-01  
Scope: Encode Steps 4/5 boundary explicitly: transfer from finite exact bridge family to RH for target `Xi`.

## Added

New module:

- `lean/Gutoe/RiemannLimitBridge.lean`

Wired into roots:

- `lean/lakefile.lean` now includes `Gutoe.RiemannLimitBridge`.

Integrated in closure:

- `lean/Gutoe/RiemannRHClosure.lean` imports/exports Layer-3 packaging theorem.

## New theorem layer

`RiemannLimitBridge` introduces explicit transfer theorems:

- `rh_of_limit_transfer`
  - assumptions:
    1. every finite level `XiN N` has exact `SpectralBridge` to `specN N`,
    2. every target-`Xi` zero appears as a zero in some finite level (`hzeroForward`).
  - conclusion: `RiemannHypothesisXi Xi`.

- `spectralBridge_of_limit_transfer`
  - adds backward transfer (`XiN` zero implies `Xi` zero),
  - conclusion: exact infinite bridge `SpectralBridge Xi (ladderSpec specN)`.

- `rh_of_exact_limit_bridge`
  - RH via the exact bridge obtained above.

And closure-level wrapper in `RiemannRHClosure`:

- `rh_from_limit_transfer`.

## Why this matters

This formalizes the exact endgame gap as assumptions, not prose:

- finite bridge family is already constructive in the finite model lane;
- remaining analytic burden is now a named transfer obligation from finite levels to target `Xi`.

In practical terms, Steps 4/5 are now represented by concrete Lean hypotheses instead of implicit narrative.

## Build verification

Executed:

```bash
cd lean
lake build Gutoe
```

Result: **passes** (`8168` jobs, warnings only).

