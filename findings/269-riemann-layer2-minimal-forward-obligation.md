# Finding 269 — RH Layer-2 Minimal Forward Obligation

Date: 2026-03-01  
Scope: Reduce Layer-2 proof target to the minimal obligation actually required for RH.

## Change

Updated:

- `lean/Gutoe/RiemannLayer2Identity.lean`
- `lean/Gutoe/RiemannRHClosure.lean`

## New theorem surfaces

In `RiemannLayer2Identity`:

- `rh_of_zero_to_finite_witness`

Statement (informal):

- If every `Xi` zero has a finite-ladder witness
  (`Xi s = 0 → ∃ N t, t ∈ specN N ∧ s = 1/2 + i t`),
  then RH-for-`Xi` follows.

In `RiemannRHClosure`:

- `rh_from_layer2_forward_only`

This packages the same minimal trigger at closure-layer API level.

## Why this matters

Layer-2 previously encoded a two-way identity (forward + backward).  
That is still useful for exact bridge equivalence, but RH itself only needs the forward side.

This splits goals cleanly:

1. **RH proof goal** (minimal): establish forward witness map.
2. **Exact Xi↔Spec(H) equivalence goal** (stronger): add backward soundness.

The remaining decisive frontier is now explicit and minimal.

## Build verification

Executed:

```bash
cd lean
lake build Gutoe
```

Result: **passes** (`8166` jobs, warnings only).

