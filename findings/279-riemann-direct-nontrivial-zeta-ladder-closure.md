# Finding 279 — RH Direct Nontrivial-ζ Ladder Closure

Date: 2026-03-01  
Scope: Peel one more layer by bypassing `XiTarget` and closing RH directly from nontrivial `riemannZeta` ladder capture.

## Updated

- `lean/Gutoe/RiemannTargetFiniteLadder.lean`

## New final-gap predicate

- `RiemannNontrivialLadderZeroCapture (specN : ℕ → Finset ℝ) : Prop`

Statement:

- every nontrivial `riemannZeta` zero (`ζ(s)=0`, not trivial negative-even, `s ≠ 1`)
  is captured as `s = criticalLinePoint t` with `t ∈ specN N` for some finite level `N`.

## New direct closure theorem

- `mathlibRH_of_riemann_nontrivial_ladder_capture`

This proves Mathlib’s `RiemannHypothesis` directly from the predicate above,
without routing through `XiTarget`/contract objects.

## Why this matters

The unresolved gap is now expressed in the most direct RH-native form available in this lane:

- a single nontrivial-`ζ` finite-ladder capture theorem.

Everything else is compiled closure.

## Build verification

Executed:

```bash
cd lean
lake build Gutoe.RiemannTargetFiniteLadder
lake build Gutoe
```

Result: both **pass** (warnings only).

