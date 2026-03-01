# Finding 278 — RH Final Gap Equivalence + Direct Closure

Date: 2026-03-01  
Scope: Make the last unresolved RH obligation maximally explicit and prove direct closure from it.

## Updated

- `lean/Gutoe/RiemannTargetFiniteLadder.lean`

## New definitions/theorems

Added:

- `XiTargetLadderZeroCapture (specN : ℕ → Finset ℝ) : Prop`
  - explicit statement: every `XiTarget` zero appears as
    `criticalLinePoint t` with `t ∈ specN N` for some finite level `N`.

- `approxZero_tolZero_iff_zeroCapture`
  - proves equivalence:
    - `ApproxZeroConvergence XiTarget (XiFiniteLadder specN) tolZero`
      iff
    - `XiTargetLadderZeroCapture specN`.

- `mathlibRH_of_target_ladder_zero_capture`
  - direct theorem:
    - from `XiTargetLadderZeroCapture specN`, derive Mathlib’s `RiemannHypothesis`.

## Why this matters

This removes contract indirection from the final unresolved gap:

- the remaining open statement is now exactly one concrete capture theorem,
- and RH follows directly from it via a compiled Lean theorem.

In other words, “what remains to prove” is no longer abstract plumbing but a single explicit proposition.

## Build verification

Executed:

```bash
cd lean
lake build Gutoe.RiemannTargetFiniteLadder
lake build Gutoe
```

Result: both **pass** (warnings only).

