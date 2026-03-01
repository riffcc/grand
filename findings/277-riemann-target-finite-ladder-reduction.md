# Finding 277 — RH Target Finite-Ladder Reduction

Date: 2026-03-01  
Scope: Collapse final RH contract to a single unresolved assumption field by discharging structural obligations automatically.

## Added

New module:

- `lean/Gutoe/RiemannTargetFiniteLadder.lean`

Wired into roots:

- `lean/lakefile.lean` now includes `Gutoe.RiemannTargetFiniteLadder`.

## New reduction layer

Defined:

- `XiFiniteLadder : (ℕ → Finset ℝ) → (ℕ → (ℂ → ℂ))`
- `tolZero : ℕ → ℝ`
- `XiTargetFiniteLadderContract` with fields:
  - `specN : ℕ → Finset ℝ`
  - `approxZero : ApproxZeroConvergence XiTarget (XiFiniteLadder specN) tolZero`

Automatically proved inside module:

- finite bridge family for ladder (`finiteBridgeFamily_XiFiniteLadder`)
- nonnegative tolerance (`zeroTol_tolZero`)
- rigidity at zero tolerance (`rigidity_tolZero`)
- builder from reduced contract to full `RHConvergenceTransferContract XiTarget`
- final closure theorem:
  - `mathlibRH_of_target_finite_ladder_contract : XiTargetFiniteLadderContract → RiemannHypothesis`

## Why this matters

The final unresolved RH gap in this lane is now a single explicit field:

- `approxZero` for `XiTarget` against the finite spectral ladder at zero tolerance.

Everything else (bridge, rigidity, tolerance plumbing, target transfer) is discharged in Lean.

## Build verification

Executed:

```bash
cd lean
lake build Gutoe.RiemannTargetFiniteLadder
lake build Gutoe
```

Result: both **pass** (warnings only).

