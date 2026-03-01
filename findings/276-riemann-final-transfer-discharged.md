# Finding 276 — RH Final Transfer Obligation Discharged

Date: 2026-03-01  
Scope: Discharge the nontrivial-zero transfer theorem in the final RH target lane.

## Updated

- `lean/Gutoe/RiemannFinalTarget.lean`

## Key correction

`XiTarget` is set to:

- `completedRiemannZeta`

instead of `completedRiemannZeta₀`.

Reason: nontrivial-zero transfer from `riemannZeta` is naturally provable for `completedRiemannZeta`
via the exact relation

- `riemannZeta s = completedRiemannZeta s / Complex.Gammaℝ s`

away from `s = 0`.

## New proven theorem

- `nontrivialZeroTransferToXiTarget : NontrivialZeroTransferToXiTarget`

Proof ingredients:

1. `riemannZeta_zero` gives `s ≠ 0` for any `riemannZeta s = 0`.
2. `Complex.Gammaℝ_eq_zero_iff` plus nontrivial-zero exclusion gives
   `Complex.Gammaℝ s ≠ 0`.
3. `riemannZeta_def_of_ne_zero` transports zero from `riemannZeta s`
   to `completedRiemannZeta s`.

## Stronger closure now available

- `mathlibRH_of_contract` (no transfer assumption argument)

Remaining assumption for full RH in this lane is now just:

- `RHConvergenceTransferContract XiTarget`

## Build verification

Executed:

```bash
cd lean
lake build Gutoe.RiemannFinalTarget
lake build Gutoe
```

Result: both **pass** (warnings only).

