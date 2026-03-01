# Finding 280 — RH Weyl + m-Function Endgame Contract

Date: 2026-03-01  
Scope: Encode the spectral-theory attack (`Weyl counts + m-function identity`) as a compiled Lean endgame interface.

## Added

New module:

- `lean/Gutoe/RiemannWeylEndgame.lean`

Wired into roots:

- `lean/lakefile.lean` now includes `Gutoe.RiemannWeylEndgame`.

## New formal objects

- `rvmMain : ℝ → ℝ`
  - Riemann–von Mangoldt main term.
- `RiemannVonMangoldtEnvelope`
  - asymptotic envelope around `rvmMain`.
- `HerglotzLike`
  - positivity condition for candidate m-functions.
- `MFunctionIdentity`
  - abstract equality marker for m-function identity.
- `RiemannWeylIdentityContract`
  - bundles:
    - Weyl envelopes for `N_H`, `N_ξ`,
    - exact count matching,
    - Herglotz conditions,
    - m-function identity,
    - explicit map from these to nontrivial-zero ordinate enumeration.

## Closure theorem

- `mathlibRH_of_weyl_identity_contract`
  - from `RiemannWeylIdentityContract`, derive Mathlib’s `RiemannHypothesis`.

## Why this matters

This exactly captures the intended final mathematical strategy in machine-checkable form:

1. match counting asymptotics / counts (`N_H = N_ξ`),
2. match m-function identity in the right analytic class,
3. conclude ordinate enumeration of nontrivial zeros,
4. conclude RH via already-compiled direct closure.

The unresolved content is now concentrated in one explicit contract field:
`ordinateEnumeration_of_weyl_and_m`.

## Build verification

Executed:

```bash
cd lean
lake build Gutoe.RiemannWeylEndgame
lake build Gutoe
```

Result: both **pass** (warnings only).

