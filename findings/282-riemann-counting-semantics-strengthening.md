# Finding 282 — RH Counting Semantics Strengthening

Date: 2026-03-01  
Scope: Strengthen Weyl endgame contract so `N_H` / `N_ξ` are explicitly counting functions of concrete point sets.

## Updated

- `lean/Gutoe/RiemannWeylEndgame.lean`

## New semantic layer

Added:

- `StepCountingSemantics (N : ℝ → ℝ) (S : Set ℝ)`
  - `cut : ℝ → Finset ℝ`
  - `cut_spec : t ∈ cut T ↔ t ∈ S ∧ 0 ≤ t ∧ t ≤ T`
  - `count_def : N T = (cut T).card`

Contract now carries:

- explicit point sets: `poleHSet`, `poleXiSet`
- explicit counting semantics: `counting_H`, `counting_xi`

instead of treating `N_H`, `N_ξ` as unconstrained black-box functions.

## Semantic bridge decomposition

`RiemannWeylIdentityContract` now requires explicit bridge fields:

1. `zero_to_poleXi`
2. `poleXi_to_poleH` (under `m_identity`)
3. `poleH_to_ordinate`

and theorem:

- `ordinateEnumeration_of_semantic_bridge`

derives ordinate enumeration from those fields.

## Why this matters

This directly addresses the abstraction-gap concern:

- counting functions are now tied to concrete finite cuts of point sets,
- bridge obligations are semantic and explicit,
- closure remains compiled to RH.

## Build verification

Executed:

```bash
cd lean
lake build Gutoe.RiemannWeylEndgame
lake build Gutoe
```

Result: both **pass** (warnings only).

