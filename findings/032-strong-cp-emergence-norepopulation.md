# 032 — Strong CP Emergence No-Repopulation Bridge (GRAND-267)

Status: theorem-level bridge added; final closure still requires model-complete emergent-SU(3) assumptions.

## New Lean module

- `lean/Gutoe/StrongCPEmergence.lean`

Key theorems:

- `no_repopulation_on_emergent_image`
- `theta_phase_unity_of_zero_charge`
- `theta_phase_unity_on_emergent_image`
- `coarse_grain_cannot_create_nontrivial_sector`
- `theta_phase_unity_of_coarse_grain_no_creation`

## What this proves

Given:

1. effective topological charge is a pullback from fundamental Z3-carrier maps,
2. fundamental constant maps carry zero charge,
3. and route-1 continuity result (fundamental continuous maps are constant),

then emergent-image topological charge is forced to zero, and therefore
`exp(i θ Q) = 1` on that image for all θ.

## Why this matters

This directly attacks the “emergent SU(3) could repopulate sectors” concern:
it provides a formal no-repopulation theorem under explicit pullback hypotheses
rather than relying on finite-lattice handwaving.

Additionally, it now includes a pointwise coarse-graining theorem:
if `CG f = φ ∘ f` and effective charge vanishes on constant fields, then
coarse-graining cannot create nontrivial effective sectors from the fundamental
Z3-carrier field space.

## Remaining gap

The still-open part is to derive the pullback/zero-normalization hypotheses from
the full GUTOE emergent-SU(3) construction itself (instead of assuming them as
bridge conditions).
