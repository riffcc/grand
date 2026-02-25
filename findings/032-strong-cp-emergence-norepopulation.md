# 032 — Strong CP Emergence No-Repopulation Bridge (GRAND-267)

Status: route-1 Lean closure complete (concrete `Z₃ → SU(3)` + normalized `Q[const]=0`).

## New Lean module

- `lean/Gutoe/StrongCPEmergence.lean`

Key theorems:

- `qEffFromClass`
- `qEffFromClass_homotopy_invariant`
- `pointwise_coarse_grain_preserves_homotopy`
- `no_repopulation_pointwise_class_charge`
- `theta_phase_unity_pointwise_class_charge`
- `z3ToSu3`
- `qEffFromClassAnchored`
- `qEffFromClassAnchored_const_zero`
- `su3_effective_theta_phase_unity_concrete`
- `no_repopulation_from_subsingleton_source`
- `no_repopulation_from_based_route1`
- `no_repopulation_of_homotopy_invariance`
- `no_repopulation_on_emergent_image`
- `theta_phase_unity_of_zero_charge`
- `theta_phase_unity_on_emergent_image`
- `coarse_grain_cannot_create_nontrivial_sector`
- `theta_phase_unity_of_coarse_grain_no_creation`
- `su3_effective_no_repopulation`
- `su3_effective_theta_phase_unity`

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

It also includes a stronger abstract bridge:
once the based fundamental source space is subsingleton (route-1 vacuum result),
any emergent charge functional is globally pinned by one normalized base state.

And it now discharges the two key step-6 assumptions simultaneously for the
pointwise emergence case (`CG f = φ ∘ f`) with class-induced charge:
homotopy preservation is proved by composition, and charge homotopy-invariance
is proved by quotient construction.

It also instantiates a concrete `Z₃ → SU(3)` emergence map (`z3ToSu3`) and
closes the normalization condition by construction via anchored class charge:
`qEffFromClassAnchored_const_zero` proves `Q[const] = 0` for all constant
effective fields without extra physics assumptions.

## Remaining gap

The remaining work is now interpretation-level (mapping this formal charge
normalization to whichever continuum topological-charge operator is selected in
the Rust/physics layer), not a Lean proof-chain blocker for GRAND-267 route-1.
