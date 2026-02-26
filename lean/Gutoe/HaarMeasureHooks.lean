/-
 * GUTOE — Haar Measure Hooks (Path-2)
 * Copyright (C) 2026  Riff Labs
 * AGPL-3.0-or-later
 *
 * GRAND-309:
 *   1) Haar existence/uniqueness hooks on the SU(3)-lane group object.
 *   2) Quotient-measure hooks for `G ⧸ Z` (center-type normal subgroup).
 *   3) Unfolding/disintegration-style integral hook through a fundamental domain.
 *
 * No `sorry`.
 -/

import Mathlib
import Gutoe.HaarBridgeScaffold
import Gutoe.GaugeGroupSU3

noncomputable section

namespace Gutoe.HaarMeasureHooks

open MeasureTheory
open MeasureTheory.Measure
open QuotientGroup
open Gutoe.GaugeGroupSU3
open Gutoe.HaarBridgeScaffold

/-- Re-exported Cl(1,3)->Z3->SU(3) structural anchor for the Haar lane. -/
theorem clifford_z3_su3_anchor :
    Finset.card quarkOrbit = 3 ∧
    Nonempty ({s // s ∈ quarkOrbit} ≃ Fin 3) ∧
    Finset.card quarkOrbit ^ 2 - 1 = 8 :=
  Gutoe.HaarBridgeScaffold.clifford_z3_su3_anchor

section HaarExistenceUniqueness

variable {G : Type*}
  [Group G] [TopologicalSpace G] [MeasurableSpace G] [BorelSpace G]
  [IsTopologicalGroup G] [LocallyCompactSpace G] [T2Space G]

/-- Standard Haar measure constructor on a compact anchor set is Haar. -/
theorem haar_measure_is_haar (K0 : TopologicalSpace.PositiveCompacts G) :
    IsHaarMeasure (haarMeasure K0) := by
  infer_instance

/-- Canonical `Measure.haar` carries the Haar-measure typeclass. -/
theorem canonical_haar_is_haar :
    IsHaarMeasure (haar : Measure G) := by
  infer_instance

variable [SecondCountableTopology G]

/-- Uniqueness up to scalar normalization for left-invariant measures, using
Mathlib's Haar uniqueness theorem. -/
theorem left_invariant_measure_eq_smul_canonical_haar
    (μ : Measure G)
    [IsFiniteMeasureOnCompacts μ] [IsMulLeftInvariant μ] :
    μ = haarScalarFactor μ (haar : Measure G) • (haar : Measure G) := by
  simpa using (isMulLeftInvariant_eq_smul μ (haar : Measure G))

end HaarExistenceUniqueness

section QuotientHooks

variable {G : Type*}
  [Group G] [TopologicalSpace G] [MeasurableSpace G] [BorelSpace G]
  [IsTopologicalGroup G] [LocallyCompactSpace G] [PolishSpace G]

variable {Γ : Subgroup G} [Subgroup.Normal Γ] [Countable Γ]

variable [T2Space (G ⧸ Γ)] [SecondCountableTopology (G ⧸ Γ)]

/-- Quotient Haar hook: if quotient measure satisfies the standard
`QuotientMeasureEqMeasurePreimage` condition against a Haar base measure, then
it is itself Haar. -/
theorem quotient_haar_of_preimage
    (ν : Measure G) (μ : Measure (G ⧸ Γ))
    [IsHaarMeasure ν] [IsMulRightInvariant ν]
    [QuotientMeasureEqMeasurePreimage ν μ]
    [HasFundamentalDomain Γ.op G ν]
    [IsFiniteMeasure μ] :
    IsHaarMeasure μ := by
  exact MeasureTheory.QuotientMeasureEqMeasurePreimage.haarMeasure_quotient
    (ν := ν) (μ := μ)

/-- Quotient measure satisfying `QuotientMeasureEqMeasurePreimage` is unique for
fixed base measure `ν`. -/
theorem quotient_measure_unique_from_preimage
    (ν : Measure G)
    (μ₁ μ₂ : Measure (G ⧸ Γ))
    [HasFundamentalDomain Γ.op G ν]
    [QuotientMeasureEqMeasurePreimage ν μ₁]
    [QuotientMeasureEqMeasurePreimage ν μ₂] :
    μ₁ = μ₂ := by
  simpa using
    (MeasureTheory.QuotientMeasureEqMeasurePreimage.unique
      (ν := ν) (μ := μ₁) (μ' := μ₂))

/-- Pushforward of a restricted measure from a fundamental domain to the
quotient space. This is the finite-fiber side of the measure decomposition hook. -/
noncomputable def quotientFiberMeasure (μ : Measure G) (𝓕 : Set G) :
    Measure (G ⧸ Γ) :=
  Measure.map (QuotientGroup.mk : G → G ⧸ Γ) (μ.restrict 𝓕)

/-- Unfolding/disintegration-style bridge: integral on `G` equals integral of
`automorphize` on the quotient with respect to the fiber pushforward measure. -/
theorem integral_unfolding_over_quotient
    {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    {μ : Measure G} [μ.IsMulRightInvariant]
    {𝓕 : Set G}
    (h𝓕 : IsFundamentalDomain Γ.op 𝓕 μ)
    {f : G → E}
    (hf₁ : Integrable f μ)
    (hf₂ : AEStronglyMeasurable (QuotientGroup.automorphize f)
      (quotientFiberMeasure (Γ := Γ) μ 𝓕)) :
    (∫ x : G, f x ∂μ) =
      ∫ x : G ⧸ Γ, QuotientGroup.automorphize f x ∂(quotientFiberMeasure (Γ := Γ) μ 𝓕) := by
  simpa [quotientFiberMeasure] using
    (QuotientGroup.integral_eq_integral_automorphize (h𝓕 := h𝓕) hf₁ hf₂)

end QuotientHooks

end Gutoe.HaarMeasureHooks
