/-
 * GUTOE — Haar Fiber Collapse (Path-2)
 * Copyright (C) 2026  Riff Labs
 * AGPL-3.0-or-later
 *
 * GRAND-311:
 *   - Gauge-invariant fiber constancy on subgroup/coset fibers.
 *   - Normalized-expectation collapse when coset factor is scalar.
 *   - Finite parity bridge to transfer-lane collapse theorem.
 *
 * No `sorry`.
 -/

import Mathlib
import Gutoe.HaarExpectationDecomposition

noncomputable section

namespace Gutoe.HaarFiberCollapse

open MeasureTheory
open QuotientGroup
open Gutoe.HaarBridgeScaffold
open Gutoe.HaarExpectationDecomposition
open Gutoe.YangMillsWilsonBridge

/-- Gauge invariance of `obs` under a group action encoded by homomorphism `ρ`.
This models "coset rotations are gauge transformations" in Path-2 form. -/
def GaugeInvariantUnder
    {H G β : Type*} [Group H] [Group G]
    (ρ : H →* G) (obs : G → β) : Prop :=
  ∀ h g, obs (g * ρ h) = obs g

/-- Observable factors through quotient projection. -/
def FactorsThroughQuotient
    {G β : Type*} [Group G] {Γ : Subgroup G}
    (obs : G → β) : Prop :=
  ∃ obsQ : G ⧸ Γ → β, obs = obsQ ∘ (QuotientGroup.mk : G → G ⧸ Γ)

/-- Fiber-constancy along subgroup right-cosets. -/
def FiberConstant
    {G β : Type*} [Group G] {Γ : Subgroup G}
    (obs : G → β) : Prop :=
  ∀ g z, z ∈ Γ → obs (g * z) = obs g

/-- Any quotient-factorized observable is fiber-constant on subgroup cosets. -/
theorem fiber_constant_of_factors_through_quotient
    {G β : Type*} [Group G] {Γ : Subgroup G}
    (obs : G → β)
    (hfact : FactorsThroughQuotient (Γ := Γ) obs) :
    FiberConstant (Γ := Γ) obs := by
  rcases hfact with ⟨obsQ, hobs⟩
  intro g z hz
  rcases hobs with rfl
  change obsQ ((QuotientGroup.mk : G → G ⧸ Γ) (g * z)) =
    obsQ ((QuotientGroup.mk : G → G ⧸ Γ) g)
  simpa using congrArg obsQ (QuotientGroup.mk_mul_of_mem g hz)

/-- Center-specialized fiber-constancy corollary. -/
theorem center_fiber_constant_of_factorization
    {G β : Type*} [Group G]
    (obs : G → β)
    (hfact : FactorsThroughQuotient (Γ := centerSubgroup (G := G)) obs) :
    FiberConstant (Γ := centerSubgroup (G := G)) obs :=
  fiber_constant_of_factors_through_quotient (Γ := centerSubgroup (G := G)) obs hfact

/-- If every subgroup fiber move is realized by a gauge transformation, then
gauge invariance implies fiber-constancy on that subgroup's cosets. -/
theorem gauge_invariant_implies_fiber_constant_of_surjective_action
    {H G β : Type*} [Group H] [Group G] {Γ : Subgroup G}
    (ρ : H →* G)
    (obs : G → β)
    (hCover : ∀ z, z ∈ Γ → ∃ h : H, ρ h = z)
    (hInv : GaugeInvariantUnder ρ obs) :
    FiberConstant (Γ := Γ) obs := by
  intro g z hz
  rcases hCover z hz with ⟨h, rfl⟩
  simpa [GaugeInvariantUnder] using hInv h g

/-- Normalized expectation helper: ratio of integral and total mass. -/
def normalizedExpectation
    {α : Type*} [MeasurableSpace α]
    (μ : Measure α) (f : α → ℝ) : ℝ :=
  (∫ x, f x ∂μ) / (μ Set.univ).toReal

/-- Scalar-factor cancellation in normalized expectations. -/
theorem normalized_expectation_scale_cancel
    (c I M : ℝ)
    (hc : c ≠ 0)
    (hM : M ≠ 0) :
    (c * I) / (c * M) = I / M := by
  field_simp [hc, hM]

/-- Coset-integral triviality in normalized form:
if both integral and total mass pick up the same scalar fiber factor `c`, that
factor cancels in normalized expectations. -/
theorem normalized_expectation_collapse_of_common_factor
    {G Q : Type*}
    [MeasurableSpace G] [MeasurableSpace Q]
    (μG : Measure G) (μQ : Measure Q)
    (fG : G → ℝ) (fQ : Q → ℝ)
    (c : ℝ)
    (hInt : (∫ x, fG x ∂μG) = c * (∫ y, fQ y ∂μQ))
    (hMass : (μG Set.univ).toReal = c * (μQ Set.univ).toReal)
    (hc : c ≠ 0)
    (hMassQ : (μQ Set.univ).toReal ≠ 0) :
    normalizedExpectation μG fG = normalizedExpectation μQ fQ := by
  unfold normalizedExpectation
  rw [hInt, hMass]
  simpa using normalized_expectation_scale_cancel c (∫ y, fQ y ∂μQ) ((μQ Set.univ).toReal) hc hMassQ

/-- Center-sector normalized expectation reduction:
if both the observable integral and total mass pick up the same coset volume
factor `c`, normalized expectations coincide. -/
theorem normalized_expectation_reduce_to_center
    {G Q : Type*}
    [MeasurableSpace G] [MeasurableSpace Q]
    (μG : Measure G) (μQ : Measure Q)
    (fG : G → ℝ) (fQ : Q → ℝ)
    (c : ℝ)
    (hInt : expectation μG fG = c * expectation μQ fQ)
    (hMass : (μG Set.univ).toReal = c * (μQ Set.univ).toReal)
    (hc : c ≠ 0)
    (hMassQ : (μQ Set.univ).toReal ≠ 0) :
    normalizedExpectation μG fG = normalizedExpectation μQ fQ := by
  exact normalized_expectation_collapse_of_common_factor
    μG μQ fG fQ c hInt hMass hc hMassQ

/-- Parity bridge theorem: the finite collapse statement used in GRAND-310 is
identical to the transfer-lane collapse theorem. -/
theorem finite_parity_bridge
    {C F : Type}
    [Fintype C] [Fintype F]
    (lift : C → F → WilsonAction)
    (f₀ : F)
    (wC : C → ℝ)
    (wF : F → ℝ)
    (hwF : ∑ f : F, wF f = 1)
    (ObsCenter : Matrix (Fin 3) (Fin 3) ℝ → ℝ)
    (hscale :
      ∀ c f,
        RowScaleEquivalent
          (wilsonWeight 1 (lift c f₀))
          (wilsonWeight 1 (lift c f))) :
    finiteFiberExpectationLHS lift wC wF ObsCenter =
      finiteFiberExpectationRHS lift f₀ wC ObsCenter := by
  exact finite_parity_with_transfer_lane lift f₀ wC wF hwF ObsCenter hscale

end Gutoe.HaarFiberCollapse
