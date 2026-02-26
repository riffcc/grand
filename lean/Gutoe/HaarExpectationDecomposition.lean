/-
 * GUTOE — Haar Expectation Decomposition (Path-2)
 * Copyright (C) 2026  Riff Labs
 * AGPL-3.0-or-later
 *
 * GRAND-310:
 *   - Prove expectation decomposition over subgroup/center fibers:
 *       E_G[f] = E_{G⧸Γ}[E_fiber[f]]
 *   - Provide finite transfer-lane parity theorem matching existing Wilson
 *     bridge collapse theorem.
 *
 * No `sorry`.
 -/

import Mathlib
import Gutoe.HaarMeasureHooks
import Gutoe.YangMillsWilsonBridge

noncomputable section

namespace Gutoe.HaarExpectationDecomposition

open MeasureTheory
open QuotientGroup
open Gutoe.HaarBridgeScaffold
open Gutoe.HaarMeasureHooks
open Gutoe.YangMillsWilsonBridge

/-- Expectation functional used in the Haar bridge lane. -/
def expectation {α : Type*} [MeasurableSpace α] (μ : Measure α) (f : α → ℝ) : ℝ :=
  ∫ x, f x ∂μ

/-- Fiber observable induced on the quotient by subgroup `Γ`. -/
def fiberExpectation
    {G : Type*} [Group G] {Γ : Subgroup G}
    (f : G → ℝ) : G ⧸ Γ → ℝ :=
  QuotientGroup.automorphize f

section Continuous

variable {G : Type*}
  [Group G] [TopologicalSpace G] [MeasurableSpace G] [BorelSpace G]
  [IsTopologicalGroup G] [LocallyCompactSpace G] [PolishSpace G]

variable {Γ : Subgroup G} [Subgroup.Normal Γ] [Countable Γ]

variable [T2Space (G ⧸ Γ)] [SecondCountableTopology (G ⧸ Γ)]

/-- Core decomposition theorem (GRAND-310): integral/expectation over `G`
factors through quotient expectation of the fiber observable. -/
theorem expectation_decomposition_over_subgroup
    {μ : Measure G} [μ.IsMulRightInvariant]
    {𝓕 : Set G}
    (h𝓕 : IsFundamentalDomain Γ.op 𝓕 μ)
    {f : G → ℝ}
    (hf₁ : Integrable f μ)
    (hf₂ : AEStronglyMeasurable (fiberExpectation (Γ := Γ) f)
      (quotientFiberMeasure (Γ := Γ) μ 𝓕)) :
    expectation μ f =
      expectation (quotientFiberMeasure (Γ := Γ) μ 𝓕)
        (fiberExpectation (Γ := Γ) f) := by
  simpa [expectation, fiberExpectation] using
    (integral_unfolding_over_quotient (Γ := Γ) (h𝓕 := h𝓕) hf₁ hf₂)

/-- Center-specialized decomposition theorem:
`E_G[f] = E_{G⧸Z(G)}[E_fiber[f]]`. -/
theorem expectation_decomposition_over_center
    [Countable ↥(centerSubgroup (G := G))]
    [T2Space (G ⧸ centerSubgroup (G := G))]
    [SecondCountableTopology (G ⧸ centerSubgroup (G := G))]
    {μ : Measure G} [μ.IsMulRightInvariant]
    {𝓕 : Set G}
    (h𝓕 : IsFundamentalDomain (centerSubgroup (G := G)).op 𝓕 μ)
    {f : G → ℝ}
    (hf₁ : Integrable f μ)
    (hf₂ : AEStronglyMeasurable
      (fiberExpectation (Γ := centerSubgroup (G := G)) f)
      (quotientFiberMeasure (Γ := centerSubgroup (G := G)) μ 𝓕)) :
    expectation μ f =
      expectation
        (quotientFiberMeasure (Γ := centerSubgroup (G := G)) μ 𝓕)
        (fiberExpectation (Γ := centerSubgroup (G := G)) f) := by
  exact expectation_decomposition_over_subgroup
    (Γ := centerSubgroup (G := G)) h𝓕 hf₁ hf₂

end Continuous

section FiniteParity

/-- Finite transfer-lane left side (center × fiber averaging form). -/
def finiteFiberExpectationLHS
    {C F : Type}
    [Fintype C] [Fintype F]
    (lift : C → F → WilsonAction)
    (wC : C → ℝ)
    (wF : F → ℝ)
    (ObsCenter : Matrix (Fin 3) (Fin 3) ℝ → ℝ) : ℝ :=
  ∑ c : C, ∑ f : F,
    (wC c) * (wF f) * ObsCenter (wilsonKernel 1 (lift c f))

/-- Finite transfer-lane right side (collapsed base expectation). -/
def finiteFiberExpectationRHS
    {C F : Type}
    [Fintype C] [Fintype F]
    (lift : C → F → WilsonAction)
    (f₀ : F)
    (wC : C → ℝ)
    (ObsCenter : Matrix (Fin 3) (Fin 3) ℝ → ℝ) : ℝ :=
  ∑ c : C, (wC c) * ObsCenter (wilsonKernel 1 (lift c f₀))

/-- Finite parity theorem (GRAND-310 acceptance): this module's finite
expectation form matches the established transfer-lane collapse theorem exactly. -/
theorem finite_parity_with_transfer_lane
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
  simpa [finiteFiberExpectationLHS, finiteFiberExpectationRHS] using
    (finite_fiber_expectation_collapse lift f₀ wC wF hwF ObsCenter hscale)

end FiniteParity

end Gutoe.HaarExpectationDecomposition
