/-
 * GUTOE — Yang-Mills Gauge Group Scope Layer
 * Copyright (C) 2026  Riff Labs
 * AGPL-3.0-or-later
 *
 * GRAND-303:
 *   Generalize Path-2 bridge statements from SU(3)-specific naming to abstract
 *   compact-group scope with finite center.
 *
 * No `sorry`.
 -/

import Mathlib
import Gutoe.HaarExpectationDecomposition
import Gutoe.HaarBridgeScaffold
import Gutoe.YangMillsContinuumSurvival

noncomputable section

namespace Gutoe.YangMillsGaugeScope

open Gutoe.HaarBridgeScaffold
open Gutoe.HaarExpectationDecomposition
open Gutoe.YangMillsContinuumSurvival
open Gutoe.YangMillsMassGap

/-- Scope proposition for the Clay-lane generalization argument:
compact-simple group with finite center. -/
def CompactSimpleGaugeScope (G : Type*) [Group G] : Prop :=
  IsSimpleGroup G ∧ Finite (Subgroup.center G)

section GeneralScope

variable {G : Type*}
  [Group G] [TopologicalSpace G] [MeasurableSpace G] [BorelSpace G]
  [IsTopologicalGroup G] [LocallyCompactSpace G] [PolishSpace G] [T2Space G]

/-- Finite center gives countable center (needed for quotient summation/measurable
interfaces in the Path-2 lane). -/
theorem finite_center_implies_countable_center
    (hFin : Finite (Subgroup.center G)) :
    Countable (Subgroup.center G) :=
  hFin.to_countable

/-- Compact groups carry Haar measure through the canonical constructor. -/
theorem compact_group_has_haar [CompactSpace G] :
    MeasureTheory.Measure.IsHaarMeasure (MeasureTheory.Measure.haar : MeasureTheory.Measure G) := by
  infer_instance

/-- Group-agnostic center decomposition theorem:
for any group `G` with finite center and standard quotient regularity hypotheses,
Path-2 expectation decomposition over `G ⧸ Z(G)` holds. -/
theorem expectation_decomposition_over_center_of_finite_center
    (hFin : Finite (Subgroup.center G))
    [SecondCountableTopology G]
    [T2Space (G ⧸ centerSubgroup (G := G))]
    [SecondCountableTopology (G ⧸ centerSubgroup (G := G))]
    {μ : MeasureTheory.Measure G} [μ.IsMulRightInvariant]
    {𝓕 : Set G}
    (h𝓕 : MeasureTheory.IsFundamentalDomain (centerSubgroup (G := G)).op 𝓕 μ)
    {f : G → ℝ}
    (hf₁ : MeasureTheory.Integrable f μ)
    (hf₂ : MeasureTheory.AEStronglyMeasurable
      (Gutoe.HaarExpectationDecomposition.fiberExpectation
        (Γ := centerSubgroup (G := G)) f)
      (Gutoe.HaarMeasureHooks.quotientFiberMeasure
        (Γ := centerSubgroup (G := G)) μ 𝓕)) :
    expectation μ f =
      expectation
        (Gutoe.HaarMeasureHooks.quotientFiberMeasure
          (Γ := centerSubgroup (G := G)) μ 𝓕)
        (Gutoe.HaarExpectationDecomposition.fiberExpectation
          (Γ := centerSubgroup (G := G)) f) := by
  haveI : Countable (Subgroup.center G) := finite_center_implies_countable_center hFin
  exact expectation_decomposition_over_center (G := G) (h𝓕 := h𝓕) hf₁ hf₂

/-- Continuum-survival mass-gap lane is group-agnostic; it can be stated in the
same compact-group scope without introducing SU(3)-specific assumptions. -/
theorem compact_scope_continuum_gap_nonvanishing
    (a_t eps : ℕ → ℝ)
    (hCont : ContinuumSurvivalHypotheses a_t eps) :
    ∃ c : ℝ, 0 < c ∧ ∀ n, c ≤ doeblinGapLowerBound (a_t n) (eps n) :=
  continuum_survival_gap_nonvanishing a_t eps hCont

/-- Scope theorem package: under compact-simple + finite-center scope, the Path-2
center decomposition interface and continuum gap interface both apply. -/
theorem compact_simple_scope_supports_path2
    (hScope : CompactSimpleGaugeScope G)
    [SecondCountableTopology G]
    [T2Space (G ⧸ centerSubgroup (G := G))]
    [SecondCountableTopology (G ⧸ centerSubgroup (G := G))]
    {μ : MeasureTheory.Measure G} [μ.IsMulRightInvariant]
    {𝓕 : Set G}
    (h𝓕 : MeasureTheory.IsFundamentalDomain (centerSubgroup (G := G)).op 𝓕 μ)
    {f : G → ℝ}
    (hf₁ : MeasureTheory.Integrable f μ)
    (hf₂ : MeasureTheory.AEStronglyMeasurable
      (Gutoe.HaarExpectationDecomposition.fiberExpectation
        (Γ := centerSubgroup (G := G)) f)
      (Gutoe.HaarMeasureHooks.quotientFiberMeasure
        (Γ := centerSubgroup (G := G)) μ 𝓕))
    (a_t eps : ℕ → ℝ)
    (hCont : ContinuumSurvivalHypotheses a_t eps) :
    (expectation μ f =
      expectation
        (Gutoe.HaarMeasureHooks.quotientFiberMeasure
          (Γ := centerSubgroup (G := G)) μ 𝓕)
        (Gutoe.HaarExpectationDecomposition.fiberExpectation
          (Γ := centerSubgroup (G := G)) f)) ∧
    (∃ c : ℝ, 0 < c ∧ ∀ n, c ≤ doeblinGapLowerBound (a_t n) (eps n)) := by
  refine ⟨?_, ?_⟩
  · exact expectation_decomposition_over_center_of_finite_center
      (G := G) hScope.2 h𝓕 hf₁ hf₂
  · exact compact_scope_continuum_gap_nonvanishing a_t eps hCont

end GeneralScope

end Gutoe.YangMillsGaugeScope
