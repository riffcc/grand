/- 
 * GUTOE — Wilson-Action Equivalence Spine (Theorem C lane)
 * Copyright (C) 2026  Riff Labs
 * AGPL-3.0-or-later
 *
 * GRAND-300:
 *   - precise equivalence domain/limits statement,
 *   - action correspondence theorem,
 *   - measure correspondence theorem (via Haar bridge + GRAND-312),
 *   - correlator correspondence theorem on the finite transfer lane,
 *   - explicit continuum-gap transfer consequence.
 *
 * No `sorry`.
 -/

import Mathlib
import Gutoe.YangMillsWilsonBridge
import Gutoe.HaarFiberCollapse

noncomputable section

namespace Gutoe.YangMillsWilsonEquivalence

open MeasureTheory
open Gutoe.YangMillsWilsonBridge
open Gutoe.YangMillsMassGap
open Gutoe.YangMillsStructuralGap
open Gutoe.HaarBridgeScaffold
open Gutoe.HaarExpectationDecomposition
open Gutoe.HaarFiberCollapse

/-- Explicit Theorem-C domain and limit package:
positive time-step schedule, bounded refinement cap, and positive Laplace floor. -/
structure WilsonEquivalenceDomain (a_t : ℕ → ℝ) (alpha : ℝ) : Prop where
  a_t_pos : ∀ n, 0 < a_t n
  a_t_cap : ∃ aCap, 0 < aCap ∧ ∀ n, a_t n ≤ aCap
  alpha_pos : 0 < alpha

/-- Action correspondence (Theorem C, action layer):
center-plaquette Wilson kernel equals the Z₃ transfer kernel at each refinement step. -/
theorem action_correspondence_of_domain
    (W : WilsonZ3Action)
    (a_t : ℕ → ℝ)
    (alpha : ℝ)
    (hdom : WilsonEquivalenceDomain a_t alpha) :
    ∀ n,
      wilsonKernel 1 (centerPlaquetteActionSchedule W alpha n) =
        smoothedTransition
          (z3NearestNeighborCounts (W.targetSchedule n))
          (rowTotalsFromCounts (z3NearestNeighborCounts (W.targetSchedule n)))
          alpha := by
  exact center_plaquette_schedule_kernel_eq_transfer W hdom.alpha_pos

/-- Gap correspondence (Theorem C, limit layer):
the Wilson lane inherits a non-vanishing continuum Doeblin mass-gap lower bound. -/
theorem c3_gap_correspondence_of_domain
    (W : WilsonZ3Action)
    (a_t : ℕ → ℝ)
    (alpha : ℝ)
    (hdom : WilsonEquivalenceDomain a_t alpha) :
    ∃ c : ℝ, 0 < c ∧
      ∀ n, c ≤ doeblinGapLowerBound (a_t n)
        (minorizationEps (wilsonRowTotalsSchedule W n) alpha) := by
  exact c3_wilson_gap_nonvanishing_from_clifford_z3
    W a_t alpha hdom.a_t_pos hdom.a_t_cap hdom.alpha_pos

/-- Two-point finite transfer-lane correlator observable from two kernel observables. -/
def twoPointKernelObservable
    (ObsA ObsB : Matrix (Fin 3) (Fin 3) ℝ → ℝ) :
    Matrix (Fin 3) (Fin 3) ℝ → ℝ :=
  fun K => ObsA K * ObsB K

/-- Correlator correspondence (finite transfer lane):
fiber averaging over the full lane collapses to center/base averaging when the
row-scale orbit hypothesis holds. -/
theorem correlator_correspondence_finite_of_row_scale_orbit
    {C F : Type}
    [Fintype C] [Fintype F]
    (lift : C → F → WilsonAction)
    (f₀ : F)
    (wC : C → ℝ)
    (wF : F → ℝ)
    (hwF : ∑ f : F, wF f = 1)
    (ObsA ObsB : Matrix (Fin 3) (Fin 3) ℝ → ℝ)
    (hscale :
      ∀ c f,
        RowScaleEquivalent
          (wilsonWeight 1 (lift c f₀))
          (wilsonWeight 1 (lift c f))) :
    (∑ c : C, ∑ f : F,
      (wC c) * (wF f) *
        (twoPointKernelObservable ObsA ObsB) (wilsonKernel 1 (lift c f))) =
      ∑ c : C, (wC c) *
        (twoPointKernelObservable ObsA ObsB) (wilsonKernel 1 (lift c f₀)) := by
  simpa [twoPointKernelObservable] using
    (finite_fiber_expectation_collapse
      lift
      f₀
      wC
      wF
      hwF
      (twoPointKernelObservable ObsA ObsB)
      hscale)

section MeasureCorrespondence

variable {G : Type*}
  [Group G] [TopologicalSpace G] [MeasurableSpace G] [BorelSpace G]
  [IsTopologicalGroup G] [LocallyCompactSpace G] [PolishSpace G]
  [Countable ↥(centerSubgroup (G := G))]
  [T2Space (G ⧸ centerSubgroup (G := G))]
  [SecondCountableTopology (G ⧸ centerSubgroup (G := G))]

/-- Measure correspondence (Theorem C, measure layer):
normalized full-lane expectation equals normalized center-quotient expectation
once quotient-normalization data is supplied. -/
theorem measure_correspondence_of_center_quotient_normalization
    (μG : Measure G) [μG.IsMulRightInvariant] [IsFiniteMeasure μG]
    (𝓕 : Set G)
    (h𝓕 : IsFundamentalDomain (centerSubgroup (G := G)).op 𝓕 μG)
    (f : G → ℝ)
    (fQ : G ⧸ centerSubgroup (G := G) → ℝ)
    (c : ℝ)
    (hfInt : Integrable f μG)
    (hfAE :
      AEStronglyMeasurable
        (fiberExpectation f : G ⧸ centerSubgroup (G := G) → ℝ)
        (Gutoe.HaarMeasureHooks.quotientFiberMeasure μG 𝓕 : Measure (G ⧸ centerSubgroup (G := G))))
    (hOneAE :
      AEStronglyMeasurable
        (fiberExpectation (fun _ : G => (1 : ℝ)) : G ⧸ centerSubgroup (G := G) → ℝ)
        (Gutoe.HaarMeasureHooks.quotientFiberMeasure μG 𝓕 : Measure (G ⧸ centerSubgroup (G := G))))
    (hFiberObs :
      (fiberExpectation f : G ⧸ centerSubgroup (G := G) → ℝ) = fun q => c * fQ q)
    (hFiberMass :
      (fiberExpectation (fun _ : G => (1 : ℝ)) : G ⧸ centerSubgroup (G := G) → ℝ) = fun _ => c)
    (hc : c ≠ 0)
    (hMassQ :
      ((Gutoe.HaarMeasureHooks.quotientFiberMeasure μG 𝓕 : Measure (G ⧸ centerSubgroup (G := G))) Set.univ).toReal ≠ 0) :
    normalizedExpectation μG f =
      normalizedExpectation
        (Gutoe.HaarMeasureHooks.quotientFiberMeasure μG 𝓕 : Measure (G ⧸ centerSubgroup (G := G))) fQ := by
  exact normalized_expectation_reduce_to_center_of_quotient_normalization
    μG 𝓕 h𝓕 f fQ c hfInt hfAE hOneAE hFiberObs hFiberMass hc hMassQ

end MeasureCorrespondence

/-- Consolidated Theorem-C spine:
under the explicit equivalence domain, action correspondence and continuum-gap
transfer both hold for the Wilson center lane. -/
theorem theorem_c_wilson_equivalence_domain_limits
    (W : WilsonZ3Action)
    (a_t : ℕ → ℝ)
    (alpha : ℝ)
    (hdom : WilsonEquivalenceDomain a_t alpha) :
    (∀ n,
      wilsonKernel 1 (centerPlaquetteActionSchedule W alpha n) =
        smoothedTransition
          (z3NearestNeighborCounts (W.targetSchedule n))
          (rowTotalsFromCounts (z3NearestNeighborCounts (W.targetSchedule n)))
          alpha) ∧
    (∃ c : ℝ, 0 < c ∧
      ∀ n, c ≤ doeblinGapLowerBound (a_t n)
        (minorizationEps (wilsonRowTotalsSchedule W n) alpha)) := by
  refine ⟨action_correspondence_of_domain W a_t alpha hdom, ?_⟩
  exact c3_gap_correspondence_of_domain W a_t alpha hdom

end Gutoe.YangMillsWilsonEquivalence

