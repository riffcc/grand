/-
 * GUTOE — Constructive Continuum QFT Interface (YM)
 * Copyright (C) 2026  Riff Labs
 * AGPL-3.0-or-later
 *
 * GRAND-302:
 *   Explicit Lean-facing interface for an OS/Wightman constructive lane and
 *   embedding of the mass-gap statement in the same framework.
 *
 * No `sorry`.
 -/

import Mathlib
import Gutoe.YangMillsContinuumSurvival
import Gutoe.YangMillsWilsonEquivalence
import Gutoe.YangMillsStructuralGap

noncomputable section

namespace Gutoe.YangMillsConstructiveQFT

open Gutoe.YangMillsContinuumSurvival
open Gutoe.YangMillsMassGap
open Gutoe.YangMillsStructuralGap
open Gutoe.YangMillsWilsonBridge
open Gutoe.YangMillsWilsonEquivalence

/-- Explicit OS-style axiom interface carried by the constructive lane. -/
structure OSAxiomInterface where
  reflectionPositivity : Prop
  euclideanInvariance : Prop
  regularity : Prop
  clusterProperty : Prop

/-- Explicit constructive milestones for the continuum YM existence lane. -/
structure ConstructiveMilestones where
  schwingerFunctionsExist : Prop
  osReconstruction : Prop
  wightmanCompatibility : Prop

/-- Lean-facing constructive object interface (non-Prop carriers) used by the
OS/Wightman lane. -/
structure ConstructiveFieldInterfaces where
  EuclideanField : Type
  MinkowskiField : Type
  SchwingerObject : Type
  WightmanObject : Type
  osReconstructMap : EuclideanField → MinkowskiField

/-- Full constructive model package used by the YM claim lane. -/
structure ConstructiveYMModel where
  os : OSAxiomInterface
  milestones : ConstructiveMilestones
  interfaces : ConstructiveFieldInterfaces

/-- Non-emptiness checks for the constructive carrier interfaces. -/
def constructiveInterfaceNonempty (M : ConstructiveYMModel) : Prop :=
  Nonempty M.interfaces.EuclideanField ∧
  Nonempty M.interfaces.MinkowskiField ∧
  Nonempty M.interfaces.SchwingerObject ∧
  Nonempty M.interfaces.WightmanObject

/-- The checklist demanded by GRAND-302 is explicit in one Lean proposition. -/
def constructiveTargetsSatisfied (M : ConstructiveYMModel) : Prop :=
  constructiveInterfaceNonempty M ∧
  M.os.reflectionPositivity ∧
  M.os.euclideanInvariance ∧
  M.os.regularity ∧
  M.os.clusterProperty ∧
  M.milestones.schwingerFunctionsExist ∧
  M.milestones.osReconstruction ∧
  M.milestones.wightmanCompatibility

/-- Unpack theorem: the explicit target list is recoverable from the model. -/
theorem constructive_targets_unpacked
    (M : ConstructiveYMModel)
    (h : constructiveTargetsSatisfied M) :
    constructiveInterfaceNonempty M ∧
    M.os.reflectionPositivity ∧
    M.os.euclideanInvariance ∧
    M.os.regularity ∧
    M.os.clusterProperty ∧
    M.milestones.schwingerFunctionsExist ∧
    M.milestones.osReconstruction ∧
    M.milestones.wightmanCompatibility := by
  exact h

/-- Mass-gap statement embedded inside the same constructive framework. -/
def massGapEmbeddedInConstructiveLane
    (M : ConstructiveYMModel)
    (a_t eps : ℕ → ℝ) : Prop :=
  constructiveTargetsSatisfied M ∧
  ContinuumSurvivalHypotheses a_t eps ∧
  ∃ c : ℝ, 0 < c ∧ ∀ n, c ≤ doeblinGapLowerBound (a_t n) (eps n)

/-- There are no hidden assumptions in the embedding proposition:
it is exactly the conjunction shown in the definition. -/
theorem no_hidden_assumptions_mass_gap_embedding
    (M : ConstructiveYMModel)
    (a_t eps : ℕ → ℝ) :
    massGapEmbeddedInConstructiveLane M a_t eps ↔
      constructiveTargetsSatisfied M ∧
      ContinuumSurvivalHypotheses a_t eps ∧
      (∃ c : ℝ, 0 < c ∧ ∀ n, c ≤ doeblinGapLowerBound (a_t n) (eps n)) := by
  rfl

/-- If constructive targets are satisfied and continuum-survival hypotheses hold,
then the mass-gap bound is embedded in the same formal lane. -/
theorem mass_gap_embedded_of_continuum_survival
    (M : ConstructiveYMModel)
    (a_t eps : ℕ → ℝ)
    (hTargets : constructiveTargetsSatisfied M)
    (hCont : ContinuumSurvivalHypotheses a_t eps) :
    massGapEmbeddedInConstructiveLane M a_t eps := by
  refine ⟨hTargets, hCont, ?_⟩
  exact continuum_survival_gap_nonvanishing a_t eps hCont

/-- Continuum-survival mass-gap lower bound can be extracted back out of the
constructive-lane package. -/
theorem embedded_mass_gap_extract
    (M : ConstructiveYMModel)
    (a_t eps : ℕ → ℝ)
    (h : massGapEmbeddedInConstructiveLane M a_t eps) :
    ∃ c : ℝ, 0 < c ∧ ∀ n, c ≤ doeblinGapLowerBound (a_t n) (eps n) := by
  exact h.2.2

/-- Targets and continuum hypotheses are each individually extractable from an
embedded mass-gap statement. -/
theorem embedded_targets_and_continuum_extract
    (M : ConstructiveYMModel)
    (a_t eps : ℕ → ℝ)
    (h : massGapEmbeddedInConstructiveLane M a_t eps) :
    constructiveTargetsSatisfied M ∧ ContinuumSurvivalHypotheses a_t eps := by
  exact ⟨h.1, h.2.1⟩

/-- If continuum-survival hypotheses fail, embedding fails. -/
theorem not_embedded_without_continuum
    (M : ConstructiveYMModel)
    (a_t eps : ℕ → ℝ)
    (hNoCont : ¬ ContinuumSurvivalHypotheses a_t eps) :
    ¬ massGapEmbeddedInConstructiveLane M a_t eps := by
  intro h
  exact hNoCont h.2.1

/-- If constructive targets fail, embedding fails. -/
theorem not_embedded_without_targets
    (M : ConstructiveYMModel)
    (a_t eps : ℕ → ℝ)
    (hNoTargets : ¬ constructiveTargetsSatisfied M) :
    ¬ massGapEmbeddedInConstructiveLane M a_t eps := by
  intro h
  exact hNoTargets h.1

/-- Embedded lane consequence: every refinement step has strictly positive
Doeblin lower-bound mass gap. -/
theorem embedded_gap_positive_each_step
    (M : ConstructiveYMModel)
    (a_t eps : ℕ → ℝ)
    (h : massGapEmbeddedInConstructiveLane M a_t eps) :
    ∀ n, 0 < doeblinGapLowerBound (a_t n) (eps n) := by
  rcases embedded_mass_gap_extract M a_t eps h with ⟨c, hcPos, hcLe⟩
  intro n
  exact lt_of_lt_of_le hcPos (hcLe n)

/-- Sequential notion used to state "vanishing gap" failures explicitly in the
constructive lane. -/
def TendsToZeroSeq (g : ℕ → ℝ) : Prop :=
  ∀ ε > 0, ∃ N, ∀ n ≥ N, |g n| < ε

/-- A sequence with a uniform strictly positive floor cannot tend to zero. -/
theorem not_tends_to_zero_of_uniform_positive_floor
    {g : ℕ → ℝ}
    (hFloor : ∃ c : ℝ, 0 < c ∧ ∀ n, c ≤ g n) :
    ¬ TendsToZeroSeq g := by
  intro hZero
  rcases hFloor with ⟨c, hcPos, hcLe⟩
  have hHalfPos : 0 < c / 2 := by nlinarith
  rcases hZero (c / 2) hHalfPos with ⟨N, hN⟩
  have hAbs : |g N| < c / 2 := hN N (le_rfl)
  have hgNonneg : 0 ≤ g N := le_trans (le_of_lt hcPos) (hcLe N)
  have hLt : g N < c / 2 := by simpa [abs_of_nonneg hgNonneg] using hAbs
  have hGe : c ≤ g N := hcLe N
  linarith

/-- Embedded constructive lane forbids a vanishing-gap sequence. -/
theorem embedded_gap_not_tends_to_zero
    (M : ConstructiveYMModel)
    (a_t eps : ℕ → ℝ)
    (h : massGapEmbeddedInConstructiveLane M a_t eps) :
    ¬ TendsToZeroSeq (fun n => doeblinGapLowerBound (a_t n) (eps n)) := by
  exact not_tends_to_zero_of_uniform_positive_floor (embedded_mass_gap_extract M a_t eps h)

/-- Combined closure statement for GRAND-302:
the constructive embedding yields both per-step positivity and non-vanishing
sequence behavior. -/
theorem constructive_lane_gap_closure
    (M : ConstructiveYMModel)
    (a_t eps : ℕ → ℝ)
    (h : massGapEmbeddedInConstructiveLane M a_t eps) :
    (∀ n, 0 < doeblinGapLowerBound (a_t n) (eps n)) ∧
    (¬ TendsToZeroSeq (fun n => doeblinGapLowerBound (a_t n) (eps n))) := by
  exact ⟨embedded_gap_positive_each_step M a_t eps h, embedded_gap_not_tends_to_zero M a_t eps h⟩

/-- Canonical constructive embedding on the Wilson-equivalence lane:
constructive targets + Wilson-equivalence domain imply a mass-gap embedding in
the same constructive framework, with `eps` induced by the Wilson row schedule. -/
theorem mass_gap_embedded_of_wilson_equivalence_domain
    (M : ConstructiveYMModel)
    (hTargets : constructiveTargetsSatisfied M)
    (W : WilsonZ3Action)
    (a_t : ℕ → ℝ)
    (alpha : ℝ)
    (hDom : WilsonEquivalenceDomain a_t alpha) :
    massGapEmbeddedInConstructiveLane M a_t
      (fun n => minorizationEps (wilsonRowTotalsSchedule W n) alpha) := by
  refine ⟨hTargets, ?_, ?_⟩
  · -- Continuum-survival hypotheses, specialized to the Wilson-induced row schedule.
    exact continuum_hypotheses_of_z3_nn_schedule
      a_t
      W.targetSchedule
      alpha
      hDom.a_t_pos
      hDom.a_t_cap
      hDom.alpha_pos
  · -- Uniform non-vanishing lower bound transferred through Theorem-C domain.
    rcases c3_gap_correspondence_of_domain W a_t alpha hDom with ⟨c, hcPos, hcLe⟩
    refine ⟨c, hcPos, ?_⟩
    intro n
    exact hcLe n

/-- Wilson-equivalence specialization of the constructive-lane closure theorem. -/
theorem constructive_lane_gap_closure_of_wilson_equivalence_domain
    (M : ConstructiveYMModel)
    (hTargets : constructiveTargetsSatisfied M)
    (W : WilsonZ3Action)
    (a_t : ℕ → ℝ)
    (alpha : ℝ)
    (hDom : WilsonEquivalenceDomain a_t alpha) :
    (∀ n, 0 < doeblinGapLowerBound
      (a_t n)
      (minorizationEps (wilsonRowTotalsSchedule W n) alpha)) ∧
    (¬ TendsToZeroSeq (fun n =>
      doeblinGapLowerBound
        (a_t n)
        (minorizationEps (wilsonRowTotalsSchedule W n) alpha))) := by
  exact constructive_lane_gap_closure
    M
    a_t
    (fun n => minorizationEps (wilsonRowTotalsSchedule W n) alpha)
    (mass_gap_embedded_of_wilson_equivalence_domain
      M hTargets W a_t alpha hDom)

end Gutoe.YangMillsConstructiveQFT
