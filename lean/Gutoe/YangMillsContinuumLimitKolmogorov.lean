/-
 * GUTOE — Kolmogorov Extension: Constructive YM Infinite Path Measure via Ionescu-Tulcea
 * Copyright (C) 2026  Riff Labs
 *
 * AGPL-3.0-or-later
 *
 * GRAND-322 (push): Full Kolmogorov extension for the Wilson Schwinger family.
 *
 * Applies Mathlib's Ionescu-Tulcea theorem (`ProbabilityTheory.Kernel.trajMeasure`)
 * to produce an actual `Measure (ℕ → Fin 3)` — the infinite path measure —
 * as a Lean probability-measure object, not a proxy carrier.
 *
 * No `sorry`.
 -/

import Mathlib
import Gutoe.YangMillsWilsonBridge
import Gutoe.YangMillsWilsonEquivalence
import Gutoe.YangMillsStructuralGap
import Gutoe.YangMillsContinuumLimit
import Mathlib.Probability.Kernel.IonescuTulcea.Traj

noncomputable section

namespace Gutoe.YangMillsContinuumLimitKolmogorov

open Gutoe.YangMillsWilsonBridge
open Gutoe.YangMillsWilsonEquivalence
open Gutoe.YangMillsStructuralGap
open Gutoe.YangMillsMassGap
open Gutoe.YangMillsContinuumLimit
open ProbabilityTheory ProbabilityTheory.Kernel
open BigOperators MeasureTheory

-- ══════════════════════════════════════════════════════════════════════════════
-- Part A: Wilson kernel as a proper ProbabilityTheory.Kernel
-- ══════════════════════════════════════════════════════════════════════════════

/-- Row PMF induced by the Wilson transition weights at refinement step `n`. -/
noncomputable def wilsonRowPMF
    (W : WilsonZ3Action) (alpha : ℝ) (n : ℕ) (i : Fin 3) : PMF (Fin 3) :=
  PMF.ofFintype
    (fun j => ENNReal.ofReal (wilsonKernel 1 (centerPlaquetteActionSchedule W alpha n) i j))
    (by
      rw [← ENNReal.ofReal_sum_of_nonneg]
      · simpa using (wilson_kernel_row_sum_one 1 (centerPlaquetteActionSchedule W alpha n) i)
      · intro j _
        unfold wilsonKernel normalizedKernelFromWeights
        exact le_of_lt (div_pos
          (wilson_weight_pos 1 (centerPlaquetteActionSchedule W alpha n) i j)
          (wilson_row_partition_pos 1 (centerPlaquetteActionSchedule W alpha n) i)))

/-- Wilson transition entries are strictly positive. -/
theorem wilson_kernel_pos
    (W : WilsonZ3Action) (alpha : ℝ) (n : ℕ) (i j : Fin 3) :
    0 < wilsonKernel 1 (centerPlaquetteActionSchedule W alpha n) i j := by
  unfold wilsonKernel normalizedKernelFromWeights
  exact div_pos
    (wilson_weight_pos 1 (centerPlaquetteActionSchedule W alpha n) i j)
    (wilson_row_partition_pos 1 (centerPlaquetteActionSchedule W alpha n) i)

/-- Lift the Wilson matrix kernel to a Markov kernel on `Fin 3`. -/
noncomputable def wilsonMarkovKernel
    (W : WilsonZ3Action) (alpha : ℝ) (n : ℕ) :
    Kernel (Fin 3) (Fin 3) :=
  Kernel.ofFunOfCountable (fun i => (wilsonRowPMF W alpha n i).toMeasure)

/-- Wilson transition kernel is Markov at every refinement step. -/
instance wilsonMarkovKernel_isMarkov
    (W : WilsonZ3Action) (alpha : ℝ) (n : ℕ) :
    IsMarkovKernel (wilsonMarkovKernel W alpha n) where
  isProbabilityMeasure i := by
    change IsProbabilityMeasure ((wilsonRowPMF W alpha n i).toMeasure)
    infer_instance

-- ══════════════════════════════════════════════════════════════════════════════
-- Part B: Transition kernel on full history via Kernel.comap
-- ══════════════════════════════════════════════════════════════════════════════

/-- Homogeneous-step history kernel: next state depends on the last state in
`Finset.Iic n` history. -/
noncomputable def wilsonHistoryKernel
    (W : WilsonZ3Action) (alpha : ℝ) (n : ℕ) :
    Kernel ((i : Finset.Iic n) → Fin 3) (Fin 3) :=
  (wilsonMarkovKernel W alpha n).comap
    (fun p => p ⟨n, Finset.mem_Iic.mpr (le_rfl : n ≤ n)⟩)
    (measurable_pi_apply _)

/-- History kernel is Markov (inherited by `Kernel.comap`). -/
instance wilsonHistoryKernel_isMarkov
    (W : WilsonZ3Action) (alpha : ℝ) (n : ℕ) :
    IsMarkovKernel (wilsonHistoryKernel W alpha n) := by
  simpa [wilsonHistoryKernel] using
    (ProbabilityTheory.Kernel.IsMarkovKernel.comap
      (κ := wilsonMarkovKernel W alpha n)
      (g := fun p : (i : Finset.Iic n) → Fin 3 =>
        p ⟨n, Finset.mem_Iic.mpr (le_rfl : n ≤ n)⟩)
      (hg := measurable_pi_apply _))

-- ══════════════════════════════════════════════════════════════════════════════
-- Part C: Infinite path measure via Ionescu-Tulcea (`trajMeasure`)
-- ══════════════════════════════════════════════════════════════════════════════

/-- Infinite Wilson path measure starting from initial state `x₀`. -/
noncomputable def wilsonPathMeasure
    (W : WilsonZ3Action) (alpha : ℝ) (x₀ : Fin 3) :
    Measure (∀ n : ℕ, Fin 3) :=
  trajMeasure (X := fun _ => Fin 3) (Measure.dirac x₀) (wilsonHistoryKernel W alpha)

/-- Path measure is a probability measure (Ionescu-Tulcea + Markov kernels). -/
instance wilsonPathMeasure_isProbability
    (W : WilsonZ3Action) (alpha : ℝ) (x₀ : Fin 3) :
    IsProbabilityMeasure (wilsonPathMeasure W alpha x₀) := by
  unfold wilsonPathMeasure
  infer_instance

-- ══════════════════════════════════════════════════════════════════════════════
-- Part D: Path expectations
-- ══════════════════════════════════════════════════════════════════════════════

/-- Expectation of a real observable under the infinite Wilson path measure. -/
noncomputable def pathExpectation
    (W : WilsonZ3Action) (alpha : ℝ) (x₀ : Fin 3)
    (f : (ℕ → Fin 3) → ℝ) : ℝ :=
  ∫ p, f p ∂(wilsonPathMeasure W alpha x₀)

/-- Normalization: expectation of constant one is one. -/
theorem pathExpectation_one
    (W : WilsonZ3Action) (alpha : ℝ) (x₀ : Fin 3) :
    pathExpectation W alpha x₀ (fun _ => 1) = 1 := by
  unfold pathExpectation
  rw [MeasureTheory.integral_const]
  have hprob : (wilsonPathMeasure W alpha x₀) Set.univ = 1 := by
    simpa using (wilsonPathMeasure_isProbability W alpha x₀).measure_univ
  simp [hprob]

-- ══════════════════════════════════════════════════════════════════════════════
-- Part E: GRAND-322 complete constructive package
-- ══════════════════════════════════════════════════════════════════════════════

/-- The complete GRAND-322 constructive package:
    (A) Wilson kernel is Markov at every step
    (B) Infinite path measure is a probability measure (Ionescu-Tulcea)
    (C) Path expectation of constant one is normalized
    (D) Finite-step Schwinger correlators are normalized
    (E) Uniform mass gap floor `c > 0` from Wilson-equivalence domain
    No proxy carriers and no standalone existential interface assumptions. -/
theorem grand322_kolmogorov_extension_complete
    (W : WilsonZ3Action) (a_t : ℕ → ℝ) (alpha : ℝ) (x₀ : Fin 3)
    (dom : WilsonEquivalenceDomain a_t alpha) :
    (∀ n, IsMarkovKernel (wilsonMarkovKernel W alpha n)) ∧
    IsProbabilityMeasure (wilsonPathMeasure W alpha x₀) ∧
    pathExpectation W alpha x₀ (fun _ => 1) = 1 ∧
    (∀ n m, wilsonSchwingerFamily W alpha n m (fun _ => 1) = 1) ∧
    (∃ c : ℝ, 0 < c ∧ ∀ n, c ≤ doeblinGapLowerBound (a_t n)
        (minorizationEps (wilsonRowTotalsSchedule W n) alpha)) := by
  refine ⟨?_, ?_, ?_, ?_, ?_⟩
  · intro n
    infer_instance
  · infer_instance
  · exact pathExpectation_one W alpha x₀
  · exact wilson_schwinger_normalized W alpha
  · exact continuum_mass_gap_from_wilson_domain W a_t alpha dom

end Gutoe.YangMillsContinuumLimitKolmogorov

end
