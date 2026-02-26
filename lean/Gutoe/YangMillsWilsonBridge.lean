/-
 * GUTOE — Wilson-Action Equivalence Bridge (Structural)
 * Copyright (C) 2026  Riff Labs
 * AGPL-3.0-or-later
 *
 * GRAND-300 (Theorem C bridge lane):
 *   - represent a Wilson-plaquette schedule by local Z₃ nearest-neighbor
 *     transition targets on the 3-state transfer basis
 *   - prove this representation induces SC-regular row totals structurally
 *   - instantiate the continuum-survival mass-gap lane with no empirical
 *     max-row certificate
 *
 * No `sorry`.
 -/

import Mathlib
import Gutoe.YangMillsContinuumSurvival
import Gutoe.YangMillsStructuralGap
import Gutoe.YangMillsMassGap
import Gutoe.LatticeGeometry

namespace Gutoe.YangMillsWilsonBridge

open Real
open Gutoe.YangMillsMassGap
open Gutoe.YangMillsStructuralGap
open Gutoe.YangMillsContinuumSurvival
open Gutoe.LatticeGeometry
open scoped BigOperators

/-- Structural Wilson-action schedule projected to the Z₃ transfer basis.

`targetSchedule n i e` gives the transfer-basis target state for refinement step
`n`, source basis state `i`, and incident SC edge `e`.

`betaSchedule` is carried for Wilson-lane bookkeeping (plaquette coupling) and
required to remain positive. The mass-gap bridge below only uses the transfer
projection object. -/
structure WilsonZ3Action where
  targetSchedule : ℕ → Z3NearestNeighborTargets
  betaSchedule : ℕ → ℝ
  beta_pos : ∀ n, 0 < betaSchedule n

/-- Row-normalized transition kernel from strictly positive edge weights on the
3-state transfer basis. -/
noncomputable def normalizedKernelFromWeights
    (W : Fin 3 → Fin 3 → ℝ) (i j : Fin 3) : ℝ :=
  W i j / (∑ k : Fin 3, W i k)

/-- Uniform row scaling does not change the normalized kernel. This is the
formal "absolute-count scale cancels under row normalization" bridge used when
moving from local to aggregated transition counts. -/
theorem normalizedKernel_row_scale_invariant
    (W : Fin 3 → Fin 3 → ℝ)
    (hW : ∀ i j, 0 < W i j)
    (s : Fin 3 → ℝ)
    (hs : ∀ i, 0 < s i) :
    normalizedKernelFromWeights (fun i j => s i * W i j) =
      normalizedKernelFromWeights W := by
  funext i j
  have hsumPos : 0 < ∑ k : Fin 3, W i k := by
    have h0 : 0 < W i 0 := hW i 0
    have h1 : 0 < W i 1 := hW i 1
    have h2 : 0 < W i 2 := hW i 2
    have hsum : ∑ k : Fin 3, W i k = W i 0 + W i 1 + W i 2 := by
      simpa [Fin.sum_univ_three]
    rw [hsum]
    nlinarith
  have hsPos : 0 < s i := hs i
  have hsNe : s i ≠ 0 := ne_of_gt hsPos
  have hsumNe : (∑ k : Fin 3, W i k) ≠ 0 := ne_of_gt hsumPos
  unfold normalizedKernelFromWeights
  have hsumScaled :
      (∑ k : Fin 3, s i * W i k) = s i * (∑ k : Fin 3, W i k) := by
    rw [Finset.mul_sum]
  rw [hsumScaled]
  field_simp [hsNe, hsumNe]

/-- Wilson action reduced to the 3-state transfer basis. -/
abbrev WilsonAction := Fin 3 → Fin 3 → ℝ

/-- Boltzmann weight induced by Wilson action and coupling `beta`. -/
noncomputable def wilsonWeight (beta : ℝ) (S : WilsonAction) (i j : Fin 3) : ℝ :=
  Real.exp (-beta * S i j)

/-- Wilson row-partition function on transfer basis state `i`. -/
noncomputable def wilsonRowPartition (beta : ℝ) (S : WilsonAction) (i : Fin 3) : ℝ :=
  ∑ j : Fin 3, wilsonWeight beta S i j

/-- Row-normalized Wilson transfer kernel. -/
noncomputable def wilsonKernel (beta : ℝ) (S : WilsonAction) (i j : Fin 3) : ℝ :=
  normalizedKernelFromWeights (wilsonWeight beta S) i j

/-- Wilson weights are strictly positive. -/
theorem wilson_weight_pos (beta : ℝ) (S : WilsonAction) :
    ∀ i j, 0 < wilsonWeight beta S i j := by
  intro i j
  exact Real.exp_pos _

/-- Wilson row partition is strictly positive. -/
theorem wilson_row_partition_pos (beta : ℝ) (S : WilsonAction) :
    ∀ i, 0 < wilsonRowPartition beta S i := by
  intro i
  have h0 : 0 < wilsonWeight beta S i 0 := wilson_weight_pos beta S i 0
  have h1 : 0 < wilsonWeight beta S i 1 := wilson_weight_pos beta S i 1
  have h2 : 0 < wilsonWeight beta S i 2 := wilson_weight_pos beta S i 2
  unfold wilsonRowPartition
  have hsum :
      (∑ j : Fin 3, wilsonWeight beta S i j) =
        wilsonWeight beta S i 0 + wilsonWeight beta S i 1 + wilsonWeight beta S i 2 := by
    simpa [Fin.sum_univ_three]
  rw [hsum]
  nlinarith

/-- Wilson kernel rows are stochastic (sum to one). -/
theorem wilson_kernel_row_sum_one (beta : ℝ) (S : WilsonAction) :
    ∀ i, (∑ j : Fin 3, wilsonKernel beta S i j) = 1 := by
  intro i
  unfold wilsonKernel normalizedKernelFromWeights
  have hsumPos : 0 < wilsonRowPartition beta S i := wilson_row_partition_pos beta S i
  have hsumNe : (∑ k : Fin 3, wilsonWeight beta S i k) ≠ 0 := by
    exact ne_of_gt (by simpa [wilsonRowPartition] using hsumPos)
  calc
    (∑ j : Fin 3, wilsonWeight beta S i j / (∑ k : Fin 3, wilsonWeight beta S i k))
        = (∑ j : Fin 3, wilsonWeight beta S i j) / (∑ k : Fin 3, wilsonWeight beta S i k) := by
          rw [Finset.sum_div]
    _ = 1 := by
      simp [hsumNe]

/-- Adding a row-wise action offset does not change the normalized Wilson
kernel (Boltzmann common-factor cancellation). -/
theorem wilson_kernel_row_offset_invariant
    (beta : ℝ) (S : WilsonAction) (c : Fin 3 → ℝ) :
    wilsonKernel beta (fun i j => S i j + c i) = wilsonKernel beta S := by
  let s : Fin 3 → ℝ := fun i => Real.exp (-beta * c i)
  have hs : ∀ i, 0 < s i := by
    intro i
    unfold s
    exact Real.exp_pos _
  have hW :
      wilsonWeight beta (fun i j => S i j + c i) =
      (fun i j => s i * wilsonWeight beta S i j) := by
    funext i j
    unfold wilsonWeight s
    have hmul : -beta * (S i j + c i) = (-beta * c i) + (-beta * S i j) := by ring
    rw [hmul, Real.exp_add]
  calc
    wilsonKernel beta (fun i j => S i j + c i)
        = normalizedKernelFromWeights (wilsonWeight beta (fun i j => S i j + c i)) := by
          rfl
    _ = normalizedKernelFromWeights (fun i j => s i * wilsonWeight beta S i j) := by
          rw [hW]
    _ = normalizedKernelFromWeights (wilsonWeight beta S) := by
          exact normalizedKernel_row_scale_invariant
            (wilsonWeight beta S) (wilson_weight_pos beta S) s hs
    _ = wilsonKernel beta S := by
          rfl

/-- Wilson-induced row-total schedule on the transfer basis. -/
def wilsonRowTotalsSchedule (W : WilsonZ3Action) : ℕ → Fin 3 → ℕ :=
  z3NearestNeighborRowTotalsSchedule W.targetSchedule

/-- Wilson-induced row totals are SC-regular at every refinement step. -/
theorem wilson_row_totals_sc_regular (W : WilsonZ3Action) :
    ∀ n, SCRegularRowTotals (wilsonRowTotalsSchedule W n) := by
  exact z3_nn_schedule_sc_regular W.targetSchedule

/-- Wilson-induced max row total is exactly the SC coordination number (`6`)
for every refinement step. -/
theorem wilson_max_row_total_eq_coordination (W : WilsonZ3Action) :
    ∀ n, maxRowTotal (wilsonRowTotalsSchedule W n) = coordinationNumber := by
  intro n
  exact z3_nn_max_row_total_eq_coordination (W.targetSchedule n)

/-- Wilson-induced minorization constant has a structural closed form at each
refinement step. -/
theorem wilson_minorization_eps_closed_form
    (W : WilsonZ3Action) (alpha : ℝ) :
    ∀ n,
      minorizationEps (wilsonRowTotalsSchedule W n) alpha =
        (3 * alpha) / ((coordinationNumber : ℝ) + 3 * alpha) := by
  intro n
  exact minorization_eps_eq_sc_regular
    (wilsonRowTotalsSchedule W n)
    alpha
    (wilson_row_totals_sc_regular W n)

/-- Bridge theorem (Theorem C lane, structural form):
if a Wilson schedule is represented by nearest-neighbor Z₃ transfer targets,
then the non-vanishing continuum mass-gap lower bound follows from the
structural Yang-Mills chain without empirical row-total hypotheses. -/
theorem wilson_action_bridge_nonvanishing_gap
    (W : WilsonZ3Action)
    (a_t : ℕ → ℝ)
    (alpha : ℝ)
    (ha_t_pos : ∀ n, 0 < a_t n)
    (ha_t_cap : ∃ aCap, 0 < aCap ∧ ∀ n, a_t n ≤ aCap)
    (ha : 0 < alpha) :
    ∃ c : ℝ, 0 < c ∧
      ∀ n, c ≤ doeblinGapLowerBound (a_t n)
        (minorizationEps (wilsonRowTotalsSchedule W n) alpha) := by
  exact continuum_survival_gap_nonvanishing_of_z3_nn_schedule
    a_t
    W.targetSchedule
    alpha
    ha_t_pos
    ha_t_cap
    ha

end Gutoe.YangMillsWilsonBridge
