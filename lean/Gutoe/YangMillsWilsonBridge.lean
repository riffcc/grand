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
import Gutoe.GaugeGroupSU3

namespace Gutoe.YangMillsWilsonBridge

open Real
open Gutoe.YangMillsMassGap
open Gutoe.YangMillsStructuralGap
open Gutoe.YangMillsContinuumSurvival
open Gutoe.LatticeGeometry
open Gutoe.GaugeGroupSU3
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

/-- Gauge-redundancy relation on positive transfer weights:
`W₂` differs from `W₁` only by a positive row-wise rescaling. -/
def RowScaleEquivalent
    (W₁ W₂ : Fin 3 → Fin 3 → ℝ) : Prop :=
  ∃ s : Fin 3 → ℝ, (∀ i, 0 < s i) ∧ ∀ i j, W₂ i j = s i * W₁ i j

/-- If two positive weight matrices induce the same normalized kernel, then they
are row-scale equivalent. This is the core "gauge redundancy" theorem for the
transfer lane: equal physical kernel implies only row-normalization gauge
freedom differs. -/
theorem kernel_eq_implies_row_scale_equivalent
    (W₁ W₂ : Fin 3 → Fin 3 → ℝ)
    (hW₁ : ∀ i j, 0 < W₁ i j)
    (hW₂ : ∀ i j, 0 < W₂ i j)
    (hK : normalizedKernelFromWeights W₁ = normalizedKernelFromWeights W₂) :
    RowScaleEquivalent W₁ W₂ := by
  let r₁ : Fin 3 → ℝ := fun i => ∑ k : Fin 3, W₁ i k
  let r₂ : Fin 3 → ℝ := fun i => ∑ k : Fin 3, W₂ i k
  have hr₁_pos : ∀ i, 0 < r₁ i := by
    intro i
    unfold r₁
    have h0 : 0 < W₁ i 0 := hW₁ i 0
    have h1 : 0 < W₁ i 1 := hW₁ i 1
    have h2 : 0 < W₁ i 2 := hW₁ i 2
    have hsum : (∑ k : Fin 3, W₁ i k) = W₁ i 0 + W₁ i 1 + W₁ i 2 := by
      simpa [Fin.sum_univ_three]
    rw [hsum]
    nlinarith
  have hr₂_pos : ∀ i, 0 < r₂ i := by
    intro i
    unfold r₂
    have h0 : 0 < W₂ i 0 := hW₂ i 0
    have h1 : 0 < W₂ i 1 := hW₂ i 1
    have h2 : 0 < W₂ i 2 := hW₂ i 2
    have hsum : (∑ k : Fin 3, W₂ i k) = W₂ i 0 + W₂ i 1 + W₂ i 2 := by
      simpa [Fin.sum_univ_three]
    rw [hsum]
    nlinarith
  refine ⟨fun i => r₂ i / r₁ i, ?_, ?_⟩
  · intro i
    exact div_pos (hr₂_pos i) (hr₁_pos i)
  · intro i j
    have hk : normalizedKernelFromWeights W₁ i j = normalizedKernelFromWeights W₂ i j := by
      simpa using congrArg (fun f => f i j) hK
    have hr₁_ne : r₁ i ≠ 0 := ne_of_gt (hr₁_pos i)
    have hr₂_ne : r₂ i ≠ 0 := ne_of_gt (hr₂_pos i)
    have hcross : W₁ i j * r₂ i = W₂ i j * r₁ i := by
      have hk' : W₁ i j / r₁ i = W₂ i j / r₂ i := by
        simpa [normalizedKernelFromWeights, r₁, r₂] using hk
      field_simp [hr₁_ne, hr₂_ne] at hk'
      nlinarith [hk']
    calc
      W₂ i j = (W₂ i j * r₁ i) / r₁ i := by field_simp [hr₁_ne]
      _ = (W₁ i j * r₂ i) / r₁ i := by rw [hcross]
      _ = (r₂ i / r₁ i) * W₁ i j := by ring

/-- Row-scale equivalent positive weights induce the same normalized kernel. -/
theorem row_scale_equivalent_implies_kernel_eq
    (W₁ W₂ : Fin 3 → Fin 3 → ℝ)
    (hW₁ : ∀ i j, 0 < W₁ i j)
    (hEq : RowScaleEquivalent W₁ W₂) :
    normalizedKernelFromWeights W₁ = normalizedKernelFromWeights W₂ := by
  rcases hEq with ⟨s, hs, hmul⟩
  have hW₂ : ∀ i j, 0 < W₂ i j := by
    intro i j
    rw [hmul i j]
    exact mul_pos (hs i) (hW₁ i j)
  have hW₂eq : W₂ = (fun i j => s i * W₁ i j) := by
    funext i j
    exact hmul i j
  calc
    normalizedKernelFromWeights W₁
        = normalizedKernelFromWeights (fun i j => s i * W₁ i j) := by
          exact (normalizedKernel_row_scale_invariant W₁ hW₁ s hs).symm
    _ = normalizedKernelFromWeights W₂ := by
          rw [hW₂eq]

/-- For strictly positive weights, kernel equality is equivalent to row-scale
gauge equivalence. -/
theorem row_scale_equivalent_iff_kernel_eq
    (W₁ W₂ : Fin 3 → Fin 3 → ℝ)
    (hW₁ : ∀ i j, 0 < W₁ i j)
    (hW₂ : ∀ i j, 0 < W₂ i j) :
    RowScaleEquivalent W₁ W₂ ↔
      normalizedKernelFromWeights W₁ = normalizedKernelFromWeights W₂ := by
  constructor
  · intro hEq
    exact row_scale_equivalent_implies_kernel_eq W₁ W₂ hW₁ hEq
  · intro hK
    exact kernel_eq_implies_row_scale_equivalent W₁ W₂ hW₁ hW₂ hK

/-- Full-SU(3)-lane completeness witness in transfer form:
if two Wilson actions have the same normalized transfer kernel, then their
Boltzmann weights differ only by row-wise gauge rescaling (no independent extra
degrees of freedom beyond the kernel). -/
theorem full_su3_kernel_completeness
    (S₁ S₂ : WilsonAction)
    (hK : wilsonKernel 1 S₁ = wilsonKernel 1 S₂) :
    RowScaleEquivalent (wilsonWeight 1 S₁) (wilsonWeight 1 S₂) := by
  exact kernel_eq_implies_row_scale_equivalent
    (wilsonWeight 1 S₁)
    (wilsonWeight 1 S₂)
    (wilson_weight_pos 1 S₁)
    (wilson_weight_pos 1 S₂)
    hK

/-- Fiber-constancy principle (Path-2 analogue of "coset integral collapse"):
any observable that factors through the normalized Wilson kernel is constant on
kernel fibers of the full Wilson-action space. -/
theorem full_su3_observable_fiber_const
    {β : Type}
    (ObsCenter : Matrix (Fin 3) (Fin 3) ℝ → β)
    (S₁ S₂ : WilsonAction)
    (hK : wilsonKernel 1 S₁ = wilsonKernel 1 S₂) :
    ObsCenter (wilsonKernel 1 S₁) = ObsCenter (wilsonKernel 1 S₂) := by
  simpa [hK]

/-- Row-scale gauge-orbit hypothesis implies kernel-level fiber constancy. -/
theorem kernel_fiber_const_of_row_scale_orbit
    {C F : Type}
    [Fintype C]
    (lift : C → F → WilsonAction)
    (f₀ : F)
    (hscale :
      ∀ c f,
        RowScaleEquivalent
          (wilsonWeight 1 (lift c f₀))
          (wilsonWeight 1 (lift c f))) :
    ∀ c f, wilsonKernel 1 (lift c f) = wilsonKernel 1 (lift c f₀) := by
  intro c f
  have hker :
      normalizedKernelFromWeights (wilsonWeight 1 (lift c f₀)) =
        normalizedKernelFromWeights (wilsonWeight 1 (lift c f)) :=
    row_scale_equivalent_implies_kernel_eq
      (wilsonWeight 1 (lift c f₀))
      (wilsonWeight 1 (lift c f))
      (wilson_weight_pos 1 (lift c f₀))
      (hscale c f)
  simpa [wilsonKernel] using hker.symm

/-- Finite-fiber expectation collapse: if the kernel is constant on the fiber
(`f`-direction), then normalized expectation over product weights collapses to
the base (`c`) expectation. This is the finite transfer-lane analogue of
Haar-fiber collapse for gauge-invariant observables. -/
theorem finite_fiber_expectation_collapse
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
    (∑ c : C, ∑ f : F,
      (wC c) * (wF f) * ObsCenter (wilsonKernel 1 (lift c f))) =
      ∑ c : C, (wC c) * ObsCenter (wilsonKernel 1 (lift c f₀)) := by
  have hconst :
      ∀ c f,
        ObsCenter (wilsonKernel 1 (lift c f)) =
          ObsCenter (wilsonKernel 1 (lift c f₀)) := by
    intro c f
    exact congrArg ObsCenter (kernel_fiber_const_of_row_scale_orbit lift f₀ hscale c f)
  calc
    (∑ c : C, ∑ f : F,
      (wC c) * (wF f) * ObsCenter (wilsonKernel 1 (lift c f)))
        = ∑ c : C, ∑ f : F, (wC c) * (wF f) * ObsCenter (wilsonKernel 1 (lift c f₀)) := by
            refine Finset.sum_congr rfl ?_
            intro c hc
            refine Finset.sum_congr rfl ?_
            intro f hf
            rw [hconst c f]
    _ = ∑ c : C, (wC c * ObsCenter (wilsonKernel 1 (lift c f₀))) * (∑ f : F, wF f) := by
          refine Finset.sum_congr rfl ?_
          intro c hc
          have hmul :
              ∑ f : F, (wC c) * (wF f) * ObsCenter (wilsonKernel 1 (lift c f₀))
                = ∑ f : F, wF f * (wC c * ObsCenter (wilsonKernel 1 (lift c f₀))) := by
                  refine Finset.sum_congr rfl ?_
                  intro f hf
                  ring
          have hsumMul :
              (∑ f : F, wF f * (wC c * ObsCenter (wilsonKernel 1 (lift c f₀)))) =
                (∑ f : F, wF f) * (wC c * ObsCenter (wilsonKernel 1 (lift c f₀))) := by
            simpa using
              ((Finset.sum_mul
                (s := (Finset.univ : Finset F))
                (f := fun f : F => wF f)
                (a := (wC c * ObsCenter (wilsonKernel 1 (lift c f₀))))).symm)
          rw [hmul, hsumMul]
          ring
    _ = ∑ c : C, (wC c * ObsCenter (wilsonKernel 1 (lift c f₀))) * 1 := by
          simp [hwF]
    _ = ∑ c : C, (wC c) * ObsCenter (wilsonKernel 1 (lift c f₀)) := by
          refine Finset.sum_congr rfl ?_
          intro c hc
          ring

/-- Effective center-sector plaquette action induced by nearest-neighbor Z₃
targets with Laplace floor `alpha`. -/
noncomputable def z3CenterPlaquetteAction
    (target : Z3NearestNeighborTargets) (alpha : ℝ) : WilsonAction :=
  fun i j => -Real.log ((z3NearestNeighborCounts target i j : ℝ) + alpha)

/-- Unit-coupling Wilson weight of the center-sector action is exactly the
Laplace-shifted transition count. -/
theorem z3_center_plaquette_weight_eq_shifted_count
    (target : Z3NearestNeighborTargets) {alpha : ℝ} (ha : 0 < alpha) :
    ∀ i j,
      wilsonWeight 1 (z3CenterPlaquetteAction target alpha) i j =
        ((z3NearestNeighborCounts target i j : ℝ) + alpha) := by
  intro i j
  unfold wilsonWeight z3CenterPlaquetteAction
  let x : ℝ := (z3NearestNeighborCounts target i j : ℝ) + alpha
  have hcount : 0 ≤ (z3NearestNeighborCounts target i j : ℝ) := by positivity
  have hpos : 0 < x := by
    unfold x
    linarith
  have hmul : Real.exp (-(1 : ℝ) * (-Real.log x)) = Real.exp (Real.log x) := by
    simp
  rw [hmul]
  simpa [x] using Real.exp_log hpos

/-- The unit-coupling Wilson kernel from the center-sector plaquette action is
exactly the Laplace-smoothed nearest-neighbor transfer kernel. -/
theorem z3_center_plaquette_kernel_eq_smoothed_transition
    (target : Z3NearestNeighborTargets) {alpha : ℝ} (ha : 0 < alpha) :
    wilsonKernel 1 (z3CenterPlaquetteAction target alpha) =
      smoothedTransition
        (z3NearestNeighborCounts target)
        (rowTotalsFromCounts (z3NearestNeighborCounts target))
        alpha := by
  funext i j
  unfold wilsonKernel normalizedKernelFromWeights smoothedTransition smoothEntry
  have hnum :
      wilsonWeight 1 (z3CenterPlaquetteAction target alpha) i j =
        ((z3NearestNeighborCounts target i j : ℝ) + alpha) :=
    z3_center_plaquette_weight_eq_shifted_count target ha i j
  have hsum :
      (∑ k : Fin 3, wilsonWeight 1 (z3CenterPlaquetteAction target alpha) i k) =
        ∑ k : Fin 3, ((z3NearestNeighborCounts target i k : ℝ) + alpha) := by
    refine Finset.sum_congr rfl ?_
    intro k hk
    exact z3_center_plaquette_weight_eq_shifted_count target ha i k
  have hsplit :
      (∑ k : Fin 3, ((z3NearestNeighborCounts target i k : ℝ) + alpha)) =
        (∑ k : Fin 3, (z3NearestNeighborCounts target i k : ℝ)) + (∑ _k : Fin 3, alpha) := by
    simpa [Finset.sum_add_distrib]
  have hconst : (∑ _k : Fin 3, alpha) = 3 * alpha := by
    norm_num [Fin.sum_univ_three]
  have hrow :
      (∑ k : Fin 3, (z3NearestNeighborCounts target i k : ℝ)) =
        (rowTotalsFromCounts (z3NearestNeighborCounts target) i : ℝ) := by
    norm_num [rowTotalsFromCounts, Fin.sum_univ_three]
  rw [hnum, hsum, hsplit, hconst, hrow]

/-- Center-sector plaquette actions over a refinement schedule. -/
noncomputable def centerPlaquetteActionSchedule
    (W : WilsonZ3Action) (alpha : ℝ) : ℕ → WilsonAction :=
  fun n => z3CenterPlaquetteAction (W.targetSchedule n) alpha

/-- The Wilson kernel from center-sector plaquette actions equals the
nearest-neighbor transfer kernel at every refinement step. -/
theorem center_plaquette_schedule_kernel_eq_transfer
    (W : WilsonZ3Action) {alpha : ℝ} (ha : 0 < alpha) :
    ∀ n,
      wilsonKernel 1 (centerPlaquetteActionSchedule W alpha n) =
        smoothedTransition
          (z3NearestNeighborCounts (W.targetSchedule n))
          (rowTotalsFromCounts (z3NearestNeighborCounts (W.targetSchedule n)))
          alpha := by
  intro n
  simpa [centerPlaquetteActionSchedule] using
    z3_center_plaquette_kernel_eq_smoothed_transition (W.targetSchedule n) ha

/-!
  Path-2 structural bridge targets

  C1. Construction: Z₃ orbit carrier gives SU(3) structural data.
  C2. Completeness/Faithfulness (center lane): count matrices are losslessly
      encoded in center-compatible Wilson actions.
  C3. Gap transfer: any Z₃ schedule gap floor transfers to the Wilson lane.
-/

/-- C1 (construction): the Z₃ quark-orbit carrier provides the SU(3) structural
data used by the transfer lane (fundamental carrier size, adjoint dimension,
and nontrivial su(3) commutator witness). -/
theorem c1_z3_to_su3_structural_construction :
    quarkOrbit.card = 3 ∧
    Nonempty ({s // s ∈ quarkOrbit} ≃ Fin 3) ∧
    quarkOrbit.card ^ 2 - 1 = 8 ∧
    gm₁ * gm₂ - gm₂ * gm₁ = (2 * Complex.I) • gm₃ := by
  exact ⟨quarkOrbit_card, quarkOrbit_equiv_fin3, quarks_predict_gluon_count, su3_comm_12⟩

/-- Center-compatible Wilson action parameterized directly by local Z₃ count
matrices. -/
noncomputable def z3CenterPlaquetteActionFromCounts
    (counts : Fin 3 → Fin 3 → ℕ) (alpha : ℝ) : WilsonAction :=
  fun i j => -Real.log ((counts i j : ℝ) + alpha)

/-- Decoder from a Wilson action back to shifted count coordinates (real-valued
form; no rounding is used in the theorem chain). -/
noncomputable def centerActionCountDecoder
    (S : WilsonAction) (alpha : ℝ) : Fin 3 → Fin 3 → ℝ :=
  fun i j => Real.exp (-S i j) - alpha

/-- The count-to-center-action map is lossless at real-valued count
coordinates. -/
theorem center_action_decoder_roundtrip
    (counts : Fin 3 → Fin 3 → ℕ) {alpha : ℝ} (ha : 0 < alpha) :
    centerActionCountDecoder (z3CenterPlaquetteActionFromCounts counts alpha) alpha =
      fun i j => (counts i j : ℝ) := by
  funext i
  funext j
  unfold centerActionCountDecoder z3CenterPlaquetteActionFromCounts
  let x : ℝ := (counts i j : ℝ) + alpha
  have hcount : 0 ≤ (counts i j : ℝ) := by positivity
  have hpos : 0 < x := by
    unfold x
    linarith
  calc
    Real.exp (-(-Real.log x)) - alpha = Real.exp (Real.log x) - alpha := by simp
    _ = x - alpha := by rw [Real.exp_log hpos]
    _ = (counts i j : ℝ) := by simp [x]

/-- C2 (completeness/fidelity, center lane): fixed-`alpha` center actions are
faithful to the underlying Z₃ local count matrix. -/
theorem center_action_from_counts_faithful
    {alpha : ℝ} (ha : 0 < alpha)
    {counts₁ counts₂ : Fin 3 → Fin 3 → ℕ}
    (hEq :
      z3CenterPlaquetteActionFromCounts counts₁ alpha =
        z3CenterPlaquetteActionFromCounts counts₂ alpha) :
    counts₁ = counts₂ := by
  have hDecoded :
      centerActionCountDecoder (z3CenterPlaquetteActionFromCounts counts₁ alpha) alpha =
        centerActionCountDecoder (z3CenterPlaquetteActionFromCounts counts₂ alpha) alpha := by
    simpa [hEq]
  rw [center_action_decoder_roundtrip counts₁ ha,
      center_action_decoder_roundtrip counts₂ ha] at hDecoded
  funext i
  funext j
  have hreal : (counts₁ i j : ℝ) = (counts₂ i j : ℝ) := by
    simpa using congrArg (fun f => f i j) hDecoded
  exact_mod_cast hreal

/-- Subtype of Wilson actions representable by Z₃ center-count data at fixed
`alpha`. -/
def CenterCompatibleAction (alpha : ℝ) :=
  {S : WilsonAction // ∃ counts : Fin 3 → Fin 3 → ℕ,
    S = z3CenterPlaquetteActionFromCounts counts alpha}

/-- Canonical encoding of local Z₃ count matrices into center-compatible
Wilson actions. -/
noncomputable def countsToCenterCompatibleAction
    (alpha : ℝ) (counts : Fin 3 → Fin 3 → ℕ) :
    CenterCompatibleAction alpha :=
  ⟨z3CenterPlaquetteActionFromCounts counts alpha, ⟨counts, rfl⟩⟩

/-- Surjectivity onto center-compatible Wilson actions is by construction. -/
theorem counts_to_center_compatible_surjective
    (alpha : ℝ) :
    Function.Surjective (countsToCenterCompatibleAction alpha) := by
  intro S
  rcases S.property with ⟨counts, hcounts⟩
  refine ⟨counts, ?_⟩
  apply Subtype.ext
  simpa [countsToCenterCompatibleAction] using hcounts.symm

/-- Injectivity of the count encoding map (fixed `alpha>0`) from C2-faithfulness. -/
theorem counts_to_center_compatible_injective
    {alpha : ℝ} (ha : 0 < alpha) :
    Function.Injective (countsToCenterCompatibleAction alpha) := by
  intro counts₁ counts₂ hEq
  apply center_action_from_counts_faithful ha
  exact congrArg Subtype.val hEq

/-- C2 as a package: for fixed `alpha>0`, Z₃ local count matrices and
center-compatible Wilson actions are in bijection. -/
theorem c2_counts_center_action_bijective
    {alpha : ℝ} (ha : 0 < alpha) :
    Function.Bijective (countsToCenterCompatibleAction alpha) := by
  exact ⟨counts_to_center_compatible_injective ha,
    counts_to_center_compatible_surjective alpha⟩

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

/-- C3 (gap transfer): any Z₃-schedule Doeblin gap floor transfers directly to
the Wilson schedule induced by the same structural target schedule. -/
theorem c3_gap_transfer_from_z3_schedule
    (W : WilsonZ3Action)
    (a_t : ℕ → ℝ)
    (alpha : ℝ)
    (hgap :
      ∃ c : ℝ, 0 < c ∧
        ∀ n, c ≤ doeblinGapLowerBound (a_t n)
          (minorizationEps (z3NearestNeighborRowTotalsSchedule W.targetSchedule n) alpha)) :
    ∃ c : ℝ, 0 < c ∧
      ∀ n, c ≤ doeblinGapLowerBound (a_t n)
        (minorizationEps (wilsonRowTotalsSchedule W n) alpha) := by
  rcases hgap with ⟨c, hc, hbound⟩
  refine ⟨c, hc, ?_⟩
  intro n
  simpa [wilsonRowTotalsSchedule, z3NearestNeighborRowTotalsSchedule] using hbound n

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

/-- C3 instantiated from the structural continuum-survival theorem:
the Wilson lane inherits a non-vanishing continuum mass-gap lower bound from
the Cl(1,3)→Z₃ nearest-neighbor schedule. -/
theorem c3_wilson_gap_nonvanishing_from_clifford_z3
    (W : WilsonZ3Action)
    (a_t : ℕ → ℝ)
    (alpha : ℝ)
    (ha_t_pos : ∀ n, 0 < a_t n)
    (ha_t_cap : ∃ aCap, 0 < aCap ∧ ∀ n, a_t n ≤ aCap)
    (ha : 0 < alpha) :
    ∃ c : ℝ, 0 < c ∧
      ∀ n, c ≤ doeblinGapLowerBound (a_t n)
        (minorizationEps (wilsonRowTotalsSchedule W n) alpha) := by
  exact c3_gap_transfer_from_z3_schedule
    W
    a_t
    alpha
    (continuum_survival_gap_nonvanishing_of_z3_nn_schedule
      a_t
      W.targetSchedule
      alpha
      ha_t_pos
      ha_t_cap
      ha)

end Gutoe.YangMillsWilsonBridge
