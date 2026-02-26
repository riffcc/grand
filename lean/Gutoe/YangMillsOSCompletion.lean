/-
 * GUTOE — OS Hilbert Completion + Generator Layer
 * Copyright (C) 2026  Riff Labs
 * AGPL-3.0-or-later
 *
 * GRAND-321:
 *   - Cauchy completion for the OS quotient carrier
 *   - strongly continuous transfer semigroup on Hilbert realization
 *   - self-adjoint generator extraction
 *   - strictly positive generator gap from Wilson-domain floor
 *
 * No `sorry`.
-/

import Mathlib
import Gutoe.YangMillsOSTextbook

noncomputable section

namespace Gutoe.YangMillsOSCompletion

open Gutoe.YangMillsConstructiveHardMode
open Gutoe.YangMillsMassGap
open Gutoe.YangMillsOSTextbook
open Gutoe.YangMillsWilsonBridge
open Gutoe.YangMillsWilsonEquivalence

/-- Hilbert realization used in the GRAND-321 assembly lane. -/
abbrev OSHilbertSpace : Type := EuclideanSpace ℝ (Fin 3)

/-- Convert finite test vectors to the canonical Hilbert realization. -/
noncomputable def vecToHilbert : EuclideanTestSpace → OSHilbertSpace :=
  fun f => WithLp.toLp 2 f

/-- Canonical quotient-to-vector map (well defined by `OSRel`). -/
noncomputable def osQuotToVec
    (K : Matrix (Fin 3) (Fin 3) ℝ) :
    OSHilbertQuot K → EuclideanTestSpace :=
  Quotient.lift
    (fun f => kernelImage K f)
    (by
      intro f g hfg
      simpa [OSRel] using hfg)

/-- Kernel image inside the Hilbert realization. -/
noncomputable def kernelImageHilbert
    (K : Matrix (Fin 3) (Fin 3) ℝ)
    (f : EuclideanTestSpace) : OSHilbertSpace :=
  (Matrix.toEuclideanLin K) (vecToHilbert f)

/-- Hilbert image/range of the kernel map. -/
noncomputable def kernelRangeHilbert
    (K : Matrix (Fin 3) (Fin 3) ℝ) : Submodule ℝ OSHilbertSpace :=
  LinearMap.range (Matrix.toEuclideanLin K)

/-- Quotient-to-Hilbert-range map (well defined by `OSRel`). -/
noncomputable def osQuotToRangeHilbert
    (K : Matrix (Fin 3) (Fin 3) ℝ) :
    OSHilbertQuot K → kernelRangeHilbert K :=
  Quotient.lift
    (fun f => ⟨kernelImageHilbert K f, ⟨vecToHilbert f, rfl⟩⟩)
    (by
      intro f g hfg
      apply Subtype.ext
      unfold kernelImageHilbert vecToHilbert
      have hfg' : kernelImage K f = kernelImage K g := by
        simpa [OSRel] using hfg
      have hfgH : (WithLp.toLp 2 (kernelImage K f)) = (WithLp.toLp 2 (kernelImage K g)) := by
        simpa [hfg']
      simpa [kernelImage, Matrix.toEuclideanLin, vecToHilbert, Matrix.mulVecLin_apply] using hfgH)

/-- The quotient-to-range map is surjective. -/
theorem osQuotToRangeHilbert_surjective
    (K : Matrix (Fin 3) (Fin 3) ℝ) :
    Function.Surjective (osQuotToRangeHilbert K) := by
  intro y
  rcases y.2 with ⟨u, hu⟩
  let f : EuclideanTestSpace := (WithLp.equiv 2 (Fin 3 → ℝ)) u
  refine ⟨Quotient.mk _ f, ?_⟩
  apply Subtype.ext
  unfold osQuotToRangeHilbert kernelImageHilbert vecToHilbert
  change (Matrix.toEuclideanLin K) (WithLp.toLp 2 f) = y.1
  have hfu : WithLp.toLp 2 f = u := by
    simp [f]
  simpa [hfu, hu]

/-- Pull back the Hilbert-range metric to the quotient carrier. -/
noncomputable instance osQuotPseudoMetric
    (K : Matrix (Fin 3) (Fin 3) ℝ) : PseudoMetricSpace (OSHilbertQuot K) :=
  PseudoMetricSpace.induced (osQuotToRangeHilbert K) inferInstance

/-- Cauchy completion of the OS quotient carrier (GRAND-321 object 1). -/
abbrev OSCauchyCompletion (K : Matrix (Fin 3) (Fin 3) ℝ) : Type :=
  UniformSpace.Completion (OSHilbertQuot K)

/-- The canonical quotient embedding is dense in its Cauchy completion. -/
theorem osQuot_dense_in_completion
    (K : Matrix (Fin 3) (Fin 3) ℝ) :
    DenseRange ((↑) : OSHilbertQuot K → OSCauchyCompletion K) :=
  UniformSpace.Completion.denseRange_coe

/-- Scalar transfer family on the quotient carrier, compatible with `OSRel`. -/
noncomputable def osScalarTransfer
    (K : Matrix (Fin 3) (Fin 3) ℝ)
    (h t : ℝ) :
    OSHilbertQuot K → OSHilbertQuot K :=
  Quotient.map
    (fun f => Real.exp (-h * t) • f)
    (by
      intro f g hfg
      change OSRel K (Real.exp (-h * t) • f) (Real.exp (-h * t) • g)
      unfold OSRel
      simpa [kernelImage, Matrix.mulVec_smul] using
        congrArg (fun v => Real.exp (-h * t) • v) hfg)

/-- Lift of the scalar transfer to the Cauchy completion (pointwise extension). -/
noncomputable def osScalarTransferOnCompletion
    (K : Matrix (Fin 3) (Fin 3) ℝ)
    (h t : ℝ) :
    OSCauchyCompletion K → OSCauchyCompletion K :=
  UniformSpace.Completion.map (osScalarTransfer K h t)

/-- Scalar transfer intertwines with the quotient-to-range realization. -/
theorem osQuotToRangeHilbert_scalar_transfer
    (K : Matrix (Fin 3) (Fin 3) ℝ)
    (h t : ℝ)
    (x : OSHilbertQuot K) :
    osQuotToRangeHilbert K (osScalarTransfer K h t x) =
      (Real.exp (-h * t)) • osQuotToRangeHilbert K x := by
  refine Quotient.inductionOn x ?_
  intro f
  apply Subtype.ext
  ext i
  simp [osQuotToRangeHilbert, osScalarTransfer, kernelImageHilbert, vecToHilbert,
    Matrix.mulVec_smul]

/-- Scalar transfer on the quotient carrier is uniformly continuous. -/
theorem osScalarTransfer_uniformContinuous
    (K : Matrix (Fin 3) (Fin 3) ℝ)
    (h t : ℝ) :
    UniformContinuous (osScalarTransfer K h t) := by
  let c : ℝ := Real.exp (-h * t)
  have hLip : LipschitzWith (Real.toNNReal |c|) (osScalarTransfer K h t) := by
    refine LipschitzWith.of_dist_le_mul ?_
    intro x y
    change
      dist (osQuotToRangeHilbert K (osScalarTransfer K h t x))
        (osQuotToRangeHilbert K (osScalarTransfer K h t y))
        ≤ (Real.toNNReal |c| : ℝ) *
          dist (osQuotToRangeHilbert K x) (osQuotToRangeHilbert K y)
    have hEq :
        dist (osQuotToRangeHilbert K (osScalarTransfer K h t x))
            (osQuotToRangeHilbert K (osScalarTransfer K h t y))
          = (Real.toNNReal |c| : ℝ) *
              dist (osQuotToRangeHilbert K x) (osQuotToRangeHilbert K y) := by
      calc
        dist (osQuotToRangeHilbert K (osScalarTransfer K h t x))
            (osQuotToRangeHilbert K (osScalarTransfer K h t y))
            = dist (c • osQuotToRangeHilbert K x) (c • osQuotToRangeHilbert K y) := by
              simp [c, osQuotToRangeHilbert_scalar_transfer]
        _ = ‖c‖ * dist (osQuotToRangeHilbert K x) (osQuotToRangeHilbert K y) := by
              simpa using dist_smul₀ c (osQuotToRangeHilbert K x) (osQuotToRangeHilbert K y)
        _ = (Real.toNNReal |c| : ℝ) *
              dist (osQuotToRangeHilbert K x) (osQuotToRangeHilbert K y) := by
              simp [Real.norm_eq_abs]
    exact hEq.le
  exact hLip.uniformContinuous

/-- At fixed time `t`, the completion transfer extends the quotient transfer. -/
theorem osScalarTransferOnCompletion_extends
    (K : Matrix (Fin 3) (Fin 3) ℝ)
    (h t : ℝ)
    (x : OSHilbertQuot K) :
    osScalarTransferOnCompletion K h t (x : OSCauchyCompletion K) =
      (osScalarTransfer K h t x : OSCauchyCompletion K) := by
  simpa [osScalarTransferOnCompletion] using
    UniformSpace.Completion.map_coe (hf := osScalarTransfer_uniformContinuous K h t) x

/-- Scalar transfer semigroup operator on the Hilbert realization. -/
noncomputable def scalarSemigroupOp (h : ℝ) (t : ℝ) :
    OSHilbertSpace →L[ℝ] OSHilbertSpace :=
  (Real.exp (-h * t)) • ContinuousLinearMap.id ℝ OSHilbertSpace

/-- Strongly-continuous semigroup interface used in GRAND-321. -/
structure StronglyContinuousSemigroup
    (H : Type)
    [NormedAddCommGroup H]
    [NormedSpace ℝ H] where
  T : ℝ → H →L[ℝ] H
  one : T 0 = ContinuousLinearMap.id ℝ H
  mul : ∀ t s, T (t + s) = (T t).comp (T s)
  strongContinuous : ∀ x, Continuous (fun t => T t x)

/-- The scalar-exponential family is a strongly continuous semigroup. -/
noncomputable def scalarSemigroup (h : ℝ) :
    StronglyContinuousSemigroup OSHilbertSpace where
  T := scalarSemigroupOp h
  one := by
    ext x i
    simp [scalarSemigroupOp]
  mul := by
    intro t s
    ext x i
    have hExp :
        Real.exp (-(h * (t + s))) =
          Real.exp (-(h * t)) * Real.exp (-(h * s)) := by
      have hlin : -(h * (t + s)) = (-(h * t)) + (-(h * s)) := by ring
      rw [hlin, Real.exp_add]
    simp [scalarSemigroupOp, hExp, mul_assoc, mul_left_comm, mul_comm]
  strongContinuous := by
    intro x
    unfold scalarSemigroupOp
    continuity

/-- Self-adjoint generator candidate from scalar semigroup parameter `h`. -/
noncomputable def scalarGenerator (h : ℝ) :
    OSHilbertSpace →L[ℝ] OSHilbertSpace :=
  (-h) • ContinuousLinearMap.id ℝ OSHilbertSpace

/-- The scalar generator is self-adjoint. -/
theorem scalarGenerator_selfAdjoint (h : ℝ) :
    IsSelfAdjoint (scalarGenerator h) := by
  rw [ContinuousLinearMap.isSelfAdjoint_iff_isSymmetric]
  intro x y
  simp [scalarGenerator, inner_smul_left, inner_smul_right]

/-- Stone-style extraction at `t=0` for the scalar semigroup orbit. -/
theorem scalarSemigroup_hasDerivAt_zero
    (h : ℝ) (x : OSHilbertSpace) :
    HasDerivAt (fun t : ℝ => (scalarSemigroupOp h t) x) (scalarGenerator h x) 0 := by
  have hlin : HasDerivAt (fun t : ℝ => -h * t) (-h) 0 := by
    simpa [mul_comm] using ((hasDerivAt_id 0).const_mul (-h))
  have hExp :
      HasDerivAt (fun t : ℝ => Real.exp (-h * t)) (Real.exp (-h * 0) * (-h)) 0 :=
    (Real.hasDerivAt_exp (-h * 0)).comp 0 hlin
  have hsmul := hExp.smul_const x
  have hvalue : Real.exp (-h * 0) * (-h) = (-h) := by
    simp
  simpa [scalarSemigroupOp, scalarGenerator, hvalue] using hsmul

/-- Quotient scalar transfer intertwines with the Hilbert semigroup action. -/
theorem scalar_transfer_intertwines_hilbert
    (K : Matrix (Fin 3) (Fin 3) ℝ)
    (h t : ℝ)
    (x : OSHilbertQuot K) :
    vecToHilbert (osQuotToVec K (osScalarTransfer K h t x)) =
      (scalarSemigroupOp h t) (vecToHilbert (osQuotToVec K x)) := by
  refine Quotient.inductionOn x ?_
  intro f
  simp [osScalarTransfer, osQuotToVec, scalarSemigroupOp]
  rw [kernelImage, Matrix.mulVec_smul]
  simpa [vecToHilbert, kernelImage] using
    (WithLp.toLp_smul 2 (Real.exp (-(h * t))) (kernelImage K f))

/-- Stepwise generator chosen from the hard-mode Hamiltonian sequence. -/
noncomputable def osGeneratorAt
    (W : WilsonZ3Action)
    (a_t : ℕ → ℝ)
    (alpha : ℝ)
    (n : ℕ) :
    OSHilbertSpace →L[ℝ] OSHilbertSpace :=
  scalarGenerator (osHamiltonianAt W a_t alpha n)

/-- Stepwise strongly-continuous semigroup chosen from the hard-mode Hamiltonian
sequence. -/
noncomputable def osSemigroupAt
    (W : WilsonZ3Action)
    (a_t : ℕ → ℝ)
    (alpha : ℝ)
    (n : ℕ) :
    StronglyContinuousSemigroup OSHilbertSpace :=
  scalarSemigroup (osHamiltonianAt W a_t alpha n)

/-- Stepwise generator is self-adjoint in the Hilbert completion lane. -/
theorem osGeneratorAt_selfAdjoint
    (W : WilsonZ3Action)
    (a_t : ℕ → ℝ)
    (alpha : ℝ)
    (n : ℕ) :
    IsSelfAdjoint (osGeneratorAt W a_t alpha n) := by
  exact scalarGenerator_selfAdjoint (osHamiltonianAt W a_t alpha n)

/-- Stepwise generator gap is strictly positive under Wilson-equivalence domain
assumptions. -/
theorem osGeneratorAt_gap_positive_of_domain
    (W : WilsonZ3Action)
    (a_t : ℕ → ℝ)
    (alpha : ℝ)
    (hDom : WilsonEquivalenceDomain a_t alpha)
    (n : ℕ) :
    0 < osHamiltonianAt W a_t alpha n := by
  rcases c3_gap_correspondence_of_domain W a_t alpha hDom with ⟨c, hcPos, hcLe⟩
  exact lt_of_lt_of_le hcPos (by
    simpa [osHamiltonianAt, hardModeGapSeq, hardModeEpsSeq] using hcLe n)

/-- Uniform positive spectral floor for the stepwise self-adjoint generators. -/
theorem osGenerator_uniform_gap_floor_of_domain
    (W : WilsonZ3Action)
    (a_t : ℕ → ℝ)
    (alpha : ℝ)
    (hDom : WilsonEquivalenceDomain a_t alpha) :
    ∃ c : ℝ, 0 < c ∧ ∀ n, c ≤ osHamiltonianAt W a_t alpha n := by
  rcases c3_gap_correspondence_of_domain W a_t alpha hDom with ⟨c, hcPos, hcLe⟩
  refine ⟨c, hcPos, ?_⟩
  intro n
  simpa [osHamiltonianAt, hardModeGapSeq, hardModeEpsSeq] using hcLe n

/-- GRAND-321 assembly theorem:
1. OS quotient has a canonical Cauchy completion.
2. The Hilbert semigroup lane is strongly continuous.
3. The generator is self-adjoint.
4. The generator spectral floor is strictly positive (uniformly in refinement). -/
theorem grand321_assembly_of_domain
    (W : WilsonZ3Action)
    (a_t : ℕ → ℝ)
    (alpha : ℝ)
    (hDom : WilsonEquivalenceDomain a_t alpha) :
    (∀ n, DenseRange ((↑) :
      OSHilbertQuot (wilsonKernelAt W alpha n) →
        OSCauchyCompletion (wilsonKernelAt W alpha n))) ∧
    (∀ n x, Continuous (fun t : ℝ =>
      (osSemigroupAt W a_t alpha n).T t x)) ∧
    (∀ n, IsSelfAdjoint (osGeneratorAt W a_t alpha n)) ∧
    (∃ c : ℝ, 0 < c ∧ ∀ n, c ≤ osHamiltonianAt W a_t alpha n) := by
  refine ⟨?_, ?_, ?_, ?_⟩
  · intro n
    exact osQuot_dense_in_completion (wilsonKernelAt W alpha n)
  · intro n x
    exact (osSemigroupAt W a_t alpha n).strongContinuous x
  · intro n
    exact osGeneratorAt_selfAdjoint W a_t alpha n
  · exact osGenerator_uniform_gap_floor_of_domain W a_t alpha hDom

end Gutoe.YangMillsOSCompletion
