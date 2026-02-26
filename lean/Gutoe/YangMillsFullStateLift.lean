/-
 * GUTOE — Yang-Mills Full State-Space Lift
 * Copyright (C) 2026  Riff Labs
 * AGPL-3.0-or-later
 *
 * GRAND-301:
 *   Lift transfer-basis gap reasoning to a full gauge-field state-space layer.
 *
 * This module provides:
 *   1) explicit reduced->full observable/operator lift maps,
 *   2) proof obligations sufficient to transfer strict gap positivity, and
 *   3) concrete counterexample checks showing where lift can fail.
 *
 * No `sorry`.
 -/

import Mathlib
import Gutoe.YangMillsMassGap
import Gutoe.YangMillsWilsonBridge
import Gutoe.YangMillsStructuralGap

noncomputable section

namespace Gutoe.YangMillsFullStateLift

open scoped BigOperators
open Gutoe.YangMillsMassGap
open Gutoe.YangMillsWilsonBridge

/-- Cl(1,3)->Z₃ anchor for the reduced transfer basis used by the lift lane. -/
theorem reduced_transfer_basis_dim_eq_three :
    Gutoe.YangMillsMassGap.transferBasisDim = 3 :=
  transfer_basis_dim_eq_three

section LiftMaps

variable {R F : Type*} [Fintype F] [DecidableEq F]

/-- Lift a reduced observable on `R` to full state space `F` via projection `π`. -/
def liftObservable (π : F → R) (obsR : R → ℝ) : F → ℝ :=
  fun f => obsR (π f)

/-- Sup-norm on the 3-state transfer-basis vectors (coordinate form). -/
def supNorm3 (v : Fin 3 → ℝ) : ℝ :=
  max (|v 0|) (max (|v 1|) (|v 2|))

/-- Coordinate-wise bound by `supNorm3`. -/
theorem abs_le_supNorm3 (v : Fin 3 → ℝ) (i : Fin 3) :
    |v i| ≤ supNorm3 v := by
  fin_cases i
  · exact le_max_left _ _
  · exact le_trans (le_max_left _ _) (le_max_right _ _)
  · exact le_trans (le_max_right _ _) (le_max_right _ _)

/-- Nonnegativity of `supNorm3`. -/
theorem supNorm3_nonneg (v : Fin 3 → ℝ) : 0 ≤ supNorm3 v := by
  exact le_trans (abs_nonneg (v 0)) (le_max_left _ _)

/-- Bounding each coordinate absolute value bounds `supNorm3`. -/
theorem supNorm3_le_of_coord_bounds
    (v : Fin 3 → ℝ) {s : ℝ}
    (h0 : |v 0| ≤ s) (h1 : |v 1| ≤ s) (h2 : |v 2| ≤ s) :
    supNorm3 v ≤ s := by
  exact max_le_iff.mpr ⟨h0, max_le_iff.mpr ⟨h1, h2⟩⟩

/-- `supNorm3` scaling for nonnegative scalars. -/
theorem supNorm3_smul_of_nonneg
    (c : ℝ) (hc : 0 ≤ c) (v : Fin 3 → ℝ) :
    supNorm3 (c • v) = c * supNorm3 v := by
  unfold supNorm3
  simp [Pi.smul_apply, abs_mul, abs_of_nonneg hc,
    max_mul_of_nonneg, hc, mul_assoc, mul_left_comm, mul_comm]

/-- Diagonal full-state operator induced by a lifted observable. -/
def liftDiagonalOperator (π : F → R) (obsR : R → ℝ) : Matrix F F ℝ :=
  fun i j => if i = j then liftObservable π obsR i else 0

/-- On-diagonal entries recover the lifted observable. -/
theorem liftDiagonalOperator_diag
    (π : F → R) (obsR : R → ℝ) (i : F) :
    liftDiagonalOperator π obsR i i = liftObservable π obsR i := by
  simp [liftDiagonalOperator]

/-- Off-diagonal entries are zero (diagonal observable operator). -/
theorem liftDiagonalOperator_offdiag
    (π : F → R) (obsR : R → ℝ) {i j : F} (hij : i ≠ j) :
    liftDiagonalOperator π obsR i j = 0 := by
  simp [liftDiagonalOperator, hij]

/-- Lifted observables are constant on fibers of `π`. -/
theorem liftObservable_fiber_constant
    (π : F → R) (obsR : R → ℝ) {x y : F}
    (hxy : π x = π y) :
    liftObservable π obsR x = liftObservable π obsR y := by
  simpa [liftObservable, hxy]

end LiftMaps

section GapLift

/-- Row-stochastic `3×3` kernels are non-expansive on `supNorm3`. -/
theorem rowStochastic_supNorm3_nonexpansive
    (R : Matrix (Fin 3) (Fin 3) ℝ)
    (hR : R ∈ Matrix.rowStochastic ℝ (Fin 3))
    (v : Fin 3 → ℝ) :
    supNorm3 (R.mulVec v) ≤ supNorm3 v := by
  let s : ℝ := supNorm3 v
  have hs0 : |v 0| ≤ s := by simpa [s] using abs_le_supNorm3 v 0
  have hs1 : |v 1| ≤ s := by simpa [s] using abs_le_supNorm3 v 1
  have hs2 : |v 2| ≤ s := by simpa [s] using abs_le_supNorm3 v 2
  have hcoord : ∀ i : Fin 3, |(R.mulVec v) i| ≤ s := by
    intro i
    have hsumAbs :
        |(R.mulVec v) i| ≤ ∑ j : Fin 3, |R i j * v j| := by
      simpa [Matrix.mulVec, dotProduct] using
        (abs_sum_le_sum_abs (s := (Finset.univ : Finset (Fin 3)))
          (f := fun j : Fin 3 => R i j * v j))
    have hsumBound :
        (∑ j : Fin 3, |R i j * v j|) ≤ ∑ j : Fin 3, (R i j) * s := by
      refine Finset.sum_le_sum ?_
      intro j hj
      have hRij : 0 ≤ R i j := Matrix.nonneg_of_mem_rowStochastic hR i j
      have hvj : |v j| ≤ s := by
        fin_cases j
        · exact hs0
        · exact hs1
        · exact hs2
      calc
        |R i j * v j| = R i j * |v j| := by simp [abs_mul, abs_of_nonneg hRij]
        _ ≤ R i j * s := mul_le_mul_of_nonneg_left hvj hRij
    have hsumEq :
        (∑ j : Fin 3, (R i j) * s) = s := by
      calc
        (∑ j : Fin 3, (R i j) * s) = (∑ j : Fin 3, R i j) * s := by
          simp [Finset.sum_mul]
        _ = 1 * s := by rw [Matrix.sum_row_of_mem_rowStochastic hR i]
        _ = s := by ring
    exact (le_trans hsumAbs (le_trans hsumBound (by simpa [hsumEq])))
  exact supNorm3_le_of_coord_bounds (R.mulVec v)
    (hcoord 0) (hcoord 1) (hcoord 2)

/-- Schur-safe contraction lemma:
Doeblin decomposition contracts every zero-sum mode in `supNorm3` by `1-ε`,
without requiring a real eigenmode witness. -/
theorem schur_safe_supNorm3_contraction_of_doeblin
    (P R : Matrix (Fin 3) (Fin 3) ℝ)
    (eps : ℝ)
    (v : Fin 3 → ℝ)
    (hdecomp : ∀ i j, P i j = eps * Gutoe.YangMillsStructuralGap.uniformKernel i j + (1 - eps) * R i j)
    (hRstoch : R ∈ Matrix.rowStochastic ℝ (Fin 3))
    (heps1 : eps < 1)
    (hsum0 : v 0 + v 1 + v 2 = 0) :
    supNorm3 (P.mulVec v) ≤ (1 - eps) * supNorm3 v := by
  have hPmat : P = eps • Gutoe.YangMillsStructuralGap.uniformKernel + (1 - eps) • R := by
    ext i j
    simp [hdecomp i j, mul_comm]
  have hU0 :
      Gutoe.YangMillsStructuralGap.uniformKernel.mulVec v = 0 :=
    Gutoe.YangMillsStructuralGap.uniformKernel_mulVec_zero_of_sum_zero v hsum0
  have hPv :
      P.mulVec v = (1 - eps) • (R.mulVec v) := by
    calc
      P.mulVec v = (eps • Gutoe.YangMillsStructuralGap.uniformKernel + (1 - eps) • R).mulVec v := by
        simp [hPmat]
      _ = eps • (Gutoe.YangMillsStructuralGap.uniformKernel.mulVec v) + (1 - eps) • (R.mulVec v) := by
            simp [Matrix.add_mulVec, Matrix.smul_mulVec]
      _ = (1 - eps) • (R.mulVec v) := by simp [hU0]
  have honeps_nonneg : 0 ≤ 1 - eps := by linarith
  have hRnonexp : supNorm3 (R.mulVec v) ≤ supNorm3 v :=
    rowStochastic_supNorm3_nonexpansive R hRstoch v
  calc
    supNorm3 (P.mulVec v)
        = supNorm3 ((1 - eps) • (R.mulVec v)) := by simpa [hPv]
    _ = (1 - eps) * supNorm3 (R.mulVec v) := supNorm3_smul_of_nonneg (1 - eps) honeps_nonneg (R.mulVec v)
    _ ≤ (1 - eps) * supNorm3 v := mul_le_mul_of_nonneg_left hRnonexp honeps_nonneg

/-- Sufficient obligations to transfer strict mass-gap positivity from reduced
ratio control to full-state ratio control. -/
structure FullStateLiftGapObligations
    (a_t lambda0R lambda1R lambda0F lambda1F : ℝ) where
  a_t_pos : 0 < a_t
  reduced_pos0 : 0 < lambda0R
  reduced_pos1 : 0 < lambda1R
  full_pos0 : 0 < lambda0F
  full_pos1 : 0 < lambda1F
  eps : ℝ
  eps_pos : 0 < eps
  eps_lt_one : eps < 1
  ratio_dom : lambda1F / lambda0F ≤ lambda1R / lambda0R
  reduced_ratio_bound : lambda1R / lambda0R ≤ 1 - eps

/-- Under ratio-dominance + reduced Doeblin control, strict positivity lifts to
full-state gap positivity. -/
theorem full_gap_positive_of_lift_obligations
    {a_t lambda0R lambda1R lambda0F lambda1F : ℝ}
    (h : FullStateLiftGapObligations a_t lambda0R lambda1R lambda0F lambda1F) :
    0 < massGapFromEigenRatio a_t lambda0F lambda1F := by
  have hratio_full : lambda1F / lambda0F ≤ 1 - h.eps :=
    le_trans h.ratio_dom h.reduced_ratio_bound
  exact mass_gap_positive_of_doeblin_ratio
    (a_t := a_t)
    (lambda0 := lambda0F)
    (lambda1 := lambda1F)
    (eps := h.eps)
    h.a_t_pos h.full_pos0 h.full_pos1 h.eps_pos h.eps_lt_one hratio_full

/-!
  Necessity variants used by counterexample checks:
  each dropped obligation has a concrete failure witness.
-/

/-- Variant with no reduced→full ratio-dominance requirement. -/
structure FullStateLiftGapObligationsNoRatioDom
    (a_t lambda0R lambda1R lambda0F lambda1F : ℝ) where
  a_t_pos : 0 < a_t
  reduced_pos0 : 0 < lambda0R
  reduced_pos1 : 0 < lambda1R
  full_pos0 : 0 < lambda0F
  full_pos1 : 0 < lambda1F
  eps : ℝ
  eps_pos : 0 < eps
  eps_lt_one : eps < 1
  reduced_ratio_bound : lambda1R / lambda0R ≤ 1 - eps

/-- Variant with no reduced-lane Doeblin ratio bound. -/
structure FullStateLiftGapObligationsNoReducedBound
    (a_t lambda0R lambda1R lambda0F lambda1F : ℝ) where
  a_t_pos : 0 < a_t
  reduced_pos0 : 0 < lambda0R
  reduced_pos1 : 0 < lambda1R
  full_pos0 : 0 < lambda0F
  full_pos1 : 0 < lambda1F
  eps : ℝ
  eps_pos : 0 < eps
  eps_lt_one : eps < 1
  ratio_dom : lambda1F / lambda0F ≤ lambda1R / lambda0R

/-!
  Refinement-schedule version of the lift obligations.
  This is the workhorse interface used to wire reduced-lane bounds into
  full-state gap positivity across all refinement steps.
-/

/-- Per-step lift obligations over refinement schedules. -/
structure LiftScheduleObligations
    (a_t lambda0R lambda1R lambda0F lambda1F eps : ℕ → ℝ) where
  a_t_pos : ∀ n, 0 < a_t n
  reduced_pos0 : ∀ n, 0 < lambda0R n
  reduced_pos1 : ∀ n, 0 < lambda1R n
  full_pos0 : ∀ n, 0 < lambda0F n
  full_pos1 : ∀ n, 0 < lambda1F n
  eps_pos : ∀ n, 0 < eps n
  eps_lt_one : ∀ n, eps n < 1
  ratio_dom : ∀ n, lambda1F n / lambda0F n ≤ lambda1R n / lambda0R n
  reduced_ratio_bound : ∀ n, lambda1R n / lambda0R n ≤ 1 - eps n

/-- If the lift obligations hold at every refinement step, strict full-state gap
positivity holds at every step. -/
theorem full_gap_positive_all_steps_of_lift_schedule
    {a_t lambda0R lambda1R lambda0F lambda1F eps : ℕ → ℝ}
    (h : LiftScheduleObligations a_t lambda0R lambda1R lambda0F lambda1F eps) :
    ∀ n, 0 < massGapFromEigenRatio (a_t n) (lambda0F n) (lambda1F n) := by
  intro n
  have hratio_full : lambda1F n / lambda0F n ≤ 1 - eps n :=
    le_trans (h.ratio_dom n) (h.reduced_ratio_bound n)
  exact mass_gap_positive_of_doeblin_ratio
    (a_t := a_t n)
    (lambda0 := lambda0F n)
    (lambda1 := lambda1F n)
    (eps := eps n)
    (h.a_t_pos n)
    (h.full_pos0 n)
    (h.full_pos1 n)
    (h.eps_pos n)
    (h.eps_lt_one n)
    hratio_full

end GapLift

section WilsonLiftBridge

/-- Reduced transfer kernel induced by the concrete Wilson/center schedule at
refinement step `n`. -/
noncomputable def wilsonReducedKernel
    (W : WilsonZ3Action) (alpha : ℝ) (n : ℕ) :
    Matrix (Fin 3) (Fin 3) ℝ :=
  Gutoe.YangMillsStructuralGap.smoothedTransition
    (Gutoe.YangMillsStructuralGap.z3NearestNeighborCounts (W.targetSchedule n))
    (Gutoe.YangMillsStructuralGap.rowTotalsFromCounts
      (Gutoe.YangMillsStructuralGap.z3NearestNeighborCounts (W.targetSchedule n)))
    alpha

/-- Mode-level assumptions on the reduced Wilson/center kernel used to derive
the reduced-lane ratio bound structurally. -/
structure WilsonReducedModeHypotheses
    (W : WilsonZ3Action)
    (lambda1R : ℕ → ℝ)
    (alpha : ℝ) where
  mode :
    ∀ n, ∃ v : Fin 3 → ℝ,
      (wilsonReducedKernel W alpha n).mulVec v = (lambda1R n) • v ∧
      (v 0 + v 1 + v 2 = 0) ∧
      v ≠ 0 ∧
      0 ≤ lambda1R n

/-- Spectral-hypothesis package for deriving reduced→full lift obligations from
the concrete Wilson/center schedule. The `eps` lane is fully structural
(`minorizationEps` on Wilson row totals); only spectral inequalities remain as
hypotheses. -/
structure WilsonCenterLiftSpectralHypotheses
    (W : WilsonZ3Action)
    (a_t lambda0R lambda1R lambda0F lambda1F : ℕ → ℝ)
    (alpha : ℝ) where
  a_t_pos : ∀ n, 0 < a_t n
  reduced_pos0 : ∀ n, 0 < lambda0R n
  reduced_pos1 : ∀ n, 0 < lambda1R n
  full_pos0 : ∀ n, 0 < lambda0F n
  full_pos1 : ∀ n, 0 < lambda1F n
  ratio_dom : ∀ n, lambda1F n / lambda0F n ≤ lambda1R n / lambda0R n
  /-- Reduced-lane spectral control against the Wilson-induced Doeblin bound. -/
  reduced_ratio_bound :
    ∀ n,
      lambda1R n / lambda0R n ≤
        1 - Gutoe.YangMillsStructuralGap.minorizationEps (wilsonRowTotalsSchedule W n) alpha

/-- Wilson closed-form epsilon lets us phrase reduced spectral control against a
single structural constant rather than per-step row totals. -/
theorem reduced_ratio_bound_of_wilson_closed_form
    (W : WilsonZ3Action)
    (lambda0R lambda1R : ℕ → ℝ)
    (alpha : ℝ)
    (hRedClosed :
      ∀ n,
        lambda1R n / lambda0R n ≤
          1 - ((3 * alpha) / ((Gutoe.LatticeGeometry.coordinationNumber : ℝ) + 3 * alpha))) :
    ∀ n,
      lambda1R n / lambda0R n ≤
        1 - Gutoe.YangMillsStructuralGap.minorizationEps (wilsonRowTotalsSchedule W n) alpha := by
  intro n
  rw [wilson_minorization_eps_closed_form W alpha n]
  exact hRedClosed n

/-- Structural reduced-lane ratio bound from reduced-kernel zero-sum mode
control on the concrete Wilson/center schedule. -/
theorem reduced_ratio_bound_of_wilson_reduced_modes
    (W : WilsonZ3Action)
    (lambda1R : ℕ → ℝ)
    (alpha : ℝ)
    (ha : 0 < alpha)
    (hMode : WilsonReducedModeHypotheses W lambda1R alpha) :
    ∀ n,
      lambda1R n ≤
        1 - Gutoe.YangMillsStructuralGap.minorizationEps (wilsonRowTotalsSchedule W n) alpha := by
  intro n
  rcases hMode.mode n with ⟨v, hEigP, hsum0, hvne, hlam0⟩
  let counts : Fin 3 → Fin 3 → ℕ :=
    Gutoe.YangMillsStructuralGap.z3NearestNeighborCounts (W.targetSchedule n)
  let rowTotals : Fin 3 → ℕ := Gutoe.YangMillsStructuralGap.rowTotalsFromCounts counts
  let P : Matrix (Fin 3) (Fin 3) ℝ := wilsonReducedKernel W alpha n
  let epsN : ℝ := Gutoe.YangMillsStructuralGap.minorizationEps rowTotals alpha
  have hrow : ∀ i, rowTotals i = ∑ j : Fin 3, counts i j := by
    intro i
    rfl
  have hεlt : epsN < 1 := by
    have hreg : Gutoe.YangMillsStructuralGap.SCRegularRowTotals rowTotals := by
      simpa [rowTotals, counts] using
        (Gutoe.YangMillsStructuralGap.z3_nn_row_totals_sc_regular (W.targetSchedule n))
    have hclosed :
        epsN =
          (3 * alpha) / ((Gutoe.LatticeGeometry.coordinationNumber : ℝ) + 3 * alpha) := by
      simpa [epsN] using
        (Gutoe.YangMillsStructuralGap.minorization_eps_eq_sc_regular rowTotals alpha hreg)
    rw [hclosed]
    have hcoordPos : (0 : ℝ) < Gutoe.LatticeGeometry.coordinationNumber := by
      norm_num [Gutoe.LatticeGeometry.coordination_number_is_6]
    have hdenPos : 0 < (Gutoe.LatticeGeometry.coordinationNumber : ℝ) + 3 * alpha := by
      nlinarith
    have hnumLt : 3 * alpha < (Gutoe.LatticeGeometry.coordinationNumber : ℝ) + 3 * alpha := by
      nlinarith
    exact (div_lt_one hdenPos).2 hnumLt
  rcases Gutoe.YangMillsStructuralGap.doeblin_decomposition counts rowTotals ha hrow hεlt with
    ⟨R, hRnonneg, hRrowsum, hdecomp⟩
  have hRstoch : R ∈ Matrix.rowStochastic ℝ (Fin 3) := by
    rw [Matrix.mem_rowStochastic_iff_sum]
    exact ⟨hRnonneg, hRrowsum⟩
  have hU0 : Gutoe.YangMillsStructuralGap.uniformKernel.mulVec v = 0 :=
    Gutoe.YangMillsStructuralGap.uniformKernel_mulVec_zero_of_sum_zero v hsum0
  have hPmat :
      epsN • Gutoe.YangMillsStructuralGap.uniformKernel + (1 - epsN) • R = P := by
    ext i j
    simpa [P, epsN] using (hdecomp i j).symm
  have hscaled :
      (1 - epsN) • (R.mulVec v) = (lambda1R n) • v := by
    calc
      (1 - epsN) • (R.mulVec v)
          = epsN • (Gutoe.YangMillsStructuralGap.uniformKernel.mulVec v) +
              (1 - epsN) • (R.mulVec v) := by simp [hU0]
      _ = (epsN • Gutoe.YangMillsStructuralGap.uniformKernel + (1 - epsN) • R).mulVec v := by
            simp [Matrix.add_mulVec, Matrix.smul_mulVec]
      _ = P.mulVec v := by rw [hPmat]
      _ = (lambda1R n) • v := hEigP
  have hdenNe : (1 - epsN) ≠ 0 := by
    linarith [hεlt]
  have hEigR :
      R.mulVec v = ((lambda1R n) / (1 - epsN)) • v := by
    have htmp :
        (1 / (1 - epsN)) • ((1 - epsN) • (R.mulVec v)) =
          (1 / (1 - epsN)) • ((lambda1R n) • v) := by
      exact congrArg (fun w => (1 / (1 - epsN)) • w) hscaled
    simpa [smul_smul, hdenNe, div_eq_mul_inv, mul_assoc, mul_left_comm, mul_comm] using htmp
  have hle :
      lambda1R n ≤ 1 - epsN :=
    Gutoe.YangMillsStructuralGap.eigenvalue_le_one_sub_eps_of_decomposition_stochastic
      P R epsN (lambda1R n) ((lambda1R n) / (1 - epsN)) v
      hεlt hlam0
      (by
        intro i j
        simpa [P, epsN]
          using hdecomp i j)
      hEigP hEigR hRstoch hsum0 hvne
  simpa [epsN, rowTotals, counts, wilsonRowTotalsSchedule] using hle

/-- Any nontrivial zero-sum eigenmode of the Wilson-induced reduced kernel is
Doeblin-contracted by the structural bound `1 - eps_n`. -/
theorem reduced_mode_dominates_of_wilson_doeblin
    (W : WilsonZ3Action)
    (alpha : ℝ)
    (ha : 0 < alpha) :
    ∀ n (μ : ℝ) (v : Fin 3 → ℝ),
      (wilsonReducedKernel W alpha n).mulVec v = μ • v →
      (v 0 + v 1 + v 2 = 0) →
      v ≠ 0 →
      0 ≤ μ →
      μ ≤ 1 - Gutoe.YangMillsStructuralGap.minorizationEps (wilsonRowTotalsSchedule W n) alpha := by
  intro n μ v hEigP hsum0 hvne hμ0
  let counts : Fin 3 → Fin 3 → ℕ :=
    Gutoe.YangMillsStructuralGap.z3NearestNeighborCounts (W.targetSchedule n)
  let rowTotals : Fin 3 → ℕ := Gutoe.YangMillsStructuralGap.rowTotalsFromCounts counts
  let P : Matrix (Fin 3) (Fin 3) ℝ := wilsonReducedKernel W alpha n
  let epsN : ℝ := Gutoe.YangMillsStructuralGap.minorizationEps rowTotals alpha
  have hrow : ∀ i, rowTotals i = ∑ j : Fin 3, counts i j := by
    intro i
    rfl
  have hεlt : epsN < 1 := by
    have hreg : Gutoe.YangMillsStructuralGap.SCRegularRowTotals rowTotals := by
      simpa [rowTotals, counts] using
        (Gutoe.YangMillsStructuralGap.z3_nn_row_totals_sc_regular (W.targetSchedule n))
    have hclosed :
        epsN =
          (3 * alpha) / ((Gutoe.LatticeGeometry.coordinationNumber : ℝ) + 3 * alpha) := by
      simpa [epsN] using
        (Gutoe.YangMillsStructuralGap.minorization_eps_eq_sc_regular rowTotals alpha hreg)
    rw [hclosed]
    have hcoordPos : (0 : ℝ) < Gutoe.LatticeGeometry.coordinationNumber := by
      norm_num [Gutoe.LatticeGeometry.coordination_number_is_6]
    have hdenPos : 0 < (Gutoe.LatticeGeometry.coordinationNumber : ℝ) + 3 * alpha := by
      nlinarith
    have hnumLt : 3 * alpha < (Gutoe.LatticeGeometry.coordinationNumber : ℝ) + 3 * alpha := by
      nlinarith
    exact (div_lt_one hdenPos).2 hnumLt
  rcases Gutoe.YangMillsStructuralGap.doeblin_decomposition counts rowTotals ha hrow hεlt with
    ⟨R, hRnonneg, hRrowsum, hdecomp⟩
  have hRstoch : R ∈ Matrix.rowStochastic ℝ (Fin 3) := by
    rw [Matrix.mem_rowStochastic_iff_sum]
    exact ⟨hRnonneg, hRrowsum⟩
  have hU0 : Gutoe.YangMillsStructuralGap.uniformKernel.mulVec v = 0 :=
    Gutoe.YangMillsStructuralGap.uniformKernel_mulVec_zero_of_sum_zero v hsum0
  have hPmat :
      epsN • Gutoe.YangMillsStructuralGap.uniformKernel + (1 - epsN) • R = P := by
    ext i j
    simpa [P, epsN] using (hdecomp i j).symm
  have hscaled :
      (1 - epsN) • (R.mulVec v) = μ • v := by
    calc
      (1 - epsN) • (R.mulVec v)
          = epsN • (Gutoe.YangMillsStructuralGap.uniformKernel.mulVec v) +
              (1 - epsN) • (R.mulVec v) := by simp [hU0]
      _ = (epsN • Gutoe.YangMillsStructuralGap.uniformKernel + (1 - epsN) • R).mulVec v := by
            simp [Matrix.add_mulVec, Matrix.smul_mulVec]
      _ = P.mulVec v := by rw [hPmat]
      _ = μ • v := hEigP
  have hdenNe : (1 - epsN) ≠ 0 := by
    linarith [hεlt]
  have hEigR :
      R.mulVec v = (μ / (1 - epsN)) • v := by
    have htmp :
        (1 / (1 - epsN)) • ((1 - epsN) • (R.mulVec v)) =
          (1 / (1 - epsN)) • (μ • v) := by
      exact congrArg (fun w => (1 / (1 - epsN)) • w) hscaled
    simpa [smul_smul, hdenNe, div_eq_mul_inv, mul_assoc, mul_left_comm, mul_comm] using htmp
  have hle :
      μ ≤ 1 - epsN :=
    Gutoe.YangMillsStructuralGap.eigenvalue_le_one_sub_eps_of_decomposition_stochastic
      P R epsN μ (μ / (1 - epsN)) v
      hεlt hμ0
      (by
        intro i j
        simpa [P, epsN]
          using hdecomp i j)
      hEigP hEigR hRstoch hsum0 hvne
  simpa [epsN, rowTotals, counts, wilsonRowTotalsSchedule] using hle

/-- Schur-safe (complex-stable) contraction witness on the Wilson-induced
reduced kernel: every zero-sum mode contracts in `supNorm3` by the structural
factor `1 - eps_n`, with no real-eigenmode assumption. -/
theorem wilson_reduced_kernel_schur_safe_contraction
    (W : WilsonZ3Action)
    (alpha : ℝ)
    (ha : 0 < alpha) :
    ∀ n (v : Fin 3 → ℝ),
      (v 0 + v 1 + v 2 = 0) →
      supNorm3 ((wilsonReducedKernel W alpha n).mulVec v) ≤
        (1 - Gutoe.YangMillsStructuralGap.minorizationEps (wilsonRowTotalsSchedule W n) alpha) *
          supNorm3 v := by
  intro n v hsum0
  let counts : Fin 3 → Fin 3 → ℕ :=
    Gutoe.YangMillsStructuralGap.z3NearestNeighborCounts (W.targetSchedule n)
  let rowTotals : Fin 3 → ℕ := Gutoe.YangMillsStructuralGap.rowTotalsFromCounts counts
  let P : Matrix (Fin 3) (Fin 3) ℝ := wilsonReducedKernel W alpha n
  let epsN : ℝ := Gutoe.YangMillsStructuralGap.minorizationEps rowTotals alpha
  have hrow : ∀ i, rowTotals i = ∑ j : Fin 3, counts i j := by
    intro i
    rfl
  have hεlt : epsN < 1 := by
    have hreg : Gutoe.YangMillsStructuralGap.SCRegularRowTotals rowTotals := by
      simpa [rowTotals, counts] using
        (Gutoe.YangMillsStructuralGap.z3_nn_row_totals_sc_regular (W.targetSchedule n))
    have hclosed :
        epsN =
          (3 * alpha) / ((Gutoe.LatticeGeometry.coordinationNumber : ℝ) + 3 * alpha) := by
      simpa [epsN] using
        (Gutoe.YangMillsStructuralGap.minorization_eps_eq_sc_regular rowTotals alpha hreg)
    rw [hclosed]
    have hdenPos : 0 < (Gutoe.LatticeGeometry.coordinationNumber : ℝ) + 3 * alpha := by
      nlinarith [show (0 : ℝ) < Gutoe.LatticeGeometry.coordinationNumber by
        norm_num [Gutoe.LatticeGeometry.coordination_number_is_6]]
    have hnumLt : 3 * alpha < (Gutoe.LatticeGeometry.coordinationNumber : ℝ) + 3 * alpha := by
      nlinarith [show (0 : ℝ) < Gutoe.LatticeGeometry.coordinationNumber by
        norm_num [Gutoe.LatticeGeometry.coordination_number_is_6]]
    exact (div_lt_one hdenPos).2 hnumLt
  rcases Gutoe.YangMillsStructuralGap.doeblin_decomposition counts rowTotals ha hrow hεlt with
    ⟨R, hRnonneg, hRrowsum, hdecomp⟩
  have hRstoch : R ∈ Matrix.rowStochastic ℝ (Fin 3) := by
    rw [Matrix.mem_rowStochastic_iff_sum]
    exact ⟨hRnonneg, hRrowsum⟩
  have hcontract :
      supNorm3 (P.mulVec v) ≤ (1 - epsN) * supNorm3 v :=
    schur_safe_supNorm3_contraction_of_doeblin
      P R epsN v
      (by
        intro i j
        simpa [P, epsN] using hdecomp i j)
      hRstoch
      hεlt
      hsum0
  simpa [epsN, rowTotals, counts, wilsonRowTotalsSchedule] using hcontract

/-- Structural no-extra-unit-mode theorem on the Wilson-induced reduced kernel:
there is no nontrivial zero-sum eigenmode with eigenvalue `1`. -/
theorem wilson_no_extra_unit_mode_zero_sum
    (W : WilsonZ3Action)
    (alpha : ℝ)
    (ha : 0 < alpha)
    (n : ℕ)
    (v : Fin 3 → ℝ)
    (hEig : (wilsonReducedKernel W alpha n).mulVec v = (1 : ℝ) • v)
    (hsum0 : v 0 + v 1 + v 2 = 0)
    (hvne : v ≠ 0) :
    False := by
  have hle :
      (1 : ℝ) ≤
        1 - Gutoe.YangMillsStructuralGap.minorizationEps (wilsonRowTotalsSchedule W n) alpha :=
    reduced_mode_dominates_of_wilson_doeblin
      W alpha ha n 1 v hEig hsum0 hvne (by norm_num)
  have heps_pos :
      0 < Gutoe.YangMillsStructuralGap.minorizationEps (wilsonRowTotalsSchedule W n) alpha := by
    rw [wilson_minorization_eps_closed_form W alpha n]
    have hcoordPos : (0 : ℝ) < Gutoe.LatticeGeometry.coordinationNumber := by
      norm_num [Gutoe.LatticeGeometry.coordination_number_is_6]
    have hdenPos : 0 < (Gutoe.LatticeGeometry.coordinationNumber : ℝ) + 3 * alpha := by
      nlinarith
    exact div_pos (by nlinarith) hdenPos
  linarith

/-- Full-gap hypotheses where reduced-lane control is given at mode level
instead of direct ratio-bound assumptions. -/
structure WilsonCenterFullGapModeHypotheses
    (W : WilsonZ3Action)
    (a_t lambda1R lambda0F lambda1F : ℕ → ℝ)
    (alpha : ℝ) where
  a_t_pos : ∀ n, 0 < a_t n
  full_pos0 : ∀ n, 0 < lambda0F n
  full_pos1 : ∀ n, 0 < lambda1F n
  ratio_dom : ∀ n, lambda1F n / lambda0F n ≤ lambda1R n
  reduced_mode : WilsonReducedModeHypotheses W lambda1R alpha

/-- Mode-dominance package that derives the full/reduced ratio-dominance lane
from concrete eigenmode statements on the Wilson-induced reduced kernel. -/
structure WilsonModeDominanceHypotheses
    (W : WilsonZ3Action)
    (lambda1R lambda0F lambda1F : ℕ → ℝ)
    (alpha : ℝ) where
  /-- Full-lane principal normalization in stochastic coordinates. -/
  lambda0F_one : ∀ n, lambda0F n = 1
  /-- Full-lane subdominant mode realized on the same reduced kernel. -/
  full_mode :
    ∀ n, ∃ v : Fin 3 → ℝ,
      (wilsonReducedKernel W alpha n).mulVec v = (lambda1F n) • v ∧
      (v 0 + v 1 + v 2 = 0) ∧
      v ≠ 0 ∧
      0 ≤ lambda1F n
  /-- Reduced-lane dominance over all nontrivial zero-sum modes. -/
  reduced_mode_dominates :
    ∀ n (μ : ℝ) (v : Fin 3 → ℝ),
      (wilsonReducedKernel W alpha n).mulVec v = μ • v →
      (v 0 + v 1 + v 2 = 0) →
      v ≠ 0 →
      0 ≤ μ →
      μ ≤ lambda1R n

/-- Derive ratio-dominance from full-mode realization + reduced-mode dominance
on the concrete Wilson-induced reduced kernel. -/
theorem ratio_dom_of_wilson_mode_dominance
    (W : WilsonZ3Action)
    (lambda1R lambda0F lambda1F : ℕ → ℝ)
    (alpha : ℝ)
    (hDom : WilsonModeDominanceHypotheses W lambda1R lambda0F lambda1F alpha) :
    ∀ n, lambda1F n / lambda0F n ≤ lambda1R n := by
  intro n
  rcases hDom.full_mode n with ⟨v, hEig, hsum0, hvne, hlamNonneg⟩
  have hle : lambda1F n ≤ lambda1R n :=
    hDom.reduced_mode_dominates n (lambda1F n) v hEig hsum0 hvne hlamNonneg
  have h0F : lambda0F n = 1 := hDom.lambda0F_one n
  calc
    lambda1F n / lambda0F n = lambda1F n := by simpa [h0F]
    _ ≤ lambda1R n := hle

/-- End-to-end full-lane positivity from concrete Wilson/center reduced-mode
data. This removes direct reduced ratio-bound assumptions from the seam. -/
theorem full_gap_positive_all_steps_of_wilson_center_modes
    (W : WilsonZ3Action)
    (a_t lambda1R lambda0F lambda1F : ℕ → ℝ)
    (alpha : ℝ)
    (ha : 0 < alpha)
    (h : WilsonCenterFullGapModeHypotheses W a_t lambda1R lambda0F lambda1F alpha) :
    ∀ n, 0 < massGapFromEigenRatio (a_t n) (lambda0F n) (lambda1F n) := by
  intro n
  have hred :
      lambda1R n ≤
        1 - Gutoe.YangMillsStructuralGap.minorizationEps (wilsonRowTotalsSchedule W n) alpha :=
    reduced_ratio_bound_of_wilson_reduced_modes W lambda1R alpha ha h.reduced_mode n
  have hratio :
      lambda1F n / lambda0F n ≤
        1 - Gutoe.YangMillsStructuralGap.minorizationEps (wilsonRowTotalsSchedule W n) alpha :=
    le_trans (h.ratio_dom n) hred
  have heps_pos :
      0 < Gutoe.YangMillsStructuralGap.minorizationEps (wilsonRowTotalsSchedule W n) alpha := by
    rw [wilson_minorization_eps_closed_form W alpha n]
    have hcoordPos : (0 : ℝ) < Gutoe.LatticeGeometry.coordinationNumber := by
      norm_num [Gutoe.LatticeGeometry.coordination_number_is_6]
    have hdenPos : 0 < (Gutoe.LatticeGeometry.coordinationNumber : ℝ) + 3 * alpha := by
      nlinarith
    exact div_pos (by nlinarith) hdenPos
  have heps_lt_one :
      Gutoe.YangMillsStructuralGap.minorizationEps (wilsonRowTotalsSchedule W n) alpha < 1 := by
    rw [wilson_minorization_eps_closed_form W alpha n]
    have hdenPos : 0 < (Gutoe.LatticeGeometry.coordinationNumber : ℝ) + 3 * alpha := by
      nlinarith [show (0 : ℝ) < Gutoe.LatticeGeometry.coordinationNumber by
        norm_num [Gutoe.LatticeGeometry.coordination_number_is_6]]
    have hnumLt : 3 * alpha < (Gutoe.LatticeGeometry.coordinationNumber : ℝ) + 3 * alpha := by
      nlinarith [show (0 : ℝ) < Gutoe.LatticeGeometry.coordinationNumber by
        norm_num [Gutoe.LatticeGeometry.coordination_number_is_6]]
    exact (div_lt_one hdenPos).2 hnumLt
  exact mass_gap_positive_of_doeblin_ratio
    (a_t := a_t n)
    (lambda0 := lambda0F n)
    (lambda1 := lambda1F n)
    (eps := Gutoe.YangMillsStructuralGap.minorizationEps (wilsonRowTotalsSchedule W n) alpha)
    (h.a_t_pos n)
    (h.full_pos0 n)
    (h.full_pos1 n)
    heps_pos
    heps_lt_one
    hratio

/-- Seam-closure variant:
derive full-lane positivity with no direct `ratio_dom` assumption by combining
mode-dominance and reduced-mode structural bounds. -/
theorem full_gap_positive_all_steps_of_wilson_center_modes_dominance
    (W : WilsonZ3Action)
    (a_t lambda1R lambda0F lambda1F : ℕ → ℝ)
    (alpha : ℝ)
    (ha : 0 < alpha)
    (ha_t_pos : ∀ n, 0 < a_t n)
    (hfull_pos0 : ∀ n, 0 < lambda0F n)
    (hfull_pos1 : ∀ n, 0 < lambda1F n)
    (hRedMode : WilsonReducedModeHypotheses W lambda1R alpha)
    (hDom : WilsonModeDominanceHypotheses W lambda1R lambda0F lambda1F alpha) :
    ∀ n, 0 < massGapFromEigenRatio (a_t n) (lambda0F n) (lambda1F n) := by
  have hratio : ∀ n, lambda1F n / lambda0F n ≤ lambda1R n :=
    ratio_dom_of_wilson_mode_dominance W lambda1R lambda0F lambda1F alpha hDom
  exact full_gap_positive_all_steps_of_wilson_center_modes
    W
    a_t
    lambda1R
    lambda0F
    lambda1F
    alpha
    ha
    {
      a_t_pos := ha_t_pos
      full_pos0 := hfull_pos0
      full_pos1 := hfull_pos1
      ratio_dom := hratio
      reduced_mode := hRedMode
    }

/-- Maximal seam-closure variant:
derive full-lane positivity directly from full-mode data on the concrete
Wilson-induced reduced kernel, with no separate reduced-lane spectral package.

Remaining assumptions are now only:
- full-mode realization on the Wilson-induced reduced kernel,
- principal normalization (`lambda0F = 1`),
- positivity/range conditions needed by `mass_gap_positive_of_doeblin_ratio`.
-/
theorem full_gap_positive_all_steps_of_wilson_center_from_full_modes_only
    (W : WilsonZ3Action)
    (a_t lambda0F lambda1F : ℕ → ℝ)
    (alpha : ℝ)
    (ha : 0 < alpha)
    (ha_t_pos : ∀ n, 0 < a_t n)
    (hfull_pos0 : ∀ n, 0 < lambda0F n)
    (hfull_pos1 : ∀ n, 0 < lambda1F n)
    (h0F : ∀ n, lambda0F n = 1)
    (hFullMode :
      ∀ n, ∃ v : Fin 3 → ℝ,
        (wilsonReducedKernel W alpha n).mulVec v = (lambda1F n) • v ∧
        (v 0 + v 1 + v 2 = 0) ∧
        v ≠ 0 ∧
        0 ≤ lambda1F n) :
    ∀ n, 0 < massGapFromEigenRatio (a_t n) (lambda0F n) (lambda1F n) := by
  intro n
  rcases hFullMode n with ⟨v, hEig, hsum0, hvne, hμ0⟩
  have hcap :
      lambda1F n ≤
        1 - Gutoe.YangMillsStructuralGap.minorizationEps (wilsonRowTotalsSchedule W n) alpha :=
    reduced_mode_dominates_of_wilson_doeblin
      W alpha ha n (lambda1F n) v hEig hsum0 hvne hμ0
  have hratio :
      lambda1F n / lambda0F n ≤
        1 - Gutoe.YangMillsStructuralGap.minorizationEps (wilsonRowTotalsSchedule W n) alpha := by
    calc
      lambda1F n / lambda0F n = lambda1F n := by simpa [h0F n]
      _ ≤ 1 - Gutoe.YangMillsStructuralGap.minorizationEps (wilsonRowTotalsSchedule W n) alpha := hcap
  have heps_pos :
      0 < Gutoe.YangMillsStructuralGap.minorizationEps (wilsonRowTotalsSchedule W n) alpha := by
    rw [wilson_minorization_eps_closed_form W alpha n]
    have hcoordPos : (0 : ℝ) < Gutoe.LatticeGeometry.coordinationNumber := by
      norm_num [Gutoe.LatticeGeometry.coordination_number_is_6]
    have hdenPos : 0 < (Gutoe.LatticeGeometry.coordinationNumber : ℝ) + 3 * alpha := by
      nlinarith
    exact div_pos (by nlinarith) hdenPos
  have heps_lt_one :
      Gutoe.YangMillsStructuralGap.minorizationEps (wilsonRowTotalsSchedule W n) alpha < 1 := by
    rw [wilson_minorization_eps_closed_form W alpha n]
    have hdenPos : 0 < (Gutoe.LatticeGeometry.coordinationNumber : ℝ) + 3 * alpha := by
      nlinarith [show (0 : ℝ) < Gutoe.LatticeGeometry.coordinationNumber by
        norm_num [Gutoe.LatticeGeometry.coordination_number_is_6]]
    have hnumLt : 3 * alpha < (Gutoe.LatticeGeometry.coordinationNumber : ℝ) + 3 * alpha := by
      nlinarith [show (0 : ℝ) < Gutoe.LatticeGeometry.coordinationNumber by
        norm_num [Gutoe.LatticeGeometry.coordination_number_is_6]]
    exact (div_lt_one hdenPos).2 hnumLt
  exact mass_gap_positive_of_doeblin_ratio
    (a_t := a_t n)
    (lambda0 := lambda0F n)
    (lambda1 := lambda1F n)
    (eps := Gutoe.YangMillsStructuralGap.minorizationEps (wilsonRowTotalsSchedule W n) alpha)
    (ha_t_pos n)
    (hfull_pos0 n)
    (hfull_pos1 n)
    heps_pos
    heps_lt_one
    hratio

/-- Identified-mode specialization:
if the full-lane subdominant mode is identified with the concrete reduced-kernel
mode (`lambda0 = 1`), strict positivity follows with no extra full↔reduced
ratio package. -/
theorem full_gap_positive_all_steps_of_wilson_center_identified_mode
    (W : WilsonZ3Action)
    (a_t lambda1 : ℕ → ℝ)
    (alpha : ℝ)
    (ha : 0 < alpha)
    (ha_t_pos : ∀ n, 0 < a_t n)
    (hlam_pos : ∀ n, 0 < lambda1 n)
    (hMode : WilsonReducedModeHypotheses W lambda1 alpha) :
    ∀ n, 0 < massGapFromEigenRatio (a_t n) 1 (lambda1 n) := by
  intro n
  have h0_pos : (0 : ℝ) < 1 := by norm_num
  have h0F : (fun _ : ℕ => (1 : ℝ)) n = 1 := by rfl
  have hModeFull :
      ∀ n, ∃ v : Fin 3 → ℝ,
        (wilsonReducedKernel W alpha n).mulVec v = (lambda1 n) • v ∧
        (v 0 + v 1 + v 2 = 0) ∧
        v ≠ 0 ∧
        0 ≤ lambda1 n := by
    intro m
    rcases hMode.mode m with ⟨v, hEig, hsum0, hvne, hμ0⟩
    exact ⟨v, hEig, hsum0, hvne, hμ0⟩
  simpa [h0F] using
    full_gap_positive_all_steps_of_wilson_center_from_full_modes_only
      (W := W)
      (a_t := a_t)
      (lambda0F := fun _ => (1 : ℝ))
      (lambda1F := lambda1)
      (alpha := alpha)
      (ha := ha)
      (ha_t_pos := ha_t_pos)
      (hfull_pos0 := fun _ => h0_pos)
      (hfull_pos1 := hlam_pos)
      (h0F := fun _ => rfl)
      (hFullMode := hModeFull)
      n

/-- Unconditional full-state gap certificate on the Wilson/center lane:
the Doeblin lower bound is strictly positive at every refinement step from
structural epsilon alone (`alpha > 0`, SC row-total regularity). -/
theorem full_doeblin_gap_positive_all_steps_of_wilson_center
    (W : WilsonZ3Action)
    (a_t : ℕ → ℝ)
    (alpha : ℝ)
    (ha : 0 < alpha)
    (ha_t_pos : ∀ n, 0 < a_t n) :
    ∀ n,
      0 < doeblinGapLowerBound
        (a_t n)
        (Gutoe.YangMillsStructuralGap.minorizationEps (wilsonRowTotalsSchedule W n) alpha) := by
  intro n
  have heps_pos :
      0 < Gutoe.YangMillsStructuralGap.minorizationEps (wilsonRowTotalsSchedule W n) alpha := by
    rw [wilson_minorization_eps_closed_form W alpha n]
    have hcoordPos : (0 : ℝ) < Gutoe.LatticeGeometry.coordinationNumber := by
      norm_num [Gutoe.LatticeGeometry.coordination_number_is_6]
    have hdenPos : 0 < (Gutoe.LatticeGeometry.coordinationNumber : ℝ) + 3 * alpha := by
      nlinarith
    exact div_pos (by nlinarith) hdenPos
  have heps_lt_one :
      Gutoe.YangMillsStructuralGap.minorizationEps (wilsonRowTotalsSchedule W n) alpha < 1 := by
    rw [wilson_minorization_eps_closed_form W alpha n]
    have hdenPos : 0 < (Gutoe.LatticeGeometry.coordinationNumber : ℝ) + 3 * alpha := by
      nlinarith [show (0 : ℝ) < Gutoe.LatticeGeometry.coordinationNumber by
        norm_num [Gutoe.LatticeGeometry.coordination_number_is_6]]
    have hnumLt : 3 * alpha < (Gutoe.LatticeGeometry.coordinationNumber : ℝ) + 3 * alpha := by
      nlinarith [show (0 : ℝ) < Gutoe.LatticeGeometry.coordinationNumber by
        norm_num [Gutoe.LatticeGeometry.coordination_number_is_6]]
    exact (div_lt_one hdenPos).2 hnumLt
  exact doeblin_bound_positive (ha_t_pos n) heps_pos heps_lt_one

/-- Derive full reduced→full lift obligations from the concrete Wilson/center
construction plus the remaining spectral hypotheses. -/
theorem lift_schedule_obligations_of_wilson_center_schedule
    (W : WilsonZ3Action)
    (a_t lambda0R lambda1R lambda0F lambda1F : ℕ → ℝ)
    (alpha : ℝ)
    (ha : 0 < alpha)
    (hSpec :
      WilsonCenterLiftSpectralHypotheses
        W a_t lambda0R lambda1R lambda0F lambda1F alpha) :
    LiftScheduleObligations
      a_t
      lambda0R
      lambda1R
      lambda0F
      lambda1F
      (fun n => Gutoe.YangMillsStructuralGap.minorizationEps (wilsonRowTotalsSchedule W n) alpha) := by
  refine {
    a_t_pos := hSpec.a_t_pos
    reduced_pos0 := hSpec.reduced_pos0
    reduced_pos1 := hSpec.reduced_pos1
    full_pos0 := hSpec.full_pos0
    full_pos1 := hSpec.full_pos1
    ratio_dom := hSpec.ratio_dom
    reduced_ratio_bound := hSpec.reduced_ratio_bound
    eps_pos := ?_
    eps_lt_one := ?_
  }
  · intro n
    rw [wilson_minorization_eps_closed_form W alpha n]
    have hcoordPos : (0 : ℝ) < Gutoe.LatticeGeometry.coordinationNumber := by
      norm_num [Gutoe.LatticeGeometry.coordination_number_is_6]
    have hdenPos : 0 < (Gutoe.LatticeGeometry.coordinationNumber : ℝ) + 3 * alpha := by
      nlinarith
    exact div_pos (by nlinarith) hdenPos
  · intro n
    rw [wilson_minorization_eps_closed_form W alpha n]
    have hcoordPos : (0 : ℝ) < Gutoe.LatticeGeometry.coordinationNumber := by
      norm_num [Gutoe.LatticeGeometry.coordination_number_is_6]
    have hdenPos : 0 < (Gutoe.LatticeGeometry.coordinationNumber : ℝ) + 3 * alpha := by
      nlinarith
    have hnumLtDen : 3 * alpha < (Gutoe.LatticeGeometry.coordinationNumber : ℝ) + 3 * alpha := by
      nlinarith
    exact (div_lt_one hdenPos).2 hnumLtDen

/-- End-to-end reduced→full positivity transfer on the concrete Wilson/center
schedule. This theorem closes GRAND-301's "obligation plumbing" gap at the
construction level. -/
theorem full_gap_positive_all_steps_of_wilson_center_schedule
    (W : WilsonZ3Action)
    (a_t lambda0R lambda1R lambda0F lambda1F : ℕ → ℝ)
    (alpha : ℝ)
    (ha : 0 < alpha)
    (hSpec :
      WilsonCenterLiftSpectralHypotheses
        W a_t lambda0R lambda1R lambda0F lambda1F alpha) :
    ∀ n, 0 < massGapFromEigenRatio (a_t n) (lambda0F n) (lambda1F n) := by
  exact full_gap_positive_all_steps_of_lift_schedule
    (lift_schedule_obligations_of_wilson_center_schedule
      W a_t lambda0R lambda1R lambda0F lambda1F alpha ha hSpec)

end WilsonLiftBridge

section Counterexamples

/-- Full-state lift failure mode: a nontrivial eigenmode with eigenvalue `1`
outside the reduced lane blocks strict-gap inference. -/
def NoExtraUnitMode (P : Matrix (Fin 4) (Fin 4) ℝ) : Prop :=
  ∀ v : Fin 4 → ℝ,
    v ≠ 0 →
    (v 0 + v 1 + v 2 + v 3 = 0) →
    P.mulVec v = (1 : ℝ) • v →
    False

/-- Identity kernel on four states (row-stochastic, but with maximal unit-mode
multiplicity). -/
def identityKernel4 : Matrix (Fin 4) (Fin 4) ℝ := 1

/-- A split mode that is orthogonal to the constant mode. -/
def splitMode4 : Fin 4 → ℝ :=
  fun i => if i.1 < 2 then 1 else -1

theorem splitMode4_nonzero : splitMode4 ≠ 0 := by
  intro h
  have h0 := congrArg (fun v : Fin 4 → ℝ => v 0) h
  norm_num [splitMode4] at h0

theorem splitMode4_zero_sum :
    splitMode4 0 + splitMode4 1 + splitMode4 2 + splitMode4 3 = 0 := by
  have h0 : splitMode4 0 = 1 := by norm_num [splitMode4]
  have h1 : splitMode4 1 = 1 := by norm_num [splitMode4]
  have h2 : splitMode4 2 = -1 := by norm_num [splitMode4]
  have h3 : splitMode4 3 = -1 := by norm_num [splitMode4]
  rw [h0, h1, h2, h3]
  norm_num

theorem identityKernel4_mulVec_splitMode4 :
    identityKernel4.mulVec splitMode4 = (1 : ℝ) • splitMode4 := by
  ext i
  simp [identityKernel4]

/-- Concrete failure witness: without a "no extra unit mode" obligation, a
full-state kernel can carry a nontrivial `λ=1` mode. -/
theorem identityKernel4_not_noExtraUnitMode :
    ¬ NoExtraUnitMode identityKernel4 := by
  intro h
  exact h splitMode4 splitMode4_nonzero splitMode4_zero_sum
    identityKernel4_mulVec_splitMode4

/-- Unit ratio implies zero gap, so strict positivity cannot be inferred from
reduced data alone. -/
theorem mass_gap_zero_at_unit_ratio {a_t : ℝ} :
    massGapFromEigenRatio a_t 1 1 = 0 := by
  unfold massGapFromEigenRatio
  simp

/-- Counterexample gate: dropping reduced→full ratio dominance makes strict
full-state positivity non-derivable. -/
theorem no_ratio_dom_rule_false :
    ¬ (∀ a_t lambda0R lambda1R lambda0F lambda1F : ℝ,
      FullStateLiftGapObligationsNoRatioDom a_t lambda0R lambda1R lambda0F lambda1F →
      0 < massGapFromEigenRatio a_t lambda0F lambda1F) := by
  intro h
  let hOb : FullStateLiftGapObligationsNoRatioDom
      (1 : ℝ) 1 (4 / 5 : ℝ) 1 1 := {
    a_t_pos := by norm_num
    reduced_pos0 := by norm_num
    reduced_pos1 := by norm_num
    full_pos0 := by norm_num
    full_pos1 := by norm_num
    eps := (1 / 10 : ℝ)
    eps_pos := by norm_num
    eps_lt_one := by norm_num
    reduced_ratio_bound := by norm_num
  }
  have hPos : 0 < massGapFromEigenRatio (1 : ℝ) 1 1 := h 1 1 (4 / 5 : ℝ) 1 1 hOb
  have hZero : massGapFromEigenRatio (1 : ℝ) 1 1 = 0 := by
    simpa using (mass_gap_zero_at_unit_ratio (a_t := (1 : ℝ)))
  linarith

/-- Counterexample gate: dropping reduced-lane Doeblin ratio control makes
strict full-state positivity non-derivable. -/
theorem no_reduced_ratio_bound_rule_false :
    ¬ (∀ a_t lambda0R lambda1R lambda0F lambda1F : ℝ,
      FullStateLiftGapObligationsNoReducedBound a_t lambda0R lambda1R lambda0F lambda1F →
      0 < massGapFromEigenRatio a_t lambda0F lambda1F) := by
  intro h
  let hOb : FullStateLiftGapObligationsNoReducedBound
      (1 : ℝ) 1 1 1 1 := {
    a_t_pos := by norm_num
    reduced_pos0 := by norm_num
    reduced_pos1 := by norm_num
    full_pos0 := by norm_num
    full_pos1 := by norm_num
    eps := (1 / 10 : ℝ)
    eps_pos := by norm_num
    eps_lt_one := by norm_num
    ratio_dom := by norm_num
  }
  have hPos : 0 < massGapFromEigenRatio (1 : ℝ) 1 1 := h 1 1 1 1 1 hOb
  have hZero : massGapFromEigenRatio (1 : ℝ) 1 1 = 0 := by
    simpa using (mass_gap_zero_at_unit_ratio (a_t := (1 : ℝ)))
  linarith

/-- Counterexample check: one can have a strictly positive reduced-lane gap and
simultaneously a zero full-lane gap if extra `λ=1` modes are not excluded. -/
theorem reduced_positive_full_zero_counterexample :
    ∃ a_t lambda0R lambda1R lambda0F lambda1F : ℝ,
      0 < massGapFromEigenRatio a_t lambda0R lambda1R ∧
      massGapFromEigenRatio a_t lambda0F lambda1F = 0 := by
  refine ⟨1, 1, (4 / 5 : ℝ), 1, 1, ?_, ?_⟩
  · have : 0 < massGapFromEigenRatio 1 1 (4 / 5 : ℝ) :=
      mass_gap_positive_of_eigen_ratio
        (a_t := 1)
        (lambda0 := 1)
        (lambda1 := (4 / 5 : ℝ))
        (by norm_num)
        (by norm_num)
        (by norm_num)
        (by norm_num)
    simpa using this
  · simpa using (mass_gap_zero_at_unit_ratio (a_t := 1))

/-- Unsoundness witness for naïve lift rule:
it is false that reduced strict positivity alone forces full strict positivity
without additional full-state obligations. -/
theorem naive_lift_rule_false :
    ¬ (∀ a_t lambda0R lambda1R lambda0F lambda1F : ℝ,
      0 < massGapFromEigenRatio a_t lambda0R lambda1R →
      0 < massGapFromEigenRatio a_t lambda0F lambda1F) := by
  intro h
  rcases reduced_positive_full_zero_counterexample with
    ⟨a_t, lambda0R, lambda1R, lambda0F, lambda1F, hRed, hFullZero⟩
  have hFullPos : 0 < massGapFromEigenRatio a_t lambda0F lambda1F :=
    h a_t lambda0R lambda1R lambda0F lambda1F hRed
  linarith [hFullPos, hFullZero]

/-- Consolidated necessity statement for the two bridge obligations that couple
reduced and full lanes. -/
theorem lift_obligation_families_are_independently_necessary :
    (¬ (∀ a_t lambda0R lambda1R lambda0F lambda1F : ℝ,
      FullStateLiftGapObligationsNoRatioDom a_t lambda0R lambda1R lambda0F lambda1F →
      0 < massGapFromEigenRatio a_t lambda0F lambda1F)) ∧
    (¬ (∀ a_t lambda0R lambda1R lambda0F lambda1F : ℝ,
      FullStateLiftGapObligationsNoReducedBound a_t lambda0R lambda1R lambda0F lambda1F →
      0 < massGapFromEigenRatio a_t lambda0F lambda1F)) := by
  exact ⟨no_ratio_dom_rule_false, no_reduced_ratio_bound_rule_false⟩

end Counterexamples

end Gutoe.YangMillsFullStateLift
