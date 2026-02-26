/-
 * GUTOE — Yang-Mills Structural Gap Preconditions
 * Copyright (C) 2026  Riff Labs
 * AGPL-3.0-or-later
 *
 * Formal structural layer toward GRAND-296/297/298 Theorem A:
 *   - Laplace-smoothed transition entries are strictly positive
 *   - strictly positive transfer matrices are primitive and irreducible
 *   - Gram/symmetric transfer proxy S·Sᵀ is strictly positive when S is
 *
 * No `sorry`.
-/

import Mathlib
import Mathlib.LinearAlgebra.Matrix.Irreducible.Defs
import Gutoe.Z3Uniqueness

namespace Gutoe.YangMillsStructuralGap

open BigOperators
open Matrix
open Gutoe.Z3Uniqueness
open scoped Matrix

/-- Transfer basis dimension from Cl(1,3) magnetic orbit cardinality. -/
def transferBasisDim : ℕ := magneticTriplet.card

/-- The Z₃ transfer basis is exactly 3-dimensional. -/
theorem transfer_basis_dim_eq_three : transferBasisDim = 3 := by
  unfold transferBasisDim
  exact su2_dim

/-- Laplace-smoothed transition entry used in empirical kernel construction. -/
noncomputable def smoothEntry (count rowTotal : ℕ) (alpha : ℝ) : ℝ :=
  ((count : ℝ) + alpha) / ((rowTotal : ℝ) + 3 * alpha)

/-- For `alpha>0`, every smoothed entry is strictly positive. -/
theorem smooth_entry_pos (count rowTotal : ℕ) {alpha : ℝ} (ha : 0 < alpha) :
    0 < smoothEntry count rowTotal alpha := by
  unfold smoothEntry
  have hcount : 0 ≤ (count : ℝ) := by positivity
  have hrow : 0 ≤ (rowTotal : ℝ) := by positivity
  have hnum : 0 < (count : ℝ) + alpha := add_pos_of_nonneg_of_pos hcount ha
  have hden : 0 < (rowTotal : ℝ) + 3 * alpha := by nlinarith
  exact div_pos hnum hden

/-- Matrix-valued smoothed transition kernel (3-state basis). -/
noncomputable def smoothedTransition
    (counts : Fin 3 → Fin 3 → ℕ)
    (rowTotals : Fin 3 → ℕ)
    (alpha : ℝ) : Matrix (Fin 3) (Fin 3) ℝ :=
  fun i j => smoothEntry (counts i j) (rowTotals i) alpha

/-- Every entry of `smoothedTransition` is strictly positive for `alpha>0`. -/
theorem smoothed_transition_entry_pos
    (counts : Fin 3 → Fin 3 → ℕ)
    (rowTotals : Fin 3 → ℕ)
    {alpha : ℝ} (ha : 0 < alpha) :
    ∀ i j, 0 < smoothedTransition counts rowTotals alpha i j := by
  intro i j
  exact smooth_entry_pos (counts i j) (rowTotals i) ha

/-- If row totals match row sums of counts, smoothed rows are normalized to `1`. -/
theorem smoothed_transition_row_sum_one
    (counts : Fin 3 → Fin 3 → ℕ)
    (rowTotals : Fin 3 → ℕ)
    {alpha : ℝ} (ha : 0 < alpha)
    (hrow : ∀ i, rowTotals i = ∑ j : Fin 3, counts i j) :
    ∀ i, (∑ j : Fin 3, smoothedTransition counts rowTotals alpha i j) = 1 := by
  intro i
  have hrowNat : rowTotals i = counts i 0 + counts i 1 + counts i 2 := by
    simpa [Fin.sum_univ_three] using hrow i
  have hrowReal : (rowTotals i : ℝ) = counts i 0 + counts i 1 + counts i 2 := by
    exact_mod_cast hrowNat
  have hden : (rowTotals i : ℝ) + 3 * alpha ≠ 0 := by
    have hdenPos : 0 < (rowTotals i : ℝ) + 3 * alpha := by nlinarith
    exact ne_of_gt hdenPos
  rw [Fin.sum_univ_three]
  unfold smoothedTransition smoothEntry
  field_simp [hden]
  nlinarith [hrowReal]

/-- Maximum row count over the fixed three-state transfer basis. -/
def maxRowTotal (rowTotals : Fin 3 → ℕ) : ℕ :=
  max (rowTotals 0) (max (rowTotals 1) (rowTotals 2))

/-- Every row total is bounded by `maxRowTotal`. -/
theorem rowTotal_le_maxRowTotal
    (rowTotals : Fin 3 → ℕ) (i : Fin 3) :
    rowTotals i ≤ maxRowTotal rowTotals := by
  fin_cases i <;> simp [maxRowTotal]

/-- Global positive floor induced by Laplace smoothing. -/
noncomputable def laplaceGlobalFloor
    (rowTotals : Fin 3 → ℕ) (alpha : ℝ) : ℝ :=
  alpha / ((maxRowTotal rowTotals : ℝ) + 3 * alpha)

/-- The global Laplace floor is strictly positive for `alpha>0`. -/
theorem laplace_global_floor_pos
    (rowTotals : Fin 3 → ℕ) {alpha : ℝ} (ha : 0 < alpha) :
    0 < laplaceGlobalFloor rowTotals alpha := by
  unfold laplaceGlobalFloor
  have hden : 0 < (maxRowTotal rowTotals : ℝ) + 3 * alpha := by nlinarith
  exact div_pos ha hden

/-- Every smoothed transition entry is bounded below by the global Laplace floor. -/
theorem smoothed_transition_entry_ge_global_floor
    (counts : Fin 3 → Fin 3 → ℕ)
    (rowTotals : Fin 3 → ℕ)
    {alpha : ℝ} (ha : 0 < alpha) :
    ∀ i j, laplaceGlobalFloor rowTotals alpha ≤ smoothedTransition counts rowTotals alpha i j := by
  intro i j
  unfold laplaceGlobalFloor smoothedTransition smoothEntry
  have hrowLeNat : rowTotals i ≤ maxRowTotal rowTotals := rowTotal_le_maxRowTotal rowTotals i
  have hrowLe : (rowTotals i : ℝ) ≤ maxRowTotal rowTotals := by exact_mod_cast hrowLeNat
  have hdenIPos : 0 < (rowTotals i : ℝ) + 3 * alpha := by nlinarith
  have hnumNonneg : 0 ≤ alpha := le_of_lt ha
  have hbase :
      alpha / (maxRowTotal rowTotals + 3 * alpha) ≤ alpha / ((rowTotals i : ℝ) + 3 * alpha) := by
    have hdenLe : (rowTotals i : ℝ) + 3 * alpha ≤ maxRowTotal rowTotals + 3 * alpha := by linarith
    exact div_le_div_of_nonneg_left hnumNonneg hdenIPos hdenLe
  have hcountNonneg : 0 ≤ (counts i j : ℝ) := by positivity
  have hnumLe : alpha ≤ (counts i j : ℝ) + alpha := by linarith
  have hcount :
      alpha / ((rowTotals i : ℝ) + 3 * alpha) ≤
        ((counts i j : ℝ) + alpha) / ((rowTotals i : ℝ) + 3 * alpha) := by
    exact div_le_div_of_nonneg_right hnumLe (le_of_lt hdenIPos)
  exact le_trans hbase hcount

/-- The global floor cannot exceed `1/3` for a normalized three-state row. -/
theorem laplace_global_floor_le_one_third
    (counts : Fin 3 → Fin 3 → ℕ)
    (rowTotals : Fin 3 → ℕ)
    {alpha : ℝ} (ha : 0 < alpha)
    (hrow : ∀ i, rowTotals i = ∑ j : Fin 3, counts i j) :
    laplaceGlobalFloor rowTotals alpha ≤ (1 : ℝ) / 3 := by
  let η : ℝ := laplaceGlobalFloor rowTotals alpha
  have hsumη :
      (∑ j : Fin 3, η) ≤
        (∑ j : Fin 3, smoothedTransition counts rowTotals alpha 0 j) := by
    refine Finset.sum_le_sum ?_
    intro j hj
    simpa [η] using smoothed_transition_entry_ge_global_floor counts rowTotals ha 0 j
  have hrow0 : (∑ j : Fin 3, smoothedTransition counts rowTotals alpha 0 j) = 1 :=
    smoothed_transition_row_sum_one counts rowTotals ha hrow 0
  have h3η : 3 * η ≤ 1 := by
    have : (∑ j : Fin 3, η) ≤ 1 := by simpa [hrow0] using hsumη
    simpa [Fin.sum_univ_three, η, add_assoc, add_left_comm, add_comm] using this
  nlinarith

/-- Doeblin-style one-step minorization constant induced by Laplace floor. -/
noncomputable def minorizationEps (rowTotals : Fin 3 → ℕ) (alpha : ℝ) : ℝ :=
  3 * laplaceGlobalFloor rowTotals alpha

/-- The induced minorization constant is strictly positive and at most `1`. -/
theorem minorization_eps_range
    (counts : Fin 3 → Fin 3 → ℕ)
    (rowTotals : Fin 3 → ℕ)
    {alpha : ℝ} (ha : 0 < alpha)
    (hrow : ∀ i, rowTotals i = ∑ j : Fin 3, counts i j) :
    0 < minorizationEps rowTotals alpha ∧ minorizationEps rowTotals alpha ≤ 1 := by
  constructor
  · unfold minorizationEps
    nlinarith [laplace_global_floor_pos rowTotals ha]
  · unfold minorizationEps
    have hη : laplaceGlobalFloor rowTotals alpha ≤ (1 : ℝ) / 3 :=
      laplace_global_floor_le_one_third counts rowTotals ha hrow
    nlinarith

/-- Uniform `3×3` kernel (all entries `1/3`). -/
noncomputable def uniformKernel : Matrix (Fin 3) (Fin 3) ℝ :=
  fun _ _ => (1 : ℝ) / 3

/-- Residual kernel after removing the uniform minorization floor. -/
noncomputable def residualKernel
    (counts : Fin 3 → Fin 3 → ℕ)
    (rowTotals : Fin 3 → ℕ)
    (alpha : ℝ) : Matrix (Fin 3) (Fin 3) ℝ :=
  let ε := minorizationEps rowTotals alpha
  fun i j =>
    (smoothedTransition counts rowTotals alpha i j - ε / 3) / (1 - ε)

/-- Doeblin decomposition: `P = εU + (1-ε)R` with `R` row-stochastic/nonnegative. -/
theorem doeblin_decomposition
    (counts : Fin 3 → Fin 3 → ℕ)
    (rowTotals : Fin 3 → ℕ)
    {alpha : ℝ} (ha : 0 < alpha)
    (hrow : ∀ i, rowTotals i = ∑ j : Fin 3, counts i j)
    (hεlt : minorizationEps rowTotals alpha < 1) :
    ∃ R : Matrix (Fin 3) (Fin 3) ℝ,
      (∀ i j, 0 ≤ R i j) ∧
      (∀ i, (∑ j : Fin 3, R i j) = 1) ∧
      (∀ i j,
        smoothedTransition counts rowTotals alpha i j =
          minorizationEps rowTotals alpha * uniformKernel i j +
            (1 - minorizationEps rowTotals alpha) * R i j) := by
  let ε : ℝ := minorizationEps rowTotals alpha
  let R : Matrix (Fin 3) (Fin 3) ℝ := residualKernel counts rowTotals alpha
  refine ⟨R, ?_, ?_, ?_⟩
  · intro i j
    have hεpos : 0 < ε := (minorization_eps_range counts rowTotals ha hrow).1
    have hdenPos : 0 < 1 - ε := by linarith
    have hεdiv : ε / 3 = laplaceGlobalFloor rowTotals alpha := by
      unfold ε minorizationEps
      ring
    have hfloor : laplaceGlobalFloor rowTotals alpha ≤ smoothedTransition counts rowTotals alpha i j :=
      smoothed_transition_entry_ge_global_floor counts rowTotals ha i j
    have hnum : 0 ≤ smoothedTransition counts rowTotals alpha i j - ε / 3 := by
      rw [hεdiv]
      linarith
    unfold R residualKernel
    exact div_nonneg hnum (le_of_lt hdenPos)
  · intro i
    have hεpos : 0 < ε := (minorization_eps_range counts rowTotals ha hrow).1
    have hdenNe : 1 - ε ≠ 0 := by linarith
    have hsumP : (∑ j : Fin 3, smoothedTransition counts rowTotals alpha i j) = 1 :=
      smoothed_transition_row_sum_one counts rowTotals ha hrow i
    have hsumP3 :
        smoothedTransition counts rowTotals alpha i 0 +
            smoothedTransition counts rowTotals alpha i 1 +
            smoothedTransition counts rowTotals alpha i 2 = 1 := by
      simpa [Fin.sum_univ_three] using hsumP
    have hsumε3 : ε / 3 + ε / 3 + ε / 3 = ε := by ring
    calc
      (∑ j : Fin 3, R i j)
          = ((smoothedTransition counts rowTotals alpha i 0 - ε / 3) +
                (smoothedTransition counts rowTotals alpha i 1 - ε / 3) +
                (smoothedTransition counts rowTotals alpha i 2 - ε / 3)) / (1 - ε) := by
              unfold R residualKernel
              rw [Fin.sum_univ_three]
              simp [ε]
              ring_nf
      _ = ((smoothedTransition counts rowTotals alpha i 0 +
              smoothedTransition counts rowTotals alpha i 1 +
              smoothedTransition counts rowTotals alpha i 2) - ε) / (1 - ε) := by
              ring_nf
      _ = (1 - ε) / (1 - ε) := by simp [hsumP3]
      _ = 1 := by field_simp [hdenNe]
  · intro i j
    have hdenNe : 1 - ε ≠ 0 := by linarith
    calc
      smoothedTransition counts rowTotals alpha i j
          = ε / 3 + (smoothedTransition counts rowTotals alpha i j - ε / 3) := by ring
      _ = ε * uniformKernel i j +
            (1 - ε) * ((smoothedTransition counts rowTotals alpha i j - ε / 3) / (1 - ε)) := by
            unfold uniformKernel
            field_simp [hdenNe]
      _ = minorizationEps rowTotals alpha * uniformKernel i j +
            (1 - minorizationEps rowTotals alpha) * R i j := by
            simp [ε, R, residualKernel]

/-- Uniform kernel sends any vector to its component-wise average. -/
theorem uniformKernel_mulVec_eq_avg
    (v : Fin 3 → ℝ) :
    uniformKernel.mulVec v = fun _ => (v 0 + v 1 + v 2) / 3 := by
  ext i
  fin_cases i <;>
    simp [uniformKernel, Matrix.mulVec, dotProduct, Fin.sum_univ_three, div_eq_mul_inv] <;>
    ring

/-- Uniform kernel annihilates zero-sum modes. -/
theorem uniformKernel_mulVec_zero_of_sum_zero
    (v : Fin 3 → ℝ)
    (hsum0 : v 0 + v 1 + v 2 = 0) :
    uniformKernel.mulVec v = 0 := by
  rw [uniformKernel_mulVec_eq_avg v]
  ext i
  simp [hsum0]

/-- Any real eigenvalue of a row-stochastic `3×3` matrix is bounded by `1` in absolute value. -/
theorem abs_eigenvalue_le_one_of_rowStochastic
    (R : Matrix (Fin 3) (Fin 3) ℝ)
    (hR : R ∈ Matrix.rowStochastic ℝ (Fin 3))
    {mu : ℝ}
    (hmu : Module.End.HasEigenvalue (Matrix.toLin' R) mu) :
    |mu| ≤ 1 := by
  rcases eigenvalue_mem_ball (A := R) hmu with ⟨k, hk⟩
  have hdiag_nonneg : 0 ≤ R k k := Matrix.nonneg_of_mem_rowStochastic hR
  have hsum_row : ∑ j, R k j = 1 := Matrix.sum_row_of_mem_rowStochastic hR k
  have hsplit : R k k + ∑ j ∈ Finset.univ.erase k, R k j = 1 := by
    have hk_univ : k ∈ (Finset.univ : Finset (Fin 3)) := Finset.mem_univ k
    have hsplit' : R k k + ∑ j ∈ Finset.univ.erase k, R k j = ∑ j, R k j := by
      simpa [add_assoc, add_left_comm, add_comm] using
        (Finset.sum_erase_add (s := Finset.univ) (a := k) (f := fun j => R k j) hk_univ)
    linarith [hsplit', hsum_row]
  have hradius_eq :
      ∑ j ∈ Finset.univ.erase k, ‖R k j‖ = 1 - R k k := by
    have habs :
        ∑ j ∈ Finset.univ.erase k, ‖R k j‖ =
          ∑ j ∈ Finset.univ.erase k, R k j := by
      refine Finset.sum_congr rfl ?_
      intro j hj
      simp [Real.norm_eq_abs, abs_of_nonneg (Matrix.nonneg_of_mem_rowStochastic hR)]
    have htail : ∑ j ∈ Finset.univ.erase k, R k j = 1 - R k k := by
      linarith [hsplit]
    calc
      ∑ j ∈ Finset.univ.erase k, ‖R k j‖ = ∑ j ∈ Finset.univ.erase k, R k j := habs
      _ = 1 - R k k := htail
  have hkraw : |mu - R k k| ≤ ∑ j ∈ Finset.univ.erase k, ‖R k j‖ := by
    simpa [Metric.mem_closedBall, Real.dist_eq] using hk
  have hkabs : |mu - R k k| ≤ 1 - R k k := hkraw.trans_eq hradius_eq
  calc
    |mu| = |(mu - R k k) + R k k| := by ring_nf
    _ ≤ |mu - R k k| + |R k k| := abs_add_le _ _
    _ ≤ (1 - R k k) + R k k := by
          refine add_le_add hkabs ?_
          simp [abs_of_nonneg hdiag_nonneg]
    _ = 1 := by ring

/-- Eigenvalue contraction bound from Doeblin decomposition on zero-sum modes. -/
theorem abs_eigenvalue_le_one_sub_eps_of_decomposition
    (P R : Matrix (Fin 3) (Fin 3) ℝ)
    (eps lam mu : ℝ)
    (v : Fin 3 → ℝ)
    (heps1 : eps < 1)
    (hdecomp : ∀ i j, P i j = eps * uniformKernel i j + (1 - eps) * R i j)
    (hEigP : P.mulVec v = lam • v)
    (hEigR : R.mulVec v = mu • v)
    (hmu : |mu| ≤ 1)
    (hsum0 : v 0 + v 1 + v 2 = 0)
    (hvne : v ≠ 0) :
    |lam| ≤ 1 - eps := by
  have hPmat : P = eps • uniformKernel + (1 - eps) • R := by
    ext i j
    simp [hdecomp i j, mul_comm]
  have hU0 : uniformKernel.mulVec v = 0 := uniformKernel_mulVec_zero_of_sum_zero v hsum0
  have hmode : (1 - eps) • (R.mulVec v) = lam • v := by
    calc
      (1 - eps) • (R.mulVec v)
          = eps • (uniformKernel.mulVec v) + (1 - eps) • (R.mulVec v) := by simp [hU0]
      _ = (eps • uniformKernel + (1 - eps) • R).mulVec v := by
            simp [Matrix.add_mulVec, Matrix.smul_mulVec]
      _ = P.mulVec v := by simp [hPmat]
      _ = lam • v := hEigP
  have hscaled : ((1 - eps) * mu) • v = lam • v := by
    calc
      ((1 - eps) * mu) • v = (1 - eps) • (mu • v) := by simp [smul_smul]
      _ = (1 - eps) • (R.mulVec v) := by simp [hEigR]
      _ = lam • v := hmode
  have hvcoord : ∃ i : Fin 3, v i ≠ 0 := by
    by_contra h
    apply hvne
    ext i
    have hi0 : v i = 0 := by
      by_contra hi
      exact h ⟨i, hi⟩
    exact hi0
  rcases hvcoord with ⟨i, hvi⟩
  have hcoord : ((1 - eps) * mu) * v i = lam * v i := by
    simpa using congrArg (fun w => w i) hscaled
  have hlam : (1 - eps) * mu = lam := by
    exact ((mul_eq_mul_right_iff).1 hcoord).resolve_right hvi
  have honeps_nonneg : 0 ≤ 1 - eps := by linarith
  calc
    |lam| = |(1 - eps) * mu| := by simp [hlam]
    _ = |1 - eps| * |mu| := by simp
    _ ≤ |1 - eps| * 1 := by
          exact mul_le_mul_of_nonneg_left hmu (abs_nonneg (1 - eps))
    _ = (1 - eps) * 1 := by simp [abs_of_nonneg honeps_nonneg]
    _ = 1 - eps := by ring

/-- Stochastic specialization of the Doeblin contraction bound (discharges `|mu| ≤ 1`). -/
theorem abs_eigenvalue_le_one_sub_eps_of_decomposition_stochastic
    (P R : Matrix (Fin 3) (Fin 3) ℝ)
    (eps lam mu : ℝ)
    (v : Fin 3 → ℝ)
    (heps1 : eps < 1)
    (hdecomp : ∀ i j, P i j = eps * uniformKernel i j + (1 - eps) * R i j)
    (hEigP : P.mulVec v = lam • v)
    (hEigR : R.mulVec v = mu • v)
    (hRstoch : R ∈ Matrix.rowStochastic ℝ (Fin 3))
    (hsum0 : v 0 + v 1 + v 2 = 0)
    (hvne : v ≠ 0) :
    |lam| ≤ 1 - eps := by
  have hmuEig : Module.End.HasEigenvalue (Matrix.toLin' R) mu := by
    have hvec : Module.End.HasEigenvector (Matrix.toLin' R) mu v := by
      refine ⟨(Module.End.mem_eigenspace_iff).2 ?_, hvne⟩
      simpa using hEigR
    exact Module.End.hasEigenvalue_of_hasEigenvector
      (f := Matrix.toLin' R) (μ := mu) (x := v) hvec
  have hmu : |mu| ≤ 1 := abs_eigenvalue_le_one_of_rowStochastic R hRstoch hmuEig
  exact abs_eigenvalue_le_one_sub_eps_of_decomposition
    P R eps lam mu v heps1 hdecomp hEigP hEigR hmu hsum0 hvne

/-- Positive-mode specialization: `lam ≤ 1-eps` on zero-sum modes under Doeblin decomposition. -/
theorem eigenvalue_le_one_sub_eps_of_decomposition
    (P R : Matrix (Fin 3) (Fin 3) ℝ)
    (eps lam mu : ℝ)
    (v : Fin 3 → ℝ)
    (heps1 : eps < 1)
    (hlam0 : 0 ≤ lam)
    (hdecomp : ∀ i j, P i j = eps * uniformKernel i j + (1 - eps) * R i j)
    (hEigP : P.mulVec v = lam • v)
    (hEigR : R.mulVec v = mu • v)
    (hmu : |mu| ≤ 1)
    (hsum0 : v 0 + v 1 + v 2 = 0)
    (hvne : v ≠ 0) :
    lam ≤ 1 - eps := by
  have habs :
      |lam| ≤ 1 - eps :=
    abs_eigenvalue_le_one_sub_eps_of_decomposition
      P R eps lam mu v heps1 hdecomp hEigP hEigR hmu hsum0 hvne
  simpa [abs_of_nonneg hlam0] using habs

/-- Positive-mode stochastic specialization (`lam ≤ 1-eps`) with no separate `|mu|` assumption. -/
theorem eigenvalue_le_one_sub_eps_of_decomposition_stochastic
    (P R : Matrix (Fin 3) (Fin 3) ℝ)
    (eps lam mu : ℝ)
    (v : Fin 3 → ℝ)
    (heps1 : eps < 1)
    (hlam0 : 0 ≤ lam)
    (hdecomp : ∀ i j, P i j = eps * uniformKernel i j + (1 - eps) * R i j)
    (hEigP : P.mulVec v = lam • v)
    (hEigR : R.mulVec v = mu • v)
    (hRstoch : R ∈ Matrix.rowStochastic ℝ (Fin 3))
    (hsum0 : v 0 + v 1 + v 2 = 0)
    (hvne : v ≠ 0) :
    lam ≤ 1 - eps := by
  have habs :
      |lam| ≤ 1 - eps :=
    abs_eigenvalue_le_one_sub_eps_of_decomposition_stochastic
      P R eps lam mu v heps1 hdecomp hEigP hEigR hRstoch hsum0 hvne
  simpa [abs_of_nonneg hlam0] using habs

/-- Any entrywise-positive `3×3` real matrix is primitive (`k=1`). -/
theorem isPrimitive_of_entrywise_pos
    (A : Matrix (Fin 3) (Fin 3) ℝ)
    (hpos : ∀ i j, 0 < A i j) :
    A.IsPrimitive := by
  refine Matrix.IsPrimitive.mk ?_ ?_
  · intro i j
    exact le_of_lt (hpos i j)
  · refine ⟨1, by norm_num, ?_⟩
    intro i j
    simpa using hpos i j

/-- Any entrywise-positive `3×3` real matrix is irreducible. -/
theorem isIrreducible_of_entrywise_pos
    (A : Matrix (Fin 3) (Fin 3) ℝ)
    (hpos : ∀ i j, 0 < A i j) :
    A.IsIrreducible := by
  exact Matrix.IsPrimitive.isIrreducible (isPrimitive_of_entrywise_pos A hpos)

/-- Structural theorem: Laplace-smoothed transition kernel is primitive (`alpha>0`). -/
theorem smoothed_transition_isPrimitive
    (counts : Fin 3 → Fin 3 → ℕ)
    (rowTotals : Fin 3 → ℕ)
    {alpha : ℝ} (ha : 0 < alpha) :
    (smoothedTransition counts rowTotals alpha).IsPrimitive := by
  refine isPrimitive_of_entrywise_pos _ ?_
  intro i j
  exact smoothed_transition_entry_pos counts rowTotals ha i j

/-- Structural theorem: Laplace-smoothed transition kernel is irreducible (`alpha>0`). -/
theorem smoothed_transition_isIrreducible
    (counts : Fin 3 → Fin 3 → ℕ)
    (rowTotals : Fin 3 → ℕ)
    {alpha : ℝ} (ha : 0 < alpha) :
    (smoothedTransition counts rowTotals alpha).IsIrreducible := by
  exact Matrix.IsPrimitive.isIrreducible (smoothed_transition_isPrimitive counts rowTotals ha)

/-- If `S` is entrywise positive, then every entry of `S * Sᵀ` is strictly positive. -/
theorem gram_entry_pos_of_entrywise_pos
    (S : Matrix (Fin 3) (Fin 3) ℝ)
    (hpos : ∀ i j, 0 < S i j) :
    ∀ i j, 0 < (S * S.transpose) i j := by
  intro i j
  have h0 : 0 < S i 0 * S j 0 := mul_pos (hpos i 0) (hpos j 0)
  have h1 : 0 ≤ S i 1 * S j 1 := le_of_lt (mul_pos (hpos i 1) (hpos j 1))
  have h2 : 0 ≤ S i 2 * S j 2 := le_of_lt (mul_pos (hpos i 2) (hpos j 2))
  have hsum_pos : 0 < ∑ k : Fin 3, S i k * S j k := by
    rw [Fin.sum_univ_three]
    nlinarith
  simpa [Matrix.mul_apply, Matrix.transpose_apply] using hsum_pos

/-- Gram/symmetric transfer proxy inherits primitivity from entrywise-positive `S`. -/
theorem gram_isPrimitive_of_entrywise_pos
    (S : Matrix (Fin 3) (Fin 3) ℝ)
    (hpos : ∀ i j, 0 < S i j) :
    (S * S.transpose).IsPrimitive := by
  refine isPrimitive_of_entrywise_pos _ ?_
  exact gram_entry_pos_of_entrywise_pos S hpos

/-- Gram/symmetric transfer proxy is irreducible when built from entrywise-positive `S`. -/
theorem gram_isIrreducible_of_entrywise_pos
    (S : Matrix (Fin 3) (Fin 3) ℝ)
    (hpos : ∀ i j, 0 < S i j) :
    (S * S.transpose).IsIrreducible := by
  exact Matrix.IsPrimitive.isIrreducible (gram_isPrimitive_of_entrywise_pos S hpos)

end Gutoe.YangMillsStructuralGap
