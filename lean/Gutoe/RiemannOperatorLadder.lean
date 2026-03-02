import Mathlib
import Gutoe.RiemannCore
import Gutoe.RiemannSelfAdjoint
import Gutoe.RiemannTargetFiniteLadder
import Gutoe.RiemannLimitBridge
import Gutoe.RiemannConvergenceTransfer
import Gutoe.RiemannFinalTarget
import Gutoe.RiemannFiniteXiModel
import Gutoe.RiemannHurwitzKernel
import Mathlib.LinearAlgebra.Matrix.Gershgorin

namespace Gutoe.RiemannOperatorLadder

open Gutoe.RiemannCore
open Gutoe.RiemannSelfAdjoint
open Gutoe.RiemannTargetFiniteLadder
open Gutoe.RiemannLimitBridge
open Gutoe.RiemannConvergenceTransfer
open Gutoe.RiemannFinalTarget
open Gutoe.RiemannFiniteXiModel
open Gutoe.RiemannHurwitzKernel

noncomputable section
open scoped Topology BigOperators

/-- Complex-lifted structural matrix lane for spectral statements. -/
def structuralRiemannMatrixC (n : ℕ) : Matrix (Fin n) (Fin n) ℂ :=
  fun i j => (structuralRiemannMatrix n i j : ℂ)

/-- Complex reversal/parity conjugation action on finite matrices. -/
def revParityConjugateMatrixC (n : ℕ) (A : Matrix (Fin n) (Fin n) ℂ) :
    Matrix (Fin n) (Fin n) ℂ :=
  fun i j => (finParitySignQ i : ℂ) * A (Fin.rev i) (Fin.rev j) * (finParitySignQ j : ℂ)

/-- Complex center matrix associated to the finite-level structural center. -/
def structuralCenterMatrixC (n : ℕ) : Matrix (Fin n) (Fin n) ℂ :=
  fun i j => if i = j then (structuralCenterQ n : ℂ) else 0

/-- Endomorphism induced by the structural matrix on coordinate vectors. -/
def structuralRiemannEnd (n : ℕ) : Module.End ℂ (Fin n → ℂ) :=
  Matrix.mulVecLin (structuralRiemannMatrixC n)

/-- A real ordinate is spectral for level `n` if it is an eigenvalue of the
complex structural endomorphism at that level. -/
def ordinateIsEigenvalue (n : ℕ) (t : ℝ) : Prop :=
  (t : ℂ) ∈ spectrum ℂ (structuralRiemannMatrixC n)

/-- Direct operator-defined spectrum lane at level `N`:
real ordinates whose complex lift lies in the matrix spectrum. -/
def operatorSpecSet (N : ℕ) : Set ℝ :=
  fun t => (t : ℂ) ∈ spectrum ℂ (structuralRiemannMatrixC (N + 1))

theorem ordinateIsEigenvalue_iff_mem_operatorSpecSet
    (N : ℕ) (t : ℝ) :
    ordinateIsEigenvalue (N + 1) t ↔ t ∈ operatorSpecSet N := by
  rfl

/-- The structural complex matrix is Hermitian (real symmetric data lifted to `ℂ`). -/
theorem structuralRiemannMatrixC_isHermitian (n : ℕ) :
    (structuralRiemannMatrixC n).IsHermitian := by
  ext i j
  simp [Matrix.conjTranspose, structuralRiemannMatrixC]
  exact (structuralRiemannMatrix_isSymm n).apply i j

/-- Complex-lifted matrix balance identity:
`A + T(A) = C`, where `A` is the structural matrix, `T` is reversal/parity
conjugation, and `C` is the center matrix. -/
theorem structuralRiemannMatrixC_matrix_balance (n : ℕ) :
    structuralRiemannMatrixC n +
      revParityConjugateMatrixC n (structuralRiemannMatrixC n)
      = structuralCenterMatrixC n := by
  ext i j
  by_cases hij : i = j
  · subst hij
    have hq :
        structuralRiemannMatrix n i i +
          finParitySignQ i * structuralRiemannMatrix n (Fin.rev i) (Fin.rev i) * finParitySignQ i
            = structuralCenterQ n := by
      simpa using structuralRiemannMatrix_rev_parity_balance n i i
    have hqC :
        (structuralRiemannMatrix n i i : ℂ) +
          (finParitySignQ i : ℂ) * (structuralRiemannMatrix n (Fin.rev i) (Fin.rev i) : ℂ) *
            (finParitySignQ i : ℂ)
            = (structuralCenterQ n : ℂ) := by
      exact_mod_cast hq
    simpa [structuralRiemannMatrixC, revParityConjugateMatrixC, structuralCenterMatrixC] using hqC
  · have hq :
      structuralRiemannMatrix n i j +
        finParitySignQ i * structuralRiemannMatrix n (Fin.rev i) (Fin.rev j) * finParitySignQ j
          = 0 := by
      simpa [hij] using structuralRiemannMatrix_rev_parity_balance n i j
    have hqC :
      (structuralRiemannMatrix n i j : ℂ) +
        (finParitySignQ i : ℂ) * (structuralRiemannMatrix n (Fin.rev i) (Fin.rev j) : ℂ) *
          (finParitySignQ j : ℂ)
          = 0 := by
      exact_mod_cast hq
    simpa [structuralRiemannMatrixC, revParityConjugateMatrixC, structuralCenterMatrixC, hij] using hqC

/-- Rearranged complex-lifted balance identity:
`T(A) = C - A` for the structural finite matrix. -/
theorem structuralRiemannMatrixC_revParity_eq_center_sub (n : ℕ) :
    revParityConjugateMatrixC n (structuralRiemannMatrixC n)
      = structuralCenterMatrixC n - structuralRiemannMatrixC n := by
  ext i j
  have hEntry := congrArg (fun M => M i j) (structuralRiemannMatrixC_matrix_balance n)
  exact (eq_sub_iff_add_eq).2 (by simpa [add_comm] using hEntry)

/-- Off-diagonal structural entries on the nearest-neighbor band have constant
complex norm `6/11`. -/
theorem norm_structuralRiemannMatrixC_offdiag_eq_hop
    (n : ℕ) (i j : Fin n) (hij : i ≠ j)
    (hband : i.1 + 1 = j.1 ∨ j.1 + 1 = i.1) :
    ‖structuralRiemannMatrixC n i j‖ = ((6 : ℝ) / 11) := by
  have hentry :
      structuralRiemannMatrixC n i j = (structuralHopQ : ℂ) := by
    simp [structuralRiemannMatrixC, structuralRiemannMatrix, hij, hband]
  rw [hentry]
  norm_num [structuralHopQ_eq_six_over_eleven]

/-- Off-diagonal structural entries away from the nearest-neighbor band vanish. -/
theorem norm_structuralRiemannMatrixC_offdiag_eq_zero
    (n : ℕ) (i j : Fin n) (hij : i ≠ j)
    (hband : ¬ (i.1 + 1 = j.1 ∨ j.1 + 1 = i.1)) :
    ‖structuralRiemannMatrixC n i j‖ = 0 := by
  have hentry : structuralRiemannMatrixC n i j = 0 := by
    simp [structuralRiemannMatrixC, structuralRiemannMatrix, hij, hband]
  simpa [hentry]

/-- Gershgorin row-radius bound for the structural tridiagonal matrix:
each row off-diagonal norm sum is at most `12/11 = 2*(6/11)`. -/
theorem structuralRiemannMatrixC_rowRadius_le_twelve_over_eleven
    (n : ℕ) (i : Fin n) :
    (∑ j ∈ Finset.univ.erase i, ‖structuralRiemannMatrixC n i j‖) ≤ ((12 : ℝ) / 11) := by
  classical
  let S := Finset.univ.erase i
  let Splus : Finset (Fin n) := S.filter (fun j => i.1 + 1 = j.1)
  let Sminus : Finset (Fin n) := S.filter (fun j => j.1 + 1 = i.1)
  have hsplit0 :
      (∑ j ∈ S, ‖structuralRiemannMatrixC n i j‖)
        = (∑ j ∈ S.filter (fun j => i.1 + 1 = j.1 ∨ j.1 + 1 = i.1),
              ‖structuralRiemannMatrixC n i j‖) := by
    have hsub :
        S.filter (fun j => i.1 + 1 = j.1 ∨ j.1 + 1 = i.1) ⊆ S := by
      intro j hj
      exact (Finset.mem_filter.mp hj).1
    have hsum_subset :
        (∑ j ∈ S.filter (fun j => i.1 + 1 = j.1 ∨ j.1 + 1 = i.1),
            ‖structuralRiemannMatrixC n i j‖)
          = (∑ j ∈ S, ‖structuralRiemannMatrixC n i j‖) := by
      refine Finset.sum_subset hsub ?_
      intro j hjS hjNot
      have hband_not : ¬ (i.1 + 1 = j.1 ∨ j.1 + 1 = i.1) := by
        intro hband
        exact hjNot (Finset.mem_filter.mpr ⟨hjS, hband⟩)
      have hij : i ≠ j := by
        have hjS' : j ∈ Finset.univ.erase i := by simpa [S] using hjS
        have hjNe : j ≠ i := (Finset.mem_erase.mp hjS').1
        exact fun h => hjNe h.symm
      have hnorm0 : ‖structuralRiemannMatrixC n i j‖ = 0 :=
        norm_structuralRiemannMatrixC_offdiag_eq_zero n i j hij hband_not
      simpa [hnorm0]
    exact hsum_subset.symm
  have hfilter_or :
      S.filter (fun j => i.1 + 1 = j.1 ∨ j.1 + 1 = i.1) = Splus ∪ Sminus := by
    simp [Splus, Sminus, S, Finset.filter_or]
  have hdisj : Disjoint Splus Sminus := by
    refine (Finset.disjoint_filter).2 ?_
    intro j hjS hjp hjm
    omega
  have hsplit :
      (∑ j ∈ S.filter (fun j => i.1 + 1 = j.1 ∨ j.1 + 1 = i.1),
          ‖structuralRiemannMatrixC n i j‖)
        = (∑ j ∈ Splus, ‖structuralRiemannMatrixC n i j‖) +
          (∑ j ∈ Sminus, ‖structuralRiemannMatrixC n i j‖) := by
    rw [hfilter_or, Finset.sum_union hdisj]
  have hcard_plus : Splus.card ≤ 1 := by
    refine (Finset.card_le_one).2 ?_
    intro a ha b hb
    have ha' : i.1 + 1 = a.1 := (Finset.mem_filter.mp ha).2
    have hb' : i.1 + 1 = b.1 := (Finset.mem_filter.mp hb).2
    apply Fin.ext
    omega
  have hcard_minus : Sminus.card ≤ 1 := by
    refine (Finset.card_le_one).2 ?_
    intro a ha b hb
    have ha' : a.1 + 1 = i.1 := (Finset.mem_filter.mp ha).2
    have hb' : b.1 + 1 = i.1 := (Finset.mem_filter.mp hb).2
    apply Fin.ext
    omega
  have hsum_plus_le : (∑ j ∈ Splus, ‖structuralRiemannMatrixC n i j‖) ≤ (6 : ℝ) / 11 := by
    have hterm :
        ∀ j ∈ Splus, ‖structuralRiemannMatrixC n i j‖ ≤ (6 : ℝ) / 11 := by
      intro j hj
      have hjS : j ∈ S := (Finset.mem_filter.mp hj).1
      have hjplus : i.1 + 1 = j.1 := (Finset.mem_filter.mp hj).2
      have hij : i ≠ j := by
        have hjS' : j ∈ Finset.univ.erase i := by simpa [S] using hjS
        have hjNe : j ≠ i := (Finset.mem_erase.mp hjS').1
        exact fun h => hjNe h.symm
      have hnorm :
          ‖structuralRiemannMatrixC n i j‖ = (6 : ℝ) / 11 :=
        norm_structuralRiemannMatrixC_offdiag_eq_hop n i j hij (Or.inl hjplus)
      simpa [hnorm]
    have hsum_le :
        (∑ j ∈ Splus, ‖structuralRiemannMatrixC n i j‖) ≤
          (∑ _j ∈ Splus, (6 : ℝ) / 11) := by
      exact Finset.sum_le_sum (fun j hj => hterm j hj)
    have hconst :
        (∑ _j ∈ Splus, (6 : ℝ) / 11) = (Splus.card : ℝ) * ((6 : ℝ) / 11) := by
      simp [Finset.sum_const, nsmul_eq_mul]
    have hmul_le : (Splus.card : ℝ) * ((6 : ℝ) / 11) ≤ (1 : ℝ) * ((6 : ℝ) / 11) := by
      exact mul_le_mul_of_nonneg_right (by exact_mod_cast hcard_plus) (by positivity)
    calc
      (∑ j ∈ Splus, ‖structuralRiemannMatrixC n i j‖)
          ≤ (∑ _j ∈ Splus, (6 : ℝ) / 11) := hsum_le
      _ = (Splus.card : ℝ) * ((6 : ℝ) / 11) := hconst
      _ ≤ (1 : ℝ) * ((6 : ℝ) / 11) := hmul_le
      _ = (6 : ℝ) / 11 := by ring
  have hsum_minus_le : (∑ j ∈ Sminus, ‖structuralRiemannMatrixC n i j‖) ≤ (6 : ℝ) / 11 := by
    have hterm :
        ∀ j ∈ Sminus, ‖structuralRiemannMatrixC n i j‖ ≤ (6 : ℝ) / 11 := by
      intro j hj
      have hjS : j ∈ S := (Finset.mem_filter.mp hj).1
      have hjminus : j.1 + 1 = i.1 := (Finset.mem_filter.mp hj).2
      have hij : i ≠ j := by
        have hjS' : j ∈ Finset.univ.erase i := by simpa [S] using hjS
        have hjNe : j ≠ i := (Finset.mem_erase.mp hjS').1
        exact fun h => hjNe h.symm
      have hnorm :
          ‖structuralRiemannMatrixC n i j‖ = (6 : ℝ) / 11 :=
        norm_structuralRiemannMatrixC_offdiag_eq_hop n i j hij (Or.inr hjminus)
      simpa [hnorm]
    have hsum_le :
        (∑ j ∈ Sminus, ‖structuralRiemannMatrixC n i j‖) ≤
          (∑ _j ∈ Sminus, (6 : ℝ) / 11) := by
      exact Finset.sum_le_sum (fun j hj => hterm j hj)
    have hconst :
        (∑ _j ∈ Sminus, (6 : ℝ) / 11) = (Sminus.card : ℝ) * ((6 : ℝ) / 11) := by
      simp [Finset.sum_const, nsmul_eq_mul]
    have hmul_le : (Sminus.card : ℝ) * ((6 : ℝ) / 11) ≤ (1 : ℝ) * ((6 : ℝ) / 11) := by
      exact mul_le_mul_of_nonneg_right (by exact_mod_cast hcard_minus) (by positivity)
    calc
      (∑ j ∈ Sminus, ‖structuralRiemannMatrixC n i j‖)
          ≤ (∑ _j ∈ Sminus, (6 : ℝ) / 11) := hsum_le
      _ = (Sminus.card : ℝ) * ((6 : ℝ) / 11) := hconst
      _ ≤ (1 : ℝ) * ((6 : ℝ) / 11) := hmul_le
      _ = (6 : ℝ) / 11 := by ring
  have hS : S = Finset.univ.erase i := rfl
  calc
    (∑ j ∈ Finset.univ.erase i, ‖structuralRiemannMatrixC n i j‖)
        = (∑ j ∈ S, ‖structuralRiemannMatrixC n i j‖) := by simp [hS]
    _ = (∑ j ∈ Splus, ‖structuralRiemannMatrixC n i j‖) +
          (∑ j ∈ Sminus, ‖structuralRiemannMatrixC n i j‖) := by
            rw [hsplit0, hsplit]
    _ ≤ (6 : ℝ) / 11 + (6 : ℝ) / 11 := add_le_add hsum_plus_le hsum_minus_le
    _ = (12 : ℝ) / 11 := by norm_num

/-- Structural Gershgorin enclosure with explicit radius `12/11` for every
complex eigenvalue of the structural finite matrix lane. -/
theorem structuralRiemannMatrixC_eigenvalue_mem_ball_twelve_over_eleven
    (n : ℕ) {μ : ℂ} (hμ : μ ∈ spectrum ℂ (structuralRiemannMatrixC n)) :
    ∃ k : Fin n, μ ∈ Metric.closedBall (structuralRiemannMatrixC n k k) ((12 : ℝ) / 11) := by
  have hμLin : μ ∈ spectrum ℂ (Matrix.toLin' (structuralRiemannMatrixC n)) := by
    simpa [Matrix.spectrum_toLin'] using hμ
  have hEig : Module.End.HasEigenvalue (Matrix.toLin' (structuralRiemannMatrixC n)) μ :=
    Module.End.HasEigenvalue.of_mem_spectrum hμLin
  rcases (eigenvalue_mem_ball (A := structuralRiemannMatrixC n) (μ := μ) hEig) with ⟨k, hk⟩
  refine ⟨k, ?_⟩
  exact Set.mem_of_subset_of_mem
    (Metric.closedBall_subset_closedBall (structuralRiemannMatrixC_rowRadius_le_twelve_over_eleven n k))
    hk

/-- The center matrix is fixed by reversal/parity conjugation. -/
theorem structuralCenterMatrixC_revParity_fixed (n : ℕ) :
    revParityConjugateMatrixC n (structuralCenterMatrixC n)
      = structuralCenterMatrixC n := by
  ext i j
  by_cases hij : i = j
  · subst hij
    have hsq : (finParitySignQ i : ℂ) * (finParitySignQ i : ℂ) = 1 := by
      exact_mod_cast finParitySignQ_sq_one i
    calc
      revParityConjugateMatrixC n (structuralCenterMatrixC n) i i
          = (finParitySignQ i : ℂ) *
              structuralCenterMatrixC n (Fin.rev i) (Fin.rev i) *
              (finParitySignQ i : ℂ) := by
                simp [revParityConjugateMatrixC]
      _ = (finParitySignQ i : ℂ) * (structuralCenterQ n : ℂ) * (finParitySignQ i : ℂ) := by
            simp [structuralCenterMatrixC]
      _ = (structuralCenterQ n : ℂ) * ((finParitySignQ i : ℂ) * (finParitySignQ i : ℂ)) := by
            ring
      _ = (structuralCenterQ n : ℂ) := by simp [hsq]
      _ = structuralCenterMatrixC n i i := by simp [structuralCenterMatrixC]
  · have hrev : Fin.rev i ≠ Fin.rev j := by
      intro h
      exact hij (Fin.rev_injective h)
    calc
      revParityConjugateMatrixC n (structuralCenterMatrixC n) i j
          = (finParitySignQ i : ℂ) *
              structuralCenterMatrixC n (Fin.rev i) (Fin.rev j) *
              (finParitySignQ j : ℂ) := by
                simp [revParityConjugateMatrixC]
      _ = 0 := by
            simp [structuralCenterMatrixC, hrev]
      _ = structuralCenterMatrixC n i j := by
            simp [structuralCenterMatrixC, hij]

/-- Reversal permutation on `Fin n` as an equivalence. -/
def revPerm (n : ℕ) : Equiv.Perm (Fin n) :=
  Function.Involutive.toPerm Fin.rev (by intro i; simpa using Fin.rev_rev i)

/-- Complex diagonal parity-sign matrix. -/
def parityDiagC (n : ℕ) : Matrix (Fin n) (Fin n) ℂ :=
  Matrix.diagonal (fun i => (finParitySignQ i : ℂ))

/-- The parity-sign diagonal matrix is an involution. -/
theorem parityDiagC_sq (n : ℕ) : parityDiagC n * parityDiagC n = 1 := by
  ext i j
  by_cases hij : i = j
  · subst hij
    have hsqC : ((finParitySignQ i * finParitySignQ i : ℚ) : ℂ) = 1 := by
      exact_mod_cast (finParitySignQ_sq_one i)
    calc
      (parityDiagC n * parityDiagC n) i i
          = (finParitySignQ i : ℂ) * (finParitySignQ i : ℂ) := by
              simp [parityDiagC]
      _ = ((finParitySignQ i * finParitySignQ i : ℚ) : ℂ) := by simp
      _ = 1 := hsqC
      _ = (1 : Matrix (Fin n) (Fin n) ℂ) i i := by simp
  · simp [parityDiagC, hij]

/-- Explicit factorization of the reversal/parity conjugation as
`D * reindex(rev) * D`. -/
theorem revParityConjugateMatrixC_eq_diag_reindex_diag
    (n : ℕ) (A : Matrix (Fin n) (Fin n) ℂ) :
    revParityConjugateMatrixC n A =
      parityDiagC n * ((Matrix.reindex (revPerm n) (revPerm n)) A) * parityDiagC n := by
  ext i j
  simp [revParityConjugateMatrixC, parityDiagC, Matrix.diagonal_mul, Matrix.mul_diagonal,
    Matrix.reindex_apply, revPerm]

/-- Characteristic polynomial invariance under the reversal/parity conjugation action. -/
theorem charpoly_revParityConjugateMatrixC
    (n : ℕ) (A : Matrix (Fin n) (Fin n) ℂ) :
    (revParityConjugateMatrixC n A).charpoly = A.charpoly := by
  let e : Equiv.Perm (Fin n) := revPerm n
  let D : Matrix (Fin n) (Fin n) ℂ := parityDiagC n
  have hrepr :
      revParityConjugateMatrixC n A = D * ((Matrix.reindex e e) A) * D := by
    simpa [e, D] using revParityConjugateMatrixC_eq_diag_reindex_diag n A
  calc
    (revParityConjugateMatrixC n A).charpoly
        = (D * ((Matrix.reindex e e) A) * D).charpoly := by simpa [hrepr]
    _ = (D * (D * ((Matrix.reindex e e) A))).charpoly := by
          simpa [Matrix.mul_assoc] using (Matrix.charpoly_mul_comm (D * (Matrix.reindex e e A)) D)
    _ = (((D * D) * ((Matrix.reindex e e) A))).charpoly := by simp [Matrix.mul_assoc]
    _ = ((Matrix.reindex e e) A).charpoly := by simp [parityDiagC_sq, D]
    _ = A.charpoly := by simpa [e] using Matrix.charpoly_reindex e A

/-- Spectrum invariance under reversal/parity conjugation. -/
theorem spectrum_revParityConjugateMatrixC_eq
    (n : ℕ) (A : Matrix (Fin n) (Fin n) ℂ) :
    spectrum ℂ (revParityConjugateMatrixC n A) = spectrum ℂ A := by
  ext l
  constructor
  · intro hl
    exact (Matrix.mem_spectrum_iff_isRoot_charpoly).2
      (by simpa [charpoly_revParityConjugateMatrixC n A] using
        (Matrix.mem_spectrum_iff_isRoot_charpoly).1 hl)
  · intro hl
    exact (Matrix.mem_spectrum_iff_isRoot_charpoly).2
      (by simpa [charpoly_revParityConjugateMatrixC n A] using
        (Matrix.mem_spectrum_iff_isRoot_charpoly).1 hl)

/-- The center matrix is exactly the matrix-algebra scalar embedding of
`structuralCenterQ n`. -/
theorem structuralCenterMatrixC_eq_algebraMap (n : ℕ) :
    structuralCenterMatrixC n =
      (algebraMap ℂ (Matrix (Fin n) (Fin n) ℂ)) (structuralCenterQ n : ℂ) := by
  ext i j
  by_cases hij : i = j
  · subst hij
    simp [structuralCenterMatrixC, Matrix.algebraMap_eq_diagonal]
  · simp [structuralCenterMatrixC, Matrix.algebraMap_eq_diagonal, hij]

/-- Finite-level spectral reflection around the structural center:
if `λ ∈ spec(A_n)`, then `c_n - λ ∈ spec(A_n)` with `c_n = structuralCenterQ n`. -/
theorem structuralRiemannMatrixC_spectrum_reflect
    (n : ℕ) (l : ℂ)
    (hl : l ∈ spectrum ℂ (structuralRiemannMatrixC n)) :
    ((structuralCenterQ n : ℂ) - l) ∈ spectrum ℂ (structuralRiemannMatrixC n) := by
  let A : Matrix (Fin n) (Fin n) ℂ := structuralRiemannMatrixC n
  let c : ℂ := (structuralCenterQ n : ℂ)
  have hT : l ∈ spectrum ℂ (revParityConjugateMatrixC n A) := by
    simpa [spectrum_revParityConjugateMatrixC_eq n A] using hl
  have hCA : l ∈ spectrum ℂ (structuralCenterMatrixC n - A) := by
    simpa [A, structuralRiemannMatrixC_revParity_eq_center_sub n] using hT
  have hNotUnitY :
      ¬ IsUnit ((algebraMap ℂ (Matrix (Fin n) (Fin n) ℂ)) l - (structuralCenterMatrixC n - A)) :=
    (spectrum.mem_iff).1 hCA
  have hEqNeg :
      ((algebraMap ℂ (Matrix (Fin n) (Fin n) ℂ)) l - (structuralCenterMatrixC n - A))
        = - (((algebraMap ℂ (Matrix (Fin n) (Fin n) ℂ)) (c - l)) - A) := by
    rw [structuralCenterMatrixC_eq_algebraMap n]
    ext i j
    by_cases hij : i = j
    · subst hij
      simp [Matrix.algebraMap_eq_diagonal]
      ring
    · simp [Matrix.algebraMap_eq_diagonal, hij]
  have hNotUnitX :
      ¬ IsUnit (((algebraMap ℂ (Matrix (Fin n) (Fin n) ℂ)) (c - l)) - A) := by
    intro hX
    apply hNotUnitY
    rw [hEqNeg]
    exact hX.neg
  have hMemX : (c - l) ∈ spectrum ℂ A := (spectrum.mem_iff).2 hNotUnitX
  simpa [A, c] using hMemX

/-- Canonical real eigenvalue enumerator for the structural matrix at level `N+1`. -/
noncomputable def operatorEigenvalues (N : ℕ) : Fin (N + 1) → ℝ :=
  (structuralRiemannMatrixC_isHermitian (N + 1)).eigenvalues

/-- Ordered (`antitone`) Hermitian eigenvalue lane for the structural matrix
at level `N+1`, taken directly from Mathlib's `eigenvalues₀`. -/
noncomputable def operatorEigenvaluesOrdered (N : ℕ) :
    Fin (Fintype.card (Fin (N + 1))) → ℝ :=
  (structuralRiemannMatrixC_isHermitian (N + 1)).eigenvalues₀

/-- The ordered structural eigenvalue lane is antitone by Mathlib. -/
theorem operatorEigenvaluesOrdered_antitone (N : ℕ) :
    Antitone (operatorEigenvaluesOrdered N) := by
  simpa [operatorEigenvaluesOrdered] using
    (Matrix.IsHermitian.eigenvalues₀_antitone
      (hA := structuralRiemannMatrixC_isHermitian (N + 1)))

/-- Canonical reindex map from the `Fin (N+1)` lane to the `eigenvalues₀` lane. -/
noncomputable def operatorEigenvaluesReindexToOrdered (N : ℕ) :
    Fin (N + 1) → Fin (Fintype.card (Fin (N + 1))) :=
  (Fintype.equivOfCardEq (Fintype.card_fin _)).symm

/-- Equivalence form of the canonical reindex map from the `Fin (N+1)` lane to
the `eigenvalues₀` lane. -/
noncomputable def operatorEigenvaluesReindexToOrderedEquiv (N : ℕ) :
    Fin (N + 1) ≃ Fin (Fintype.card (Fin (N + 1))) :=
  (Fintype.equivOfCardEq (Fintype.card_fin _)).symm

/-- Current operator eigenvalue enumeration is exactly the ordered lane composed
with the canonical reindex map used by Mathlib's `eigenvalues` definition. -/
theorem operatorEigenvalues_eq_ordered_reindex (N : ℕ) (i : Fin (N + 1)) :
    operatorEigenvalues N i =
      operatorEigenvaluesOrdered N (operatorEigenvaluesReindexToOrdered N i) := by
  rfl

/-- Every concrete indexed operator eigenvalue is Gershgorin-close to some structural
diagonal center with explicit radius `12/11`. -/
theorem operatorEigenvalue_exists_center_gap_le_twelve_over_eleven
    (M : ℕ) (i : Fin (M + 1)) :
    ∃ k : Fin (M + 1),
      ‖((operatorEigenvalues M i : ℂ) - structuralRiemannMatrixC (M + 1) k k)‖ ≤ ((12 : ℝ) / 11) := by
  have hiReal : (operatorEigenvalues M i) ∈ spectrum ℝ (structuralRiemannMatrixC (M + 1)) := by
    exact Matrix.IsHermitian.eigenvalues_mem_spectrum_real
      (hA := structuralRiemannMatrixC_isHermitian (M + 1)) i
  have hiComplex : ((operatorEigenvalues M i : ℝ) : ℂ) ∈
      spectrum ℂ (structuralRiemannMatrixC (M + 1)) := by
    exact (spectrum.algebraMap_mem_iff ℂ).2 hiReal
  rcases structuralRiemannMatrixC_eigenvalue_mem_ball_twelve_over_eleven
      (n := M + 1) (μ := (operatorEigenvalues M i : ℂ)) hiComplex with ⟨k, hk⟩
  refine ⟨k, ?_⟩
  simpa [Metric.mem_closedBall, dist_eq_norm, sub_eq_add_neg, add_comm, add_left_comm, add_assoc]
    using hk

/-- Sharp uniform explicit lower bound on all indexed operator eigenvalues, obtained
from the structural Gershgorin enclosure and diagonal formula. -/
theorem operatorEigenvalue_lower_uniform_sharp
    (M : ℕ) (i : Fin (M + 1)) :
    ((127 : ℝ) / 176) ≤ operatorEigenvalues M i := by
  rcases operatorEigenvalue_exists_center_gap_le_twelve_over_eleven M i with ⟨k, hgap⟩
  let c : ℝ := (k.1 : ℝ) + (29 : ℝ) / 16
  have hdiag :
      structuralRiemannMatrixC (M + 1) k k = (c : ℂ) := by
    simp [structuralRiemannMatrixC, structuralRiemannMatrix_diag, c, timelikeOffsetQ]
    ring_nf
  have habs :
      |operatorEigenvalues M i - c| ≤ ((12 : ℝ) / 11) := by
    have hgap' : ‖((operatorEigenvalues M i : ℂ) - (c : ℂ))‖ ≤ ((12 : ℝ) / 11) := by
      simpa [hdiag, sub_eq_add_neg, add_comm, add_left_comm, add_assoc] using hgap
    have hgap'' : ‖(((operatorEigenvalues M i - c : ℝ) : ℂ))‖ ≤ ((12 : ℝ) / 11) := by
      convert hgap' using 1
      simp
    calc
      |operatorEigenvalues M i - c| = ‖(((operatorEigenvalues M i - c : ℝ) : ℂ))‖ := by
        simpa using (Complex.norm_real (operatorEigenvalues M i - c)).symm
      _ ≤ ((12 : ℝ) / 11) := hgap''
  have hc_ge : (29 : ℝ) / 16 ≤ c := by
    have hk_nonneg : (0 : ℝ) ≤ k.1 := by positivity
    dsimp [c]
    linarith
  have hlower1 : c - ((12 : ℝ) / 11) ≤ operatorEigenvalues M i := by
    have hpair := abs_le.mp habs
    linarith [hpair.1]
  have hlower_const : (127 : ℝ) / 176 ≤ c - (12 : ℝ) / 11 := by
    have h127_eq : (127 : ℝ) / 176 = ((29 : ℝ) / 16) - ((12 : ℝ) / 11) := by norm_num
    linarith [hc_ge, h127_eq]
  calc
    (127 : ℝ) / 176 ≤ c - (12 : ℝ) / 11 := hlower_const
    _ ≤ operatorEigenvalues M i := hlower1

/-- Backward-compatible coarse uniform lower bound, derived from the sharp bound. -/
theorem operatorEigenvalue_lower_uniform
    (M : ℕ) (i : Fin (M + 1)) :
    ((117 : ℝ) / 176) ≤ operatorEigenvalues M i := by
  have hsharp : ((127 : ℝ) / 176) ≤ operatorEigenvalues M i :=
    operatorEigenvalue_lower_uniform_sharp M i
  have hcoarse : ((117 : ℝ) / 176) ≤ ((127 : ℝ) / 176) := by norm_num
  exact le_trans hcoarse hsharp

/-- Weyl-style per-index center-gap contract for the structural operator lane.
If true, this is the single perturbative bridge from `D + C` to indexed growth. -/
def OperatorWeylCenterGap : Prop :=
  ∀ M : ℕ, ∀ i : Fin (M + 1),
    |operatorEigenvalues M i - ((i.1 : ℝ) + (29 : ℝ) / 16)| ≤ (12 : ℝ) / 11

/-- Permutation-invariant center-gap contract for the structural operator lane.
At each finite level, eigenvalues can be bijectively paired with diagonal centers
`k + 29/16` with uniform gap `≤ 12/11`; no index ordering assumption is used. -/
def OperatorCenterGapPermutationInvariant : Prop :=
  ∀ M : ℕ, ∃ σ : Fin (M + 1) ≃ Fin (M + 1), ∀ i : Fin (M + 1),
    |operatorEigenvalues M i - (((σ i).1 : ℝ) + (29 : ℝ) / 16)| ≤ (12 : ℝ) / 11

/-- Candidate center indices for the `i`-th finite-level eigenvalue:
all structural centers within radius `12/11` in the real gap metric. -/
def operatorCenterCandidates (M : ℕ) (i : Fin (M + 1)) : Finset (Fin (M + 1)) :=
  Finset.univ.filter (fun k =>
    |operatorEigenvalues M i - ((k.1 : ℝ) + (29 : ℝ) / 16)| ≤ (12 : ℝ) / 11)

/-- Real-gap form of the concrete Gershgorin enclosure:
every indexed operator eigenvalue is within `12/11` of some structural center
`k + 29/16`. -/
theorem operatorEigenvalue_exists_center_gap_real_le_twelve_over_eleven
    (M : ℕ) (i : Fin (M + 1)) :
    ∃ k : Fin (M + 1),
      |operatorEigenvalues M i - ((k.1 : ℝ) + (29 : ℝ) / 16)| ≤ (12 : ℝ) / 11 := by
  rcases operatorEigenvalue_exists_center_gap_le_twelve_over_eleven M i with ⟨k, hgap⟩
  let c : ℝ := (k.1 : ℝ) + (29 : ℝ) / 16
  have hdiag :
      structuralRiemannMatrixC (M + 1) k k = (c : ℂ) := by
    simp [structuralRiemannMatrixC, structuralRiemannMatrix_diag, c, timelikeOffsetQ]
    ring_nf
  have hgap' : ‖((operatorEigenvalues M i : ℂ) - (c : ℂ))‖ ≤ ((12 : ℝ) / 11) := by
    simpa [hdiag, sub_eq_add_neg, add_comm, add_left_comm, add_assoc] using hgap
  have hgap'' : ‖(((operatorEigenvalues M i - c : ℝ) : ℂ))‖ ≤ ((12 : ℝ) / 11) := by
    convert hgap' using 1
    simp
  have habs : |operatorEigenvalues M i - c| ≤ ((12 : ℝ) / 11) := by
    calc
      |operatorEigenvalues M i - c| = ‖(((operatorEigenvalues M i - c : ℝ) : ℂ))‖ := by
        simpa using (Complex.norm_real (operatorEigenvalues M i - c)).symm
      _ ≤ ((12 : ℝ) / 11) := hgap''
  refine ⟨k, ?_⟩
  simpa [c] using habs

/-- Ordered-lane real Gershgorin gap witness:
each ordered eigenvalue is within `12/11` of some structural center. -/
theorem operatorEigenvalueOrdered_exists_center_gap_real_le_twelve_over_eleven
    (M : ℕ) (j : Fin (Fintype.card (Fin (M + 1)))) :
    ∃ k : Fin (M + 1),
      |operatorEigenvaluesOrdered M j - ((k.1 : ℝ) + (29 : ℝ) / 16)| ≤ (12 : ℝ) / 11 := by
  let e : Fin (M + 1) ≃ Fin (Fintype.card (Fin (M + 1))) :=
    operatorEigenvaluesReindexToOrderedEquiv M
  let i : Fin (M + 1) := e.symm j
  rcases operatorEigenvalue_exists_center_gap_real_le_twelve_over_eleven M i with ⟨k, hk⟩
  refine ⟨k, ?_⟩
  have hi :
      operatorEigenvalues M i = operatorEigenvaluesOrdered M j := by
    calc
      operatorEigenvalues M i
          = operatorEigenvaluesOrdered M (operatorEigenvaluesReindexToOrdered M i) := by
            simpa using operatorEigenvalues_eq_ordered_reindex M i
      _ = operatorEigenvaluesOrdered M (e i) := by rfl
      _ = operatorEigenvaluesOrdered M j := by simp [i]
  simpa [hi] using hk

/-- Chosen center index for each ordered eigenvalue on level `M`, extracted from
the ordered-lane Gershgorin real-gap witness. -/
noncomputable def operatorEigenvalueOrderedCenterChoice
    (M : ℕ) (j : Fin (Fintype.card (Fin (M + 1)))) : Fin (M + 1) :=
  Classical.choose (operatorEigenvalueOrdered_exists_center_gap_real_le_twelve_over_eleven M j)

/-- Ordered-lane candidate center set at index `j`. -/
def operatorCenterCandidatesOrdered
    (M : ℕ) (j : Fin (Fintype.card (Fin (M + 1)))) : Finset (Fin (M + 1)) :=
  (Finset.univ : Finset (Fin (M + 1))).filter
    (fun k =>
      |operatorEigenvaluesOrdered M j - ((k.1 : ℝ) + (29 : ℝ) / 16)| ≤ (12 : ℝ) / 11)

/-- The ordered-lane candidate set is nonempty (Gershgorin witness). -/
theorem operatorCenterCandidatesOrdered_nonempty
    (M : ℕ) (j : Fin (Fintype.card (Fin (M + 1)))) :
    (operatorCenterCandidatesOrdered M j).Nonempty := by
  rcases operatorEigenvalueOrdered_exists_center_gap_real_le_twelve_over_eleven M j with ⟨k, hk⟩
  refine ⟨k, ?_⟩
  exact Finset.mem_filter.mpr ⟨Finset.mem_univ k, hk⟩

/-- Canonical deterministic ordered-lane center selector: pick the maximal
admissible candidate center index. -/
noncomputable def operatorEigenvalueOrderedCenterChoiceMax
    (M : ℕ) (j : Fin (Fintype.card (Fin (M + 1)))) : Fin (M + 1) :=
  (operatorCenterCandidatesOrdered M j).max' (operatorCenterCandidatesOrdered_nonempty M j)

/-- Membership certificate for the maximal ordered-lane center choice. -/
theorem operatorEigenvalueOrderedCenterChoiceMax_mem
    (M : ℕ) (j : Fin (Fintype.card (Fin (M + 1)))) :
    operatorEigenvalueOrderedCenterChoiceMax M j ∈ operatorCenterCandidatesOrdered M j := by
  exact Finset.max'_mem (operatorCenterCandidatesOrdered M j)
    (operatorCenterCandidatesOrdered_nonempty M j)

/-- Gap certificate for the maximal ordered-lane center choice. -/
theorem operatorEigenvalueOrderedCenterChoiceMax_spec
    (M : ℕ) (j : Fin (Fintype.card (Fin (M + 1)))) :
    |operatorEigenvaluesOrdered M j
      - (((operatorEigenvalueOrderedCenterChoiceMax M j).1 : ℝ) + (29 : ℝ) / 16)| ≤ (12 : ℝ) / 11 := by
  have hmem : operatorEigenvalueOrderedCenterChoiceMax M j ∈ operatorCenterCandidatesOrdered M j :=
    operatorEigenvalueOrderedCenterChoiceMax_mem M j
  exact (Finset.mem_filter.mp hmem).2

/-- Canonical deterministic ordered-lane lower endpoint selector: pick the minimal
admissible candidate center index. -/
noncomputable def operatorEigenvalueOrderedCenterChoiceMin
    (M : ℕ) (j : Fin (Fintype.card (Fin (M + 1)))) : Fin (M + 1) :=
  (operatorCenterCandidatesOrdered M j).min' (operatorCenterCandidatesOrdered_nonempty M j)

/-- Membership certificate for the minimal ordered-lane center choice. -/
theorem operatorEigenvalueOrderedCenterChoiceMin_mem
    (M : ℕ) (j : Fin (Fintype.card (Fin (M + 1)))) :
    operatorEigenvalueOrderedCenterChoiceMin M j ∈ operatorCenterCandidatesOrdered M j := by
  exact Finset.min'_mem (operatorCenterCandidatesOrdered M j)
    (operatorCenterCandidatesOrdered_nonempty M j)

/-- Gap certificate for the minimal ordered-lane center choice. -/
theorem operatorEigenvalueOrderedCenterChoiceMin_spec
    (M : ℕ) (j : Fin (Fintype.card (Fin (M + 1)))) :
    |operatorEigenvaluesOrdered M j
      - (((operatorEigenvalueOrderedCenterChoiceMin M j).1 : ℝ) + (29 : ℝ) / 16)| ≤ (12 : ℝ) / 11 := by
  have hmem : operatorEigenvalueOrderedCenterChoiceMin M j ∈ operatorCenterCandidatesOrdered M j :=
    operatorEigenvalueOrderedCenterChoiceMin_mem M j
  exact (Finset.mem_filter.mp hmem).2

/-- Gap certificate for the chosen ordered-lane center index. -/
theorem operatorEigenvalueOrderedCenterChoice_spec
    (M : ℕ) (j : Fin (Fintype.card (Fin (M + 1)))) :
    |operatorEigenvaluesOrdered M j
      - (((operatorEigenvalueOrderedCenterChoice M j).1 : ℝ) + (29 : ℝ) / 16)| ≤ (12 : ℝ) / 11 := by
  exact Classical.choose_spec
    (operatorEigenvalueOrdered_exists_center_gap_real_le_twelve_over_eleven M j)

/-- Generic ordered-lane jump exclusion:
for any center selector satisfying the `12/11` gap spec, indices cannot jump up
by `3` or more along the ordered eigenvalue lane. -/
theorem operatorEigenvalueOrdered_no_up_jump_three_of_spec
    (M : ℕ)
    (f : Fin (Fintype.card (Fin (M + 1))) → Fin (M + 1))
    (hSpec : ∀ j : Fin (Fintype.card (Fin (M + 1))),
      |operatorEigenvaluesOrdered M j - (((f j).1 : ℝ) + (29 : ℝ) / 16)| ≤ (12 : ℝ) / 11)
    {i j : Fin (Fintype.card (Fin (M + 1)))} (hij : i < j)
    (hjump : (f i).1 + 3 ≤ (f j).1) :
    False := by
  let ki : ℕ := (f i).1
  let kj : ℕ := (f j).1
  let ci : ℝ := (ki : ℝ) + (29 : ℝ) / 16
  let cj : ℝ := (kj : ℝ) + (29 : ℝ) / 16
  have hgi :
      |operatorEigenvaluesOrdered M i - ci| ≤ (12 : ℝ) / 11 := by
    simpa [ki, ci] using hSpec i
  have hgj :
      |operatorEigenvaluesOrdered M j - cj| ≤ (12 : ℝ) / 11 := by
    simpa [kj, cj] using hSpec j
  have hi_upper : operatorEigenvaluesOrdered M i ≤ ci + (12 : ℝ) / 11 := by
    linarith [(abs_le.mp hgi).2]
  have hj_lower : cj - (12 : ℝ) / 11 ≤ operatorEigenvaluesOrdered M j := by
    linarith [(abs_le.mp hgj).1]
  have hkj_ge : (ki : ℝ) + 3 ≤ (kj : ℝ) := by exact_mod_cast hjump
  have hcj_ge : ci + 3 ≤ cj := by
    dsimp [ci, cj]
    linarith
  have hstrict : ci + (12 : ℝ) / 11 < cj - (12 : ℝ) / 11 := by
    linarith [hcj_ge]
  have hlt : operatorEigenvaluesOrdered M i < operatorEigenvaluesOrdered M j := by
    linarith [hi_upper, hj_lower, hstrict]
  have hanti := operatorEigenvaluesOrdered_antitone M
  have hle : operatorEigenvaluesOrdered M j ≤ operatorEigenvaluesOrdered M i := hanti hij.le
  exact not_lt_of_ge hle hlt

/-- Ordered-lane jump exclusion:
for `i < j`, the chosen center index cannot increase by `3` or more. -/
theorem operatorEigenvalueOrderedCenterChoice_no_up_jump_three
    (M : ℕ)
    {i j : Fin (Fintype.card (Fin (M + 1)))} (hij : i < j)
    (hjump :
      (operatorEigenvalueOrderedCenterChoice M i).1 + 3
        ≤ (operatorEigenvalueOrderedCenterChoice M j).1) :
    False := by
  exact operatorEigenvalueOrdered_no_up_jump_three_of_spec M
    (f := operatorEigenvalueOrderedCenterChoice M)
    (hSpec := operatorEigenvalueOrderedCenterChoice_spec M)
    hij hjump

/-- Ordered-lane jump exclusion for the canonical maximal center selector. -/
theorem operatorEigenvalueOrderedCenterChoiceMax_no_up_jump_three
    (M : ℕ)
    {i j : Fin (Fintype.card (Fin (M + 1)))} (hij : i < j)
    (hjump :
      (operatorEigenvalueOrderedCenterChoiceMax M i).1 + 3
        ≤ (operatorEigenvalueOrderedCenterChoiceMax M j).1) :
    False := by
  exact operatorEigenvalueOrdered_no_up_jump_three_of_spec M
    (f := operatorEigenvalueOrderedCenterChoiceMax M)
    (hSpec := operatorEigenvalueOrderedCenterChoiceMax_spec M)
    hij hjump

/-- Ordered-lane jump exclusion for the canonical minimal center selector. -/
theorem operatorEigenvalueOrderedCenterChoiceMin_no_up_jump_three
    (M : ℕ)
    {i j : Fin (Fintype.card (Fin (M + 1)))} (hij : i < j)
    (hjump :
      (operatorEigenvalueOrderedCenterChoiceMin M i).1 + 3
        ≤ (operatorEigenvalueOrderedCenterChoiceMin M j).1) :
    False := by
  exact operatorEigenvalueOrdered_no_up_jump_three_of_spec M
    (f := operatorEigenvalueOrderedCenterChoiceMin M)
    (hSpec := operatorEigenvalueOrderedCenterChoiceMin_spec M)
    hij hjump

/-- Monotonicity of the canonical maximal ordered-lane center selector:
as ordered eigenvalue index increases, the selected center index cannot increase. -/
theorem operatorEigenvalueOrderedCenterChoiceMax_antitone
    (M : ℕ) :
    Antitone (fun j : Fin (Fintype.card (Fin (M + 1))) =>
      (operatorEigenvalueOrderedCenterChoiceMax M j).1) := by
  intro i j hij
  let Ki : Fin (M + 1) := operatorEigenvalueOrderedCenterChoiceMax M i
  let Kj : Fin (M + 1) := operatorEigenvalueOrderedCenterChoiceMax M j
  let xi : ℝ := operatorEigenvaluesOrdered M i
  let xj : ℝ := operatorEigenvaluesOrdered M j
  let ci : ℝ := (Ki.1 : ℝ) + (29 : ℝ) / 16
  let cj : ℝ := (Kj.1 : ℝ) + (29 : ℝ) / 16
  have hi_abs : |xi - ci| ≤ (12 : ℝ) / 11 := by
    simpa [Ki, xi, ci] using operatorEigenvalueOrderedCenterChoiceMax_spec M i
  have hj_abs : |xj - cj| ≤ (12 : ℝ) / 11 := by
    simpa [Kj, xj, cj] using operatorEigenvalueOrderedCenterChoiceMax_spec M j
  have hanti := operatorEigenvaluesOrdered_antitone M
  have hxj_le_xi : xj ≤ xi := hanti hij
  have hcj_le_xj_plus : cj ≤ xj + (12 : ℝ) / 11 := by
    linarith [(abs_le.mp hj_abs).1]
  have hcj_le_xi_plus : cj ≤ xi + (12 : ℝ) / 11 := by
    linarith [hcj_le_xj_plus, hxj_le_xi]
  by_cases hxi_minus_le_cj : xi - (12 : ℝ) / 11 ≤ cj
  · have habs_xi_cj : |xi - cj| ≤ (12 : ℝ) / 11 := by
      refine abs_le.mpr ?_
      constructor
      · linarith [hcj_le_xi_plus]
      · linarith [hxi_minus_le_cj]
    have hmem :
        Kj ∈ operatorCenterCandidatesOrdered M i := by
      refine Finset.mem_filter.mpr ?_
      refine ⟨Finset.mem_univ Kj, ?_⟩
      simpa [xi, cj] using habs_xi_cj
    have hle_fin : Kj ≤ Ki := Finset.le_max' (operatorCenterCandidatesOrdered M i) Kj hmem
    exact hle_fin
  · have hcj_lt_xi_minus : cj < xi - (12 : ℝ) / 11 := lt_of_not_ge hxi_minus_le_cj
    have hxi_minus_le_ci : xi - (12 : ℝ) / 11 ≤ ci := by
      linarith [(abs_le.mp hi_abs).2]
    have hcj_lt_ci : cj < ci := by
      linarith [hcj_lt_xi_minus, hxi_minus_le_ci]
    have hkj_lt_ki_real : ((Kj.1 : ℝ)) < (Ki.1 : ℝ) := by
      dsimp [ci, cj] at hcj_lt_ci ⊢
      linarith
    have hkj_lt_ki : Kj.1 < Ki.1 := by
      exact_mod_cast hkj_lt_ki_real
    exact le_of_lt hkj_lt_ki

/-- Any ordered-lane candidate center lies within the last three indices
ending at the maximal ordered-lane center choice. -/
theorem operatorCenterCandidatesOrdered_index_bounds_of_max
    (M : ℕ) (j : Fin (Fintype.card (Fin (M + 1))))
    {i : Fin (M + 1)}
    (hi : i ∈ operatorCenterCandidatesOrdered M j) :
    (operatorEigenvalueOrderedCenterChoiceMax M j).1 - 2 ≤ i.1 ∧
      i.1 ≤ (operatorEigenvalueOrderedCenterChoiceMax M j).1 := by
  let k : Fin (M + 1) := operatorEigenvalueOrderedCenterChoiceMax M j
  let x : ℝ := operatorEigenvaluesOrdered M j
  let ci : ℝ := (i.1 : ℝ) + (29 : ℝ) / 16
  let ck : ℝ := (k.1 : ℝ) + (29 : ℝ) / 16
  have hi_abs : |x - ci| ≤ (12 : ℝ) / 11 := by
    exact (Finset.mem_filter.mp hi).2
  have hk_abs : |x - ck| ≤ (12 : ℝ) / 11 := by
    simpa [k, x, ck] using operatorEigenvalueOrderedCenterChoiceMax_spec M j
  have hi_le_k_fin : i ≤ k := Finset.le_max' (operatorCenterCandidatesOrdered M j) i hi
  have hi_le_k : i.1 ≤ k.1 := hi_le_k_fin
  have hci_lower : x - (12 : ℝ) / 11 ≤ ci := by
    linarith [(abs_le.mp hi_abs).2]
  have hck_upper : ck ≤ x + (12 : ℝ) / 11 := by
    linarith [(abs_le.mp hk_abs).1]
  have hdiff_le : ck - ci ≤ (24 : ℝ) / 11 := by
    linarith [hck_upper, hci_lower]
  have hdiff_lt_three : ck - ci < (3 : ℝ) := by
    linarith [hdiff_le]
  have hk_lt_i_plus_three_real : (k.1 : ℝ) < (i.1 : ℝ) + 3 := by
    dsimp [ck, ci] at hdiff_lt_three ⊢
    linarith
  have hk_lt_i_plus_three : k.1 < i.1 + 3 := by
    exact_mod_cast hk_lt_i_plus_three_real
  have hk_le_i_plus_two : k.1 ≤ i.1 + 2 := by
    have htmp : k.1 < i.1 + 2 + 1 := by
      simpa [Nat.add_assoc, Nat.add_comm, Nat.add_left_comm] using hk_lt_i_plus_three
    exact Nat.lt_succ_iff.mp htmp
  have hk_sub_two_le_i : k.1 - 2 ≤ i.1 := by
    omega
  exact ⟨hk_sub_two_le_i, hi_le_k⟩

/-- Cardinality bound for ordered-lane candidate centers: each candidate set has
at most three indices, concentrated near its maximal candidate center. -/
theorem operatorCenterCandidatesOrdered_card_le_three
    (M : ℕ) (j : Fin (Fintype.card (Fin (M + 1)))) :
    (operatorCenterCandidatesOrdered M j).card ≤ 3 := by
  classical
  let k : Fin (M + 1) := operatorEigenvalueOrderedCenterChoiceMax M j
  let Sfin : Finset (Fin (M + 1)) := operatorCenterCandidatesOrdered M j
  let Snat : Finset ℕ := Sfin.image (fun i : Fin (M + 1) => i.1)
  let T : Finset ℕ := {k.1 - 2, k.1 - 1, k.1}
  have hcard_img : Snat.card = Sfin.card := by
    simpa [Snat, Sfin] using
      (Finset.card_image_of_injOn (s := Sfin) (f := fun i : Fin (M + 1) => i.1)
        (by intro a ha b hb hab; exact Fin.ext hab))
  have hsubset : Snat ⊆ T := by
    intro n hn
    rcases Finset.mem_image.mp hn with ⟨i, hiS, rfl⟩
    have hbounds :=
      operatorCenterCandidatesOrdered_index_bounds_of_max M j (i := i) hiS
    have hcases : i.1 = k.1 - 2 ∨ i.1 = k.1 - 1 ∨ i.1 = k.1 := by
      omega
    rcases hcases with hEq | hEq | hEq
    · simpa [T, hEq]
    · simpa [T, hEq]
    · simpa [T, hEq]
  have hT_card_le_three : T.card ≤ 3 := by
    simpa [T] using (Finset.card_le_three (a := k.1 - 2) (b := k.1 - 1) (c := k.1))
  calc
    (operatorCenterCandidatesOrdered M j).card = Sfin.card := by rfl
    _ = Snat.card := hcard_img.symm
    _ ≤ T.card := Finset.card_le_card hsubset
    _ ≤ 3 := hT_card_le_three

/-- Unordered-lane candidate-center cardinality bound, transferred from the
ordered eigenvalue lane via the reindex equivalence. -/
theorem operatorCenterCandidates_card_le_three
    (M : ℕ) (i : Fin (M + 1)) :
    (operatorCenterCandidates M i).card ≤ 3 := by
  let j : Fin (Fintype.card (Fin (M + 1))) := operatorEigenvaluesReindexToOrdered M i
  have hij :
      operatorEigenvalues M i = operatorEigenvaluesOrdered M j := by
    simpa [j] using operatorEigenvalues_eq_ordered_reindex M i
  have hsetEq :
      operatorCenterCandidates M i = operatorCenterCandidatesOrdered M j := by
    ext k
    constructor <;> intro hk
    · have hk' := Finset.mem_filter.mp hk
      refine Finset.mem_filter.mpr ?_
      refine ⟨hk'.1, ?_⟩
      simpa [operatorCenterCandidates, operatorCenterCandidatesOrdered, hij] using hk'.2
    · have hk' := Finset.mem_filter.mp hk
      refine Finset.mem_filter.mpr ?_
      refine ⟨hk'.1, ?_⟩
      simpa [operatorCenterCandidates, operatorCenterCandidatesOrdered, hij] using hk'.2
  calc
    (operatorCenterCandidates M i).card = (operatorCenterCandidatesOrdered M j).card := by
      simpa [hsetEq]
    _ ≤ 3 := operatorCenterCandidatesOrdered_card_le_three M j

/-- Ordered candidate sets are index-convex:
if two center indices are admissible at level `j`, every index between them is
also admissible. -/
theorem operatorCenterCandidatesOrdered_mem_of_between
    (M : ℕ) (j : Fin (Fintype.card (Fin (M + 1))))
    {i k m : Fin (M + 1)}
    (hi : i ∈ operatorCenterCandidatesOrdered M j)
    (hk : k ∈ operatorCenterCandidatesOrdered M j)
    (him : i.1 ≤ m.1) (hmk : m.1 ≤ k.1) :
    m ∈ operatorCenterCandidatesOrdered M j := by
  have hi_abs : |operatorEigenvaluesOrdered M j - ((i.1 : ℝ) + (29 : ℝ) / 16)| ≤ (12 : ℝ) / 11 :=
    (Finset.mem_filter.mp hi).2
  have hk_abs : |operatorEigenvaluesOrdered M j - ((k.1 : ℝ) + (29 : ℝ) / 16)| ≤ (12 : ℝ) / 11 :=
    (Finset.mem_filter.mp hk).2
  have hi_lo : operatorEigenvaluesOrdered M j - (12 : ℝ) / 11 ≤ (i.1 : ℝ) + (29 : ℝ) / 16 := by
    linarith [(abs_le.mp hi_abs).2]
  have hk_hi : (k.1 : ℝ) + (29 : ℝ) / 16 ≤ operatorEigenvaluesOrdered M j + (12 : ℝ) / 11 := by
    linarith [(abs_le.mp hk_abs).1]
  have himR : (i.1 : ℝ) ≤ (m.1 : ℝ) := by exact_mod_cast him
  have hmkR : (m.1 : ℝ) ≤ (k.1 : ℝ) := by exact_mod_cast hmk
  have hm_lo :
      operatorEigenvaluesOrdered M j - (12 : ℝ) / 11 ≤ (m.1 : ℝ) + (29 : ℝ) / 16 := by
    linarith [hi_lo, himR]
  have hm_hi :
      (m.1 : ℝ) + (29 : ℝ) / 16 ≤ operatorEigenvaluesOrdered M j + (12 : ℝ) / 11 := by
    linarith [hk_hi, hmkR]
  have hm_abs : |operatorEigenvaluesOrdered M j - ((m.1 : ℝ) + (29 : ℝ) / 16)| ≤ (12 : ℝ) / 11 := by
    refine abs_le.mpr ?_
    constructor <;> linarith [hm_lo, hm_hi]
  refine Finset.mem_filter.mpr ?_
  exact ⟨Finset.mem_univ m, hm_abs⟩

/-- Exact interval characterization of ordered candidate sets:
membership is equivalent to lying between the canonical min and max selectors. -/
theorem operatorCenterCandidatesOrdered_mem_iff_between_min_max
    (M : ℕ) (j : Fin (Fintype.card (Fin (M + 1))))
    (i : Fin (M + 1)) :
    i ∈ operatorCenterCandidatesOrdered M j ↔
      (operatorEigenvalueOrderedCenterChoiceMin M j).1 ≤ i.1 ∧
        i.1 ≤ (operatorEigenvalueOrderedCenterChoiceMax M j).1 := by
  constructor
  · intro hi
    have hmin_mem : operatorEigenvalueOrderedCenterChoiceMin M j ∈ operatorCenterCandidatesOrdered M j :=
      operatorEigenvalueOrderedCenterChoiceMin_mem M j
    have hmax_mem : operatorEigenvalueOrderedCenterChoiceMax M j ∈ operatorCenterCandidatesOrdered M j :=
      operatorEigenvalueOrderedCenterChoiceMax_mem M j
    have hmin_le_i_fin :
        operatorEigenvalueOrderedCenterChoiceMin M j ≤ i :=
      Finset.min'_le (operatorCenterCandidatesOrdered M j) i hi
    have hi_le_max_fin :
        i ≤ operatorEigenvalueOrderedCenterChoiceMax M j :=
      Finset.le_max' (operatorCenterCandidatesOrdered M j) i hi
    exact ⟨hmin_le_i_fin, hi_le_max_fin⟩
  · intro hi
    have hmin_mem : operatorEigenvalueOrderedCenterChoiceMin M j ∈ operatorCenterCandidatesOrdered M j :=
      operatorEigenvalueOrderedCenterChoiceMin_mem M j
    have hmax_mem : operatorEigenvalueOrderedCenterChoiceMax M j ∈ operatorCenterCandidatesOrdered M j :=
      operatorEigenvalueOrderedCenterChoiceMax_mem M j
    exact operatorCenterCandidatesOrdered_mem_of_between M j
      (i := operatorEigenvalueOrderedCenterChoiceMin M j)
      (k := operatorEigenvalueOrderedCenterChoiceMax M j)
      (m := i) hmin_mem hmax_mem hi.1 hi.2

/-- The ordered candidate interval has width at most `2` in index units. -/
theorem operatorCenterCandidatesOrdered_width_le_two
    (M : ℕ) (j : Fin (Fintype.card (Fin (M + 1)))) :
    (operatorEigenvalueOrderedCenterChoiceMax M j).1
      ≤ (operatorEigenvalueOrderedCenterChoiceMin M j).1 + 2 := by
  have hmin_mem : operatorEigenvalueOrderedCenterChoiceMin M j ∈ operatorCenterCandidatesOrdered M j :=
    operatorEigenvalueOrderedCenterChoiceMin_mem M j
  have hbounds :=
    operatorCenterCandidatesOrdered_index_bounds_of_max M j
      (i := operatorEigenvalueOrderedCenterChoiceMin M j) hmin_mem
  omega

/-- Arithmetic overlap kernel for the structural radius `12/11` and unit center spacing:
if `|y| ≤ 12/11`, then at least one shifted neighbor is also within `12/11`. -/
theorem abs_shift_neighbor_le_twelve_over_eleven
    (y : ℝ) (hy : |y| ≤ (12 : ℝ) / 11) :
    |y + 1| ≤ (12 : ℝ) / 11 ∨ |y - 1| ≤ (12 : ℝ) / 11 := by
  by_contra h
  have hplus : ¬ |y + 1| ≤ (12 : ℝ) / 11 := by
    exact (not_or.mp h).1
  have hminus : ¬ |y - 1| ≤ (12 : ℝ) / 11 := by
    exact (not_or.mp h).2
  have hylo : -((12 : ℝ) / 11) ≤ y := by
    linarith [(abs_le.mp hy).1]
  have hyhi : y ≤ (12 : ℝ) / 11 := by
    linarith [(abs_le.mp hy).2]
  have hy_gt_one_over_eleven : (1 : ℝ) / 11 < y := by
    by_cases hnonneg : 0 ≤ y + 1
    · have hAbsEq : |y + 1| = y + 1 := abs_of_nonneg hnonneg
      have hgt : (12 : ℝ) / 11 < y + 1 := by
        simpa [hAbsEq] using (lt_of_not_ge hplus)
      linarith [hgt]
    · have hy1lt : y + 1 < 0 := lt_of_not_ge hnonneg
      have hAbsEq : |y + 1| = -(y + 1) := abs_of_neg hy1lt
      have hgt : (12 : ℝ) / 11 < -(y + 1) := by simpa [hAbsEq] using (lt_of_not_ge hplus)
      linarith [hgt, hylo]
  have hy_lt_neg_one_over_eleven : y < -((1 : ℝ) / 11) := by
    by_cases hnonneg : 0 ≤ y - 1
    · have hAbsEq : |y - 1| = y - 1 := abs_of_nonneg hnonneg
      have hgt : (12 : ℝ) / 11 < y - 1 := by
        simpa [hAbsEq] using (lt_of_not_ge hminus)
      linarith [hgt, hyhi]
    · have hy1lt : y - 1 < 0 := lt_of_not_ge hnonneg
      have hAbsEq : |y - 1| = -(y - 1) := abs_of_neg hy1lt
      have hgt : (12 : ℝ) / 11 < -(y - 1) := by simpa [hAbsEq] using (lt_of_not_ge hminus)
      linarith [hgt]
  linarith

/-- Interior ordered candidate intervals cannot collapse to a singleton:
if the minimum candidate index is strictly inside `0..M`, then the maximum is
strictly larger than the minimum. -/
theorem operatorCenterCandidatesOrdered_min_lt_max_of_interior
    (M : ℕ) (j : Fin (Fintype.card (Fin (M + 1))))
    (hmin_pos : 0 < (operatorEigenvalueOrderedCenterChoiceMin M j).1)
    (hmax_lt : (operatorEigenvalueOrderedCenterChoiceMax M j).1 < M) :
    (operatorEigenvalueOrderedCenterChoiceMin M j).1 <
      (operatorEigenvalueOrderedCenterChoiceMax M j).1 := by
  let m : Fin (M + 1) := operatorEigenvalueOrderedCenterChoiceMin M j
  let kmax : Fin (M + 1) := operatorEigenvalueOrderedCenterChoiceMax M j
  let x : ℝ := operatorEigenvaluesOrdered M j
  have hm_mem : m ∈ operatorCenterCandidatesOrdered M j :=
    operatorEigenvalueOrderedCenterChoiceMin_mem M j
  have hkmax_mem : kmax ∈ operatorCenterCandidatesOrdered M j :=
    operatorEigenvalueOrderedCenterChoiceMax_mem M j
  have hmin_le_max_fin : m ≤ kmax := by
    exact Finset.min'_le (operatorCenterCandidatesOrdered M j) kmax hkmax_mem
  have hmin_le_max : m.1 ≤ kmax.1 := hmin_le_max_fin
  have hm_lt_M : m.1 < M := lt_of_le_of_lt hmin_le_max hmax_lt
  let kPrev : Fin (M + 1) := ⟨m.1 - 1, by omega⟩
  let kNext : Fin (M + 1) := ⟨m.1 + 1, by omega⟩
  have hy_center :
      |x - ((m.1 : ℝ) + (29 : ℝ) / 16)| ≤ (12 : ℝ) / 11 := by
    simpa [x, m] using (Finset.mem_filter.mp hm_mem).2
  let y : ℝ := x - ((m.1 : ℝ) + (29 : ℝ) / 16)
  have hy : |y| ≤ (12 : ℝ) / 11 := by
    simpa [y] using hy_center
  have hneighbor :
      |x - ((kPrev.1 : ℝ) + (29 : ℝ) / 16)| ≤ (12 : ℝ) / 11 ∨
        |x - ((kNext.1 : ℝ) + (29 : ℝ) / 16)| ≤ (12 : ℝ) / 11 := by
    have hshift := abs_shift_neighbor_le_twelve_over_eleven y hy
    have hshift' :
        |(x - ((m.1 : ℝ) + (29 : ℝ) / 16)) + 1| ≤ (12 : ℝ) / 11 ∨
          |(x - ((m.1 : ℝ) + (29 : ℝ) / 16)) - 1| ≤ (12 : ℝ) / 11 := by
      simpa [y] using hshift
    rcases hshift' with hprev | hnext
    · left
      have hprev_rewrite :
          x - ((kPrev.1 : ℝ) + (29 : ℝ) / 16) =
            (x - ((m.1 : ℝ) + (29 : ℝ) / 16)) + 1 := by
        have hkPrev_val : (kPrev.1 : ℝ) = (m.1 : ℝ) - 1 := by
          have hm_pos : 0 < m.1 := by simpa [m] using hmin_pos
          have hm1_le : 1 ≤ m.1 := Nat.succ_le_of_lt hm_pos
          have hkPrev_nat : kPrev.1 = m.1 - 1 := by simp [kPrev]
          rw [hkPrev_nat, Nat.cast_sub hm1_le]
          norm_num
        linarith [hkPrev_val]
      simpa [hprev_rewrite] using hprev
    · right
      simpa [kNext, sub_eq_add_neg, add_assoc, add_left_comm, add_comm] using hnext
  have hsingleton_impossible : kmax.1 ≠ m.1 := by
    intro hEq
    have hprev_not_mem : kPrev ∉ operatorCenterCandidatesOrdered M j := by
      intro hkPrev_mem
      have hle : m ≤ kPrev := Finset.min'_le (operatorCenterCandidatesOrdered M j) kPrev hkPrev_mem
      have hkPrev_lt : kPrev.1 < m.1 := by
        have hm_pos : 0 < m.1 := by simpa [m] using hmin_pos
        have hkPrev_val : kPrev.1 = m.1 - 1 := by
          simp [kPrev]
        rw [hkPrev_val]
        omega
      exact (not_lt_of_ge hle) hkPrev_lt
    have hnext_not_mem : kNext ∉ operatorCenterCandidatesOrdered M j := by
      intro hkNext_mem
      have hle : kNext ≤ kmax := Finset.le_max' (operatorCenterCandidatesOrdered M j) kNext hkNext_mem
      have hkNext_gt : kmax.1 < kNext.1 := by
        have hm_lt : m.1 < m.1 + 1 := Nat.lt_succ_self m.1
        simpa [hEq, kNext] using hm_lt
      exact (not_lt_of_ge hle) hkNext_gt
    rcases hneighbor with hprev_abs | hnext_abs
    · have hkPrev_mem : kPrev ∈ operatorCenterCandidatesOrdered M j := by
        exact Finset.mem_filter.mpr ⟨by simp [kPrev], by simpa [x, kPrev] using hprev_abs⟩
      exact hprev_not_mem hkPrev_mem
    · have hkNext_mem : kNext ∈ operatorCenterCandidatesOrdered M j := by
        exact Finset.mem_filter.mpr ⟨by simp [kNext], by simpa [x, kNext] using hnext_abs⟩
      exact hnext_not_mem hkNext_mem
  have hkmax_not_le_m : ¬ kmax.1 ≤ m.1 := by
    intro hkmax_le_m
    exact hsingleton_impossible (le_antisymm hkmax_le_m hmin_le_max)
  exact lt_of_not_ge hkmax_not_le_m

/-- Monotonicity of the canonical minimal ordered-lane center selector:
as ordered eigenvalue index increases, the selected center index cannot increase. -/
theorem operatorEigenvalueOrderedCenterChoiceMin_antitone
    (M : ℕ) :
    Antitone (fun j : Fin (Fintype.card (Fin (M + 1))) =>
      (operatorEigenvalueOrderedCenterChoiceMin M j).1) := by
  intro i j hij
  let Ki : Fin (M + 1) := operatorEigenvalueOrderedCenterChoiceMin M i
  let Kj : Fin (M + 1) := operatorEigenvalueOrderedCenterChoiceMin M j
  let xi : ℝ := operatorEigenvaluesOrdered M i
  let xj : ℝ := operatorEigenvaluesOrdered M j
  let ci : ℝ := (Ki.1 : ℝ) + (29 : ℝ) / 16
  let cj : ℝ := (Kj.1 : ℝ) + (29 : ℝ) / 16
  have hi_abs : |xi - ci| ≤ (12 : ℝ) / 11 := by
    simpa [Ki, xi, ci] using operatorEigenvalueOrderedCenterChoiceMin_spec M i
  have hj_abs : |xj - cj| ≤ (12 : ℝ) / 11 := by
    simpa [Kj, xj, cj] using operatorEigenvalueOrderedCenterChoiceMin_spec M j
  have hanti := operatorEigenvaluesOrdered_antitone M
  have hxj_le_xi : xj ≤ xi := hanti hij
  have hci_ge_xi_minus : xi - (12 : ℝ) / 11 ≤ ci := by
    linarith [(abs_le.mp hi_abs).2]
  have hci_ge_xj_minus : xj - (12 : ℝ) / 11 ≤ ci := by
    linarith [hci_ge_xi_minus, hxj_le_xi]
  by_cases hci_le_xj_plus : ci ≤ xj + (12 : ℝ) / 11
  · have habs_xj_ci : |xj - ci| ≤ (12 : ℝ) / 11 := by
      refine abs_le.mpr ?_
      constructor
      · linarith [hci_le_xj_plus]
      · linarith [hci_ge_xj_minus]
    have hmem :
        Ki ∈ operatorCenterCandidatesOrdered M j := by
      refine Finset.mem_filter.mpr ?_
      refine ⟨Finset.mem_univ Ki, ?_⟩
      simpa [xj, ci] using habs_xj_ci
    have hle_fin : Kj ≤ Ki := Finset.min'_le (operatorCenterCandidatesOrdered M j) Ki hmem
    exact hle_fin
  · have hcj_le_xj_plus : cj ≤ xj + (12 : ℝ) / 11 := by
      linarith [(abs_le.mp hj_abs).1]
    have hci_gt_xj_plus : xj + (12 : ℝ) / 11 < ci := lt_of_not_ge hci_le_xj_plus
    have hcj_lt_ci : cj < ci := by
      linarith [hcj_le_xj_plus, hci_gt_xj_plus]
    have hkj_lt_ki_real : ((Kj.1 : ℝ)) < (Ki.1 : ℝ) := by
      dsimp [ci, cj] at hcj_lt_ci ⊢
      linarith
    have hkj_lt_ki : Kj.1 < Ki.1 := by
      exact_mod_cast hkj_lt_ki_real
    exact le_of_lt hkj_lt_ki

/-- Deterministic ordered-lane tie-break center index: reverse index in the
canonical `Fin (M+1) ↔ Fin(card)` reindexing. -/
noncomputable def operatorOrderedTieBreakCenter
    (M : ℕ) (j : Fin (Fintype.card (Fin (M + 1)))) : Fin (M + 1) :=
  let e : Fin (M + 1) ≃ Fin (Fintype.card (Fin (M + 1))) :=
    operatorEigenvaluesReindexToOrderedEquiv M
  ⟨M - (e.symm j).1, by
    exact Nat.lt_succ_of_le (Nat.sub_le M (e.symm j).1)⟩

/-- The ordered tie-break center selector is injective. -/
theorem operatorOrderedTieBreakCenter_injective
    (M : ℕ) :
    Function.Injective (operatorOrderedTieBreakCenter M) := by
  intro j1 j2 hEq
  let e : Fin (M + 1) ≃ Fin (Fintype.card (Fin (M + 1))) :=
    operatorEigenvaluesReindexToOrderedEquiv M
  have hNat :
      M - (e.symm j1).1 = M - (e.symm j2).1 := by
    simpa [operatorOrderedTieBreakCenter, e] using congrArg Fin.val hEq
  have hj1 : (e.symm j1).1 ≤ M := Nat.le_of_lt_succ (e.symm j1).2
  have hj2 : (e.symm j2).1 ≤ M := Nat.le_of_lt_succ (e.symm j2).2
  have hidx : (e.symm j1).1 = (e.symm j2).1 := by
    omega
  have hsymm : e.symm j1 = e.symm j2 := Fin.ext hidx
  exact e.symm.injective hsymm

/-- Offset form of the deterministic ordered tie-break center:
`tieBreak = maxChoice - offset`. -/
def operatorOrderedTieBreakOffset
    (M : ℕ) (j : Fin (Fintype.card (Fin (M + 1)))) : ℕ :=
  (operatorEigenvalueOrderedCenterChoiceMax M j).1
    - (operatorOrderedTieBreakCenter M j).1

/-- Three-line admissibility bridge:
if the tie-break offset is bounded by the candidate interval width, then the
tie-break index lies in the ordered candidate interval. -/
theorem operatorOrderedTieBreakCenter_mem_of_offset_le_width
    (M : ℕ) (j : Fin (Fintype.card (Fin (M + 1))))
    (hTieLeMax :
      (operatorOrderedTieBreakCenter M j).1
        ≤ (operatorEigenvalueOrderedCenterChoiceMax M j).1)
    (hOffsetWidth :
      operatorOrderedTieBreakOffset M j
        ≤ (operatorEigenvalueOrderedCenterChoiceMax M j).1
            - (operatorEigenvalueOrderedCenterChoiceMin M j).1) :
    operatorOrderedTieBreakCenter M j ∈ operatorCenterCandidatesOrdered M j := by
  have hTie_le_max :
      (operatorOrderedTieBreakCenter M j).1
        ≤ (operatorEigenvalueOrderedCenterChoiceMax M j).1 := hTieLeMax
  have hMin_mem :
      operatorEigenvalueOrderedCenterChoiceMin M j ∈ operatorCenterCandidatesOrdered M j :=
    operatorEigenvalueOrderedCenterChoiceMin_mem M j
  have hMax_mem :
      operatorEigenvalueOrderedCenterChoiceMax M j ∈ operatorCenterCandidatesOrdered M j :=
    operatorEigenvalueOrderedCenterChoiceMax_mem M j
  have hMin_le_max :
      (operatorEigenvalueOrderedCenterChoiceMin M j).1
        ≤ (operatorEigenvalueOrderedCenterChoiceMax M j).1 := by
    exact Finset.min'_le (operatorCenterCandidatesOrdered M j)
      (operatorEigenvalueOrderedCenterChoiceMax M j) hMax_mem
  have hMin_le_tie :
      (operatorEigenvalueOrderedCenterChoiceMin M j).1
        ≤ (operatorOrderedTieBreakCenter M j).1 := by
    have hOffsetWidth' :
        (operatorEigenvalueOrderedCenterChoiceMax M j).1
            - (operatorOrderedTieBreakCenter M j).1
          ≤ (operatorEigenvalueOrderedCenterChoiceMax M j).1
              - (operatorEigenvalueOrderedCenterChoiceMin M j).1 := by
      simpa [operatorOrderedTieBreakOffset] using hOffsetWidth
    omega
  exact (operatorCenterCandidatesOrdered_mem_iff_between_min_max M j
    (operatorOrderedTieBreakCenter M j)).2 ⟨hMin_le_tie, hTie_le_max⟩

/-- Clamped tie-break offset:
start from reverse-index raw offset and clamp to individual candidate width. -/
def operatorOrderedTieBreakOffsetClamped
    (M : ℕ) (j : Fin (Fintype.card (Fin (M + 1)))) : ℕ :=
  let e : Fin (M + 1) ≃ Fin (Fintype.card (Fin (M + 1))) :=
    operatorEigenvaluesReindexToOrderedEquiv M
  let kmax := (operatorEigenvalueOrderedCenterChoiceMax M j).1
  let kmin := (operatorEigenvalueOrderedCenterChoiceMin M j).1
  let width := kmax - kmin
  let rawOffset := kmax - (M - (e.symm j).1)
  Nat.min rawOffset width

/-- Clamped ordered tie-break center:
`max(j) - offset(j)` with clamped per-index offset. -/
noncomputable def operatorOrderedTieBreakCenterClamped
    (M : ℕ) (j : Fin (Fintype.card (Fin (M + 1)))) : Fin (M + 1) :=
  let kmax : ℕ := (operatorEigenvalueOrderedCenterChoiceMax M j).1
  let off : ℕ := operatorOrderedTieBreakOffsetClamped M j
  ⟨kmax - off, by
    have hkmax_lt : kmax < M + 1 := (operatorEigenvalueOrderedCenterChoiceMax M j).2
    exact lt_of_le_of_lt (Nat.sub_le _ _) hkmax_lt⟩

/-- `hTieLeMax` for the clamped tie-break is immediate by construction. -/
theorem operatorOrderedTieBreakCenterClamped_le_max
    (M : ℕ) (j : Fin (Fintype.card (Fin (M + 1)))) :
    (operatorOrderedTieBreakCenterClamped M j).1
      ≤ (operatorEigenvalueOrderedCenterChoiceMax M j).1 := by
  dsimp [operatorOrderedTieBreakCenterClamped]
  exact Nat.sub_le _ _

/-- `hOffsetWidth` for the clamped tie-break is immediate by construction. -/
theorem operatorOrderedTieBreakOffsetClamped_le_width
    (M : ℕ) (j : Fin (Fintype.card (Fin (M + 1)))) :
    operatorOrderedTieBreakOffsetClamped M j
      ≤ (operatorEigenvalueOrderedCenterChoiceMax M j).1
          - (operatorEigenvalueOrderedCenterChoiceMin M j).1 := by
  unfold operatorOrderedTieBreakOffsetClamped
  exact Nat.min_le_right _ _

/-- The clamped tie-break center is always admissible in the ordered candidate interval. -/
theorem operatorOrderedTieBreakCenterClamped_mem
    (M : ℕ) (j : Fin (Fintype.card (Fin (M + 1)))) :
    operatorOrderedTieBreakCenterClamped M j ∈ operatorCenterCandidatesOrdered M j := by
  have hTie_le_max :
      (operatorOrderedTieBreakCenterClamped M j).1
        ≤ (operatorEigenvalueOrderedCenterChoiceMax M j).1 :=
    operatorOrderedTieBreakCenterClamped_le_max M j
  have hMax_mem :
      operatorEigenvalueOrderedCenterChoiceMax M j ∈ operatorCenterCandidatesOrdered M j :=
    operatorEigenvalueOrderedCenterChoiceMax_mem M j
  have hMin_le_max :
      (operatorEigenvalueOrderedCenterChoiceMin M j).1
        ≤ (operatorEigenvalueOrderedCenterChoiceMax M j).1 := by
    exact Finset.min'_le (operatorCenterCandidatesOrdered M j)
      (operatorEigenvalueOrderedCenterChoiceMax M j) hMax_mem
  have hOffset_le_width :
      operatorOrderedTieBreakOffsetClamped M j
        ≤ (operatorEigenvalueOrderedCenterChoiceMax M j).1
            - (operatorEigenvalueOrderedCenterChoiceMin M j).1 :=
    operatorOrderedTieBreakOffsetClamped_le_width M j
  have hMin_le_tie :
      (operatorEigenvalueOrderedCenterChoiceMin M j).1
        ≤ (operatorOrderedTieBreakCenterClamped M j).1 := by
    unfold operatorOrderedTieBreakCenterClamped
    dsimp
    omega
  exact (operatorCenterCandidatesOrdered_mem_iff_between_min_max M j
    (operatorOrderedTieBreakCenterClamped M j)).2 ⟨hMin_le_tie, hTie_le_max⟩

/-- Clamped tie-break closure:
if the clamped ordered tie-break is injective at every finite level, then the
permutation-invariant center-gap contract follows immediately from admissibility. -/
theorem operatorCenterGapPermutationInvariant_of_orderedTieBreakClamped
    (hInj :
      ∀ M : ℕ,
        Function.Injective
          (fun j : Fin (Fintype.card (Fin (M + 1))) =>
            operatorOrderedTieBreakCenterClamped M j)) :
    OperatorCenterGapPermutationInvariant := by
  intro M
  classical
  let e : Fin (M + 1) ≃ Fin (Fintype.card (Fin (M + 1))) :=
    operatorEigenvaluesReindexToOrderedEquiv M
  let f : Fin (M + 1) → Fin (M + 1) := fun i =>
    operatorOrderedTieBreakCenterClamped M (e i)
  have hf_injective : Function.Injective f := by
    intro i1 i2 h
    exact e.injective ((hInj M) h)
  have hf_surjective : Function.Surjective f := (Finite.injective_iff_surjective).1 hf_injective
  let σ : Fin (M + 1) ≃ Fin (M + 1) := Equiv.ofBijective f ⟨hf_injective, hf_surjective⟩
  refine ⟨σ, ?_⟩
  intro i
  have hTie_i :
      operatorOrderedTieBreakCenterClamped M (e i) ∈
        operatorCenterCandidatesOrdered M (e i) :=
    operatorOrderedTieBreakCenterClamped_mem M (e i)
  have hTie_abs :
      |operatorEigenvaluesOrdered M (e i) -
        (((operatorOrderedTieBreakCenterClamped M (e i)).1 : ℝ) + (29 : ℝ) / 16)| ≤ (12 : ℝ) / 11 := by
    exact (Finset.mem_filter.mp hTie_i).2
  have hEqEig :
      operatorEigenvalues M i = operatorEigenvaluesOrdered M (e i) := by
    simpa [e] using operatorEigenvalues_eq_ordered_reindex M i
  have hσEq : σ i = operatorOrderedTieBreakCenterClamped M (e i) := rfl
  have hAbs :
      |operatorEigenvalues M i - (((σ i).1 : ℝ) + (29 : ℝ) / 16)| ≤ (12 : ℝ) / 11 := by
    simpa [hEqEig, hσEq] using hTie_abs
  exact hAbs

/-- Finite ordered-lane cardinality at level `M`. -/
private def operatorGreedyCard (M : ℕ) : ℕ := Fintype.card (Fin (M + 1))

mutual

/-- Recursive available-center set for greedy assignment:
start from all centers, then erase each previously chosen greedy center. -/
noncomputable def operatorGreedyAvailableNat (M : ℕ) : ℕ → Finset (Fin (M + 1))
  | 0 => Finset.univ
  | n + 1 =>
      let A := operatorGreedyAvailableNat M n
      if hn : n < operatorGreedyCard M then
        A.erase (operatorGreedyChoiceNat M n)
      else
        A

/-- Recursive greedy choice at ordered index `n`:
pick the maximal available candidate-center if possible; otherwise fallback to
the maximal currently available center (and finally `0` in the impossible empty
available-set branch). -/
noncomputable def operatorGreedyChoiceNat (M : ℕ) (n : ℕ) : Fin (M + 1) :=
  let A := operatorGreedyAvailableNat M n
  if hn : n < operatorGreedyCard M then
    let j : Fin (operatorGreedyCard M) := ⟨n, hn⟩
    let CandAvail := (operatorCenterCandidatesOrdered M j).filter (fun k => k ∈ A)
    if hCA : CandAvail.Nonempty then
      CandAvail.max' hCA
    else if hA : A.Nonempty then
      A.max' hA
    else
      ⟨0, Nat.succ_pos M⟩
  else
    ⟨0, Nat.succ_pos M⟩

end

/-- True greedy center selector on ordered indices (recursive available-set lane). -/
noncomputable def operatorGreedyCenter
    (M : ℕ) : Fin (operatorGreedyCard M) → Fin (M + 1) :=
  fun j => operatorGreedyChoiceNat M j.1

theorem operatorGreedyChoiceNat_mem_available_of_nonempty
    (M : ℕ) (n : ℕ)
    (hn : n < operatorGreedyCard M)
    (hA : (operatorGreedyAvailableNat M n).Nonempty) :
    operatorGreedyChoiceNat M n ∈ operatorGreedyAvailableNat M n := by
  classical
  have hn' : n < operatorGreedyCard M := hn
  by_cases hCA :
      (((operatorCenterCandidatesOrdered M ⟨n, hn'⟩).filter
        (fun k => k ∈ operatorGreedyAvailableNat M n)).Nonempty)
  · have hmem :
        ((operatorCenterCandidatesOrdered M ⟨n, hn'⟩).filter
          (fun k => k ∈ operatorGreedyAvailableNat M n)).max' hCA ∈
        ((operatorCenterCandidatesOrdered M ⟨n, hn'⟩).filter
          (fun k => k ∈ operatorGreedyAvailableNat M n)) :=
      Finset.max'_mem _ hCA
    have hchoice :
        operatorGreedyChoiceNat M n =
          ((operatorCenterCandidatesOrdered M ⟨n, hn'⟩).filter
            (fun k => k ∈ operatorGreedyAvailableNat M n)).max' hCA := by
      let P : Prop :=
        (((operatorCenterCandidatesOrdered M ⟨n, hn'⟩).filter
          (fun k => k ∈ operatorGreedyAvailableNat M n)).Nonempty)
      let F : P → Fin (M + 1) := fun h =>
        ((operatorCenterCandidatesOrdered M ⟨n, hn'⟩).filter
          (fun k => k ∈ operatorGreedyAvailableNat M n)).max' h
      let G : ¬P → Fin (M + 1) := fun _ =>
        if hA : (operatorGreedyAvailableNat M n).Nonempty then
          (operatorGreedyAvailableNat M n).max' hA
        else
          ⟨0, Nat.succ_pos M⟩
      have hdif : dite P F G = F hCA := by
        exact dif_pos hCA
      simpa [operatorGreedyChoiceNat, hn', P, F, G] using hdif
    exact hchoice ▸ (Finset.mem_filter.mp hmem).2
  · by_cases hA' : (operatorGreedyAvailableNat M n).Nonempty
    · have hmem : (operatorGreedyAvailableNat M n).max' hA' ∈ operatorGreedyAvailableNat M n :=
        Finset.max'_mem _ hA'
      have hchoice :
          operatorGreedyChoiceNat M n = (operatorGreedyAvailableNat M n).max' hA' := by
        let P : Prop :=
          (((operatorCenterCandidatesOrdered M ⟨n, hn'⟩).filter
            (fun k => k ∈ operatorGreedyAvailableNat M n)).Nonempty)
        let F : P → Fin (M + 1) := fun h =>
          ((operatorCenterCandidatesOrdered M ⟨n, hn'⟩).filter
            (fun k => k ∈ operatorGreedyAvailableNat M n)).max' h
        let G : ¬P → Fin (M + 1) := fun _ =>
          if hA : (operatorGreedyAvailableNat M n).Nonempty then
            (operatorGreedyAvailableNat M n).max' hA
          else
            ⟨0, Nat.succ_pos M⟩
        have hPneg : ¬P := by
          simpa [P] using hCA
        have hdif : dite P F G = G hPneg := by
          exact dif_neg hPneg
        have hG : G hPneg = (operatorGreedyAvailableNat M n).max' hA' := by
          simp [G, hA']
        simpa [operatorGreedyChoiceNat, hn', P, F, G] using hdif.trans hG
      exact hchoice ▸ hmem
    · exact False.elim (hA' hA)

theorem operatorGreedyAvailableNat_card
    (M : ℕ) :
    ∀ n : ℕ, n ≤ operatorGreedyCard M →
      (operatorGreedyAvailableNat M n).card = operatorGreedyCard M - n := by
  intro n
  induction' n with n ih
  · intro h0
    simp [operatorGreedyAvailableNat, operatorGreedyCard]
  · intro hn1
    have hnlt : n < operatorGreedyCard M := by omega
    have hnle : n ≤ operatorGreedyCard M := Nat.le_of_lt hnlt
    have hcardPrev : (operatorGreedyAvailableNat M n).card = operatorGreedyCard M - n := ih hnle
    have hApos : 0 < (operatorGreedyAvailableNat M n).card := by
      rw [hcardPrev]
      exact Nat.sub_pos_of_lt hnlt
    have hAne : (operatorGreedyAvailableNat M n).Nonempty := Finset.card_pos.mp hApos
    have hmem :
        operatorGreedyChoiceNat M n ∈ operatorGreedyAvailableNat M n :=
      operatorGreedyChoiceNat_mem_available_of_nonempty M n hnlt hAne
    have hcardErase :
        (operatorGreedyAvailableNat M (n + 1)).card
          = (operatorGreedyAvailableNat M n).card - 1 := by
      simp [operatorGreedyAvailableNat, hnlt, Finset.card_erase_of_mem, hmem]
    rw [hcardErase, hcardPrev]
    omega

theorem operatorGreedyAvailableNat_nonempty_of_lt
    (M : ℕ) (n : ℕ) (hn : n < operatorGreedyCard M) :
    (operatorGreedyAvailableNat M n).Nonempty := by
  have hcard := operatorGreedyAvailableNat_card M n (Nat.le_of_lt hn)
  have hpos : 0 < (operatorGreedyCard M - n) := Nat.sub_pos_of_lt hn
  have hcardPos : 0 < (operatorGreedyAvailableNat M n).card := by simpa [hcard] using hpos
  exact Finset.card_pos.mp hcardPos

theorem operatorGreedyChoiceNat_mem_available_of_lt
    (M : ℕ) (n : ℕ) (hn : n < operatorGreedyCard M) :
    operatorGreedyChoiceNat M n ∈ operatorGreedyAvailableNat M n := by
  exact operatorGreedyChoiceNat_mem_available_of_nonempty M n hn
    (operatorGreedyAvailableNat_nonempty_of_lt M n hn)

theorem operatorGreedyAvailableNat_succ_subset
    (M : ℕ) (n : ℕ) (hn : n < operatorGreedyCard M) :
    operatorGreedyAvailableNat M (n + 1) ⊆ operatorGreedyAvailableNat M n := by
  intro x hx
  have hx' : x ≠ operatorGreedyChoiceNat M n ∧ x ∈ operatorGreedyAvailableNat M n := by
    simpa [operatorGreedyAvailableNat, hn] using hx
  exact hx'.2

theorem operatorGreedyChoiceNat_not_mem_available_succ
    (M : ℕ) (n : ℕ) (hn : n < operatorGreedyCard M) :
    operatorGreedyChoiceNat M n ∉ operatorGreedyAvailableNat M (n + 1) := by
  simp [operatorGreedyAvailableNat, hn]

theorem operatorGreedyAvailableNat_antitone
    (M : ℕ) :
    ∀ {m n : ℕ}, m ≤ n → n ≤ operatorGreedyCard M →
      operatorGreedyAvailableNat M n ⊆ operatorGreedyAvailableNat M m := by
  intro m n hmn
  induction' hmn with n hmn ih
  · intro hn
    exact subset_rfl
  · intro hn
    have hnlt : n < operatorGreedyCard M := by
      exact lt_of_lt_of_le (Nat.lt_succ_self n) hn
    exact Set.Subset.trans
      (operatorGreedyAvailableNat_succ_subset M n hnlt)
      (ih (Nat.le_trans (Nat.le_succ n) hn))

theorem operatorGreedyChoiceNat_not_mem_available_of_lt
    (M : ℕ) (i j : ℕ)
    (hij : i < j) (hj : j ≤ operatorGreedyCard M) :
    operatorGreedyChoiceNat M i ∉ operatorGreedyAvailableNat M j := by
  have hi : i < operatorGreedyCard M := lt_of_lt_of_le hij hj
  have hnot1 :
      operatorGreedyChoiceNat M i ∉ operatorGreedyAvailableNat M (i + 1) := by
    simpa using operatorGreedyChoiceNat_not_mem_available_succ M i hi
  have hsub :
      operatorGreedyAvailableNat M j ⊆ operatorGreedyAvailableNat M (i + 1) :=
    operatorGreedyAvailableNat_antitone M (Nat.succ_le_of_lt hij) hj
  exact fun hmem => hnot1 (hsub hmem)

/-- If `x` is not available at step `j`, then some earlier greedy step chose `x`. -/
theorem operatorGreedyChoiceNat_exists_prev_of_not_mem_available
    (M : ℕ) :
    ∀ j : ℕ, ∀ x : Fin (M + 1),
      x ∉ operatorGreedyAvailableNat M j →
        ∃ i : ℕ, i < j ∧ operatorGreedyChoiceNat M i = x := by
  intro j
  induction' j with j ih
  · intro x hx
    have hxU : x ∉ (Finset.univ : Finset (Fin (M + 1))) := by
      simpa [operatorGreedyAvailableNat] using hx
    exact False.elim (hxU (Finset.mem_univ x))
  · intro x hx
    by_cases hjlt : j < operatorGreedyCard M
    · have hx' :
        x = operatorGreedyChoiceNat M j ∨ x ∉ operatorGreedyAvailableNat M j := by
        have hxImp :
            x ≠ operatorGreedyChoiceNat M j →
              x ∉ operatorGreedyAvailableNat M j := by
          simpa [operatorGreedyAvailableNat, hjlt, Finset.mem_erase] using hx
        by_cases hEq : x = operatorGreedyChoiceNat M j
        · exact Or.inl hEq
        · exact Or.inr (hxImp hEq)
      rcases hx' with hEq | hNotPrev
      · refine ⟨j, Nat.lt_succ_self j, hEq.symm⟩
      · rcases ih x hNotPrev with ⟨i, hi, hix⟩
        exact ⟨i, Nat.lt_trans hi (Nat.lt_succ_self j), hix⟩
    · have hNotPrev : x ∉ operatorGreedyAvailableNat M j := by
        simpa [operatorGreedyAvailableNat, hjlt] using hx
      rcases ih x hNotPrev with ⟨i, hi, hix⟩
      exact ⟨i, Nat.lt_trans hi (Nat.lt_succ_self j), hix⟩

/-- Generic availability/no-previous-choice equivalence:
an arbitrary center `x` is available at step `j` iff no earlier greedy step
chose `x`. -/
theorem operatorAvailable_mem_iff_no_prev
    (M : ℕ) (j : ℕ) (hj : j ≤ operatorGreedyCard M) (x : Fin (M + 1)) :
    x ∈ operatorGreedyAvailableNat M j ↔
      ∀ i : ℕ, i < j → operatorGreedyChoiceNat M i ≠ x := by
  constructor
  · intro hAvail i hij hEq
    have hnot :=
      operatorGreedyChoiceNat_not_mem_available_of_lt M i j hij hj
    exact hnot (hEq ▸ hAvail)
  · intro hNoPrev
    by_contra hnot
    rcases operatorGreedyChoiceNat_exists_prev_of_not_mem_available M j x hnot with
      ⟨i, hij, hEq⟩
    exact (hNoPrev i hij) hEq

/-- Single-gap reduction:
if no earlier step can choose the future ordered minimum center, then that
minimum center is available at its own step. -/
theorem operatorEigenvalueOrderedCenterChoiceMin_mem_available_of_noPremature
    (M : ℕ)
    (hNoPremature :
      ∀ i j : ℕ, ∀ hij : i < j, ∀ hj : j < operatorGreedyCard M,
        operatorGreedyChoiceNat M i ≠
          operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩) :
    ∀ j : ℕ, ∀ hj : j < operatorGreedyCard M,
      operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩ ∈
        operatorGreedyAvailableNat M j := by
  intro j hj
  by_contra hnot
  rcases operatorGreedyChoiceNat_exists_prev_of_not_mem_available M j
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩) hnot with
    ⟨i, hij, hEq⟩
  exact (hNoPremature i j hij hj) hEq

/-- Max-version of the single-gap reduction:
if no earlier step can choose the future ordered maximal center, then that
maximal center is available at its own step. -/
theorem operatorEigenvalueOrderedCenterChoiceMax_mem_available_of_noPremature
    (M : ℕ)
    (hNoPremature :
      ∀ i j : ℕ, ∀ hij : i < j, ∀ hj : j < operatorGreedyCard M,
        operatorGreedyChoiceNat M i ≠
          operatorEigenvalueOrderedCenterChoiceMax M ⟨j, hj⟩) :
    ∀ j : ℕ, ∀ hj : j < operatorGreedyCard M,
      operatorEigenvalueOrderedCenterChoiceMax M ⟨j, hj⟩ ∈
        operatorGreedyAvailableNat M j := by
  intro j hj
  by_contra hnot
  rcases operatorGreedyChoiceNat_exists_prev_of_not_mem_available M j
      (operatorEigenvalueOrderedCenterChoiceMax M ⟨j, hj⟩) hnot with
    ⟨i, hij, hEq⟩
  exact (hNoPremature i j hij hj) hEq

/-- Availability/no-premature equivalence for the future ordered minimum center. -/
theorem operatorMin_mem_available_iff_no_prev
    (M : ℕ) (j : ℕ) (hj : j < operatorGreedyCard M) :
    operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩ ∈ operatorGreedyAvailableNat M j ↔
      ∀ i : ℕ, i < j →
        operatorGreedyChoiceNat M i ≠ operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩ := by
  constructor
  · intro hAvail i hij hEq
    have hjle : j ≤ operatorGreedyCard M := Nat.le_of_lt hj
    have hnot := operatorGreedyChoiceNat_not_mem_available_of_lt M i j hij hjle
    exact hnot (hEq ▸ hAvail)
  · intro hNoPrev
    by_contra hnot
    rcases operatorGreedyChoiceNat_exists_prev_of_not_mem_available M j
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩) hnot with
      ⟨i, hij, hEq⟩
    exact (hNoPrev i hij) hEq

/-- No-premature-use follows directly from minimum availability. -/
theorem operatorGreedyChoiceNat_ne_future_min_of_min_available
    (M : ℕ)
    (hminAvail :
      ∀ j : ℕ, ∀ hj : j < operatorGreedyCard M,
        operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩ ∈
          operatorGreedyAvailableNat M j) :
    ∀ i j : ℕ, ∀ hij : i < j, ∀ hj : j < operatorGreedyCard M,
      operatorGreedyChoiceNat M i ≠
        operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩ := by
  intro i j hij hj
  exact (operatorMin_mem_available_iff_no_prev M j hj).1 (hminAvail j hj) i hij

/-- If the canonical ordered-lane minimum center is available at step `n`,
then the candidate∩available set at step `n` is nonempty. -/
theorem operatorGreedyCandAvail_nonempty_of_min_available
    (M : ℕ) (n : ℕ) (hn : n < operatorGreedyCard M)
    (hminAvail :
      operatorEigenvalueOrderedCenterChoiceMin M ⟨n, hn⟩ ∈ operatorGreedyAvailableNat M n) :
    (((operatorCenterCandidatesOrdered M ⟨n, hn⟩).filter
      (fun k => k ∈ operatorGreedyAvailableNat M n)).Nonempty) := by
  refine ⟨operatorEigenvalueOrderedCenterChoiceMin M ⟨n, hn⟩, ?_⟩
  exact Finset.mem_filter.mpr
    ⟨operatorEigenvalueOrderedCenterChoiceMin_mem M ⟨n, hn⟩, hminAvail⟩

/-- If the ordered-lane minimum center is available at step `n`, the greedy
choice index is at least that minimum index. -/
theorem operatorGreedyChoiceNat_ge_min_of_min_available
    (M : ℕ) (n : ℕ) (hn : n < operatorGreedyCard M)
    (hminAvail :
      operatorEigenvalueOrderedCenterChoiceMin M ⟨n, hn⟩ ∈ operatorGreedyAvailableNat M n) :
    (operatorEigenvalueOrderedCenterChoiceMin M ⟨n, hn⟩).1
      ≤ (operatorGreedyChoiceNat M n).1 := by
  let CandAvail : Finset (Fin (M + 1)) :=
    ((operatorCenterCandidatesOrdered M ⟨n, hn⟩).filter
      (fun k => k ∈ operatorGreedyAvailableNat M n))
  have hCA :
      CandAvail.Nonempty :=
    operatorGreedyCandAvail_nonempty_of_min_available M n hn hminAvail
  have hchoice :
      operatorGreedyChoiceNat M n = CandAvail.max' hCA := by
    let P : Prop := CandAvail.Nonempty
    let F : P → Fin (M + 1) := fun h => CandAvail.max' h
    let G : ¬P → Fin (M + 1) := fun _ =>
      if hA : (operatorGreedyAvailableNat M n).Nonempty then
        (operatorGreedyAvailableNat M n).max' hA
      else
        ⟨0, Nat.succ_pos M⟩
    have hdif : dite P F G = F hCA := dif_pos hCA
    simpa [operatorGreedyChoiceNat, hn, CandAvail, P, F, G] using hdif
  have hmemCandAvail : CandAvail.max' hCA ∈ CandAvail := Finset.max'_mem CandAvail hCA
  have hmemCand :
      CandAvail.max' hCA ∈ operatorCenterCandidatesOrdered M ⟨n, hn⟩ := by
    exact (Finset.mem_filter.mp hmemCandAvail).1
  have hbetween :
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨n, hn⟩).1
        ≤ (CandAvail.max' hCA).1 := by
    exact (operatorCenterCandidatesOrdered_mem_iff_between_min_max M ⟨n, hn⟩
      (CandAvail.max' hCA)).1 hmemCand |>.1
  have hchoice_val : (operatorGreedyChoiceNat M n).1 = (CandAvail.max' hCA).1 := by
    simpa [hchoice]
  simpa [hchoice_val] using hbetween

/-- Local-prefix reduction:
under prefix minimum-availability up to `j-1`, any equality
`greedy(i) = min(j)` must occur at the adjacent index `i = j-1`. -/
theorem operatorGreedyChoiceNat_eq_future_min_implies_adjacent_of_min_available_prefix
    (M : ℕ) (j : ℕ) (hj : j < operatorGreedyCard M)
    (hminAvailPrefix :
      ∀ k : ℕ, ∀ hk : k < j,
        operatorEigenvalueOrderedCenterChoiceMin M ⟨k, lt_trans hk hj⟩ ∈
          operatorGreedyAvailableNat M k) :
    ∀ i : ℕ, i < j →
      operatorGreedyChoiceNat M i =
        operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩ →
      i + 1 = j := by
  intro i hij hEq
  have hiCard : i < operatorGreedyCard M := lt_trans hij hj
  have hiSuccCard : i + 1 < operatorGreedyCard M := by omega
  have hmin_i_avail0 :=
    hminAvailPrefix i hij
  have hmin_i_avail :
      operatorEigenvalueOrderedCenterChoiceMin M ⟨i, hiCard⟩ ∈
        operatorGreedyAvailableNat M i := by
    simpa using hmin_i_avail0
  have hmin_i_le_greedy :
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨i, hiCard⟩).1 ≤
        (operatorGreedyChoiceNat M i).1 :=
    operatorGreedyChoiceNat_ge_min_of_min_available M i hiCard hmin_i_avail
  have hmin_antitone := operatorEigenvalueOrderedCenterChoiceMin_antitone M
  have hmin_j_le_i :
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1 ≤
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨i, hiCard⟩).1 := by
    exact hmin_antitone (show (⟨i, hiCard⟩ : Fin (operatorGreedyCard M)) ≤ ⟨j, hj⟩ by
      exact Nat.le_of_lt hij)
  have hmin_i_le_j :
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨i, hiCard⟩).1 ≤
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1 := by
    calc
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨i, hiCard⟩).1
          ≤ (operatorGreedyChoiceNat M i).1 := hmin_i_le_greedy
      _ = (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1 := by
          simpa [hEq]
  have hmin_eq :
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨i, hiCard⟩).1 =
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1 := by
    exact le_antisymm hmin_i_le_j hmin_j_le_i
  by_contra hNotAdj
  have hi1ltj : i + 1 < j := by omega
  have hmin_i1_avail0 :=
    hminAvailPrefix (i + 1) hi1ltj
  have hmin_i1_avail :
      operatorEigenvalueOrderedCenterChoiceMin M ⟨i + 1, hiSuccCard⟩ ∈
        operatorGreedyAvailableNat M (i + 1) := by
    simpa using hmin_i1_avail0
  have hmin_i1_le_i :
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨i + 1, hiSuccCard⟩).1 ≤
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨i, hiCard⟩).1 := by
    exact hmin_antitone
      (show (⟨i, hiCard⟩ : Fin (operatorGreedyCard M)) ≤ ⟨i + 1, hiSuccCard⟩ by
        exact Nat.le_succ i)
  have hmin_j_le_i1 :
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1 ≤
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨i + 1, hiSuccCard⟩).1 := by
    exact hmin_antitone
      (show (⟨i + 1, hiSuccCard⟩ : Fin (operatorGreedyCard M)) ≤ ⟨j, hj⟩ by
        exact Nat.le_of_lt hi1ltj)
  have hmin_i1_eq_j :
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨i + 1, hiSuccCard⟩).1 =
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1 := by
    apply le_antisymm
    · calc
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨i + 1, hiSuccCard⟩).1
            ≤ (operatorEigenvalueOrderedCenterChoiceMin M ⟨i, hiCard⟩).1 := hmin_i1_le_i
        _ = (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1 := hmin_eq
    · exact hmin_j_le_i1
  have hnot_i1 :
      operatorGreedyChoiceNat M i ∉ operatorGreedyAvailableNat M (i + 1) := by
    simpa using operatorGreedyChoiceNat_not_mem_available_succ M i hiCard
  have hEq_i1 :
      operatorGreedyChoiceNat M i =
        operatorEigenvalueOrderedCenterChoiceMin M ⟨i + 1, hiSuccCard⟩ := by
    apply Fin.ext
    simpa [hEq, hmin_i1_eq_j]
  exact hnot_i1 (hEq_i1 ▸ hmin_i1_avail)

/-- Prefix no-premature form:
under prefix minimum-availability up to `j-1`, any index strictly before `j-1`
cannot choose `min(j)`. -/
theorem operatorGreedyChoiceNat_ne_future_min_of_min_available_prefix
    (M : ℕ) (j : ℕ) (hj : j < operatorGreedyCard M)
    (hminAvailPrefix :
      ∀ k : ℕ, ∀ hk : k < j,
        operatorEigenvalueOrderedCenterChoiceMin M ⟨k, lt_trans hk hj⟩ ∈
          operatorGreedyAvailableNat M k) :
    ∀ i : ℕ, i + 1 < j →
      operatorGreedyChoiceNat M i ≠
        operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩ := by
  intro i hij hEq
  have hadj :
      i + 1 = j :=
    operatorGreedyChoiceNat_eq_future_min_implies_adjacent_of_min_available_prefix
      M j hj hminAvailPrefix i (lt_of_lt_of_le (Nat.lt_succ_self i) (Nat.le_of_lt hij)) hEq
  omega

/-- If `min(j)` is unavailable under prefix minimum-availability up to `j-1`,
then the immediate predecessor `j-1` must be the step that chose it. -/
theorem operatorMin_not_mem_available_implies_pred_choice_of_min_available_prefix
    (M : ℕ) (j : ℕ) (hj : j < operatorGreedyCard M)
    (hjpos : 0 < j)
    (hminAvailPrefix :
      ∀ k : ℕ, ∀ hk : k < j,
        operatorEigenvalueOrderedCenterChoiceMin M ⟨k, lt_trans hk hj⟩ ∈
          operatorGreedyAvailableNat M k)
    (hnot :
      operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩ ∉
        operatorGreedyAvailableNat M j) :
    operatorGreedyChoiceNat M (j - 1) =
      operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩ := by
  rcases operatorGreedyChoiceNat_exists_prev_of_not_mem_available M j
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩) hnot with
    ⟨i, hij, hEq⟩
  have hle : i + 1 ≤ j := Nat.succ_le_of_lt hij
  have hge : j ≤ i + 1 := by
    by_contra hlt
    have hlt' : i + 1 < j := Nat.lt_of_not_ge hlt
    exact (operatorGreedyChoiceNat_ne_future_min_of_min_available_prefix M j hj
      hminAvailPrefix i hlt' hEq).elim
  have hi : i = j - 1 := by omega
  simpa [hi] using hEq

/-- Same-min branch micro-kernel:
if `min(j-1)=min(j)` and there is any available candidate above that minimum at
step `j-1`, then the greedy predecessor choice cannot equal `min(j)`. -/
theorem operatorPredChoice_ne_futureMin_of_sameMin_and_aboveAvailable
    (M : ℕ) (j : ℕ) (hj : j < operatorGreedyCard M)
    (hjm1 : j - 1 < operatorGreedyCard M)
    (hSameMin :
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 =
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1)
    (hAboveAvail :
      ∃ c : Fin (M + 1),
        c ∈ operatorCenterCandidatesOrdered M ⟨j - 1, hjm1⟩ ∧
        c ∈ operatorGreedyAvailableNat M (j - 1) ∧
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 < c.1) :
    operatorGreedyChoiceNat M (j - 1) ≠
      operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩ := by
  rcases hAboveAvail with ⟨c, hcCand, hcAvail, hmin_lt_c⟩
  let CandAvail : Finset (Fin (M + 1)) :=
    ((operatorCenterCandidatesOrdered M ⟨j - 1, hjm1⟩).filter
      (fun k => k ∈ operatorGreedyAvailableNat M (j - 1)))
  have hcIn : c ∈ CandAvail := by
    exact Finset.mem_filter.mpr ⟨by simpa [hjm1] using hcCand, hcAvail⟩
  have hCA : CandAvail.Nonempty := ⟨c, hcIn⟩
  have hchoice :
      operatorGreedyChoiceNat M (j - 1) = CandAvail.max' hCA := by
    let P : Prop := CandAvail.Nonempty
    let F : P → Fin (M + 1) := fun h => CandAvail.max' h
    let G : ¬P → Fin (M + 1) := fun _ =>
      if hA : (operatorGreedyAvailableNat M (j - 1)).Nonempty then
        (operatorGreedyAvailableNat M (j - 1)).max' hA
      else
        ⟨0, Nat.succ_pos M⟩
    have hdif : dite P F G = F hCA := dif_pos hCA
    simpa [operatorGreedyChoiceNat, hjm1, CandAvail, P, F, G] using hdif
  have hcle :
      c.1 ≤ (CandAvail.max' hCA).1 := by
    exact Finset.le_max' CandAvail c hcIn
  have hgt_choice :
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 <
        (operatorGreedyChoiceNat M (j - 1)).1 := by
    have hchoiceVal : (operatorGreedyChoiceNat M (j - 1)).1 = (CandAvail.max' hCA).1 := by
      simpa [hchoice]
    have hgtMax :
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 <
          (CandAvail.max' hCA).1 := lt_of_lt_of_le hmin_lt_c hcle
    simpa [hchoiceVal] using hgtMax
  intro hEq
  have hleEq :
      (operatorGreedyChoiceNat M (j - 1)).1 =
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 := by
    calc
      (operatorGreedyChoiceNat M (j - 1)).1
          = (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1 := by simpa [hEq]
      _ = (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 := by
          simpa [hSameMin]
  exact (not_lt_of_ge (le_of_eq hleEq)) hgt_choice

/-- Two-case predecessor exclusion reduction:
the predecessor non-collision obligation follows from a single same-min branch
obligation; the strict-min branch is discharged structurally from
`greedy ≥ min` + antitone minima. -/
theorem operatorPredExcl_of_sameMinBranch
    (M : ℕ)
    (hSameMinBranch :
      ∀ j : ℕ, ∀ hj : j < operatorGreedyCard M, ∀ hjpos : 0 < j,
      ∀ hjm1 : j - 1 < operatorGreedyCard M,
        (∀ k : ℕ, ∀ hk : k < j,
          operatorEigenvalueOrderedCenterChoiceMin M ⟨k, lt_trans hk hj⟩ ∈
            operatorGreedyAvailableNat M k) →
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 =
          (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1 →
        operatorGreedyChoiceNat M (j - 1) ≠
          operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩) :
    ∀ j : ℕ, ∀ hj : j < operatorGreedyCard M, ∀ hjpos : 0 < j,
      (∀ k : ℕ, ∀ hk : k < j,
        operatorEigenvalueOrderedCenterChoiceMin M ⟨k, lt_trans hk hj⟩ ∈
          operatorGreedyAvailableNat M k) →
      operatorGreedyChoiceNat M (j - 1) ≠
        operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩ := by
  intro j hj hjpos hPrefix
  have hjm1Card : j - 1 < operatorGreedyCard M := by omega
  have hminPrevAvail :
      operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1Card⟩ ∈
        operatorGreedyAvailableNat M (j - 1) := by
    exact hPrefix (j - 1) (by omega)
  have hgreedy_ge_prevmin :
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1Card⟩).1 ≤
        (operatorGreedyChoiceNat M (j - 1)).1 :=
    operatorGreedyChoiceNat_ge_min_of_min_available M (j - 1) hjm1Card hminPrevAvail
  have hantiMin := operatorEigenvalueOrderedCenterChoiceMin_antitone M
  have hmin_j_le_prev :
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1 ≤
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1Card⟩).1 := by
    exact hantiMin (show (⟨j - 1, hjm1Card⟩ : Fin (operatorGreedyCard M)) ≤ ⟨j, hj⟩ by
      exact Nat.sub_le j 1)
  by_cases hEqMin :
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1Card⟩).1 =
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1
  · exact hSameMinBranch j hj hjpos hjm1Card hPrefix hEqMin
  · have hltMin :
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1 <
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1Card⟩).1 := by
      exact lt_of_le_of_ne hmin_j_le_prev (Ne.symm hEqMin)
    intro hEq
    have hgeViaEq :
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1Card⟩).1 ≤
          (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1 := by
      calc
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1Card⟩).1
            ≤ (operatorGreedyChoiceNat M (j - 1)).1 := hgreedy_ge_prevmin
        _ = (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1 := by
            simpa [hEq]
    exact (not_le_of_gt hltMin) hgeViaEq

/-- Predecessor exclusion from an explicit same-min witness condition:
if every same-min branch provides an available candidate strictly above that
minimum at step `j-1`, then predecessor exclusion follows globally. -/
theorem operatorPredExcl_of_sameMin_aboveAvailable
    (M : ℕ)
    (hAboveInSameMin :
      ∀ j : ℕ, ∀ hj : j < operatorGreedyCard M, ∀ hjpos : 0 < j,
      ∀ hjm1 : j - 1 < operatorGreedyCard M,
        (∀ k : ℕ, ∀ hk : k < j,
          operatorEigenvalueOrderedCenterChoiceMin M ⟨k, lt_trans hk hj⟩ ∈
            operatorGreedyAvailableNat M k) →
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 =
          (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1 →
        ∃ c : Fin (M + 1),
          c ∈ operatorCenterCandidatesOrdered M ⟨j - 1, hjm1⟩ ∧
          c ∈ operatorGreedyAvailableNat M (j - 1) ∧
          (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 < c.1) :
    ∀ j : ℕ, ∀ hj : j < operatorGreedyCard M, ∀ hjpos : 0 < j,
      (∀ k : ℕ, ∀ hk : k < j,
        operatorEigenvalueOrderedCenterChoiceMin M ⟨k, lt_trans hk hj⟩ ∈
          operatorGreedyAvailableNat M k) →
      operatorGreedyChoiceNat M (j - 1) ≠
        operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩ := by
  refine operatorPredExcl_of_sameMinBranch M ?_
  intro j hj hjpos hjm1 hPrefix hSameMin
  have hAbove :=
    hAboveInSameMin j hj hjpos hjm1 hPrefix hSameMin
  exact operatorPredChoice_ne_futureMin_of_sameMin_and_aboveAvailable
    M j hj hjm1 hSameMin hAbove

/-- Same-min witness from available predecessor-maximum:
if `min(j-1)=min(j)`, predecessor maximum is strictly above predecessor minimum,
and that maximum is available at step `j-1`, then the required above-min witness
exists for the same-min branch. -/
theorem operatorSameMin_aboveAvailable_of_maxAvailable
    (M : ℕ) (j : ℕ) (hj : j < operatorGreedyCard M)
    (hjm1 : j - 1 < operatorGreedyCard M)
    (hSameMin :
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 =
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1)
    (hMaxAbove :
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 <
        (operatorEigenvalueOrderedCenterChoiceMax M ⟨j - 1, hjm1⟩).1)
    (hMaxAvail :
      operatorEigenvalueOrderedCenterChoiceMax M ⟨j - 1, hjm1⟩ ∈
        operatorGreedyAvailableNat M (j - 1)) :
    ∃ c : Fin (M + 1),
      c ∈ operatorCenterCandidatesOrdered M ⟨j - 1, hjm1⟩ ∧
      c ∈ operatorGreedyAvailableNat M (j - 1) ∧
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 < c.1 := by
  refine ⟨operatorEigenvalueOrderedCenterChoiceMax M ⟨j - 1, hjm1⟩, ?_⟩
  refine ⟨operatorEigenvalueOrderedCenterChoiceMax_mem M ⟨j - 1, hjm1⟩, ?_⟩
  refine ⟨hMaxAvail, ?_⟩
  exact hMaxAbove

/-- Same-min witness from `min+1` availability:
if `min(j-1)=min(j)`, the predecessor candidate interval is non-singleton
(`max(j-1)>min(j-1)`), and `min(j-1)+1` is available at step `j-1`, then the
required above-min witness exists with `c := min(j-1)+1`. -/
theorem operatorPlusOne_mem_candidates_of_min_lt_max
    (M : ℕ) (j : ℕ) (hjm1 : j - 1 < operatorGreedyCard M)
    (hMaxAbove :
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 <
        (operatorEigenvalueOrderedCenterChoiceMax M ⟨j - 1, hjm1⟩).1) :
    ∃ c : Fin (M + 1),
      c.1 =
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 + 1 ∧
      c ∈ operatorCenterCandidatesOrdered M ⟨j - 1, hjm1⟩ := by
  let m := (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1
  let k := (operatorEigenvalueOrderedCenterChoiceMax M ⟨j - 1, hjm1⟩).1
  have hm1_le_k : m + 1 ≤ k := by omega
  have hk_lt : k < M + 1 :=
    (operatorEigenvalueOrderedCenterChoiceMax M ⟨j - 1, hjm1⟩).2
  let c : Fin (M + 1) := ⟨m + 1, lt_of_le_of_lt hm1_le_k hk_lt⟩
  refine ⟨c, ?_, ?_⟩
  · simp [c, m]
  · refine (operatorCenterCandidatesOrdered_mem_iff_between_min_max M ⟨j - 1, hjm1⟩ c).2 ?_
    constructor
    · simp [c, m]
    · simpa [c, m, k] using hm1_le_k

theorem operatorSameMin_aboveAvailable_of_plusOneAvail_and_maxAbove
    (M : ℕ) (j : ℕ) (hj : j < operatorGreedyCard M)
    (hjm1 : j - 1 < operatorGreedyCard M)
    (hSameMin :
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 =
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1)
    (hMaxAbove :
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 <
        (operatorEigenvalueOrderedCenterChoiceMax M ⟨j - 1, hjm1⟩).1)
    (hPlusOneAvail :
      ∃ c : Fin (M + 1),
        c.1 =
          (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 + 1 ∧
        c ∈ operatorGreedyAvailableNat M (j - 1)) :
    ∃ c : Fin (M + 1),
      c ∈ operatorCenterCandidatesOrdered M ⟨j - 1, hjm1⟩ ∧
      c ∈ operatorGreedyAvailableNat M (j - 1) ∧
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 < c.1 := by
  rcases hPlusOneAvail with ⟨c, hcEq, hcAvail⟩
  have hCand :
      c ∈ operatorCenterCandidatesOrdered M ⟨j - 1, hjm1⟩ := by
    refine (operatorCenterCandidatesOrdered_mem_iff_between_min_max M ⟨j - 1, hjm1⟩ c).2 ?_
    constructor
    · omega
    · have hm1_le_k :
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 + 1
          ≤ (operatorEigenvalueOrderedCenterChoiceMax M ⟨j - 1, hjm1⟩).1 := by
        omega
      omega
  have hMinLt :
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 < c.1 := by
    omega
  exact ⟨c, hCand, hcAvail, hMinLt⟩

theorem operatorSameMin_aboveAvailable_of_plusOneAvailable
    (M : ℕ) (j : ℕ) (hj : j < operatorGreedyCard M)
    (hjm1 : j - 1 < operatorGreedyCard M)
    (hSameMin :
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 =
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1)
    (hPlusOneCandAvail :
      ∃ c : Fin (M + 1),
        c.1 =
          (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 + 1 ∧
        c ∈ operatorCenterCandidatesOrdered M ⟨j - 1, hjm1⟩ ∧
        c ∈ operatorGreedyAvailableNat M (j - 1)) :
    ∃ c : Fin (M + 1),
      c ∈ operatorCenterCandidatesOrdered M ⟨j - 1, hjm1⟩ ∧
      c ∈ operatorGreedyAvailableNat M (j - 1) ∧
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 < c.1 := by
  rcases hPlusOneCandAvail with ⟨c, hcEq, hcCand, hcAvail⟩
  have hmin_lt_c :
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 < c.1 := by
    omega
  refine ⟨c, hcCand, hcAvail, hmin_lt_c⟩

/-- Same-min predecessor exclusion from `min+1` availability:
this is the direct 5-step route
`min+1 ∈ candidates ∩ available ⇒ greedy(j-1) > min(j)`. -/
theorem operatorPredExcl_of_sameMin_plusOneAvailable
    (M : ℕ)
    (hPlusOneCandAvailInSameMin :
      ∀ j : ℕ, ∀ hj : j < operatorGreedyCard M, ∀ hjpos : 0 < j,
      ∀ hjm1 : j - 1 < operatorGreedyCard M,
        (∀ k : ℕ, ∀ hk : k < j,
          operatorEigenvalueOrderedCenterChoiceMin M ⟨k, lt_trans hk hj⟩ ∈
            operatorGreedyAvailableNat M k) →
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 =
          (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1 →
        ∃ c : Fin (M + 1),
          c.1 =
            (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 + 1 ∧
          c ∈ operatorCenterCandidatesOrdered M ⟨j - 1, hjm1⟩ ∧
          c ∈ operatorGreedyAvailableNat M (j - 1)) :
    ∀ j : ℕ, ∀ hj : j < operatorGreedyCard M, ∀ hjpos : 0 < j,
      (∀ k : ℕ, ∀ hk : k < j,
        operatorEigenvalueOrderedCenterChoiceMin M ⟨k, lt_trans hk hj⟩ ∈
          operatorGreedyAvailableNat M k) →
      operatorGreedyChoiceNat M (j - 1) ≠
        operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩ := by
  refine operatorPredExcl_of_sameMin_aboveAvailable M ?_
  intro j hj hjpos hjm1 hPrefix hSameMin
  have hPlusOneCandAvail :
      ∃ c : Fin (M + 1),
        c.1 =
          (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 + 1 ∧
        c ∈ operatorCenterCandidatesOrdered M ⟨j - 1, hjm1⟩ ∧
        c ∈ operatorGreedyAvailableNat M (j - 1) :=
    hPlusOneCandAvailInSameMin j hj hjpos hjm1 hPrefix hSameMin
  exact operatorSameMin_aboveAvailable_of_plusOneAvailable
    M j hj hjm1 hSameMin hPlusOneCandAvail

/-- Same-min predecessor exclusion from two reduced obligations:
1) structural non-singleton branch (`max(j-1) > min(j-1)`), and
2) availability of `min(j-1)+1` at step `j-1`.

Candidate admissibility of `min+1` is derived internally from (1). -/
theorem operatorPredExcl_of_sameMin_plusOneAvail_and_maxAbove
    (M : ℕ)
    (hMaxAboveInSameMin :
      ∀ j : ℕ, ∀ hj : j < operatorGreedyCard M, ∀ hjpos : 0 < j,
      ∀ hjm1 : j - 1 < operatorGreedyCard M,
        (∀ k : ℕ, ∀ hk : k < j,
          operatorEigenvalueOrderedCenterChoiceMin M ⟨k, lt_trans hk hj⟩ ∈
            operatorGreedyAvailableNat M k) →
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 =
          (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1 →
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 <
          (operatorEigenvalueOrderedCenterChoiceMax M ⟨j - 1, hjm1⟩).1)
    (hPlusOneAvailInSameMin :
      ∀ j : ℕ, ∀ hj : j < operatorGreedyCard M, ∀ hjpos : 0 < j,
      ∀ hjm1 : j - 1 < operatorGreedyCard M,
        (∀ k : ℕ, ∀ hk : k < j,
          operatorEigenvalueOrderedCenterChoiceMin M ⟨k, lt_trans hk hj⟩ ∈
            operatorGreedyAvailableNat M k) →
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 =
          (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1 →
        ∃ c : Fin (M + 1),
          c.1 =
            (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 + 1 ∧
          c ∈ operatorGreedyAvailableNat M (j - 1)) :
    ∀ j : ℕ, ∀ hj : j < operatorGreedyCard M, ∀ hjpos : 0 < j,
      (∀ k : ℕ, ∀ hk : k < j,
        operatorEigenvalueOrderedCenterChoiceMin M ⟨k, lt_trans hk hj⟩ ∈
          operatorGreedyAvailableNat M k) →
      operatorGreedyChoiceNat M (j - 1) ≠
        operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩ := by
  refine operatorPredExcl_of_sameMin_aboveAvailable M ?_
  intro j hj hjpos hjm1 hPrefix hSameMin
  have hMaxAbove :
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 <
        (operatorEigenvalueOrderedCenterChoiceMax M ⟨j - 1, hjm1⟩).1 :=
    hMaxAboveInSameMin j hj hjpos hjm1 hPrefix hSameMin
  have hPlusOneAvail :
      ∃ c : Fin (M + 1),
        c.1 =
          (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 + 1 ∧
        c ∈ operatorGreedyAvailableNat M (j - 1) :=
    hPlusOneAvailInSameMin j hj hjpos hjm1 hPrefix hSameMin
  exact operatorSameMin_aboveAvailable_of_plusOneAvail_and_maxAbove
    M j hj hjm1 hSameMin hMaxAbove hPlusOneAvail

/-- From a `min+1` candidate witness in the same-min branch, `max(j-1) > min(j-1)`
follows immediately by interval bounds. -/
theorem operatorMaxAboveInSameMin_of_plusOneCandidate
    (M : ℕ) (j : ℕ) (hj : j < operatorGreedyCard M)
    (hjm1 : j - 1 < operatorGreedyCard M)
    (hPlusOneCand :
      ∃ c : Fin (M + 1),
        c.1 =
          (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 + 1 ∧
        c ∈ operatorCenterCandidatesOrdered M ⟨j - 1, hjm1⟩) :
    (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 <
      (operatorEigenvalueOrderedCenterChoiceMax M ⟨j - 1, hjm1⟩).1 := by
  rcases hPlusOneCand with ⟨c, hcEq, hcCand⟩
  have hcLeMax :
      c.1 ≤ (operatorEigenvalueOrderedCenterChoiceMax M ⟨j - 1, hjm1⟩).1 :=
    (operatorCenterCandidatesOrdered_mem_iff_between_min_max M ⟨j - 1, hjm1⟩ c).1 hcCand |>.2
  omega

/-- Canonical `min+1` witness reduction:
under `max(j-1)>min(j-1)`, there is a canonical center index `c = min+1`, and
its availability at step `j-1` is equivalent to a pure no-previous-choice
condition. This isolates the remaining seam to predecessor exclusion for that
single explicit center. -/
theorem operatorPlusOneAvail_iff_no_prev_choice_of_min_lt_max
    (M : ℕ) (j : ℕ) (hjm1 : j - 1 < operatorGreedyCard M)
    (hMaxAbove :
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 <
        (operatorEigenvalueOrderedCenterChoiceMax M ⟨j - 1, hjm1⟩).1) :
    ∃ c : Fin (M + 1),
      c.1 =
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 + 1 ∧
      (c ∈ operatorGreedyAvailableNat M (j - 1) ↔
        ∀ i : ℕ, i < (j - 1) → operatorGreedyChoiceNat M i ≠ c) := by
  rcases operatorPlusOne_mem_candidates_of_min_lt_max M j hjm1 hMaxAbove with
    ⟨c, hcEq, _hcCand⟩
  refine ⟨c, hcEq, ?_⟩
  exact operatorAvailable_mem_iff_no_prev M (j - 1) (Nat.le_of_lt hjm1) c

/-- Clean same-min predecessor exclusion from predecessor-max availability.
This isolates the core six-line argument and leaves endpoint/interior arithmetic
to upstream lemmas. -/
theorem operatorPredExcl_of_sameMin_maxAvailable
    (M : ℕ)
    (hMaxAboveInSameMin :
      ∀ j : ℕ, ∀ hj : j < operatorGreedyCard M, ∀ hjpos : 0 < j,
      ∀ hjm1 : j - 1 < operatorGreedyCard M,
        (∀ k : ℕ, ∀ hk : k < j,
          operatorEigenvalueOrderedCenterChoiceMin M ⟨k, lt_trans hk hj⟩ ∈
            operatorGreedyAvailableNat M k) →
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 =
          (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1 →
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 <
          (operatorEigenvalueOrderedCenterChoiceMax M ⟨j - 1, hjm1⟩).1)
    (hMaxAvailInSameMin :
      ∀ j : ℕ, ∀ hj : j < operatorGreedyCard M, ∀ hjpos : 0 < j,
      ∀ hjm1 : j - 1 < operatorGreedyCard M,
        (∀ k : ℕ, ∀ hk : k < j,
          operatorEigenvalueOrderedCenterChoiceMin M ⟨k, lt_trans hk hj⟩ ∈
            operatorGreedyAvailableNat M k) →
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 =
          (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1 →
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 <
          (operatorEigenvalueOrderedCenterChoiceMax M ⟨j - 1, hjm1⟩).1 →
        operatorEigenvalueOrderedCenterChoiceMax M ⟨j - 1, hjm1⟩ ∈
          operatorGreedyAvailableNat M (j - 1)) :
    ∀ j : ℕ, ∀ hj : j < operatorGreedyCard M, ∀ hjpos : 0 < j,
      (∀ k : ℕ, ∀ hk : k < j,
        operatorEigenvalueOrderedCenterChoiceMin M ⟨k, lt_trans hk hj⟩ ∈
          operatorGreedyAvailableNat M k) →
      operatorGreedyChoiceNat M (j - 1) ≠
        operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩ := by
  refine operatorPredExcl_of_sameMin_aboveAvailable M ?_
  intro j hj hjpos hjm1 hPrefix hSameMin
  have hMaxAbove :
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 <
        (operatorEigenvalueOrderedCenterChoiceMax M ⟨j - 1, hjm1⟩).1 :=
    hMaxAboveInSameMin j hj hjpos hjm1 hPrefix hSameMin
  have hMaxAvail :
      operatorEigenvalueOrderedCenterChoiceMax M ⟨j - 1, hjm1⟩ ∈
        operatorGreedyAvailableNat M (j - 1) :=
    hMaxAvailInSameMin j hj hjpos hjm1 hPrefix hSameMin hMaxAbove
  exact operatorSameMin_aboveAvailable_of_maxAvailable
    M j hj hjm1 hSameMin hMaxAbove hMaxAvail

/-- Same-min predecessor exclusion with max-availability discharged from a
global no-premature-max hypothesis. This reuses the max-availability reduction
pattern directly and leaves only the structural max-above-min branch. -/
theorem operatorPredExcl_of_sameMin_maxAvailable_of_noPrematureMax
    (M : ℕ)
    (hMaxAboveInSameMin :
      ∀ j : ℕ, ∀ hj : j < operatorGreedyCard M, ∀ hjpos : 0 < j,
      ∀ hjm1 : j - 1 < operatorGreedyCard M,
        (∀ k : ℕ, ∀ hk : k < j,
          operatorEigenvalueOrderedCenterChoiceMin M ⟨k, lt_trans hk hj⟩ ∈
            operatorGreedyAvailableNat M k) →
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 =
          (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1 →
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 <
          (operatorEigenvalueOrderedCenterChoiceMax M ⟨j - 1, hjm1⟩).1)
    (hNoPrematureMax :
      ∀ i j : ℕ, ∀ hij : i < j, ∀ hj : j < operatorGreedyCard M,
        operatorGreedyChoiceNat M i ≠
          operatorEigenvalueOrderedCenterChoiceMax M ⟨j, hj⟩) :
    ∀ j : ℕ, ∀ hj : j < operatorGreedyCard M, ∀ hjpos : 0 < j,
      (∀ k : ℕ, ∀ hk : k < j,
        operatorEigenvalueOrderedCenterChoiceMin M ⟨k, lt_trans hk hj⟩ ∈
          operatorGreedyAvailableNat M k) →
      operatorGreedyChoiceNat M (j - 1) ≠
        operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩ := by
  refine operatorPredExcl_of_sameMin_maxAvailable M hMaxAboveInSameMin ?_
  intro j hj hjpos hjm1 _ hSameMin hMaxAbove
  exact operatorEigenvalueOrderedCenterChoiceMax_mem_available_of_noPremature
    M hNoPrematureMax (j - 1) hjm1

/-- Same-min branch max-above-min from the existing interior non-collapse theorem. -/
theorem operatorMaxAboveInSameMin_of_interior
    (M : ℕ) (j : ℕ) (hj : j < operatorGreedyCard M) (hjpos : 0 < j)
    (hjm1 : j - 1 < operatorGreedyCard M)
    (hminPos :
      0 < (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1)
    (hmaxLt :
      (operatorEigenvalueOrderedCenterChoiceMax M ⟨j - 1, hjm1⟩).1 < M) :
    (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 <
      (operatorEigenvalueOrderedCenterChoiceMax M ⟨j - 1, hjm1⟩).1 := by
  have h := operatorCenterCandidatesOrdered_min_lt_max_of_interior
    M ⟨j - 1, hjm1⟩ hminPos hmaxLt
  simpa using h


/-- Strong-induction closure:
if the predecessor step never chooses `min(j)` under the prefix min-availability
hypothesis, then `min(j)` is available at every step. -/
theorem operatorEigenvalueOrderedCenterChoiceMin_mem_available_of_pred_exclusion
    (M : ℕ)
    (hPredExcl :
      ∀ j : ℕ, ∀ hj : j < operatorGreedyCard M, ∀ hjpos : 0 < j,
        (∀ k : ℕ, ∀ hk : k < j,
          operatorEigenvalueOrderedCenterChoiceMin M ⟨k, lt_trans hk hj⟩ ∈
            operatorGreedyAvailableNat M k) →
        operatorGreedyChoiceNat M (j - 1) ≠
          operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩) :
    ∀ j : ℕ, ∀ hj : j < operatorGreedyCard M,
      operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩ ∈
        operatorGreedyAvailableNat M j := by
  intro j
  refine Nat.strongRecOn j ?_
  intro j ih
  intro hj
  by_cases hj0 : j = 0
  · subst hj0
    simpa [operatorGreedyAvailableNat] using
      (Finset.mem_univ (operatorEigenvalueOrderedCenterChoiceMin M ⟨0, hj⟩))
  · have hjpos : 0 < j := Nat.pos_of_ne_zero hj0
    have hPrefix :
        ∀ k : ℕ, ∀ hk : k < j,
          operatorEigenvalueOrderedCenterChoiceMin M ⟨k, lt_trans hk hj⟩ ∈
            operatorGreedyAvailableNat M k := by
      intro k hk
      exact ih k hk (lt_trans hk hj)
    by_contra hnot
    have hPredChoice :
        operatorGreedyChoiceNat M (j - 1) =
          operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩ :=
      operatorMin_not_mem_available_implies_pred_choice_of_min_available_prefix
        M j hj hjpos hPrefix hnot
    exact (hPredExcl j hj hjpos hPrefix) hPredChoice

/-- Packaged `hCandAvail` from the single target property
`min(j) ∈ available(j)` at each ordered step. -/
theorem operatorGreedy_hCandAvail_of_min_available
    (M : ℕ)
    (hminAvail :
      ∀ n : ℕ, ∀ hn : n < operatorGreedyCard M,
        operatorEigenvalueOrderedCenterChoiceMin M ⟨n, hn⟩ ∈
          operatorGreedyAvailableNat M n) :
    ∀ n : ℕ, ∀ hn : n < operatorGreedyCard M,
      (((operatorCenterCandidatesOrdered M ⟨n, hn⟩).filter
        (fun k => k ∈ operatorGreedyAvailableNat M n)).Nonempty) := by
  intro n hn
  exact operatorGreedyCandAvail_nonempty_of_min_available M n hn (hminAvail n hn)

/-- Global greedy closure from predecessor exclusion:
the candidate-available branch is nonempty at every ordered step. -/
theorem operatorGreedy_hCandAvail_of_pred_exclusion
    (M : ℕ)
    (hPredExcl :
      ∀ j : ℕ, ∀ hj : j < operatorGreedyCard M, ∀ hjpos : 0 < j,
        (∀ k : ℕ, ∀ hk : k < j,
          operatorEigenvalueOrderedCenterChoiceMin M ⟨k, lt_trans hk hj⟩ ∈
            operatorGreedyAvailableNat M k) →
        operatorGreedyChoiceNat M (j - 1) ≠
          operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩) :
    ∀ n : ℕ, ∀ hn : n < operatorGreedyCard M,
      (((operatorCenterCandidatesOrdered M ⟨n, hn⟩).filter
        (fun k => k ∈ operatorGreedyAvailableNat M n)).Nonempty) := by
  exact operatorGreedy_hCandAvail_of_min_available M
    (operatorEigenvalueOrderedCenterChoiceMin_mem_available_of_pred_exclusion M hPredExcl)

/-- Two-stage reduction:
if no earlier step can choose the future ordered minimum center, the greedy
candidate∩available branch is nonempty at every ordered step. -/
theorem operatorGreedy_hCandAvail_of_noPremature
    (M : ℕ)
    (hNoPremature :
      ∀ i j : ℕ, ∀ hij : i < j, ∀ hj : j < operatorGreedyCard M,
        operatorGreedyChoiceNat M i ≠
          operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩) :
    ∀ n : ℕ, ∀ hn : n < operatorGreedyCard M,
      (((operatorCenterCandidatesOrdered M ⟨n, hn⟩).filter
        (fun k => k ∈ operatorGreedyAvailableNat M n)).Nonempty) := by
  exact operatorGreedy_hCandAvail_of_min_available M
    (operatorEigenvalueOrderedCenterChoiceMin_mem_available_of_noPremature M hNoPremature)

/-- Greedy-center injectivity is by construction:
each step chooses from the current available set and is erased before the next step. -/
theorem operatorGreedyCenter_injective
    (M : ℕ) :
    Function.Injective (operatorGreedyCenter M) := by
  intro j1 j2 hEq
  by_cases hlt : j1.1 < j2.1
  · have hj2le : j2.1 ≤ operatorGreedyCard M := Nat.le_of_lt j2.2
    have hnot :
        operatorGreedyChoiceNat M j1.1 ∉ operatorGreedyAvailableNat M j2.1 :=
      operatorGreedyChoiceNat_not_mem_available_of_lt M j1.1 j2.1 hlt hj2le
    have hmem :
        operatorGreedyChoiceNat M j2.1 ∈ operatorGreedyAvailableNat M j2.1 :=
      operatorGreedyChoiceNat_mem_available_of_lt M j2.1 j2.2
    have hEqNat : operatorGreedyChoiceNat M j1.1 = operatorGreedyChoiceNat M j2.1 := by
      simpa [operatorGreedyCenter] using hEq
    have hmem' :
        operatorGreedyChoiceNat M j1.1 ∈ operatorGreedyAvailableNat M j2.1 := by
      simpa [hEqNat] using hmem
    exact (hnot hmem').elim
  · have hge : j2.1 ≤ j1.1 := Nat.le_of_not_lt hlt
    by_cases hlt' : j2.1 < j1.1
    · have hj1le : j1.1 ≤ operatorGreedyCard M := Nat.le_of_lt j1.2
      have hnot :
          operatorGreedyChoiceNat M j2.1 ∉ operatorGreedyAvailableNat M j1.1 :=
        operatorGreedyChoiceNat_not_mem_available_of_lt M j2.1 j1.1 hlt' hj1le
      have hmem :
          operatorGreedyChoiceNat M j1.1 ∈ operatorGreedyAvailableNat M j1.1 :=
        operatorGreedyChoiceNat_mem_available_of_lt M j1.1 j1.2
      have hEqNat : operatorGreedyChoiceNat M j1.1 = operatorGreedyChoiceNat M j2.1 := by
        simpa [operatorGreedyCenter] using hEq
      have hmem' :
          operatorGreedyChoiceNat M j2.1 ∈ operatorGreedyAvailableNat M j1.1 := by
        simpa [hEqNat] using hmem
      exact (hnot hmem').elim
    · have hle12 : j1.1 ≤ j2.1 := Nat.le_of_not_lt hlt'
      have hNat : j1.1 = j2.1 := Nat.le_antisymm hle12 hge
      exact Fin.ext hNat

/-- Greedy-center admissibility in ordered candidate sets, provided the
candidate∩available branch is nonempty at each ordered step. -/
theorem operatorGreedyCenter_mem_candidates
    (M : ℕ)
    (hCandAvail :
      ∀ n : ℕ, ∀ hn : n < operatorGreedyCard M,
        (((operatorCenterCandidatesOrdered M ⟨n, hn⟩).filter
          (fun k => k ∈ operatorGreedyAvailableNat M n)).Nonempty))
    (j : Fin (operatorGreedyCard M)) :
    operatorGreedyCenter M j ∈ operatorCenterCandidatesOrdered M j := by
  classical
  have hn : j.1 < operatorGreedyCard M := j.2
  unfold operatorGreedyCenter
  unfold operatorGreedyChoiceNat
  simp [hn]
  let A : Finset (Fin (M + 1)) := operatorGreedyAvailableNat M j.1
  let CandAvail : Finset (Fin (M + 1)) := (operatorCenterCandidatesOrdered M j).filter (fun k => k ∈ A)
  have hCA : CandAvail.Nonempty := by
    simpa [A, CandAvail] using hCandAvail j.1 hn
  have hmem : CandAvail.max' hCA ∈ CandAvail := Finset.max'_mem CandAvail hCA
  have hmemC : CandAvail.max' hCA ∈ operatorCenterCandidatesOrdered M j :=
    (Finset.mem_filter.mp hmem).1
  simpa [A, CandAvail, hn, hCA] using hmemC

/-- Greedy-center closure to the permutation-invariant center-gap contract,
assuming candidate∩available nonempty at each greedy step. -/
theorem operatorCenterGapPermutationInvariant_of_operatorGreedyCenter
    (hCandAvail :
      ∀ M : ℕ, ∀ n : ℕ, ∀ hn : n < operatorGreedyCard M,
        (((operatorCenterCandidatesOrdered M ⟨n, hn⟩).filter
          (fun k => k ∈ operatorGreedyAvailableNat M n)).Nonempty)) :
    OperatorCenterGapPermutationInvariant := by
  intro M
  classical
  let e : Fin (M + 1) ≃ Fin (operatorGreedyCard M) :=
    operatorEigenvaluesReindexToOrderedEquiv M
  let f : Fin (M + 1) → Fin (M + 1) := fun i => operatorGreedyCenter M (e i)
  have hf_injective : Function.Injective f := by
    intro i1 i2 h
    exact e.injective (operatorGreedyCenter_injective M h)
  have hf_surjective : Function.Surjective f := (Finite.injective_iff_surjective).1 hf_injective
  let σ : Fin (M + 1) ≃ Fin (M + 1) := Equiv.ofBijective f ⟨hf_injective, hf_surjective⟩
  refine ⟨σ, ?_⟩
  intro i
  have hTie_i :
      operatorGreedyCenter M (e i) ∈ operatorCenterCandidatesOrdered M (e i) :=
    operatorGreedyCenter_mem_candidates M (hCandAvail M) (e i)
  have hTie_abs :
      |operatorEigenvaluesOrdered M (e i) -
        (((operatorGreedyCenter M (e i)).1 : ℝ) + (29 : ℝ) / 16)| ≤ (12 : ℝ) / 11 := by
    exact (Finset.mem_filter.mp hTie_i).2
  have hEqEig :
      operatorEigenvalues M i = operatorEigenvaluesOrdered M (e i) := by
    simpa [e] using operatorEigenvalues_eq_ordered_reindex M i
  have hσEq : σ i = operatorGreedyCenter M (e i) := rfl
  have hAbs :
      |operatorEigenvalues M i - (((σ i).1 : ℝ) + (29 : ℝ) / 16)| ≤ (12 : ℝ) / 11 := by
    simpa [hEqEig, hσEq] using hTie_abs
  exact hAbs

/-- Reduction: if we prove `min(j)` is always available at step `j`, the greedy
route closes `OperatorCenterGapPermutationInvariant` directly. -/
theorem operatorCenterGapPermutationInvariant_of_operatorGreedyMinAvailable
    (hminAvail :
      ∀ M : ℕ, ∀ n : ℕ, ∀ hn : n < operatorGreedyCard M,
        operatorEigenvalueOrderedCenterChoiceMin M ⟨n, hn⟩ ∈
          operatorGreedyAvailableNat M n) :
    OperatorCenterGapPermutationInvariant := by
  refine operatorCenterGapPermutationInvariant_of_operatorGreedyCenter ?_
  intro M n hn
  exact operatorGreedy_hCandAvail_of_min_available M (hminAvail M) n hn

/-- Global closure to permutation-invariant center-gap from predecessor exclusion:
if each level satisfies the predecessor non-collision condition under prefix
minimum-availability, the greedy route closes the full center-gap contract. -/
theorem operatorCenterGapPermutationInvariant_of_operatorGreedyPredExclusion
    (hPredExcl :
      ∀ M : ℕ, ∀ j : ℕ, ∀ hj : j < operatorGreedyCard M, ∀ hjpos : 0 < j,
        (∀ k : ℕ, ∀ hk : k < j,
          operatorEigenvalueOrderedCenterChoiceMin M ⟨k, lt_trans hk hj⟩ ∈
            operatorGreedyAvailableNat M k) →
        operatorGreedyChoiceNat M (j - 1) ≠
          operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩) :
    OperatorCenterGapPermutationInvariant := by
  refine operatorCenterGapPermutationInvariant_of_operatorGreedyCenter ?_
  intro M n hn
  exact operatorGreedy_hCandAvail_of_pred_exclusion M (hPredExcl M) n hn

/-- Global closure reduction:
if the same-min branch always supplies a `min+1` candidate that is both
admissible and available at step `j-1`, then the greedy predecessor-exclusion
route closes `OperatorCenterGapPermutationInvariant`. -/
theorem operatorCenterGapPermutationInvariant_of_operatorGreedyPredExclusion_sameMin_plusOne
    (hPlusOneCandAvailInSameMin :
      ∀ M : ℕ, ∀ j : ℕ, ∀ hj : j < operatorGreedyCard M, ∀ hjpos : 0 < j,
      ∀ hjm1 : j - 1 < operatorGreedyCard M,
        (∀ k : ℕ, ∀ hk : k < j,
          operatorEigenvalueOrderedCenterChoiceMin M ⟨k, lt_trans hk hj⟩ ∈
            operatorGreedyAvailableNat M k) →
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 =
          (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1 →
        ∃ c : Fin (M + 1),
          c.1 =
            (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 + 1 ∧
          c ∈ operatorCenterCandidatesOrdered M ⟨j - 1, hjm1⟩ ∧
          c ∈ operatorGreedyAvailableNat M (j - 1)) :
    OperatorCenterGapPermutationInvariant := by
  refine operatorCenterGapPermutationInvariant_of_operatorGreedyPredExclusion ?_
  intro M j hj hjpos hPrefix
  exact operatorPredExcl_of_sameMin_plusOneAvailable M
    (hPlusOneCandAvailInSameMin M) j hj hjpos hPrefix

/-- Global closure reduction (reduced same-min seam):
if same-min branches provide
1) `max(j-1) > min(j-1)` and
2) availability of `min(j-1)+1`,
then the greedy predecessor-exclusion route closes
`OperatorCenterGapPermutationInvariant`. -/
theorem operatorCenterGapPermutationInvariant_of_operatorGreedyPredExclusion_sameMin_plusOneAvail_and_maxAbove
    (hMaxAboveInSameMin :
      ∀ M : ℕ, ∀ j : ℕ, ∀ hj : j < operatorGreedyCard M, ∀ hjpos : 0 < j,
      ∀ hjm1 : j - 1 < operatorGreedyCard M,
        (∀ k : ℕ, ∀ hk : k < j,
          operatorEigenvalueOrderedCenterChoiceMin M ⟨k, lt_trans hk hj⟩ ∈
            operatorGreedyAvailableNat M k) →
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 =
          (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1 →
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 <
          (operatorEigenvalueOrderedCenterChoiceMax M ⟨j - 1, hjm1⟩).1)
    (hPlusOneAvailInSameMin :
      ∀ M : ℕ, ∀ j : ℕ, ∀ hj : j < operatorGreedyCard M, ∀ hjpos : 0 < j,
      ∀ hjm1 : j - 1 < operatorGreedyCard M,
        (∀ k : ℕ, ∀ hk : k < j,
          operatorEigenvalueOrderedCenterChoiceMin M ⟨k, lt_trans hk hj⟩ ∈
            operatorGreedyAvailableNat M k) →
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 =
          (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1 →
        ∃ c : Fin (M + 1),
          c.1 =
            (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 + 1 ∧
          c ∈ operatorGreedyAvailableNat M (j - 1)) :
    OperatorCenterGapPermutationInvariant := by
  refine operatorCenterGapPermutationInvariant_of_operatorGreedyPredExclusion ?_
  intro M j hj hjpos hPrefix
  exact operatorPredExcl_of_sameMin_plusOneAvail_and_maxAbove M
    (hMaxAboveInSameMin M) (hPlusOneAvailInSameMin M) j hj hjpos hPrefix

/-- Compatibility bridge: the older same-min witness format
`min+1 ∈ candidates ∩ available` implies the reduced pair
`max>min` plus `min+1 ∈ available`. -/
theorem operatorCenterGapPermutationInvariant_of_operatorGreedyPredExclusion_sameMin_plusOne_via_reduced
    (hPlusOneCandAvailInSameMin :
      ∀ M : ℕ, ∀ j : ℕ, ∀ hj : j < operatorGreedyCard M, ∀ hjpos : 0 < j,
      ∀ hjm1 : j - 1 < operatorGreedyCard M,
        (∀ k : ℕ, ∀ hk : k < j,
          operatorEigenvalueOrderedCenterChoiceMin M ⟨k, lt_trans hk hj⟩ ∈
            operatorGreedyAvailableNat M k) →
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 =
          (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1 →
        ∃ c : Fin (M + 1),
          c.1 =
            (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 + 1 ∧
          c ∈ operatorCenterCandidatesOrdered M ⟨j - 1, hjm1⟩ ∧
          c ∈ operatorGreedyAvailableNat M (j - 1)) :
    OperatorCenterGapPermutationInvariant := by
  refine operatorCenterGapPermutationInvariant_of_operatorGreedyPredExclusion_sameMin_plusOneAvail_and_maxAbove ?_ ?_
  · intro M j hj hjpos hjm1 hPrefix hSameMin
    have hPlusOneCandAvail :
        ∃ c : Fin (M + 1),
          c.1 =
            (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 + 1 ∧
          c ∈ operatorCenterCandidatesOrdered M ⟨j - 1, hjm1⟩ ∧
          c ∈ operatorGreedyAvailableNat M (j - 1) :=
      hPlusOneCandAvailInSameMin M j hj hjpos hjm1 hPrefix hSameMin
    have hPlusOneCand :
        ∃ c : Fin (M + 1),
          c.1 =
            (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 + 1 ∧
          c ∈ operatorCenterCandidatesOrdered M ⟨j - 1, hjm1⟩ := by
      rcases hPlusOneCandAvail with ⟨c, hcEq, hcCand, _⟩
      exact ⟨c, hcEq, hcCand⟩
    exact operatorMaxAboveInSameMin_of_plusOneCandidate M j hj hjm1 hPlusOneCand
  · intro M j hj hjpos hjm1 hPrefix hSameMin
    rcases hPlusOneCandAvailInSameMin M j hj hjpos hjm1 hPrefix hSameMin with
      ⟨c, hcEq, _hcCand, hcAvail⟩
    exact ⟨c, hcEq, hcAvail⟩

/-- Final constructive seam (explicit form):
in each same-min branch, no earlier greedy step picks the canonical center
`min(j-1)+1`. -/
def OperatorSameMinPlusOneNoPrevObligation : Prop :=
  ∀ M : ℕ, ∀ j : ℕ, ∀ hj : j < operatorGreedyCard M, ∀ hjpos : 0 < j,
  ∀ hjm1 : j - 1 < operatorGreedyCard M,
    (∀ k : ℕ, ∀ hk : k < j,
      operatorEigenvalueOrderedCenterChoiceMin M ⟨k, lt_trans hk hj⟩ ∈
        operatorGreedyAvailableNat M k) →
    (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 =
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1 →
    ∀ c : Fin (M + 1),
      c.1 =
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 + 1 →
      ∀ i : ℕ, i < (j - 1) → operatorGreedyChoiceNat M i ≠ c

/-- Final same-min structural branch obligation:
in each same-min branch, predecessor candidate interval is non-singleton
(`max(j-1) > min(j-1)`). -/
def OperatorSameMinMaxAboveObligation : Prop :=
  ∀ M : ℕ, ∀ j : ℕ, ∀ hj : j < operatorGreedyCard M, ∀ hjpos : 0 < j,
  ∀ hjm1 : j - 1 < operatorGreedyCard M,
    (∀ k : ℕ, ∀ hk : k < j,
      operatorEigenvalueOrderedCenterChoiceMin M ⟨k, lt_trans hk hj⟩ ∈
        operatorGreedyAvailableNat M k) →
    (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 =
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1 →
    (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 <
      (operatorEigenvalueOrderedCenterChoiceMax M ⟨j - 1, hjm1⟩).1

/-- Bridge: same-min `min+1` availability implies the explicit no-previous-pick
obligation for canonical `min+1`. -/
theorem operatorSameMinPlusOneNoPrevObligation_of_plusOneAvailInSameMin
    (hPlusOneAvailInSameMin :
      ∀ M : ℕ, ∀ j : ℕ, ∀ hj : j < operatorGreedyCard M, ∀ hjpos : 0 < j,
      ∀ hjm1 : j - 1 < operatorGreedyCard M,
        (∀ k : ℕ, ∀ hk : k < j,
          operatorEigenvalueOrderedCenterChoiceMin M ⟨k, lt_trans hk hj⟩ ∈
            operatorGreedyAvailableNat M k) →
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 =
          (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1 →
        ∃ c : Fin (M + 1),
          c.1 =
            (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 + 1 ∧
          c ∈ operatorGreedyAvailableNat M (j - 1)) :
    OperatorSameMinPlusOneNoPrevObligation := by
  intro M j hj hjpos hjm1 hPrefix hSameMin c hcEq i hi
  rcases hPlusOneAvailInSameMin M j hj hjpos hjm1 hPrefix hSameMin with
    ⟨c', hc'Eq, hc'Avail⟩
  have hcc' : c = c' := by
    apply Fin.ext
    omega
  have hAvail : c ∈ operatorGreedyAvailableNat M (j - 1) := by
    simpa [hcc'] using hc'Avail
  exact (operatorAvailable_mem_iff_no_prev M (j - 1) (Nat.le_of_lt hjm1) c).1 hAvail i hi

/-- Reverse bridge: under same-min `max(j-1)>min(j-1)`, the explicit no-prev
obligation implies constructive availability of canonical `min+1`. -/
theorem operatorPlusOneAvailInSameMin_of_noPrev_and_maxAbove
    (hNoPrevPlusOneInSameMin : OperatorSameMinPlusOneNoPrevObligation)
    (hMaxAboveInSameMin :
      ∀ M : ℕ, ∀ j : ℕ, ∀ hj : j < operatorGreedyCard M, ∀ hjpos : 0 < j,
      ∀ hjm1 : j - 1 < operatorGreedyCard M,
        (∀ k : ℕ, ∀ hk : k < j,
          operatorEigenvalueOrderedCenterChoiceMin M ⟨k, lt_trans hk hj⟩ ∈
            operatorGreedyAvailableNat M k) →
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 =
          (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1 →
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 <
          (operatorEigenvalueOrderedCenterChoiceMax M ⟨j - 1, hjm1⟩).1) :
    ∀ M : ℕ, ∀ j : ℕ, ∀ hj : j < operatorGreedyCard M, ∀ hjpos : 0 < j,
    ∀ hjm1 : j - 1 < operatorGreedyCard M,
      (∀ k : ℕ, ∀ hk : k < j,
        operatorEigenvalueOrderedCenterChoiceMin M ⟨k, lt_trans hk hj⟩ ∈
          operatorGreedyAvailableNat M k) →
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 =
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1 →
      ∃ c : Fin (M + 1),
        c.1 =
          (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 + 1 ∧
        c ∈ operatorGreedyAvailableNat M (j - 1) := by
  intro M j hj hjpos hjm1 hPrefix hSameMin
  have hMaxAbove :
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 <
        (operatorEigenvalueOrderedCenterChoiceMax M ⟨j - 1, hjm1⟩).1 :=
    hMaxAboveInSameMin M j hj hjpos hjm1 hPrefix hSameMin
  rcases operatorPlusOneAvail_iff_no_prev_choice_of_min_lt_max M j hjm1 hMaxAbove with
    ⟨c, hcEq, hAvailIff⟩
  have hNoPrev :
      ∀ i : ℕ, i < (j - 1) → operatorGreedyChoiceNat M i ≠ c :=
    hNoPrevPlusOneInSameMin M j hj hjpos hjm1 hPrefix hSameMin c hcEq
  exact ⟨c, hcEq, (hAvailIff).2 hNoPrev⟩

/-- Endgame closure from the explicit final seam:
if same-min branches satisfy `max(j-1)>min(j-1)` and the no-previous-pick
obligation for canonical `min+1`, then the greedy route closes
`OperatorCenterGapPermutationInvariant`. -/
theorem operatorCenterGapPermutationInvariant_of_operatorGreedyPredExclusion_sameMin_plusOneNoPrev_and_maxAbove
    (hMaxAboveInSameMin :
      ∀ M : ℕ, ∀ j : ℕ, ∀ hj : j < operatorGreedyCard M, ∀ hjpos : 0 < j,
      ∀ hjm1 : j - 1 < operatorGreedyCard M,
        (∀ k : ℕ, ∀ hk : k < j,
          operatorEigenvalueOrderedCenterChoiceMin M ⟨k, lt_trans hk hj⟩ ∈
            operatorGreedyAvailableNat M k) →
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 =
          (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1 →
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 <
          (operatorEigenvalueOrderedCenterChoiceMax M ⟨j - 1, hjm1⟩).1)
    (hNoPrevPlusOneInSameMin : OperatorSameMinPlusOneNoPrevObligation) :
    OperatorCenterGapPermutationInvariant := by
  refine operatorCenterGapPermutationInvariant_of_operatorGreedyPredExclusion_sameMin_plusOneAvail_and_maxAbove
    hMaxAboveInSameMin ?_
  intro M j hj hjpos hjm1 hPrefix hSameMin
  have hMaxAbove :
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 <
        (operatorEigenvalueOrderedCenterChoiceMax M ⟨j - 1, hjm1⟩).1 :=
    hMaxAboveInSameMin M j hj hjpos hjm1 hPrefix hSameMin
  rcases operatorPlusOneAvail_iff_no_prev_choice_of_min_lt_max M j hjm1 hMaxAbove with
    ⟨c, hcEq, hAvailIff⟩
  have hNoPrev :
      ∀ i : ℕ, i < (j - 1) → operatorGreedyChoiceNat M i ≠ c :=
    hNoPrevPlusOneInSameMin M j hj hjpos hjm1 hPrefix hSameMin c hcEq
  have hAvail : c ∈ operatorGreedyAvailableNat M (j - 1) :=
    (hAvailIff).2 hNoPrev
  exact ⟨c, hcEq, hAvail⟩

/-- Packaged center-gap closure from the two explicit same-min endgame
obligations: non-singleton same-min branch plus canonical `min+1` no-prev. -/
theorem operatorCenterGapPermutationInvariant_of_sameMinEndgameObligations
    (hMaxAboveInSameMin : OperatorSameMinMaxAboveObligation)
    (hNoPrevPlusOneInSameMin : OperatorSameMinPlusOneNoPrevObligation) :
    OperatorCenterGapPermutationInvariant := by
  exact
    operatorCenterGapPermutationInvariant_of_operatorGreedyPredExclusion_sameMin_plusOneNoPrev_and_maxAbove
      hMaxAboveInSameMin hNoPrevPlusOneInSameMin

/-- Geometric-algebraic same-min symmetry contract (Clifford-lane interface):
in every same-min branch, the predecessor candidate interval is interior
(`0 < min` and `max < M`) and the canonical center `min+1` is available.

This is the minimal structural package that feeds the two explicit seams:
`maxAbove` (via interior non-collapse) and `noPrev` (via availability bridge). -/
def OperatorCliffordSameMinSymmetryContract : Prop :=
  (∀ M : ℕ, ∀ j : ℕ, ∀ hj : j < operatorGreedyCard M, ∀ hjpos : 0 < j,
    ∀ hjm1 : j - 1 < operatorGreedyCard M,
      (∀ k : ℕ, ∀ hk : k < j,
        operatorEigenvalueOrderedCenterChoiceMin M ⟨k, lt_trans hk hj⟩ ∈
          operatorGreedyAvailableNat M k) →
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 =
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1 →
      0 < (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 ∧
      (operatorEigenvalueOrderedCenterChoiceMax M ⟨j - 1, hjm1⟩).1 < M)
  ∧
  (∀ M : ℕ, ∀ j : ℕ, ∀ hj : j < operatorGreedyCard M, ∀ hjpos : 0 < j,
    ∀ hjm1 : j - 1 < operatorGreedyCard M,
      (∀ k : ℕ, ∀ hk : k < j,
        operatorEigenvalueOrderedCenterChoiceMin M ⟨k, lt_trans hk hj⟩ ∈
          operatorGreedyAvailableNat M k) →
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 =
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1 →
      ∃ c : Fin (M + 1),
        c.1 =
          (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 + 1 ∧
        c ∈ operatorGreedyAvailableNat M (j - 1))

/-- Clifford same-min symmetry contract closes the operator center-gap geometry
lane by discharging both explicit same-min seams (`maxAbove`, `noPrev`). -/
theorem operatorCenterGapPermutationInvariant_of_cliffordSameMinSymmetryContract
    (hCliff : OperatorCliffordSameMinSymmetryContract) :
    OperatorCenterGapPermutationInvariant := by
  rcases hCliff with ⟨hInteriorInSameMin, hPlusOneAvailInSameMin⟩
  have hMaxAboveInSameMin :
      ∀ M : ℕ, ∀ j : ℕ, ∀ hj : j < operatorGreedyCard M, ∀ hjpos : 0 < j,
      ∀ hjm1 : j - 1 < operatorGreedyCard M,
        (∀ k : ℕ, ∀ hk : k < j,
          operatorEigenvalueOrderedCenterChoiceMin M ⟨k, lt_trans hk hj⟩ ∈
            operatorGreedyAvailableNat M k) →
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 =
          (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1 →
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 <
          (operatorEigenvalueOrderedCenterChoiceMax M ⟨j - 1, hjm1⟩).1 := by
    intro M j hj hjpos hjm1 hPrefix hSameMin
    rcases hInteriorInSameMin M j hj hjpos hjm1 hPrefix hSameMin with
      ⟨hminPos, hmaxLt⟩
    exact operatorMaxAboveInSameMin_of_interior M j hj hjpos hjm1 hminPos hmaxLt
  have hNoPrevPlusOneInSameMin : OperatorSameMinPlusOneNoPrevObligation :=
    operatorSameMinPlusOneNoPrevObligation_of_plusOneAvailInSameMin
      hPlusOneAvailInSameMin
  exact
    operatorCenterGapPermutationInvariant_of_operatorGreedyPredExclusion_sameMin_plusOneNoPrev_and_maxAbove
      hMaxAboveInSameMin hNoPrevPlusOneInSameMin

/-- Obligation instantiator: extract same-min `maxAbove` from the Clifford
same-min symmetry contract. -/
theorem operatorMaxAboveInSameMin_of_cliffordSameMinSymmetryContract
    (hCliff : OperatorCliffordSameMinSymmetryContract) :
    ∀ M : ℕ, ∀ j : ℕ, ∀ hj : j < operatorGreedyCard M, ∀ hjpos : 0 < j,
    ∀ hjm1 : j - 1 < operatorGreedyCard M,
      (∀ k : ℕ, ∀ hk : k < j,
        operatorEigenvalueOrderedCenterChoiceMin M ⟨k, lt_trans hk hj⟩ ∈
          operatorGreedyAvailableNat M k) →
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 =
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1 →
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 <
        (operatorEigenvalueOrderedCenterChoiceMax M ⟨j - 1, hjm1⟩).1 := by
  rcases hCliff with ⟨hInteriorInSameMin, _hPlusOneAvailInSameMin⟩
  intro M j hj hjpos hjm1 hPrefix hSameMin
  rcases hInteriorInSameMin M j hj hjpos hjm1 hPrefix hSameMin with
    ⟨hminPos, hmaxLt⟩
  exact operatorMaxAboveInSameMin_of_interior M j hj hjpos hjm1 hminPos hmaxLt

/-- Obligation instantiator: extract same-min `min+1` availability from the
Clifford same-min symmetry contract. -/
theorem operatorPlusOneAvailInSameMin_of_cliffordSameMinSymmetryContract
    (hCliff : OperatorCliffordSameMinSymmetryContract) :
    ∀ M : ℕ, ∀ j : ℕ, ∀ hj : j < operatorGreedyCard M, ∀ hjpos : 0 < j,
    ∀ hjm1 : j - 1 < operatorGreedyCard M,
      (∀ k : ℕ, ∀ hk : k < j,
        operatorEigenvalueOrderedCenterChoiceMin M ⟨k, lt_trans hk hj⟩ ∈
          operatorGreedyAvailableNat M k) →
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 =
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1 →
      ∃ c : Fin (M + 1),
        c.1 =
          (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 + 1 ∧
        c ∈ operatorGreedyAvailableNat M (j - 1) := by
  exact (And.right hCliff)

/-- Obligation instantiator: extract explicit same-min no-prev seam from the
Clifford same-min symmetry contract. -/
theorem operatorSameMinPlusOneNoPrevObligation_of_cliffordSameMinSymmetryContract
    (hCliff : OperatorCliffordSameMinSymmetryContract) :
    OperatorSameMinPlusOneNoPrevObligation := by
  exact operatorSameMinPlusOneNoPrevObligation_of_plusOneAvailInSameMin
    (operatorPlusOneAvailInSameMin_of_cliffordSameMinSymmetryContract hCliff)

/-- Step-1 style instantiation:
Clifford same-min symmetry contract implies the explicit same-min `max>min`
obligation. -/
theorem operatorSameMinMaxAboveObligation_clifford
    (hCliff : OperatorCliffordSameMinSymmetryContract) :
    OperatorSameMinMaxAboveObligation := by
  exact operatorMaxAboveInSameMin_of_cliffordSameMinSymmetryContract hCliff

/-- Step-2 style instantiation:
Clifford same-min symmetry contract implies the explicit same-min
canonical-`min+1` no-prev obligation. -/
theorem operatorSameMinPlusOneNoPrevObligation_clifford
    (hCliff : OperatorCliffordSameMinSymmetryContract) :
    OperatorSameMinPlusOneNoPrevObligation := by
  exact operatorSameMinPlusOneNoPrevObligation_of_cliffordSameMinSymmetryContract
    hCliff

/-- Hodge/even-odd contract for the same-min branch.

This names the minimal constructive payload expected from the Cl(1,3) parity
lane: interior non-collapse (`0 < min`, `max < M`) and canonical `min+1`
availability in every same-min branch, grounded in concrete operator symmetry
(`T(A)=C-A` and centered spectral reflection). -/
structure OperatorHodgeParityContract : Prop where
  revParityCenterSub :
    ∀ n : ℕ,
      revParityConjugateMatrixC n (structuralRiemannMatrixC n)
        = structuralCenterMatrixC n - structuralRiemannMatrixC n
  spectrumReflectCentered :
    ∀ n : ℕ, ∀ l : ℂ,
      l ∈ spectrum ℂ (structuralRiemannMatrixC n) →
      ((structuralCenterQ n : ℂ) - l) ∈ spectrum ℂ (structuralRiemannMatrixC n)
  interiorInSameMin :
    ∀ M : ℕ, ∀ j : ℕ, ∀ hj : j < operatorGreedyCard M, ∀ hjpos : 0 < j,
    ∀ hjm1 : j - 1 < operatorGreedyCard M,
      (∀ k : ℕ, ∀ hk : k < j,
        operatorEigenvalueOrderedCenterChoiceMin M ⟨k, lt_trans hk hj⟩ ∈
          operatorGreedyAvailableNat M k) →
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 =
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1 →
      0 < (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 ∧
      (operatorEigenvalueOrderedCenterChoiceMax M ⟨j - 1, hjm1⟩).1 < M
  plusOneAvailInSameMin :
    ∀ M : ℕ, ∀ j : ℕ, ∀ hj : j < operatorGreedyCard M, ∀ hjpos : 0 < j,
    ∀ hjm1 : j - 1 < operatorGreedyCard M,
      (∀ k : ℕ, ∀ hk : k < j,
        operatorEigenvalueOrderedCenterChoiceMin M ⟨k, lt_trans hk hj⟩ ∈
          operatorGreedyAvailableNat M k) →
      (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 =
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1 →
      ∃ c : Fin (M + 1),
        c.1 =
          (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 + 1 ∧
        c ∈ operatorGreedyAvailableNat M (j - 1)

/-- Concrete operator-symmetry core required by the Hodge parity lane:
the finite structural matrix satisfies centered reversal/parity balance and the
induced centered spectral reflection. -/
def OperatorHodgeParityCore : Prop :=
  (∀ n : ℕ,
    revParityConjugateMatrixC n (structuralRiemannMatrixC n)
      = structuralCenterMatrixC n - structuralRiemannMatrixC n)
  ∧
  (∀ n : ℕ, ∀ l : ℂ,
    l ∈ spectrum ℂ (structuralRiemannMatrixC n) →
    ((structuralCenterQ n : ℂ) - l) ∈ spectrum ℂ (structuralRiemannMatrixC n))

/-- The operator already provides the Hodge parity core by construction. -/
theorem operatorHodgeParityCore_default : OperatorHodgeParityCore := by
  refine ⟨?_, ?_⟩
  · intro n
    exact structuralRiemannMatrixC_revParity_eq_center_sub n
  · intro n l hl
    exact structuralRiemannMatrixC_spectrum_reflect n l hl

/-- Constructor: the Hodge parity contract is exactly the concrete operator
symmetry core plus the same-min branch payload. -/
theorem operatorHodgeParityContract_of_core_and_sameMinPayload
    (hCore : OperatorHodgeParityCore)
    (hInteriorInSameMin :
      ∀ M : ℕ, ∀ j : ℕ, ∀ hj : j < operatorGreedyCard M, ∀ hjpos : 0 < j,
      ∀ hjm1 : j - 1 < operatorGreedyCard M,
        (∀ k : ℕ, ∀ hk : k < j,
          operatorEigenvalueOrderedCenterChoiceMin M ⟨k, lt_trans hk hj⟩ ∈
            operatorGreedyAvailableNat M k) →
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 =
          (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1 →
        0 < (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 ∧
        (operatorEigenvalueOrderedCenterChoiceMax M ⟨j - 1, hjm1⟩).1 < M)
    (hPlusOneAvailInSameMin :
      ∀ M : ℕ, ∀ j : ℕ, ∀ hj : j < operatorGreedyCard M, ∀ hjpos : 0 < j,
      ∀ hjm1 : j - 1 < operatorGreedyCard M,
        (∀ k : ℕ, ∀ hk : k < j,
          operatorEigenvalueOrderedCenterChoiceMin M ⟨k, lt_trans hk hj⟩ ∈
            operatorGreedyAvailableNat M k) →
        (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 =
          (operatorEigenvalueOrderedCenterChoiceMin M ⟨j, hj⟩).1 →
        ∃ c : Fin (M + 1),
          c.1 =
            (operatorEigenvalueOrderedCenterChoiceMin M ⟨j - 1, hjm1⟩).1 + 1 ∧
          c ∈ operatorGreedyAvailableNat M (j - 1)) :
    OperatorHodgeParityContract := by
  rcases hCore with ⟨hRev, hSpec⟩
  exact
    { revParityCenterSub := hRev
      spectrumReflectCentered := hSpec
      interiorInSameMin := hInteriorInSameMin
      plusOneAvailInSameMin := hPlusOneAvailInSameMin }

/-- Hodge parity contract instantiates the existing Clifford same-min symmetry
contract surface. -/
theorem operatorCliffordSameMinSymmetryContract_of_hodgeParity
    (hHP : OperatorHodgeParityContract) :
    OperatorCliffordSameMinSymmetryContract := by
  refine ⟨hHP.interiorInSameMin, hHP.plusOneAvailInSameMin⟩

/-- Step-1 explicit obligation directly from the Hodge parity contract. -/
theorem operatorSameMinMaxAboveObligation_of_hodgeParity
    (hHP : OperatorHodgeParityContract) :
    OperatorSameMinMaxAboveObligation := by
  exact operatorSameMinMaxAboveObligation_clifford
    (operatorCliffordSameMinSymmetryContract_of_hodgeParity hHP)

/-- Step-2 explicit obligation directly from the Hodge parity contract. -/
theorem operatorSameMinPlusOneNoPrevObligation_of_hodgeParity
    (hHP : OperatorHodgeParityContract) :
    OperatorSameMinPlusOneNoPrevObligation := by
  exact operatorSameMinPlusOneNoPrevObligation_clifford
    (operatorCliffordSameMinSymmetryContract_of_hodgeParity hHP)

/-- Tie-break closure theorem: if the deterministic ordered tie-break center is
admissible in every ordered candidate set, then the full permutation-invariant
center-gap contract holds. -/
theorem operatorCenterGapPermutationInvariant_of_orderedTieBreak
    (hTie :
      ∀ M : ℕ, ∀ j : Fin (Fintype.card (Fin (M + 1))),
        operatorOrderedTieBreakCenter M j ∈ operatorCenterCandidatesOrdered M j) :
    OperatorCenterGapPermutationInvariant := by
  intro M
  classical
  let e : Fin (M + 1) ≃ Fin (Fintype.card (Fin (M + 1))) :=
    operatorEigenvaluesReindexToOrderedEquiv M
  let f : Fin (M + 1) → Fin (M + 1) := fun i => operatorOrderedTieBreakCenter M (e i)
  have hf_injective : Function.Injective f := by
    intro i1 i2 h
    exact e.injective (operatorOrderedTieBreakCenter_injective M h)
  have hf_surjective : Function.Surjective f := (Finite.injective_iff_surjective).1 hf_injective
  let σ : Fin (M + 1) ≃ Fin (M + 1) := Equiv.ofBijective f ⟨hf_injective, hf_surjective⟩
  refine ⟨σ, ?_⟩
  intro i
  have hTie_i :
      operatorOrderedTieBreakCenter M (e i) ∈ operatorCenterCandidatesOrdered M (e i) :=
    hTie M (e i)
  have hTie_abs :
      |operatorEigenvaluesOrdered M (e i) -
        (((operatorOrderedTieBreakCenter M (e i)).1 : ℝ) + (29 : ℝ) / 16)| ≤ (12 : ℝ) / 11 := by
    exact (Finset.mem_filter.mp hTie_i).2
  have hEqEig :
      operatorEigenvalues M i = operatorEigenvaluesOrdered M (e i) := by
    simpa [e] using operatorEigenvalues_eq_ordered_reindex M i
  have hσEq : σ i = operatorOrderedTieBreakCenter M (e i) := rfl
  have hAbs :
      |operatorEigenvalues M i - (((σ i).1 : ℝ) + (29 : ℝ) / 16)| ≤ (12 : ℝ) / 11 := by
    simpa [hEqEig, hσEq] using hTie_abs
  exact hAbs

/-- Hall-style finite-level center-gap condition:
every finite subfamily of eigenvalue candidate-center sets has union cardinality
at least the subfamily cardinality. -/
def OperatorCenterGapHallCondition : Prop :=
  ∀ M : ℕ, ∀ s : Finset (Fin (M + 1)),
    s.card ≤ (s.biUnion (operatorCenterCandidates M)).card

/-- Hall condition yields permutation-invariant center-gap pairing. -/
theorem operatorCenterGapPermutationInvariant_of_hallCondition
    (hHall : OperatorCenterGapHallCondition) :
    OperatorCenterGapPermutationInvariant := by
  intro M
  classical
  let t : Fin (M + 1) → Finset (Fin (M + 1)) := operatorCenterCandidates M
  have hcard : ∀ s : Finset (Fin (M + 1)), s.card ≤ (s.biUnion t).card := by
    intro s
    simpa [t] using hHall M s
  rcases (Finset.all_card_le_biUnion_card_iff_existsInjective' t).1 hcard with
      ⟨f, hf, hmem⟩
  have hsurj : Function.Surjective f := (Finite.injective_iff_surjective).1 hf
  let σ : Fin (M + 1) ≃ Fin (M + 1) := Equiv.ofBijective f ⟨hf, hsurj⟩
  refine ⟨σ, ?_⟩
  intro i
  have hfi : f i ∈ t i := hmem i
  simpa [t, operatorCenterCandidates, σ] using hfi

/-- Ordered Weyl center-gap implies Hall condition by taking the diagonal index
as a witness in each candidate-center set. -/
theorem operatorCenterGapHallCondition_of_weylCenterGap
    (hW : OperatorWeylCenterGap) :
    OperatorCenterGapHallCondition := by
  intro M s
  have hsubset : s ⊆ s.biUnion (operatorCenterCandidates M) := by
    intro i hi
    refine Finset.mem_biUnion.mpr ?_
    refine ⟨i, hi, ?_⟩
    refine Finset.mem_filter.mpr ?_
    exact ⟨Finset.mem_univ i, by simpa using hW M i⟩
  exact Finset.card_le_card hsubset

/-- Ordered Weyl center-gap implies the permutation-invariant center-gap
via Hall SDR on the candidate-center family. -/
theorem operatorCenterGapPermutationInvariant_of_weylCenterGap
    (hW : OperatorWeylCenterGap) :
    OperatorCenterGapPermutationInvariant := by
  exact operatorCenterGapPermutationInvariant_of_hallCondition
    (operatorCenterGapHallCondition_of_weylCenterGap hW)

/-- Structural center on the real axis for the `k`-th Gershgorin lane. -/
def operatorCenterAt (k : ℕ) : ℝ := (k : ℝ) + (29 : ℝ) / 16

/-- Concrete tridiagonal diagonal coefficient `d_k` for the structural matrix:
`d_k = (k+1) + 13/16 = k + 29/16`. -/
def operatorSturmDiag (k : ℕ) : ℝ := ((k + 1 : ℕ) : ℝ) + (13 : ℝ) / 16

/-- Concrete nearest-neighbor coefficient `a_k` for the structural matrix. -/
def operatorSturmOff : ℝ := (6 : ℝ) / 11

/-- Concrete characteristic-polynomial recursion for the structural tridiagonal lane:
`p₀(x)=1`, `p₁(x)=d₀-x`, and
`p_{k+1}(x) = (d_k-x)p_k(x) - a² p_{k-1}(x)` for `k≥1`. -/
def operatorSturmP : ℕ → ℝ → ℝ
  | 0, _ => 1
  | 1, x => operatorSturmDiag 0 - x
  | n + 2, x =>
      (operatorSturmDiag (n + 1) - x) * operatorSturmP (n + 1) x
        - (operatorSturmOff ^ (2 : ℕ)) * operatorSturmP n x

theorem operatorSturmDiag_eq_centerAt (k : ℕ) :
    operatorSturmDiag k = operatorCenterAt k := by
  unfold operatorSturmDiag operatorCenterAt
  norm_num
  ring_nf

@[simp] theorem operatorSturmP_zero (x : ℝ) :
    operatorSturmP 0 x = 1 := by
  rfl

@[simp] theorem operatorSturmP_one (x : ℝ) :
    operatorSturmP 1 x = operatorSturmDiag 0 - x := by
  rfl

theorem operatorSturmP_step (n : ℕ) (x : ℝ) :
    operatorSturmP (n + 2) x =
      (operatorSturmDiag (n + 1) - x) * operatorSturmP (n + 1) x
        - (operatorSturmOff ^ (2 : ℕ)) * operatorSturmP n x := by
  rfl

theorem operatorSturmP_step_center (n : ℕ) (x : ℝ) :
    operatorSturmP (n + 2) x =
      (operatorCenterAt (n + 1) - x) * operatorSturmP (n + 1) x
        - (operatorSturmOff ^ (2 : ℕ)) * operatorSturmP n x := by
  simpa [operatorSturmDiag_eq_centerAt] using operatorSturmP_step n x

theorem operatorSturmOff_sq_pos : 0 < (operatorSturmOff ^ (2 : ℕ)) := by
  unfold operatorSturmOff
  positivity

@[simp] theorem operatorSturmP_one_eq_center_sub (x : ℝ) :
    operatorSturmP 1 x = operatorCenterAt 0 - x := by
  simpa [operatorSturmP, operatorSturmDiag_eq_centerAt]

/-- Ternary sign marker used for finite Sturm sign-variation counting. -/
def operatorSturmSign (r : ℝ) : Int :=
  if r > 0 then 1 else if r < 0 then -1 else 0

theorem operatorSturmSign_eq_one_of_pos {r : ℝ} (hr : 0 < r) :
    operatorSturmSign r = 1 := by
  unfold operatorSturmSign
  have hnotlt : ¬ r < 0 := by linarith
  simp [hr, hnotlt]

theorem operatorSturmSign_eq_neg_one_of_neg {r : ℝ} (hr : r < 0) :
    operatorSturmSign r = -1 := by
  unfold operatorSturmSign
  have hnotgt : ¬ r > 0 := by linarith
  simp [hnotgt, hr]

theorem operatorSturmSign_ne_one_of_nonpos {r : ℝ} (hr : r ≤ 0) :
    operatorSturmSign r ≠ 1 := by
  intro hs
  by_cases hgt : r > 0
  · linarith
  · by_cases hlt : r < 0
    · have hcontra : (-1 : Int) = 1 := by
        simpa [operatorSturmSign, hgt, hlt] using hs
      exact (by decide : (-1 : Int) ≠ 1) hcontra
    · have hz : r = 0 := by linarith
      have hcontra : (0 : Int) = 1 := by
        simpa [operatorSturmSign, hgt, hlt, hz] using hs
      exact (by decide : (0 : Int) ≠ 1) hcontra

theorem operatorSturmSign_ne_neg_one_of_nonneg {r : ℝ} (hr : 0 ≤ r) :
    operatorSturmSign r ≠ -1 := by
  intro hs
  by_cases hgt : r > 0
  · have hcontra : (1 : Int) = -1 := by
      simpa [operatorSturmSign, hgt] using hs
    exact (by decide : (1 : Int) ≠ -1) hcontra
  · by_cases hlt : r < 0
    · linarith
    · have hz : r = 0 := by linarith
      have hcontra : (0 : Int) = -1 := by
        simpa [operatorSturmSign, hgt, hlt, hz] using hs
      exact (by decide : (0 : Int) ≠ -1) hcontra

/-- Recurrence-local sign dominance (positive branch):
if `x` is below the next structural center, the terminal value is positive, and
the last edge flips sign, then the previous recurrence value is positive. -/
theorem operatorSturm_prev_pos_of_center_gt_and_flip_pos
    (M : ℕ) (x : ℝ)
    (hcx : x < operatorCenterAt (M + 1))
    (hp1 : 0 < operatorSturmP (M + 1) x)
    (hflip : operatorSturmSign (operatorSturmP (M + 1) x) ≠
      operatorSturmSign (operatorSturmP (M + 2) x)) :
    0 < operatorSturmP M x := by
  have hs1 : operatorSturmSign (operatorSturmP (M + 1) x) = 1 :=
    operatorSturmSign_eq_one_of_pos hp1
  have hp2_nonpos : operatorSturmP (M + 2) x ≤ 0 := by
    by_contra hp2_pos
    have hs2 : operatorSturmSign (operatorSturmP (M + 2) x) = 1 :=
      operatorSturmSign_eq_one_of_pos (lt_of_not_ge hp2_pos)
    exact hflip (by simpa [hs1, hs2])
  have hcpos : 0 < operatorCenterAt (M + 1) - x := sub_pos.mpr hcx
  have hdpos : 0 < (operatorSturmOff ^ (2 : ℕ)) := operatorSturmOff_sq_pos
  have hrec := operatorSturmP_step_center M x
  have hrhs_pos :
      0 < (operatorCenterAt (M + 1) - x) * operatorSturmP (M + 1) x
            - operatorSturmP (M + 2) x := by
    have hmulpos : 0 < (operatorCenterAt (M + 1) - x) * operatorSturmP (M + 1) x :=
      mul_pos hcpos hp1
    linarith
  have hmul_pos : 0 < (operatorSturmOff ^ (2 : ℕ)) * operatorSturmP M x := by
    linarith [hrec, hrhs_pos]
  nlinarith [hdpos, hmul_pos]

/-- Recurrence-local sign dominance (negative branch):
if `x` is below the next structural center, the terminal value is negative, and
the last edge flips sign, then the previous recurrence value is negative. -/
theorem operatorSturm_prev_neg_of_center_gt_and_flip_neg
    (M : ℕ) (x : ℝ)
    (hcx : x < operatorCenterAt (M + 1))
    (hp1 : operatorSturmP (M + 1) x < 0)
    (hflip : operatorSturmSign (operatorSturmP (M + 1) x) ≠
      operatorSturmSign (operatorSturmP (M + 2) x)) :
    operatorSturmP M x < 0 := by
  have hs1 : operatorSturmSign (operatorSturmP (M + 1) x) = -1 :=
    operatorSturmSign_eq_neg_one_of_neg hp1
  have hp2_nonneg : 0 ≤ operatorSturmP (M + 2) x := by
    by_contra hp2_neg
    have hs2 : operatorSturmSign (operatorSturmP (M + 2) x) = -1 :=
      operatorSturmSign_eq_neg_one_of_neg (lt_of_not_ge hp2_neg)
    exact hflip (by simpa [hs1, hs2])
  have hcpos : 0 < operatorCenterAt (M + 1) - x := sub_pos.mpr hcx
  have hdpos : 0 < (operatorSturmOff ^ (2 : ℕ)) := operatorSturmOff_sq_pos
  have hrec := operatorSturmP_step_center M x
  have hrhs_neg :
      (operatorCenterAt (M + 1) - x) * operatorSturmP (M + 1) x
        - operatorSturmP (M + 2) x < 0 := by
    have hmulneg : (operatorCenterAt (M + 1) - x) * operatorSturmP (M + 1) x < 0 :=
      mul_neg_of_pos_of_neg hcpos hp1
    linarith
  have hmul_neg : (operatorSturmOff ^ (2 : ℕ)) * operatorSturmP M x < 0 := by
    linarith [hrec, hrhs_neg]
  nlinarith [hdpos, hmul_neg]

/-- Recurrence-local sign lock below the next center:
if the terminal edge flips and `p_{M+1}(x) ≠ 0`, then `p_M(x)` and `p_{M+1}(x)`
have the same sign. -/
theorem operatorSturm_prev_sign_eq_of_center_gt_and_flip
    (M : ℕ) (x : ℝ)
    (hcx : x < operatorCenterAt (M + 1))
    (hflip : operatorSturmSign (operatorSturmP (M + 1) x) ≠
      operatorSturmSign (operatorSturmP (M + 2) x))
    (hp1nz : operatorSturmP (M + 1) x ≠ 0) :
    operatorSturmSign (operatorSturmP M x) =
      operatorSturmSign (operatorSturmP (M + 1) x) := by
  have hp1_cases :
      operatorSturmP (M + 1) x < 0 ∨ 0 < operatorSturmP (M + 1) x :=
    lt_or_gt_of_ne hp1nz
  rcases hp1_cases with hp1_neg | hp1_pos
  · have hp0_neg : operatorSturmP M x < 0 :=
      operatorSturm_prev_neg_of_center_gt_and_flip_neg M x hcx hp1_neg hflip
    have hs0 : operatorSturmSign (operatorSturmP M x) = -1 :=
      operatorSturmSign_eq_neg_one_of_neg hp0_neg
    have hs1 : operatorSturmSign (operatorSturmP (M + 1) x) = -1 :=
      operatorSturmSign_eq_neg_one_of_neg hp1_neg
    simpa [hs0, hs1]
  · have hp0_pos : 0 < operatorSturmP M x :=
      operatorSturm_prev_pos_of_center_gt_and_flip_pos M x hcx hp1_pos hflip
    have hs0 : operatorSturmSign (operatorSturmP M x) = 1 :=
      operatorSturmSign_eq_one_of_pos hp0_pos
    have hs1 : operatorSturmSign (operatorSturmP (M + 1) x) = 1 :=
      operatorSturmSign_eq_one_of_pos hp1_pos
    simpa [hs0, hs1]

/-- Indicator form of the previous-sign lock:
below the next center, if the terminal edge flips and `p_{M+1}(x) ≠ 0`, then
the previous edge indicator is forced to `0`. -/
theorem operatorSturm_prev_edge_indicator_zero_of_center_gt_and_flip
    (M : ℕ) (x : ℝ)
    (hcx : x < operatorCenterAt (M + 1))
    (hflip : operatorSturmSign (operatorSturmP (M + 1) x) ≠
      operatorSturmSign (operatorSturmP (M + 2) x))
    (hp1nz : operatorSturmP (M + 1) x ≠ 0) :
    (if operatorSturmSign (operatorSturmP M x) ≠
          operatorSturmSign (operatorSturmP (M + 1) x) then 1 else 0) = 0 := by
  have hsame :
      operatorSturmSign (operatorSturmP M x) =
        operatorSturmSign (operatorSturmP (M + 1) x) :=
    operatorSturm_prev_sign_eq_of_center_gt_and_flip M x hcx hflip hp1nz
  simp [hsame]

/-- Upper-boundary pattern connector:
if the local increment indicators satisfy `(a,b)=(1,0)` and `p_{M+1}(x) ≠ 0`,
then the previous edge indicator is forced to `0` by the recurrence-local
sign lock. -/
theorem operatorSturm_prev_edge_indicator_zero_of_upper_pattern_and_nonzero
    (M : ℕ) (x : ℝ)
    (ha1 : (if operatorSturmSign (operatorSturmP (M + 1) x) ≠
                operatorSturmSign (operatorSturmP (M + 2) x) then 1 else 0) = 1)
    (hb0 : (if operatorCenterAt (M + 1) ≤ x then 1 else 0) = 0)
    (hp1nz : operatorSturmP (M + 1) x ≠ 0) :
    (if operatorSturmSign (operatorSturmP M x) ≠
          operatorSturmSign (operatorSturmP (M + 1) x) then 1 else 0) = 0 := by
  have hflip : operatorSturmSign (operatorSturmP (M + 1) x) ≠
      operatorSturmSign (operatorSturmP (M + 2) x) := by
    by_cases h :
        operatorSturmSign (operatorSturmP (M + 1) x) =
          operatorSturmSign (operatorSturmP (M + 2) x)
    · simp [h] at ha1
    · exact h
  have hcx : x < operatorCenterAt (M + 1) := by
    by_cases hle : operatorCenterAt (M + 1) ≤ x
    · simp [hle] at hb0
    · exact lt_of_not_ge hle
  exact operatorSturm_prev_edge_indicator_zero_of_center_gt_and_flip M x hcx hflip hp1nz

/-- Recurrence-local sign separation above (or at) the next center:
if the terminal edge does not flip and `p_{M+1}(x) ≠ 0`, then `p_M(x)` and
`p_{M+1}(x)` have opposite signs. -/
theorem operatorSturm_prev_sign_ne_of_center_le_and_noflip_nonzero
    (M : ℕ) (x : ℝ)
    (hcx : operatorCenterAt (M + 1) ≤ x)
    (hnoflip : operatorSturmSign (operatorSturmP (M + 1) x) =
      operatorSturmSign (operatorSturmP (M + 2) x))
    (hp1nz : operatorSturmP (M + 1) x ≠ 0) :
    operatorSturmSign (operatorSturmP M x) ≠
      operatorSturmSign (operatorSturmP (M + 1) x) := by
  have hc_nonpos : operatorCenterAt (M + 1) - x ≤ 0 := sub_nonpos.mpr hcx
  have hd_pos : 0 < (operatorSturmOff ^ (2 : ℕ)) := operatorSturmOff_sq_pos
  have hrec := operatorSturmP_step_center M x
  have hp1_cases :
      operatorSturmP (M + 1) x < 0 ∨ 0 < operatorSturmP (M + 1) x :=
    lt_or_gt_of_ne hp1nz
  rcases hp1_cases with hp1_neg | hp1_pos
  · have hs1 : operatorSturmSign (operatorSturmP (M + 1) x) = -1 :=
      operatorSturmSign_eq_neg_one_of_neg hp1_neg
    have hs2 : operatorSturmSign (operatorSturmP (M + 2) x) = -1 := by
      simpa [hs1] using hnoflip.symm
    have hp2_neg : operatorSturmP (M + 2) x < 0 := by
      by_contra hp2_nonneg
      exact (operatorSturmSign_ne_neg_one_of_nonneg (not_lt.mp hp2_nonneg)) hs2
    have hterm_nonneg :
        0 ≤ (operatorCenterAt (M + 1) - x) * operatorSturmP (M + 1) x :=
      mul_nonneg_of_nonpos_of_nonpos hc_nonpos (le_of_lt hp1_neg)
    have hmul_pos : 0 < (operatorSturmOff ^ (2 : ℕ)) * operatorSturmP M x := by
      linarith [hrec, hp2_neg, hterm_nonneg]
    have hp0_pos : 0 < operatorSturmP M x := by
      nlinarith [hd_pos, hmul_pos]
    have hs0 : operatorSturmSign (operatorSturmP M x) = 1 :=
      operatorSturmSign_eq_one_of_pos hp0_pos
    intro hEq
    have : (1 : Int) = -1 := by simpa [hs0, hs1] using hEq
    exact (by decide : (1 : Int) ≠ -1) this
  · have hs1 : operatorSturmSign (operatorSturmP (M + 1) x) = 1 :=
      operatorSturmSign_eq_one_of_pos hp1_pos
    have hs2 : operatorSturmSign (operatorSturmP (M + 2) x) = 1 := by
      simpa [hs1] using hnoflip.symm
    have hp2_pos : 0 < operatorSturmP (M + 2) x := by
      by_contra hp2_nonpos
      exact (operatorSturmSign_ne_one_of_nonpos (le_of_not_gt hp2_nonpos)) hs2
    have hterm_nonpos :
        (operatorCenterAt (M + 1) - x) * operatorSturmP (M + 1) x ≤ 0 :=
      mul_nonpos_of_nonpos_of_nonneg hc_nonpos (le_of_lt hp1_pos)
    have hmul_neg : (operatorSturmOff ^ (2 : ℕ)) * operatorSturmP M x < 0 := by
      linarith [hrec, hp2_pos, hterm_nonpos]
    have hp0_neg : operatorSturmP M x < 0 := by
      nlinarith [hd_pos, hmul_neg]
    have hs0 : operatorSturmSign (operatorSturmP M x) = -1 :=
      operatorSturmSign_eq_neg_one_of_neg hp0_neg
    intro hEq
    have : (-1 : Int) = 1 := by simpa [hs0, hs1] using hEq
    exact (by decide : (-1 : Int) ≠ 1) this

/-- Indicator form of the previous-sign separation above (or at) center:
if the terminal edge does not flip and `p_{M+1}(x) ≠ 0`, the previous edge
indicator is forced to `1`. -/
theorem operatorSturm_prev_edge_indicator_one_of_center_le_and_noflip_nonzero
    (M : ℕ) (x : ℝ)
    (hcx : operatorCenterAt (M + 1) ≤ x)
    (hnoflip : operatorSturmSign (operatorSturmP (M + 1) x) =
      operatorSturmSign (operatorSturmP (M + 2) x))
    (hp1nz : operatorSturmP (M + 1) x ≠ 0) :
    (if operatorSturmSign (operatorSturmP M x) ≠
          operatorSturmSign (operatorSturmP (M + 1) x) then 1 else 0) = 1 := by
  have hneq :
      operatorSturmSign (operatorSturmP M x) ≠
        operatorSturmSign (operatorSturmP (M + 1) x) :=
    operatorSturm_prev_sign_ne_of_center_le_and_noflip_nonzero M x hcx hnoflip hp1nz
  simp [hneq]

/-- Lower-boundary pattern connector:
if the local increment indicators satisfy `(b,a)=(1,0)` and `p_{M+1}(x) ≠ 0`,
then the previous edge indicator is forced to `1` by recurrence-local sign
separation. -/
theorem operatorSturm_prev_edge_indicator_one_of_lower_pattern_and_nonzero
    (M : ℕ) (x : ℝ)
    (hb1 : (if operatorCenterAt (M + 1) ≤ x then 1 else 0) = 1)
    (ha0 : (if operatorSturmSign (operatorSturmP (M + 1) x) ≠
                operatorSturmSign (operatorSturmP (M + 2) x) then 1 else 0) = 0)
    (hp1nz : operatorSturmP (M + 1) x ≠ 0) :
    (if operatorSturmSign (operatorSturmP M x) ≠
          operatorSturmSign (operatorSturmP (M + 1) x) then 1 else 0) = 1 := by
  have hcx : operatorCenterAt (M + 1) ≤ x := by
    by_cases hle : operatorCenterAt (M + 1) ≤ x
    · exact hle
    · simp [hle] at hb1
  have hnoflip : operatorSturmSign (operatorSturmP (M + 1) x) =
      operatorSturmSign (operatorSturmP (M + 2) x) := by
    by_cases hneq : operatorSturmSign (operatorSturmP (M + 1) x) ≠
        operatorSturmSign (operatorSturmP (M + 2) x)
    · simp [hneq] at ha0
    · exact not_not.mp hneq
  exact operatorSturm_prev_edge_indicator_one_of_center_le_and_noflip_nonzero
    M x hcx hnoflip hp1nz

/-- Finite sign-variation count on consecutive Sturm recurrence values
`p₀(x), p₁(x), ..., p_{M+1}(x)`. -/
def operatorSturmSignVariationCount (M : ℕ) (x : ℝ) : ℕ :=
  (((Finset.range (M + 1)).filter
      (fun k =>
        operatorSturmSign (operatorSturmP k x) ≠
          operatorSturmSign (operatorSturmP (k + 1) x))).card)

@[simp] theorem operatorSturmSign_p0_eq_one (x : ℝ) :
    operatorSturmSign (operatorSturmP 0 x) = 1 := by
  change operatorSturmSign 1 = 1
  simp [operatorSturmSign]

/-- Level-`0` Sturm sign-variation count is exactly the single-center indicator. -/
theorem operatorSturmSignVariationCount_zero_eq_center_indicator (x : ℝ) :
    operatorSturmSignVariationCount 0 x =
      (if operatorCenterAt 0 ≤ x then 1 else 0) := by
  have hbase :
      operatorSturmSignVariationCount 0 x
        = (if operatorSturmSign 1 = operatorSturmSign (operatorSturmDiag 0 - x)
            then 0 else 1) := by
    unfold operatorSturmSignVariationCount
    by_cases h : operatorSturmSign 1 = operatorSturmSign (operatorSturmDiag 0 - x)
    · simp [operatorSturmP, h]
    · let F : Finset ℕ :=
        (Finset.range (0 + 1)).filter
          (fun k =>
            operatorSturmSign (operatorSturmP k x) ≠
              operatorSturmSign (operatorSturmP (k + 1) x))
      have hmem0 : 0 ∈ F := by
        refine Finset.mem_filter.mpr ?_
        constructor
        · simp
        · simpa [operatorSturmP] using h
      have hcard_le : F.card ≤ 1 := by
        unfold F
        calc
          ((Finset.range (0 + 1)).filter
              (fun k =>
                operatorSturmSign (operatorSturmP k x) ≠
                  operatorSturmSign (operatorSturmP (k + 1) x))).card
              ≤ (Finset.range (0 + 1)).card := by
                simpa using Finset.card_filter_le
                  (s := Finset.range (0 + 1))
                  (p := fun k =>
                    operatorSturmSign (operatorSturmP k x) ≠
                      operatorSturmSign (operatorSturmP (k + 1) x))
          _ = 1 := by simp
      have hcard_pos : 0 < F.card := Finset.card_pos.mpr ⟨0, hmem0⟩
      have hcard_eq : F.card = 1 := by omega
      have hcard_eq_raw :
          ((Finset.range (0 + 1)).filter
              (fun k =>
                operatorSturmSign (operatorSturmP k x) ≠
                  operatorSturmSign (operatorSturmP (k + 1) x))).card = 1 := by
        simpa [F] using hcard_eq
      calc
        ((Finset.range (0 + 1)).filter
            (fun k =>
              operatorSturmSign (operatorSturmP k x) ≠
                operatorSturmSign (operatorSturmP (k + 1) x))).card = 1 := hcard_eq_raw
        _ = (if operatorSturmSign 1 = operatorSturmSign (operatorSturmDiag 0 - x)
              then 0 else 1) := by simp [h]
  rw [hbase]
  have hs1 : operatorSturmSign 1 = 1 := by simp [operatorSturmSign]
  by_cases hx : operatorCenterAt 0 ≤ x
  · have hdiag_le : operatorSturmDiag 0 ≤ x := by
      simpa [operatorSturmDiag_eq_centerAt] using hx
    have hnonpos : operatorSturmDiag 0 - x ≤ 0 := by linarith
    have hsDiag_ne_one : operatorSturmSign (operatorSturmDiag 0 - x) ≠ 1 :=
      operatorSturmSign_ne_one_of_nonpos hnonpos
    have hEqFalse : ¬ (operatorSturmSign 1 = operatorSturmSign (operatorSturmDiag 0 - x)) := by
      intro hEq
      have hsDiag_eq_one : operatorSturmSign (operatorSturmDiag 0 - x) = 1 := by
        calc
          operatorSturmSign (operatorSturmDiag 0 - x) = operatorSturmSign 1 := hEq.symm
          _ = 1 := hs1
      exact hsDiag_ne_one hsDiag_eq_one
    simp [hx, hEqFalse]
  · have hdiag_gt : x < operatorSturmDiag 0 := by
      have : x < operatorCenterAt 0 := lt_of_not_ge hx
      simpa [operatorSturmDiag_eq_centerAt] using this
    have hpos : 0 < operatorSturmDiag 0 - x := by linarith
    have hsDiag_eq_one : operatorSturmSign (operatorSturmDiag 0 - x) = 1 :=
      operatorSturmSign_eq_one_of_pos hpos
    have hEqTrue : operatorSturmSign 1 = operatorSturmSign (operatorSturmDiag 0 - x) := by
      calc
        operatorSturmSign 1 = 1 := hs1
        _ = operatorSturmSign (operatorSturmDiag 0 - x) := hsDiag_eq_one.symm
    simp [hx, hEqTrue]

theorem operatorSturmSignVariationCount_le (M : ℕ) (x : ℝ) :
    operatorSturmSignVariationCount M x ≤ M + 1 := by
  unfold operatorSturmSignVariationCount
  calc
    (((Finset.range (M + 1)).filter
        (fun k =>
          operatorSturmSign (operatorSturmP k x) ≠
            operatorSturmSign (operatorSturmP (k + 1) x))).card)
        ≤ (Finset.range (M + 1)).card := by
            simpa using (Finset.card_filter_le
              (s := Finset.range (M + 1))
              (p := fun k =>
                operatorSturmSign (operatorSturmP k x) ≠
                  operatorSturmSign (operatorSturmP (k + 1) x)))
    _ = M + 1 := by simp

/-- Appending one recurrence step can change the finite sign-variation count
by at most one. -/
theorem operatorSturmSignVariationCount_succ_le (M : ℕ) (x : ℝ) :
    operatorSturmSignVariationCount (M + 1) x ≤
      operatorSturmSignVariationCount M x + 1 := by
  classical
  let P : ℕ → Prop := fun k =>
    operatorSturmSign (operatorSturmP k x) ≠
      operatorSturmSign (operatorSturmP (k + 1) x)
  have hsplit :
      ((Finset.range (M + 2)).filter P).card
        ≤ ((Finset.range (M + 1)).filter P).card + 1 := by
    have hsub :
        ((Finset.range (M + 1)).filter P) ⊆ ((Finset.range (M + 2)).filter P) := by
      intro k hk
      have hk' := Finset.mem_filter.mp hk
      have hklt : k < M + 1 := Finset.mem_range.mp hk'.1
      have hklt' : k < M + 2 := by omega
      exact Finset.mem_filter.mpr ⟨Finset.mem_range.mpr hklt', hk'.2⟩
    have hadd :
        ((Finset.range (M + 2)).filter P)
          ⊆ insert (M + 1) ((Finset.range (M + 1)).filter P) := by
      intro k hk
      have hk' := Finset.mem_filter.mp hk
      have hklt : k < M + 2 := Finset.mem_range.mp hk'.1
      by_cases hEq : k = M + 1
      · exact Finset.mem_insert.mpr (Or.inl hEq)
      · have hlt : k < M + 1 := by omega
        exact Finset.mem_insert.mpr (Or.inr (Finset.mem_filter.mpr ⟨Finset.mem_range.mpr hlt, hk'.2⟩))
    calc
      ((Finset.range (M + 2)).filter P).card
          ≤ (insert (M + 1) ((Finset.range (M + 1)).filter P)).card := Finset.card_le_card hadd
      _ ≤ (((Finset.range (M + 1)).filter P).card) + 1 := by
            simpa [Nat.add_comm, Nat.add_left_comm, Nat.add_assoc] using Finset.card_insert_le (a := M + 1) (s := ((Finset.range (M + 1)).filter P))
  simpa [operatorSturmSignVariationCount, P]
    using hsplit

/-- Exact successor decomposition for Sturm sign-variation count:
appending one recurrence edge contributes exactly one extra variation iff the
new terminal edge flips sign. -/
theorem operatorSturmSignVariationCount_succ_eq_add_indicator (M : ℕ) (x : ℝ) :
    operatorSturmSignVariationCount (M + 1) x =
      operatorSturmSignVariationCount M x +
        (if operatorSturmSign (operatorSturmP (M + 1) x) ≠
              operatorSturmSign (operatorSturmP (M + 2) x) then 1 else 0) := by
  classical
  let P : ℕ → Prop := fun k =>
    operatorSturmSign (operatorSturmP k x) ≠
      operatorSturmSign (operatorSturmP (k + 1) x)
  have hnot_mem : M + 1 ∉ (Finset.range (M + 1)).filter P := by
    intro hmem
    exact (Nat.lt_irrefl (M + 1)) ((Finset.mem_filter.mp hmem).1
      |> Finset.mem_range.mp)
  rw [operatorSturmSignVariationCount, operatorSturmSignVariationCount]
  rw [Finset.range_add_one, Finset.filter_insert]
  by_cases hP : P (M + 1)
  · have hcard : (insert (M + 1) ((Finset.range (M + 1)).filter P)).card
        = ((Finset.range (M + 1)).filter P).card + 1 := by
      simpa [Nat.add_comm, Nat.add_left_comm, Nat.add_assoc] using
        Finset.card_insert_of_not_mem hnot_mem
    simpa [P, hP] using hcard
  · simp [P, hP]

/-- Recurrence-explicit form of the successor decomposition:
the new terminal-edge indicator is written using the concrete three-term recurrence
for `p_{M+2}`. -/
theorem operatorSturmSignVariationCount_succ_eq_add_indicator_recurrence (M : ℕ) (x : ℝ) :
    operatorSturmSignVariationCount (M + 1) x =
      operatorSturmSignVariationCount M x +
        (if operatorSturmSign (operatorSturmP (M + 1) x) ≠
              operatorSturmSign
                ((operatorSturmDiag (M + 1) - x) * operatorSturmP (M + 1) x
                  - (operatorSturmOff ^ (2 : ℕ)) * operatorSturmP M x)
         then 1 else 0) := by
  simpa [operatorSturmP_step] using
    operatorSturmSignVariationCount_succ_eq_add_indicator M x

/-- Gershgorin-window separation geometry used by the Sturm route:
if indices are at least `3` apart (forward), radius-`12/11` windows are strictly separated. -/
theorem operatorCenterWindow_separated_of_add_three_le
    {k l : ℕ} (hkl : k + 3 ≤ l) :
    operatorCenterAt k + (12 : ℝ) / 11 < operatorCenterAt l - (12 : ℝ) / 11 := by
  have hklR : (k : ℝ) + 3 ≤ (l : ℝ) := by
    exact_mod_cast hkl
  unfold operatorCenterAt
  linarith

/-- Finite-level eigenvalue counting function below threshold `x` on the concrete lane. -/
def operatorEigenvalueCountLE (M : ℕ) (x : ℝ) : ℕ :=
  Finset.card
    ((Finset.univ : Finset (Fin (M + 1))).filter
      (fun i => operatorEigenvalues M i ≤ x))

/-- Finite-level center counting function below threshold `x` for structural centers
`k + 29/16`. -/
def operatorCenterCountLE (M : ℕ) (x : ℝ) : ℕ :=
  Finset.card
    ((Finset.univ : Finset (Fin (M + 1))).filter
      (fun i => operatorCenterAt i.1 ≤ x))

/-- `operatorCenterCountLE` as a concrete filtered-cardinality over `range (M+1)`. -/
theorem operatorCenterCountLE_eq_range_filter_card (M : ℕ) (x : ℝ) :
    operatorCenterCountLE M x =
      ((Finset.range (M + 1)).filter (fun k => operatorCenterAt k ≤ x)).card := by
  unfold operatorCenterCountLE
  classical
  let Sfin : Finset (Fin (M + 1)) :=
    (Finset.univ : Finset (Fin (M + 1))).filter (fun i => operatorCenterAt i.1 ≤ x)
  let Snat : Finset ℕ := (Finset.range (M + 1)).filter (fun k => operatorCenterAt k ≤ x)
  have himage_eq : Sfin.image (fun i : Fin (M + 1) => i.1) = Snat := by
    ext k
    constructor
    · intro hk
      rcases Finset.mem_image.mp hk with ⟨i, hi, rfl⟩
      have hi' := Finset.mem_filter.mp hi
      exact Finset.mem_filter.mpr ⟨Finset.mem_range.mpr i.2, hi'.2⟩
    · intro hk
      have hk' := Finset.mem_filter.mp hk
      refine Finset.mem_image.mpr ?_
      refine ⟨⟨k, Finset.mem_range.mp hk'.1⟩, ?_, rfl⟩
      exact Finset.mem_filter.mpr ⟨by simp, hk'.2⟩
  have hcard_img :
      (Sfin.image (fun i : Fin (M + 1) => i.1)).card = Sfin.card := by
    simpa using Finset.card_image_of_injOn (s := Sfin) (f := fun i : Fin (M + 1) => i.1)
      (by intro a ha b hb hab; exact Fin.ext hab)
  calc
    Sfin.card = (Sfin.image (fun i : Fin (M + 1) => i.1)).card := hcard_img.symm
    _ = Snat.card := by simpa [himage_eq]

theorem operatorCenterCountLE_zero_eq_center_indicator (x : ℝ) :
    operatorCenterCountLE 0 x = (if operatorCenterAt 0 ≤ x then 1 else 0) := by
  rw [operatorCenterCountLE_eq_range_filter_card]
  let F : Finset ℕ := (Finset.range (0 + 1)).filter (fun k => operatorCenterAt k ≤ x)
  by_cases hx : operatorCenterAt 0 ≤ x
  · have hmem0 : 0 ∈ F := by
      unfold F
      simp [hx]
    have hcard_le : F.card ≤ 1 := by
      unfold F
      calc
        ((Finset.range (0 + 1)).filter (fun k => operatorCenterAt k ≤ x)).card
            ≤ (Finset.range (0 + 1)).card := by
              simpa using Finset.card_filter_le
                (s := Finset.range (0 + 1))
                (p := fun k => operatorCenterAt k ≤ x)
        _ = 1 := by simp
    have hcard_pos : 0 < F.card := Finset.card_pos.mpr ⟨0, hmem0⟩
    have hcard_eq : F.card = 1 := by omega
    have hraw :
        ((Finset.range (0 + 1)).filter (fun k => operatorCenterAt k ≤ x)).card = 1 := by
      simpa [F] using hcard_eq
    calc
      ((Finset.range (0 + 1)).filter (fun k => operatorCenterAt k ≤ x)).card = 1 := hraw
      _ = (if operatorCenterAt 0 ≤ x then 1 else 0) := by simp [hx]
  · have hcard_eq : F.card = 0 := by
      have hFempty : F = ∅ := by
        ext k
        constructor
        · intro hk
          have hk0 : k = 0 := by
            have hkRange : k ∈ Finset.range (0 + 1) := (Finset.mem_filter.mp hk).1
            have hklt : k < 1 := Finset.mem_range.mp hkRange
            omega
          subst hk0
          exact False.elim (hx ((Finset.mem_filter.mp hk).2))
        · intro hk
          exfalso
          simpa using hk
      simpa [hFempty]
    have hraw :
        ((Finset.range (0 + 1)).filter (fun k => operatorCenterAt k ≤ x)).card = 0 := by
      simpa [F] using hcard_eq
    calc
      ((Finset.range (0 + 1)).filter (fun k => operatorCenterAt k ≤ x)).card = 0 := hraw
      _ = (if operatorCenterAt 0 ≤ x then 1 else 0) := by simp [hx]

theorem operatorSturmSignVariationCount_centerCount_gap_zero (x : ℝ) :
    operatorSturmSignVariationCount 0 x ≤ operatorCenterCountLE 0 x + 1 ∧
    operatorCenterCountLE 0 x ≤ operatorSturmSignVariationCount 0 x + 1 := by
  rw [operatorSturmSignVariationCount_zero_eq_center_indicator,
    operatorCenterCountLE_zero_eq_center_indicator]
  by_cases hx : operatorCenterAt 0 ≤ x <;> simp [hx]

/-- Appending one structural center index changes center-count at most by one. -/
theorem operatorCenterCountLE_succ_le_add_one (M : ℕ) (x : ℝ) :
    operatorCenterCountLE (M + 1) x ≤ operatorCenterCountLE M x + 1 := by
  classical
  let P : ℕ → Prop := fun k => operatorCenterAt k ≤ x
  have hstep :
      ((Finset.range (M + 2)).filter P).card
        ≤ ((Finset.range (M + 1)).filter P).card + 1 := by
    have hadd :
        ((Finset.range (M + 2)).filter P)
          ⊆ insert (M + 1) ((Finset.range (M + 1)).filter P) := by
      intro k hk
      have hk' := Finset.mem_filter.mp hk
      have hklt : k < M + 2 := Finset.mem_range.mp hk'.1
      by_cases hEq : k = M + 1
      · exact Finset.mem_insert.mpr (Or.inl hEq)
      · have hlt : k < M + 1 := by omega
        exact Finset.mem_insert.mpr (Or.inr (Finset.mem_filter.mpr ⟨Finset.mem_range.mpr hlt, hk'.2⟩))
    calc
      ((Finset.range (M + 2)).filter P).card
          ≤ (insert (M + 1) ((Finset.range (M + 1)).filter P)).card := Finset.card_le_card hadd
      _ ≤ (((Finset.range (M + 1)).filter P).card) + 1 := by
            simpa [Nat.add_comm, Nat.add_left_comm, Nat.add_assoc]
              using Finset.card_insert_le (a := M + 1) (s := ((Finset.range (M + 1)).filter P))
  rw [operatorCenterCountLE_eq_range_filter_card, operatorCenterCountLE_eq_range_filter_card]
  simpa [P] using hstep

/-- Center-count monotonicity under level extension `M ↦ M+1`. -/
theorem operatorCenterCountLE_mono_succ (M : ℕ) (x : ℝ) :
    operatorCenterCountLE M x ≤ operatorCenterCountLE (M + 1) x := by
  classical
  rw [operatorCenterCountLE_eq_range_filter_card, operatorCenterCountLE_eq_range_filter_card]
  refine Finset.card_le_card ?_
  intro k hk
  have hk' := Finset.mem_filter.mp hk
  have hklt : k < M + 1 := Finset.mem_range.mp hk'.1
  have hklt' : k < M + 2 := by omega
  exact Finset.mem_filter.mpr ⟨Finset.mem_range.mpr hklt', hk'.2⟩

/-- One-step center-count increment is bounded by `1`. -/
theorem operatorCenterCountLE_succ_sub_le_one (M : ℕ) (x : ℝ) :
    operatorCenterCountLE (M + 1) x - operatorCenterCountLE M x ≤ 1 := by
  have hle := operatorCenterCountLE_succ_le_add_one M x
  omega

/-- Exact successor decomposition for center counting:
appending one center index contributes exactly one iff that center is below
threshold. -/
theorem operatorCenterCountLE_succ_eq_add_indicator (M : ℕ) (x : ℝ) :
    operatorCenterCountLE (M + 1) x =
      operatorCenterCountLE M x +
        (if operatorCenterAt (M + 1) ≤ x then 1 else 0) := by
  classical
  let P : ℕ → Prop := fun k => operatorCenterAt k ≤ x
  have hnot_mem : M + 1 ∉ (Finset.range (M + 1)).filter P := by
    intro hmem
    exact (Nat.lt_irrefl (M + 1)) ((Finset.mem_filter.mp hmem).1
      |> Finset.mem_range.mp)
  rw [operatorCenterCountLE_eq_range_filter_card, operatorCenterCountLE_eq_range_filter_card]
  rw [Finset.range_add_one, Finset.filter_insert]
  by_cases hP : P (M + 1)
  · have hcard : (insert (M + 1) ((Finset.range (M + 1)).filter P)).card
        = ((Finset.range (M + 1)).filter P).card + 1 := by
      simpa [Nat.add_comm, Nat.add_left_comm, Nat.add_assoc] using
        Finset.card_insert_of_not_mem hnot_mem
    simpa [P, hP] using hcard
  · simp [P, hP]

/-- Center count as a finite sum of 0/1 indicators over structural centers. -/
theorem operatorCenterCountLE_eq_sum_indicator (M : ℕ) (x : ℝ) :
    operatorCenterCountLE M x =
      Finset.sum (Finset.range (M + 1))
        (fun k => if operatorCenterAt k ≤ x then 1 else 0) := by
  rw [operatorCenterCountLE_eq_range_filter_card]
  rw [Finset.card_eq_sum_ones, Finset.sum_filter]

/-- Sturm sign-variation count as a finite sum of 0/1 edge-flip indicators. -/
theorem operatorSturmSignVariationCount_eq_sum_indicator (M : ℕ) (x : ℝ) :
    operatorSturmSignVariationCount M x =
      Finset.sum (Finset.range (M + 1))
        (fun k =>
          if operatorSturmSign (operatorSturmP k x) ≠
                operatorSturmSign (operatorSturmP (k + 1) x)
          then 1 else 0) := by
  rw [operatorSturmSignVariationCount]
  rw [Finset.card_eq_sum_ones, Finset.sum_filter]

/-- Sturm-route contract (corrected):
eigenvalue and center counting functions differ by at most one at every threshold
(`|Δcount| ≤ 1` in natural-number form). -/
def OperatorSturmCountContract : Prop :=
  ∀ M : ℕ, ∀ x : ℝ,
    operatorEigenvalueCountLE M x ≤ operatorCenterCountLE M x + 1 ∧
    operatorCenterCountLE M x ≤ operatorEigenvalueCountLE M x + 1

/-- Finite-step preservation of the `±1` gap between Sturm sign-variation count
and center count, provided the edge-lock implications hold at level `M`. -/
theorem sturmCenter_gap_step_preserved
    (M : ℕ) (x : ℝ)
    (hGapM :
      operatorSturmSignVariationCount M x ≤ operatorCenterCountLE M x + 1 ∧
      operatorCenterCountLE M x ≤ operatorSturmSignVariationCount M x + 1)
    (hEdgeLock :
      (operatorSturmSignVariationCount M x = operatorCenterCountLE M x + 1 →
        operatorCenterAt (M + 1) ≤ x) ∧
      (operatorCenterCountLE M x = operatorSturmSignVariationCount M x + 1 →
        operatorSturmSign (operatorSturmP (M + 1) x) ≠
          operatorSturmSign (operatorSturmP (M + 2) x))) :
    operatorSturmSignVariationCount (M + 1) x ≤ operatorCenterCountLE (M + 1) x + 1 ∧
    operatorCenterCountLE (M + 1) x ≤ operatorSturmSignVariationCount (M + 1) x + 1 := by
  let A : ℕ := operatorSturmSignVariationCount M x
  let B : ℕ := operatorCenterCountLE M x
  let a : ℕ :=
    if operatorSturmSign (operatorSturmP (M + 1) x) ≠
          operatorSturmSign (operatorSturmP (M + 2) x) then 1 else 0
  let b : ℕ := if operatorCenterAt (M + 1) ≤ x then 1 else 0
  have hA' : operatorSturmSignVariationCount (M + 1) x = A + a := by
    unfold A a
    simpa [operatorSturmP_step] using
      operatorSturmSignVariationCount_succ_eq_add_indicator_recurrence M x
  have hB' : operatorCenterCountLE (M + 1) x = B + b := by
    unfold B b
    simpa using operatorCenterCountLE_succ_eq_add_indicator M x
  have ha_le_one : a ≤ 1 := by
    unfold a
    split_ifs <;> omega
  have hb_le_one : b ≤ 1 := by
    unfold b
    split_ifs <;> omega
  have hA_le_B1 : A ≤ B + 1 := hGapM.1
  have hB_le_A1 : B ≤ A + 1 := hGapM.2
  have hA1_to_b1 :
      A = B + 1 → b = 1 := by
    intro hEq
    have hx : operatorCenterAt (M + 1) ≤ x := hEdgeLock.1 (by simpa [A, B] using hEq)
    unfold b
    simp [hx]
  have hB1_to_a1 :
      B = A + 1 → a = 1 := by
    intro hEq
    have hx :
        operatorSturmSign (operatorSturmP (M + 1) x) ≠
          operatorSturmSign (operatorSturmP (M + 2) x) :=
      hEdgeLock.2 (by simpa [A, B] using hEq)
    unfold a
    simp [hx]
  have hA'_le_B'1 : A + a ≤ B + b + 1 := by
    by_cases hAB : A ≤ B
    · have hA_add : A + a ≤ B + a := Nat.add_le_add_right hAB a
      have hA_add' : A + a ≤ B + 1 := le_trans hA_add (by omega)
      exact le_trans hA_add' (by omega)
    · have hBA_lt : B < A := Nat.lt_of_not_ge hAB
      have hA_eq : A = B + 1 := by omega
      have hb1 : b = 1 := hA1_to_b1 hA_eq
      have hA_add' : A + a ≤ B + 2 := by omega
      have hA_add'' : A + a ≤ B + b + 1 := by
        omega
      exact hA_add''
  have hB'_le_A'1 : B + b ≤ A + a + 1 := by
    by_cases hBA : B ≤ A
    · have hB_add : B + b ≤ A + b := Nat.add_le_add_right hBA b
      have hB_add' : B + b ≤ A + 1 := le_trans hB_add (by omega)
      exact le_trans hB_add' (by omega)
    · have hAB_lt : A < B := Nat.lt_of_not_ge hBA
      have hB_eq : B = A + 1 := by omega
      have ha1 : a = 1 := hB1_to_a1 hB_eq
      have hB_add' : B + b ≤ A + 2 := by omega
      have hB_add'' : B + b ≤ A + a + 1 := by
        omega
      exact hB_add''
  constructor
  · simpa [hA', hB', Nat.add_assoc, Nat.add_left_comm, Nat.add_comm] using hA'_le_B'1
  · simpa [hA', hB', Nat.add_assoc, Nat.add_left_comm, Nat.add_comm] using hB'_le_A'1

/-- Minimal step-compatibility form needed for `±1`-gap preservation:
at the upper boundary (`A = B + 1`), the Sturm increment cannot exceed the center
increment; symmetrically at the lower boundary (`B = A + 1`). -/
def OperatorSturmStepCompatibility : Prop :=
  ∀ M : ℕ, ∀ x : ℝ,
    let A := operatorSturmSignVariationCount M x
    let B := operatorCenterCountLE M x
    let a := (if operatorSturmSign (operatorSturmP (M + 1) x) ≠
                  operatorSturmSign (operatorSturmP (M + 2) x) then 1 else 0)
    let b := (if operatorCenterAt (M + 1) ≤ x then 1 else 0)
    (A = B + 1 → a ≤ b) ∧ (B = A + 1 → b ≤ a)

/-- Pointwise unpacking of `OperatorSturmStepCompatibility`:
at each `(M,x)`, the boundary conditions are exactly the two forbidden outward
increment patterns `(a,b) = (1,0)` and `(b,a) = (1,0)`. -/
theorem operatorSturmStepCompatibility_at_iff_forbidden_boundary_patterns
    (M : ℕ) (x : ℝ) :
    (let A := operatorSturmSignVariationCount M x
     let B := operatorCenterCountLE M x
     let a := (if operatorSturmSign (operatorSturmP (M + 1) x) ≠
                   operatorSturmSign (operatorSturmP (M + 2) x) then 1 else 0)
     let b := (if operatorCenterAt (M + 1) ≤ x then 1 else 0)
     (A = B + 1 → a ≤ b) ∧ (B = A + 1 → b ≤ a))
    ↔
    (let A := operatorSturmSignVariationCount M x
     let B := operatorCenterCountLE M x
     let a := (if operatorSturmSign (operatorSturmP (M + 1) x) ≠
                   operatorSturmSign (operatorSturmP (M + 2) x) then 1 else 0)
     let b := (if operatorCenterAt (M + 1) ≤ x then 1 else 0)
     (A = B + 1 → ¬ (a = 1 ∧ b = 0)) ∧
     (B = A + 1 → ¬ (b = 1 ∧ a = 0))) := by
  dsimp
  let A : ℕ := operatorSturmSignVariationCount M x
  let B : ℕ := operatorCenterCountLE M x
  let a : ℕ :=
    if operatorSturmSign (operatorSturmP (M + 1) x) ≠
          operatorSturmSign (operatorSturmP (M + 2) x) then 1 else 0
  let b : ℕ := if operatorCenterAt (M + 1) ≤ x then 1 else 0
  have ha_cases : a = 0 ∨ a = 1 := by
    unfold a
    split_ifs <;> simp
  have hb_cases : b = 0 ∨ b = 1 := by
    unfold b
    split_ifs <;> simp
  constructor
  · intro h
    constructor
    · intro hAB hab
      rcases hab with ⟨ha1, hb0⟩
      have hle : a ≤ b := h.1 hAB
      omega
    · intro hBA hba
      rcases hba with ⟨hb1, ha0⟩
      have hle : b ≤ a := h.2 hBA
      omega
  · intro h
    constructor
    · intro hAB
      rcases ha_cases with ha0 | ha1
      · have hle : a ≤ b := by
          omega
        simpa [a, b] using hle
      · have hb1 : b = 1 := by
          by_contra hb1ne
          have hb0 : b = 0 := by
            rcases hb_cases with hb0 | hb1
            · exact hb0
            · exfalso
              exact hb1ne hb1
          exact (h.1 hAB) ⟨ha1, hb0⟩
        have hle : a ≤ b := by
          omega
        simpa [a, b] using hle
    · intro hBA
      rcases hb_cases with hb0 | hb1
      · have hle : b ≤ a := by
          omega
        simpa [a, b] using hle
      · have ha1 : a = 1 := by
          by_contra ha1ne
          have ha0 : a = 0 := by
            rcases ha_cases with ha0 | ha1
            · exact ha0
            · exfalso
              exact ha1ne ha1
          exact (h.2 hBA) ⟨hb1, ha0⟩
        have hle : b ≤ a := by
          omega
        simpa [a, b] using hle

/-- Global reformulation of `OperatorSturmStepCompatibility` as forbidden boundary
patterns at every level and threshold. -/
theorem operatorSturmStepCompatibility_iff_forbidden_boundary_patterns :
    OperatorSturmStepCompatibility
      ↔
    (∀ M : ℕ, ∀ x : ℝ,
      let A := operatorSturmSignVariationCount M x
      let B := operatorCenterCountLE M x
      let a := (if operatorSturmSign (operatorSturmP (M + 1) x) ≠
                    operatorSturmSign (operatorSturmP (M + 2) x) then 1 else 0)
      let b := (if operatorCenterAt (M + 1) ≤ x then 1 else 0)
      (A = B + 1 → ¬ (a = 1 ∧ b = 0)) ∧
      (B = A + 1 → ¬ (b = 1 ∧ a = 0))) := by
  constructor
  · intro h M x
    exact (operatorSturmStepCompatibility_at_iff_forbidden_boundary_patterns M x).1 (h M x)
  · intro h M x
    exact (operatorSturmStepCompatibility_at_iff_forbidden_boundary_patterns M x).2 (h M x)

/-- Upper-boundary exclusion extracted from `OperatorSturmStepCompatibility`:
when `A = B + 1`, the outward pattern `(a,b) = (1,0)` is impossible. -/
theorem operatorSturm_forbid_upper_boundary_pattern_of_stepCompatibility
    (hCompat : OperatorSturmStepCompatibility) :
    ∀ M : ℕ, ∀ x : ℝ,
      let A := operatorSturmSignVariationCount M x
      let B := operatorCenterCountLE M x
      let a := (if operatorSturmSign (operatorSturmP (M + 1) x) ≠
                    operatorSturmSign (operatorSturmP (M + 2) x) then 1 else 0)
      let b := (if operatorCenterAt (M + 1) ≤ x then 1 else 0)
      A = B + 1 → ¬ (a = 1 ∧ b = 0) := by
  intro M x
  exact (operatorSturmStepCompatibility_iff_forbidden_boundary_patterns.mp hCompat M x).1

/-- Lower-boundary exclusion extracted from `OperatorSturmStepCompatibility`:
when `B = A + 1`, the outward pattern `(b,a) = (1,0)` is impossible. -/
theorem operatorSturm_forbid_lower_boundary_pattern_of_stepCompatibility
    (hCompat : OperatorSturmStepCompatibility) :
    ∀ M : ℕ, ∀ x : ℝ,
      let A := operatorSturmSignVariationCount M x
      let B := operatorCenterCountLE M x
      let a := (if operatorSturmSign (operatorSturmP (M + 1) x) ≠
                    operatorSturmSign (operatorSturmP (M + 2) x) then 1 else 0)
      let b := (if operatorCenterAt (M + 1) ≤ x then 1 else 0)
      B = A + 1 → ¬ (b = 1 ∧ a = 0) := by
  intro M x
  exact (operatorSturmStepCompatibility_iff_forbidden_boundary_patterns.mp hCompat M x).2

/-- Constructor: the two boundary exclusion lemmas are exactly enough to recover
`OperatorSturmStepCompatibility`. -/
theorem operatorSturmStepCompatibility_of_boundary_exclusions
    (hUpper :
      ∀ M : ℕ, ∀ x : ℝ,
        let A := operatorSturmSignVariationCount M x
        let B := operatorCenterCountLE M x
        let a := (if operatorSturmSign (operatorSturmP (M + 1) x) ≠
                      operatorSturmSign (operatorSturmP (M + 2) x) then 1 else 0)
        let b := (if operatorCenterAt (M + 1) ≤ x then 1 else 0)
        A = B + 1 → ¬ (a = 1 ∧ b = 0))
    (hLower :
      ∀ M : ℕ, ∀ x : ℝ,
        let A := operatorSturmSignVariationCount M x
        let B := operatorCenterCountLE M x
        let a := (if operatorSturmSign (operatorSturmP (M + 1) x) ≠
                      operatorSturmSign (operatorSturmP (M + 2) x) then 1 else 0)
        let b := (if operatorCenterAt (M + 1) ≤ x then 1 else 0)
        B = A + 1 → ¬ (b = 1 ∧ a = 0)) :
    OperatorSturmStepCompatibility := by
  apply (operatorSturmStepCompatibility_iff_forbidden_boundary_patterns).2
  intro M x
  exact ⟨hUpper M x, hLower M x⟩

/-- Upper-boundary exclusion from the global Sturm-route counting contract and
the eigenvalue↔Sturm sign-variation bridge. -/
theorem operatorSturm_forbid_upper_boundary_pattern_of_signVariationBridge_and_sturmContract
    (hEigSturm :
      ∀ M : ℕ, ∀ x : ℝ,
        operatorEigenvalueCountLE M x = operatorSturmSignVariationCount M x)
    (hS : OperatorSturmCountContract) :
    ∀ M : ℕ, ∀ x : ℝ,
      let A := operatorSturmSignVariationCount M x
      let B := operatorCenterCountLE M x
      let a := (if operatorSturmSign (operatorSturmP (M + 1) x) ≠
                    operatorSturmSign (operatorSturmP (M + 2) x) then 1 else 0)
      let b := (if operatorCenterAt (M + 1) ≤ x then 1 else 0)
      A = B + 1 → ¬ (a = 1 ∧ b = 0) := by
  intro M x
  dsimp
  let A : ℕ := operatorSturmSignVariationCount M x
  let B : ℕ := operatorCenterCountLE M x
  let a : ℕ :=
    if operatorSturmSign (operatorSturmP (M + 1) x) ≠
          operatorSturmSign (operatorSturmP (M + 2) x) then 1 else 0
  let b : ℕ := if operatorCenterAt (M + 1) ≤ x then 1 else 0
  have hA' : operatorSturmSignVariationCount (M + 1) x = A + a := by
    unfold A a
    simpa [operatorSturmP_step] using
      operatorSturmSignVariationCount_succ_eq_add_indicator_recurrence M x
  have hB' : operatorCenterCountLE (M + 1) x = B + b := by
    unfold B b
    simpa using operatorCenterCountLE_succ_eq_add_indicator M x
  intro hAB hab
  rcases hab with ⟨ha1, hb0⟩
  have ha1' : a = 1 := by simpa [a] using ha1
  have hb0' : b = 0 := by simpa [b] using hb0
  have hGapSucc :
      operatorSturmSignVariationCount (M + 1) x
        ≤ operatorCenterCountLE (M + 1) x + 1 := by
    have hEigBound : operatorEigenvalueCountLE (M + 1) x
        ≤ operatorCenterCountLE (M + 1) x + 1 := (hS (M + 1) x).1
    simpa [hEigSturm (M + 1) x] using hEigBound
  have hGapAB : A + a ≤ B + b + 1 := by
    simpa [hA', hB', Nat.add_assoc, Nat.add_left_comm, Nat.add_comm] using hGapSucc
  have hA_le_B : A ≤ B := by
    have hA_le_B' : A + 1 ≤ B + 1 := by
      calc
        A + 1 = A + a := by simpa [ha1']
        _ ≤ B + b + 1 := hGapAB
        _ = B + 1 := by
              have : b = 0 := hb0'
              omega
    omega
  omega

/-- Lower-boundary exclusion from the global Sturm-route counting contract and
the eigenvalue↔Sturm sign-variation bridge. -/
theorem operatorSturm_forbid_lower_boundary_pattern_of_signVariationBridge_and_sturmContract
    (hEigSturm :
      ∀ M : ℕ, ∀ x : ℝ,
        operatorEigenvalueCountLE M x = operatorSturmSignVariationCount M x)
    (hS : OperatorSturmCountContract) :
    ∀ M : ℕ, ∀ x : ℝ,
      let A := operatorSturmSignVariationCount M x
      let B := operatorCenterCountLE M x
      let a := (if operatorSturmSign (operatorSturmP (M + 1) x) ≠
                    operatorSturmSign (operatorSturmP (M + 2) x) then 1 else 0)
      let b := (if operatorCenterAt (M + 1) ≤ x then 1 else 0)
      B = A + 1 → ¬ (b = 1 ∧ a = 0) := by
  intro M x
  dsimp
  let A : ℕ := operatorSturmSignVariationCount M x
  let B : ℕ := operatorCenterCountLE M x
  let a : ℕ :=
    if operatorSturmSign (operatorSturmP (M + 1) x) ≠
          operatorSturmSign (operatorSturmP (M + 2) x) then 1 else 0
  let b : ℕ := if operatorCenterAt (M + 1) ≤ x then 1 else 0
  have hA' : operatorSturmSignVariationCount (M + 1) x = A + a := by
    unfold A a
    simpa [operatorSturmP_step] using
      operatorSturmSignVariationCount_succ_eq_add_indicator_recurrence M x
  have hB' : operatorCenterCountLE (M + 1) x = B + b := by
    unfold B b
    simpa using operatorCenterCountLE_succ_eq_add_indicator M x
  intro hBA hba
  rcases hba with ⟨hb1, ha0⟩
  have hb1' : b = 1 := by simpa [b] using hb1
  have ha0' : a = 0 := by simpa [a] using ha0
  have hGapSucc :
      operatorCenterCountLE (M + 1) x
        ≤ operatorSturmSignVariationCount (M + 1) x + 1 := by
    have hCenterEigBound : operatorCenterCountLE (M + 1) x
        ≤ operatorEigenvalueCountLE (M + 1) x + 1 := (hS (M + 1) x).2
    simpa [hEigSturm (M + 1) x] using hCenterEigBound
  have hGapBA : B + b ≤ A + a + 1 := by
    simpa [hA', hB', Nat.add_assoc, Nat.add_left_comm, Nat.add_comm] using hGapSucc
  have hB_le_A : B ≤ A := by
    have hB_le_A' : B + 1 ≤ A + 1 := by
      calc
        B + 1 = B + b := by simpa [hb1']
        _ ≤ A + a + 1 := hGapBA
        _ = A + 1 := by
              have : a = 0 := ha0'
              omega
    omega
  omega

/-- `OperatorSturmStepCompatibility` follows from the global corrected Sturm
counting contract once eigenvalue counting is bridged to Sturm sign variation. -/
theorem operatorSturmStepCompatibility_of_signVariationBridge_and_sturmContract
    (hEigSturm :
      ∀ M : ℕ, ∀ x : ℝ,
        operatorEigenvalueCountLE M x = operatorSturmSignVariationCount M x)
    (hS : OperatorSturmCountContract) :
    OperatorSturmStepCompatibility := by
  apply operatorSturmStepCompatibility_of_boundary_exclusions
  · exact operatorSturm_forbid_upper_boundary_pattern_of_signVariationBridge_and_sturmContract
      hEigSturm hS
  · exact operatorSturm_forbid_lower_boundary_pattern_of_signVariationBridge_and_sturmContract
      hEigSturm hS

/-- Strong edge-lock implies the minimal step-compatibility condition. -/
theorem operatorSturmStepCompatibility_of_edgeLock
    (hEdgeLock :
      ∀ M : ℕ, ∀ x : ℝ,
        (operatorSturmSignVariationCount M x = operatorCenterCountLE M x + 1 →
          operatorCenterAt (M + 1) ≤ x) ∧
        (operatorCenterCountLE M x = operatorSturmSignVariationCount M x + 1 →
          operatorSturmSign (operatorSturmP (M + 1) x) ≠
            operatorSturmSign (operatorSturmP (M + 2) x))) :
    OperatorSturmStepCompatibility := by
  intro M x
  dsimp
  constructor
  · intro hAB
    have hx : operatorCenterAt (M + 1) ≤ x := (hEdgeLock M x).1 hAB
    split_ifs <;> omega
  · intro hBA
    have hx :
        operatorSturmSign (operatorSturmP (M + 1) x) ≠
          operatorSturmSign (operatorSturmP (M + 2) x) := (hEdgeLock M x).2 hBA
    split_ifs <;> omega

/-- One-step `±1` gap preservation from local boundary exclusions.
This is the induction-step core: assuming `|A-B| ≤ 1` at level `M`, the next
gap stays within `±1` provided the two outward boundary patterns are excluded
at this `(M,x)`. -/
theorem sturmCenter_gap_step_preserved_of_boundary_exclusions_split_upper
    (M : ℕ) (x : ℝ)
    (hGapM :
      operatorSturmSignVariationCount M x ≤ operatorCenterCountLE M x + 1 ∧
      operatorCenterCountLE M x ≤ operatorSturmSignVariationCount M x + 1)
    (hUpperPrev :
      let A := operatorSturmSignVariationCount M x
      let B := operatorCenterCountLE M x
      let a := (if operatorSturmSign (operatorSturmP (M + 1) x) ≠
                    operatorSturmSign (operatorSturmP (M + 2) x) then 1 else 0)
      let b := (if operatorCenterAt (M + 1) ≤ x then 1 else 0)
      let prev := (if operatorSturmSign (operatorSturmP M x) ≠
                     operatorSturmSign (operatorSturmP (M + 1) x) then 1 else 0)
      A = B + 1 → prev = 0 → ¬ (a = 1 ∧ b = 0))
    (hUpperZero :
      let A := operatorSturmSignVariationCount M x
      let B := operatorCenterCountLE M x
      let a := (if operatorSturmSign (operatorSturmP (M + 1) x) ≠
                    operatorSturmSign (operatorSturmP (M + 2) x) then 1 else 0)
      let b := (if operatorCenterAt (M + 1) ≤ x then 1 else 0)
      A = B + 1 → operatorSturmP (M + 1) x = 0 → ¬ (a = 1 ∧ b = 0))
    (hLowerPrev :
      let A := operatorSturmSignVariationCount M x
      let B := operatorCenterCountLE M x
      let a := (if operatorSturmSign (operatorSturmP (M + 1) x) ≠
                    operatorSturmSign (operatorSturmP (M + 2) x) then 1 else 0)
      let b := (if operatorCenterAt (M + 1) ≤ x then 1 else 0)
      let prev := (if operatorSturmSign (operatorSturmP M x) ≠
                     operatorSturmSign (operatorSturmP (M + 1) x) then 1 else 0)
      B = A + 1 → prev = 1 → ¬ (b = 1 ∧ a = 0))
    (hLowerZero :
      let A := operatorSturmSignVariationCount M x
      let B := operatorCenterCountLE M x
      let a := (if operatorSturmSign (operatorSturmP (M + 1) x) ≠
                    operatorSturmSign (operatorSturmP (M + 2) x) then 1 else 0)
      let b := (if operatorCenterAt (M + 1) ≤ x then 1 else 0)
      B = A + 1 → operatorSturmP (M + 1) x = 0 → ¬ (b = 1 ∧ a = 0)) :
    operatorSturmSignVariationCount (M + 1) x ≤ operatorCenterCountLE (M + 1) x + 1 ∧
    operatorCenterCountLE (M + 1) x ≤ operatorSturmSignVariationCount (M + 1) x + 1 := by
  let A : ℕ := operatorSturmSignVariationCount M x
  let B : ℕ := operatorCenterCountLE M x
  let a : ℕ :=
    if operatorSturmSign (operatorSturmP (M + 1) x) ≠
          operatorSturmSign (operatorSturmP (M + 2) x) then 1 else 0
  let b : ℕ := if operatorCenterAt (M + 1) ≤ x then 1 else 0
  have hA' : operatorSturmSignVariationCount (M + 1) x = A + a := by
    unfold A a
    simpa [operatorSturmP_step] using
      operatorSturmSignVariationCount_succ_eq_add_indicator_recurrence M x
  have hB' : operatorCenterCountLE (M + 1) x = B + b := by
    unfold B b
    simpa using operatorCenterCountLE_succ_eq_add_indicator M x
  have ha_cases : a = 0 ∨ a = 1 := by
    unfold a
    split_ifs <;> simp
  have hb_cases : b = 0 ∨ b = 1 := by
    unfold b
    split_ifs <;> simp
  have hA'_le_B'1 : A + a ≤ B + b + 1 := by
    by_cases hAB : A ≤ B
    · have hA_add : A + a ≤ B + a := Nat.add_le_add_right hAB a
      have hA_add' : A + a ≤ B + 1 := le_trans hA_add (by omega)
      exact le_trans hA_add' (by omega)
    · have hA_eq : A = B + 1 := by omega
      have hab : a ≤ b := by
        rcases ha_cases with ha0 | ha1
        · omega
        · rcases hb_cases with hb0 | hb1
          · have hPrev0_of_nonzero :
                operatorSturmP (M + 1) x ≠ 0 →
                  (if operatorSturmSign (operatorSturmP M x) ≠
                        operatorSturmSign (operatorSturmP (M + 1) x) then 1 else 0) = 0 := by
              intro hp1nz
              exact operatorSturm_prev_edge_indicator_zero_of_upper_pattern_and_nonzero
                M x (by simpa [a] using ha1) (by simpa [b] using hb0) hp1nz
            by_cases hp1nz : operatorSturmP (M + 1) x ≠ 0
            · exfalso
              exact (hUpperPrev hA_eq (hPrev0_of_nonzero hp1nz)) ⟨ha1, hb0⟩
            · exfalso
              have hp1z : operatorSturmP (M + 1) x = 0 := by
                by_contra hp1ne
                exact hp1nz hp1ne
              exact (hUpperZero hA_eq hp1z) ⟨ha1, hb0⟩
          · omega
      have hA_add' : A + a ≤ B + b + 1 := by omega
      exact hA_add'
  have hB'_le_A'1 : B + b ≤ A + a + 1 := by
    by_cases hBA : B ≤ A
    · have hB_add : B + b ≤ A + b := Nat.add_le_add_right hBA b
      have hB_add' : B + b ≤ A + 1 := le_trans hB_add (by omega)
      exact le_trans hB_add' (by omega)
    · have hB_eq : B = A + 1 := by omega
      have hba : b ≤ a := by
        rcases hb_cases with hb0 | hb1
        · omega
        · rcases ha_cases with ha0 | ha1
          · exfalso
            have hPrev1_of_nonzero :
                  operatorSturmP (M + 1) x ≠ 0 →
                    (if operatorSturmSign (operatorSturmP M x) ≠
                          operatorSturmSign (operatorSturmP (M + 1) x) then 1 else 0) = 1 := by
                intro hp1nz
                exact operatorSturm_prev_edge_indicator_one_of_lower_pattern_and_nonzero
                  M x (by simpa [b] using hb1) (by simpa [a] using ha0) hp1nz
            by_cases hp1nz : operatorSturmP (M + 1) x ≠ 0
            · exfalso
              exact (hLowerPrev hB_eq (hPrev1_of_nonzero hp1nz)) ⟨hb1, ha0⟩
            · exfalso
              have hp1z : operatorSturmP (M + 1) x = 0 := by
                by_contra hp1ne
                exact hp1nz hp1ne
              exact (hLowerZero hB_eq hp1z) ⟨hb1, ha0⟩
          · omega
      have hB_add' : B + b ≤ A + a + 1 := by omega
      exact hB_add'
  constructor
  · simpa [hA', hB', Nat.add_assoc, Nat.add_left_comm, Nat.add_comm] using hA'_le_B'1
  · simpa [hA', hB', Nat.add_assoc, Nat.add_left_comm, Nat.add_comm] using hB'_le_A'1

/-- One-step `±1` gap preservation from local boundary exclusions.
This is the induction-step core: assuming `|A-B| ≤ 1` at level `M`, the next
gap stays within `±1` provided the two outward boundary patterns are excluded
at this `(M,x)`. -/
theorem sturmCenter_gap_step_preserved_of_boundary_exclusions
    (M : ℕ) (x : ℝ)
    (hGapM :
      operatorSturmSignVariationCount M x ≤ operatorCenterCountLE M x + 1 ∧
      operatorCenterCountLE M x ≤ operatorSturmSignVariationCount M x + 1)
    (hUpper :
      let A := operatorSturmSignVariationCount M x
      let B := operatorCenterCountLE M x
      let a := (if operatorSturmSign (operatorSturmP (M + 1) x) ≠
                    operatorSturmSign (operatorSturmP (M + 2) x) then 1 else 0)
      let b := (if operatorCenterAt (M + 1) ≤ x then 1 else 0)
      A = B + 1 → ¬ (a = 1 ∧ b = 0))
    (hLower :
      let A := operatorSturmSignVariationCount M x
      let B := operatorCenterCountLE M x
      let a := (if operatorSturmSign (operatorSturmP (M + 1) x) ≠
                    operatorSturmSign (operatorSturmP (M + 2) x) then 1 else 0)
      let b := (if operatorCenterAt (M + 1) ≤ x then 1 else 0)
      B = A + 1 → ¬ (b = 1 ∧ a = 0)) :
    operatorSturmSignVariationCount (M + 1) x ≤ operatorCenterCountLE (M + 1) x + 1 ∧
    operatorCenterCountLE (M + 1) x ≤ operatorSturmSignVariationCount (M + 1) x + 1 := by
  exact sturmCenter_gap_step_preserved_of_boundary_exclusions_split_upper M x hGapM
    (fun Aeq prevEq hab => (hUpper Aeq) hab)
    (fun Aeq hp1z hab => (hUpper Aeq) hab)
    (fun Beq prevEq hba => (hLower Beq) hba)
    (fun Beq hp1z hba => (hLower Beq) hba)

/-- Finite-step preservation of the `±1` gap under minimal step-compatibility. -/
theorem sturmCenter_gap_step_preserved_of_stepCompatibility
    (M : ℕ) (x : ℝ)
    (hGapM :
      operatorSturmSignVariationCount M x ≤ operatorCenterCountLE M x + 1 ∧
      operatorCenterCountLE M x ≤ operatorSturmSignVariationCount M x + 1)
    (hCompat : OperatorSturmStepCompatibility) :
    operatorSturmSignVariationCount (M + 1) x ≤ operatorCenterCountLE (M + 1) x + 1 ∧
    operatorCenterCountLE (M + 1) x ≤ operatorSturmSignVariationCount (M + 1) x + 1 := by
  exact sturmCenter_gap_step_preserved_of_boundary_exclusions M x hGapM
    (operatorSturm_forbid_upper_boundary_pattern_of_stepCompatibility hCompat M x)
    (operatorSturm_forbid_lower_boundary_pattern_of_stepCompatibility hCompat M x)

/-- Reduction theorem: if eigenvalue count equals Sturm sign-variation count and
the edge-lock implications hold at every level, then the Sturm-route contract
(`|Δcount| ≤ 1`) follows for all thresholds. -/
theorem operatorSturmCountContract_of_signVariationBridge_and_edgeLock
    (hEigSturm :
      ∀ M : ℕ, ∀ x : ℝ,
        operatorEigenvalueCountLE M x = operatorSturmSignVariationCount M x)
    (hEdgeLock :
      ∀ M : ℕ, ∀ x : ℝ,
        (operatorSturmSignVariationCount M x = operatorCenterCountLE M x + 1 →
          operatorCenterAt (M + 1) ≤ x) ∧
        (operatorCenterCountLE M x = operatorSturmSignVariationCount M x + 1 →
          operatorSturmSign (operatorSturmP (M + 1) x) ≠
            operatorSturmSign (operatorSturmP (M + 2) x))) :
    OperatorSturmCountContract := by
  intro M x
  have hGapSC : operatorSturmSignVariationCount M x ≤ operatorCenterCountLE M x + 1 ∧
      operatorCenterCountLE M x ≤ operatorSturmSignVariationCount M x + 1 := by
    induction' M with M ih
    · simpa using operatorSturmSignVariationCount_centerCount_gap_zero x
    · exact sturmCenter_gap_step_preserved M x ih (hEdgeLock M x)
  constructor
  · simpa [hEigSturm M x] using hGapSC.1
  · simpa [hEigSturm M x] using hGapSC.2

/-- Same reduction as above, but with the minimal step-compatibility obligation. -/
theorem operatorSturmCountContract_of_signVariationBridge_and_stepCompatibility
    (hEigSturm :
      ∀ M : ℕ, ∀ x : ℝ,
        operatorEigenvalueCountLE M x = operatorSturmSignVariationCount M x)
    (hCompat : OperatorSturmStepCompatibility) :
    OperatorSturmCountContract := by
  intro M x
  have hGapSC : operatorSturmSignVariationCount M x ≤ operatorCenterCountLE M x + 1 ∧
      operatorCenterCountLE M x ≤ operatorSturmSignVariationCount M x + 1 := by
    induction' M with M ih
    · simpa using operatorSturmSignVariationCount_centerCount_gap_zero x
    · exact sturmCenter_gap_step_preserved_of_stepCompatibility M x ih hCompat
  constructor
  · simpa [hEigSturm M x] using hGapSC.1
  · simpa [hEigSturm M x] using hGapSC.2

/-- Direct induction route from local boundary exclusions:
if eigenvalue counting is bridged to Sturm sign-variation counting, and at each
`(M,x)` the two outward boundary patterns are excluded, then the corrected
Sturm-route counting contract follows globally. -/
theorem operatorSturmCountContract_of_signVariationBridge_and_boundary_exclusions
    (hEigSturm :
      ∀ M : ℕ, ∀ x : ℝ,
        operatorEigenvalueCountLE M x = operatorSturmSignVariationCount M x)
    (hBoundary :
      ∀ M : ℕ, ∀ x : ℝ,
        (let A := operatorSturmSignVariationCount M x
         let B := operatorCenterCountLE M x
         let a := (if operatorSturmSign (operatorSturmP (M + 1) x) ≠
                       operatorSturmSign (operatorSturmP (M + 2) x) then 1 else 0)
         let b := (if operatorCenterAt (M + 1) ≤ x then 1 else 0)
         (A = B + 1 → ¬ (a = 1 ∧ b = 0)) ∧
         (B = A + 1 → ¬ (b = 1 ∧ a = 0)))) :
    OperatorSturmCountContract := by
  intro M x
  have hGapSC : operatorSturmSignVariationCount M x ≤ operatorCenterCountLE M x + 1 ∧
      operatorCenterCountLE M x ≤ operatorSturmSignVariationCount M x + 1 := by
    induction' M with M ih
    · simpa using operatorSturmSignVariationCount_centerCount_gap_zero x
    · exact sturmCenter_gap_step_preserved_of_boundary_exclusions M x ih
        (hBoundary M x).1 (hBoundary M x).2
  constructor
  · simpa [hEigSturm M x] using hGapSC.1
  · simpa [hEigSturm M x] using hGapSC.2

/-- Under the eigenvalue↔Sturm bridge, the corrected Sturm-route contract and
minimal step-compatibility are equivalent. -/
theorem operatorSturmCountContract_iff_stepCompatibility_of_signVariationBridge
    (hEigSturm :
      ∀ M : ℕ, ∀ x : ℝ,
        operatorEigenvalueCountLE M x = operatorSturmSignVariationCount M x) :
    OperatorSturmCountContract ↔ OperatorSturmStepCompatibility := by
  constructor
  · intro hS
    exact operatorSturmStepCompatibility_of_signVariationBridge_and_sturmContract hEigSturm hS
  · intro hCompat
    exact operatorSturmCountContract_of_signVariationBridge_and_stepCompatibility hEigSturm hCompat

/-- Structural center-window occupancy bound:
for fixed center index `k`, at most three structural centers can lie in the
radius-`12/11` window around `operatorCenterAt k` (namely indices `k-1,k,k+1`
after clipping to the finite level). -/
theorem operatorCenterWindow_card_le_three
    (M : ℕ) (k : Fin (M + 1)) :
    ((Finset.univ : Finset (Fin (M + 1))).filter
      (fun i => operatorCenterAt k.1 - (12 : ℝ) / 11 < operatorCenterAt i.1 ∧
        operatorCenterAt i.1 ≤ operatorCenterAt k.1 + (12 : ℝ) / 11)).card ≤ 3 := by
  classical
  let Sfin : Finset (Fin (M + 1)) :=
    (Finset.univ : Finset (Fin (M + 1))).filter
      (fun i => operatorCenterAt k.1 - (12 : ℝ) / 11 < operatorCenterAt i.1 ∧
        operatorCenterAt i.1 ≤ operatorCenterAt k.1 + (12 : ℝ) / 11)
  let Snat : Finset ℕ := (Finset.range (M + 1)).filter
    (fun n => k.1 - 1 ≤ n ∧ n ≤ k.1 + 1)
  have hcard_img : (Sfin.image (fun i : Fin (M + 1) => i.1)).card = Sfin.card := by
    simpa using Finset.card_image_of_injOn (s := Sfin) (f := fun i : Fin (M + 1) => i.1)
      (by intro a ha b hb hab; exact Fin.ext hab)
  have hsubset_img : (Sfin.image (fun i : Fin (M + 1) => i.1)) ⊆ Snat := by
    intro n hn
    rcases Finset.mem_image.mp hn with ⟨i, hiS, rfl⟩
    have hiS' := Finset.mem_filter.mp hiS
    have hiRange : i.1 < M + 1 := i.2
    have hlow : operatorCenterAt k.1 - (12 : ℝ) / 11 < operatorCenterAt i.1 := hiS'.2.1
    have hup : operatorCenterAt i.1 ≤ operatorCenterAt k.1 + (12 : ℝ) / 11 := hiS'.2.2
    have hlt2 : (i.1 : ℝ) < (k.1 : ℝ) + 2 := by
      unfold operatorCenterAt at hup
      linarith
    have hi_upper : i.1 ≤ k.1 + 1 := Nat.lt_succ_iff.mp (by exact_mod_cast hlt2)
    have hklt2 : (k.1 : ℝ) < (i.1 : ℝ) + 2 := by
      unfold operatorCenterAt at hlow
      linarith
    have hk_le : k.1 ≤ i.1 + 1 := Nat.lt_succ_iff.mp (by exact_mod_cast hklt2)
    have hi_lower : k.1 - 1 ≤ i.1 := by omega
    exact Finset.mem_filter.mpr ⟨Finset.mem_range.mpr hiRange, ⟨hi_lower, hi_upper⟩⟩
  have hSnat_card_le_three : Snat.card ≤ 3 := by
    let T : Finset ℕ := {k.1 - 1, k.1, k.1 + 1}
    have hsubset_T : Snat ⊆ T := by
      intro n hn
      have hn' := Finset.mem_filter.mp hn
      have hkn : k.1 - 1 ≤ n := hn'.2.1
      have hnk : n ≤ k.1 + 1 := hn'.2.2
      have hcases : n = k.1 - 1 ∨ n = k.1 ∨ n = k.1 + 1 := by
        omega
      rcases hcases with rfl | rfl | rfl <;> simp [T]
    exact le_trans (Finset.card_le_card hsubset_T)
      (by simpa [T] using (Finset.card_le_three (a := k.1 - 1) (b := k.1) (c := k.1 + 1)))
  calc
    ((Finset.univ : Finset (Fin (M + 1))).filter
        (fun i => operatorCenterAt k.1 - (12 : ℝ) / 11 < operatorCenterAt i.1 ∧
          operatorCenterAt i.1 ≤ operatorCenterAt k.1 + (12 : ℝ) / 11)).card
        = Sfin.card := by rfl
    _ = (Sfin.image (fun i : Fin (M + 1) => i.1)).card := hcard_img.symm
    _ ≤ Snat.card := Finset.card_le_card hsubset_img
    _ ≤ 3 := hSnat_card_le_three

/-- Center counting-function window bound:
the finite-level center counting function can increase by at most `3` across the
`±12/11` window around any structural center. -/
theorem operatorCenterCountLE_window_sub_le_three
    (M : ℕ) (k : Fin (M + 1)) :
    operatorCenterCountLE M (operatorCenterAt k.1 + (12 : ℝ) / 11) -
      operatorCenterCountLE M (operatorCenterAt k.1 - (12 : ℝ) / 11) ≤ 3 := by
  classical
  let upper : ℝ := operatorCenterAt k.1 + (12 : ℝ) / 11
  let lower : ℝ := operatorCenterAt k.1 - (12 : ℝ) / 11
  let U : Finset (Fin (M + 1)) :=
    (Finset.univ : Finset (Fin (M + 1))).filter (fun i => operatorCenterAt i.1 ≤ upper)
  let L : Finset (Fin (M + 1)) :=
    (Finset.univ : Finset (Fin (M + 1))).filter (fun i => operatorCenterAt i.1 ≤ lower)
  have hlower_le_upper : lower ≤ upper := by
    unfold lower upper
    linarith
  have hsubset : L ⊆ U := by
    intro i hiL
    have hiL' := Finset.mem_filter.mp hiL
    exact Finset.mem_filter.mpr ⟨hiL'.1, le_trans hiL'.2 hlower_le_upper⟩
  have hU_sdiff_L :
      (U \ L) =
        ((Finset.univ : Finset (Fin (M + 1))).filter
          (fun i => lower < operatorCenterAt i.1 ∧ operatorCenterAt i.1 ≤ upper)) := by
    ext i
    constructor
    · intro hi
      have hiU : i ∈ U := (Finset.mem_sdiff.mp hi).1
      have hiNotL : i ∉ L := (Finset.mem_sdiff.mp hi).2
      have hiU' := Finset.mem_filter.mp hiU
      have hiUpper : operatorCenterAt i.1 ≤ upper := hiU'.2
      have hiLowerNot : ¬ operatorCenterAt i.1 ≤ lower := by
        intro hle
        exact hiNotL (Finset.mem_filter.mpr ⟨by simpa using hiU'.1, hle⟩)
      have hiLower : lower < operatorCenterAt i.1 := lt_of_not_ge hiLowerNot
      exact Finset.mem_filter.mpr ⟨by simpa using hiU'.1, ⟨hiLower, hiUpper⟩⟩
    · intro hi
      have hi' := Finset.mem_filter.mp hi
      have hiLower : lower < operatorCenterAt i.1 := hi'.2.1
      have hiUpper : operatorCenterAt i.1 ≤ upper := hi'.2.2
      have hiU : i ∈ U := Finset.mem_filter.mpr ⟨by simpa using hi'.1, hiUpper⟩
      have hiNotL : i ∉ L := by
        intro hiL
        have hiL' := Finset.mem_filter.mp hiL
        exact (not_le_of_gt hiLower) hiL'.2
      exact Finset.mem_sdiff.mpr ⟨hiU, hiNotL⟩
  have hcard_sub :
      operatorCenterCountLE M upper - operatorCenterCountLE M lower = (U \ L).card := by
    unfold operatorCenterCountLE U L
    have hcard : (U \ L).card = U.card - (L ∩ U).card := by
      simpa [Finset.inter_comm] using (Finset.card_sdiff (s := L) (t := U))
    have hinter : L ∩ U = L := by
      exact Finset.inter_eq_left.mpr hsubset
    rw [hcard, hinter]
  have hwindow_card :
      ((Finset.univ : Finset (Fin (M + 1))).filter
        (fun i => lower < operatorCenterAt i.1 ∧ operatorCenterAt i.1 ≤ upper)).card ≤ 3 := by
    simpa [lower, upper] using operatorCenterWindow_card_le_three M k
  calc
    operatorCenterCountLE M upper - operatorCenterCountLE M lower
        = (U \ L).card := hcard_sub
    _ = ((Finset.univ : Finset (Fin (M + 1))).filter
          (fun i => lower < operatorCenterAt i.1 ∧ operatorCenterAt i.1 ≤ upper)).card := by
          simpa [hU_sdiff_L]
    _ ≤ 3 := hwindow_card

/-- Exact half-window center count jump:
across the structural window of radius `1/2` around center `k+29/16`,
the center counting function increases by exactly `1`. -/
theorem operatorCenterCountLE_halfWindow_sub_eq_one
    (M : ℕ) (k : Fin (M + 1)) :
    operatorCenterCountLE M (operatorCenterAt k.1 + (1 : ℝ) / 2) -
      operatorCenterCountLE M (operatorCenterAt k.1 - (1 : ℝ) / 2) = 1 := by
  classical
  let upper : ℝ := operatorCenterAt k.1 + (1 : ℝ) / 2
  let lower : ℝ := operatorCenterAt k.1 - (1 : ℝ) / 2
  let U : Finset (Fin (M + 1)) :=
    (Finset.univ : Finset (Fin (M + 1))).filter (fun i => operatorCenterAt i.1 ≤ upper)
  let L : Finset (Fin (M + 1)) :=
    (Finset.univ : Finset (Fin (M + 1))).filter (fun i => operatorCenterAt i.1 ≤ lower)
  have hlower_le_upper : lower ≤ upper := by
    unfold lower upper
    linarith
  have hsubset : L ⊆ U := by
    intro i hiL
    have hiL' := Finset.mem_filter.mp hiL
    exact Finset.mem_filter.mpr ⟨hiL'.1, le_trans hiL'.2 hlower_le_upper⟩
  have hU_sdiff_L :
      (U \ L) =
        ((Finset.univ : Finset (Fin (M + 1))).filter
          (fun i => lower < operatorCenterAt i.1 ∧ operatorCenterAt i.1 ≤ upper)) := by
    ext i
    constructor
    · intro hi
      have hiU : i ∈ U := (Finset.mem_sdiff.mp hi).1
      have hiNotL : i ∉ L := (Finset.mem_sdiff.mp hi).2
      have hiU' := Finset.mem_filter.mp hiU
      have hiUpper : operatorCenterAt i.1 ≤ upper := hiU'.2
      have hiLowerNot : ¬ operatorCenterAt i.1 ≤ lower := by
        intro hle
        exact hiNotL (Finset.mem_filter.mpr ⟨by simpa using hiU'.1, hle⟩)
      have hiLower : lower < operatorCenterAt i.1 := lt_of_not_ge hiLowerNot
      exact Finset.mem_filter.mpr ⟨by simpa using hiU'.1, ⟨hiLower, hiUpper⟩⟩
    · intro hi
      have hi' := Finset.mem_filter.mp hi
      have hiLower : lower < operatorCenterAt i.1 := hi'.2.1
      have hiUpper : operatorCenterAt i.1 ≤ upper := hi'.2.2
      have hiU : i ∈ U := Finset.mem_filter.mpr ⟨by simpa using hi'.1, hiUpper⟩
      have hiNotL : i ∉ L := by
        intro hiL
        have hiL' := Finset.mem_filter.mp hiL
        exact (not_le_of_gt hiLower) hiL'.2
      exact Finset.mem_sdiff.mpr ⟨hiU, hiNotL⟩
  have hcard_sub :
      operatorCenterCountLE M upper - operatorCenterCountLE M lower = (U \ L).card := by
    unfold operatorCenterCountLE U L
    have hcard : (U \ L).card = U.card - (L ∩ U).card := by
      simpa [Finset.inter_comm] using (Finset.card_sdiff (s := L) (t := U))
    have hinter : L ∩ U = L := by
      exact Finset.inter_eq_left.mpr hsubset
    rw [hcard, hinter]
  have hwindow_singleton :
      ((Finset.univ : Finset (Fin (M + 1))).filter
        (fun i => lower < operatorCenterAt i.1 ∧ operatorCenterAt i.1 ≤ upper)) = ({k} : Finset (Fin (M + 1))) := by
    ext i
    constructor
    · intro hi
      have hi' := Finset.mem_filter.mp hi
      have hlow : lower < operatorCenterAt i.1 := hi'.2.1
      have hup : operatorCenterAt i.1 ≤ upper := hi'.2.2
      have hi_le_k : i.1 ≤ k.1 := by
        have hlt : (i.1 : ℝ) < (k.1 : ℝ) + 1 := by
          dsimp [lower, upper, operatorCenterAt] at hlow hup ⊢
          linarith
        exact Nat.lt_succ_iff.mp (by exact_mod_cast hlt)
      have hk_le_i : k.1 ≤ i.1 := by
        have hlt : (k.1 : ℝ) < (i.1 : ℝ) + 1 := by
          dsimp [lower, upper, operatorCenterAt] at hlow hup ⊢
          linarith
        exact Nat.lt_succ_iff.mp (by exact_mod_cast hlt)
      have hik : i = k := Fin.ext (Nat.le_antisymm hi_le_k hk_le_i)
      simpa [hik]
    · intro hi
      have hik : i = k := by simpa using hi
      subst i
      have hlow : lower < operatorCenterAt k.1 := by
        dsimp [lower, operatorCenterAt]
        linarith
      have hup : operatorCenterAt k.1 ≤ upper := by
        dsimp [upper, operatorCenterAt]
        linarith
      exact Finset.mem_filter.mpr ⟨Finset.mem_univ k, ⟨hlow, hup⟩⟩
  calc
    operatorCenterCountLE M (operatorCenterAt k.1 + (1 : ℝ) / 2) -
      operatorCenterCountLE M (operatorCenterAt k.1 - (1 : ℝ) / 2)
        = operatorCenterCountLE M upper - operatorCenterCountLE M lower := by rfl
    _ = (U \ L).card := hcard_sub
    _ = ((Finset.univ : Finset (Fin (M + 1))).filter
          (fun i => lower < operatorCenterAt i.1 ∧ operatorCenterAt i.1 ≤ upper)).card := by
          simpa [hU_sdiff_L]
    _ = ({k} : Finset (Fin (M + 1))).card := by simpa [hwindow_singleton]
    _ = 1 := by simp

/-- Sturm-route half-window jump bound (corrected contract):
if `|Δcount| ≤ 1` at every threshold, then each structural half-window captures
at most three eigenvalue count increments. -/
theorem operatorEigenvalueCountLE_halfWindow_sub_le_three_of_sturm
    (hS : OperatorSturmCountContract) (M : ℕ) (k : Fin (M + 1)) :
    operatorEigenvalueCountLE M (operatorCenterAt k.1 + (1 : ℝ) / 2) -
      operatorEigenvalueCountLE M (operatorCenterAt k.1 - (1 : ℝ) / 2) ≤ 3 := by
  let upper : ℝ := operatorCenterAt k.1 + (1 : ℝ) / 2
  let lower : ℝ := operatorCenterAt k.1 - (1 : ℝ) / 2
  let eu : ℕ := operatorEigenvalueCountLE M upper
  let el : ℕ := operatorEigenvalueCountLE M lower
  let cu : ℕ := operatorCenterCountLE M upper
  let cl : ℕ := operatorCenterCountLE M lower
  change eu - el ≤ 3
  have hSupper : eu ≤ cu + 1 := (hS M upper).1
  have hSlower : cl ≤ el + 1 := (hS M lower).2
  have hcenter1 : cu - cl = 1 := by
    unfold cu cl upper lower
    simpa using operatorCenterCountLE_halfWindow_sub_eq_one M k
  have hcu_eq : cu = cl + 1 := by
    omega
  have heu_le_el3 : eu ≤ el + 3 := by
    calc
      eu ≤ cu + 1 := hSupper
      _ = cl + 2 := by simpa [hcu_eq, Nat.add_assoc]
      _ ≤ el + 3 := by
            calc
              cl + 2 ≤ (el + 1) + 2 := Nat.add_le_add_right hSlower 2
              _ = el + 3 := by omega
  exact (Nat.sub_le_iff_le_add).2 (by
    simpa [Nat.add_assoc, Nat.add_left_comm, Nat.add_comm] using heu_le_el3)

/-- Sturm-route counting consequence (corrected contract):
if `|Δcount| ≤ 1` at every threshold, then the eigenvalue counting function
increases by at most `5` across the structural `±12/11` window around each center. -/
theorem operatorEigenvalueCountLE_window_sub_le_three_of_sturm
    (hS : OperatorSturmCountContract) (M : ℕ) (k : Fin (M + 1)) :
    operatorEigenvalueCountLE M (operatorCenterAt k.1 + (12 : ℝ) / 11) -
      operatorEigenvalueCountLE M (operatorCenterAt k.1 - (12 : ℝ) / 11) ≤ 5 := by
  let upper : ℝ := operatorCenterAt k.1 + (12 : ℝ) / 11
  let lower : ℝ := operatorCenterAt k.1 - (12 : ℝ) / 11
  let eu : ℕ := operatorEigenvalueCountLE M upper
  let el : ℕ := operatorEigenvalueCountLE M lower
  let cu : ℕ := operatorCenterCountLE M upper
  let cl : ℕ := operatorCenterCountLE M lower
  change eu - el ≤ 5
  have hSupper : eu ≤ cu + 1 := (hS M upper).1
  have hSlower : cl ≤ el + 1 := (hS M lower).2
  have hcenter0 : cu - cl ≤ 3 := by
    unfold cu cl upper lower
    simpa using operatorCenterCountLE_window_sub_le_three M k
  have hcu_le : cu ≤ cl + 3 := by
    omega
  have heu_le_el5 : eu ≤ el + 5 := by
    calc
      eu ≤ cu + 1 := hSupper
      _ ≤ (cl + 3) + 1 := Nat.add_le_add_right hcu_le 1
      _ = cl + 4 := by omega
      _ ≤ el + 5 := by
            calc
              cl + 4 ≤ (el + 1) + 4 := Nat.add_le_add_right hSlower 4
              _ = el + 5 := by omega
  exact (Nat.sub_le_iff_le_add).2 (by
    simpa [Nat.add_assoc, Nat.add_left_comm, Nat.add_comm] using heu_le_el5)

/-- The Weyl center-gap contract implies indexed linear growth of operator eigenvalues. -/
theorem indexed_linear_growth_of_operatorWeylCenterGap
    (hW : OperatorWeylCenterGap) :
    ∀ M : ℕ, ∀ i : Fin (M + 1),
      (((i.1 + 1 : ℕ) : ℝ) / 2) ≤ operatorEigenvalues M i := by
  intro M i
  have hgap : |operatorEigenvalues M i - ((i.1 : ℝ) + (29 : ℝ) / 16)| ≤ (12 : ℝ) / 11 := hW M i
  have hpair := abs_le.mp hgap
  have hlower :
      (i.1 : ℝ) + (127 : ℝ) / 176 ≤ operatorEigenvalues M i := by
    have hconst : (i.1 : ℝ) + (29 : ℝ) / 16 - (12 : ℝ) / 11
        = (i.1 : ℝ) + (127 : ℝ) / 176 := by
      ring_nf
    linarith [hpair.1, hconst]
  have hhalf_le :
      (((i.1 + 1 : ℕ) : ℝ) / 2) ≤ (i.1 : ℝ) + (127 : ℝ) / 176 := by
    have hi_nonneg : (0 : ℝ) ≤ i.1 := by positivity
    have hle_half : (((i.1 + 1 : ℕ) : ℝ) / 2) ≤ (i.1 : ℝ) + (1 : ℝ) / 2 := by
      calc
        (((i.1 + 1 : ℕ) : ℝ) / 2)
            = ((i.1 : ℝ) + 1) / 2 := by norm_num [Nat.cast_add]
        _ = (i.1 : ℝ) / 2 + (1 : ℝ) / 2 := by ring
        _ ≤ (i.1 : ℝ) + (1 : ℝ) / 2 := by nlinarith
    have hconst : (1 : ℝ) / 2 ≤ (127 : ℝ) / 176 := by norm_num
    have hle_const : (i.1 : ℝ) + (1 : ℝ) / 2 ≤ (i.1 : ℝ) + (127 : ℝ) / 176 := by
      simpa [add_comm, add_left_comm, add_assoc] using add_le_add_left hconst (i.1 : ℝ)
    exact le_trans hle_half hle_const
  exact le_trans hhalf_le hlower

/-- Concrete operator-derived finite ladder:
real eigenvalues of the structural Hermitian matrix at each level `N+1`. -/
noncomputable def operatorSpecN (N : ℕ) : Finset ℝ :=
  Finset.univ.image (operatorEigenvalues N)

/-- Every ordinate listed in `operatorSpecN N` lies in the real spectrum of
the structural matrix at level `N+1`. -/
theorem mem_operatorSpecN_implies_mem_real_spectrum
    {N : ℕ} {t : ℝ} (ht : t ∈ operatorSpecN N) :
    t ∈ spectrum ℝ (structuralRiemannMatrixC (N + 1)) := by
  rcases Finset.mem_image.mp ht with ⟨i, _hi, rfl⟩
  exact Matrix.IsHermitian.eigenvalues_mem_spectrum_real
    (hA := structuralRiemannMatrixC_isHermitian (N + 1)) i

/-- Exact membership characterization for the concrete operator ladder:
listed ordinates are exactly spectral ordinates of the structural matrix. -/
theorem mem_operatorSpecN_iff_ordinateIsEigenvalue
    (N : ℕ) (t : ℝ) :
    t ∈ operatorSpecN N ↔ ordinateIsEigenvalue (N + 1) t := by
  constructor
  · intro ht
    have hreal : t ∈ spectrum ℝ (structuralRiemannMatrixC (N + 1)) :=
      mem_operatorSpecN_implies_mem_real_spectrum ht
    exact (spectrum.algebraMap_mem_iff ℂ).2 hreal
  · intro ht
    have hreal : t ∈ spectrum ℝ (structuralRiemannMatrixC (N + 1)) :=
      (spectrum.algebraMap_mem_iff ℂ).1 ht
    have hrange :
        t ∈ Set.range ((structuralRiemannMatrixC_isHermitian (N + 1)).eigenvalues) := by
      simpa [Matrix.IsHermitian.spectrum_real_eq_range_eigenvalues
        (hA := structuralRiemannMatrixC_isHermitian (N + 1))] using hreal
    rcases hrange with ⟨i, rfl⟩
    exact Finset.mem_image.mpr ⟨i, Finset.mem_univ i, rfl⟩

/-- Concrete operator ladder membership is exactly operator-spectrum set membership. -/
theorem mem_operatorSpecN_iff_mem_operatorSpecSet
    (N : ℕ) (t : ℝ) :
    t ∈ operatorSpecN N ↔ t ∈ operatorSpecSet N := by
  exact (mem_operatorSpecN_iff_ordinateIsEigenvalue N t).trans
    (ordinateIsEigenvalue_iff_mem_operatorSpecSet N t)

/-- Finite-level operator ladder symmetry:
membership is closed under reflection about the structural center. -/
theorem mem_operatorSpecN_reflect_structuralCenter
    {N : ℕ} {t : ℝ} (ht : t ∈ operatorSpecN N) :
    ((structuralCenterQ (N + 1) : ℝ) - t) ∈ operatorSpecN N := by
  have htReal : t ∈ spectrum ℝ (structuralRiemannMatrixC (N + 1)) :=
    mem_operatorSpecN_implies_mem_real_spectrum ht
  have htComplex : (t : ℂ) ∈ spectrum ℂ (structuralRiemannMatrixC (N + 1)) :=
    (spectrum.algebraMap_mem_iff ℂ).2 htReal
  have hReflectComplex :
      ((structuralCenterQ (N + 1) : ℂ) - (t : ℂ)) ∈
        spectrum ℂ (structuralRiemannMatrixC (N + 1)) :=
    structuralRiemannMatrixC_spectrum_reflect (N + 1) (t : ℂ) htComplex
  have hReflectEigen :
      ordinateIsEigenvalue (N + 1) ((structuralCenterQ (N + 1) : ℝ) - t) :=
    (spectrum.algebraMap_mem_iff ℂ).1 (by simpa using hReflectComplex)
  exact (mem_operatorSpecN_iff_ordinateIsEigenvalue N
    ((structuralCenterQ (N + 1) : ℝ) - t)).2 hReflectEigen

/-- Set-level version of finite reflection symmetry on the concrete operator lane. -/
theorem mem_operatorSpecSet_reflect_structuralCenter
    {N : ℕ} {t : ℝ} (ht : t ∈ operatorSpecSet N) :
    ((structuralCenterQ (N + 1) : ℝ) - t) ∈ operatorSpecSet N := by
  have htN : t ∈ operatorSpecN N :=
    (mem_operatorSpecN_iff_mem_operatorSpecSet N t).2 ht
  have hRefN :
      ((structuralCenterQ (N + 1) : ℝ) - t) ∈ operatorSpecN N :=
    mem_operatorSpecN_reflect_structuralCenter htN
  exact (mem_operatorSpecN_iff_mem_operatorSpecSet N
    ((structuralCenterQ (N + 1) : ℝ) - t)).1 hRefN

/-- The finite operator ladder is invariant under reflection about the
structural center at each level `N`. -/
theorem operatorSpecN_image_reflect_structuralCenter (N : ℕ) :
    (operatorSpecN N).image (fun t => (structuralCenterQ (N + 1) : ℝ) - t) =
      operatorSpecN N := by
  ext t
  constructor
  · intro ht
    rcases Finset.mem_image.mp ht with ⟨u, hu, rfl⟩
    exact mem_operatorSpecN_reflect_structuralCenter hu
  · intro ht
    let u : ℝ := (structuralCenterQ (N + 1) : ℝ) - t
    have hu : u ∈ operatorSpecN N := by
      exact mem_operatorSpecN_reflect_structuralCenter ht
    refine Finset.mem_image.mpr ?_
    refine ⟨u, hu, ?_⟩
    dsimp [u]
    ring

/-- Set-level invariance of `operatorSpecSet` under center reflection. -/
theorem operatorSpecSet_image_reflect_structuralCenter (N : ℕ) :
    (operatorSpecSet N).image (fun t => (structuralCenterQ (N + 1) : ℝ) - t) =
      operatorSpecSet N := by
  ext t
  constructor
  · intro ht
    rcases ht with ⟨u, hu, rfl⟩
    exact mem_operatorSpecSet_reflect_structuralCenter hu
  · intro ht
    let u : ℝ := (structuralCenterQ (N + 1) : ℝ) - t
    have hu : u ∈ operatorSpecSet N := by
      exact mem_operatorSpecSet_reflect_structuralCenter ht
    refine ⟨u, hu, ?_⟩
    dsimp [u]
    ring

/-- Reflection around a fixed real center is injective. -/
theorem reflectAround_injective (c : ℝ) :
    Function.Injective (fun t : ℝ => c - t) := by
  intro a b hab
  linarith

/-- Finite product identity on the concrete operator ladder: replacing the
ordinate set by its reflected image rewrites each factor by reflected ordinates. -/
theorem XiFinite_operatorSpecN_reflectedFactors
    (N : ℕ) (s : ℂ) :
    XiFinite (operatorSpecN N) s =
      Finset.prod (operatorSpecN N)
        (fun t => s - criticalLinePoint ((structuralCenterQ (N + 1) : ℝ) - t)) := by
  classical
  let c : ℝ := structuralCenterQ (N + 1)
  have hEq : operatorSpecN N = (operatorSpecN N).image (fun t => c - t) := by
    simpa [c] using (operatorSpecN_image_reflect_structuralCenter N).symm
  calc
    XiFinite (operatorSpecN N) s
        = XiFinite ((operatorSpecN N).image (fun t => c - t)) s := by
            exact congrArg (fun spec => XiFinite spec s) hEq
    _ = Finset.prod ((operatorSpecN N).image (fun t => c - t))
          (fun u => s - criticalLinePoint u) := rfl
    _ = Finset.prod (operatorSpecN N) (fun t => s - criticalLinePoint (c - t)) := by
          simpa using
            (Finset.prod_image (s := operatorSpecN N)
              (g := fun t : ℝ => c - t)
              (f := fun u : ℝ => s - criticalLinePoint u)
              (reflectAround_injective c).injOn)
    _ = Finset.prod (operatorSpecN N) (fun t => s - criticalLinePoint (c - t)) := rfl
    _ = Finset.prod (operatorSpecN N)
          (fun t => s - criticalLinePoint ((structuralCenterQ (N + 1) : ℝ) - t)) := by
          simp [c]

/-- Centered affine involution in the `s`-plane induced by the finite structural
center at level `N`. -/
def operatorCenterInvolutionArg (N : ℕ) (s : ℂ) : ℂ :=
  ((1 : ℂ) + (structuralCenterQ (N + 1) : ℂ) * Complex.I) - s

/-- The centered affine map is an involution. -/
theorem operatorCenterInvolutionArg_involutive (N : ℕ) :
    Function.Involutive (operatorCenterInvolutionArg N) := by
  intro s
  simp [operatorCenterInvolutionArg, sub_eq_add_neg, add_assoc, add_left_comm, add_comm]

/-- Finite centered involution identity on the concrete operator ladder:
reflecting `s` about the structural center line multiplies the finite product by
the parity phase `(-1)^|spec_N|`. -/
theorem XiFinite_operatorSpecN_centered_involution
    (N : ℕ) (s : ℂ) :
    XiFinite (operatorSpecN N) (operatorCenterInvolutionArg N s) =
      (-1 : ℂ) ^ (operatorSpecN N).card * XiFinite (operatorSpecN N) s := by
  classical
  let c : ℝ := structuralCenterQ (N + 1)
  calc
    XiFinite (operatorSpecN N) (operatorCenterInvolutionArg N s)
        = Finset.prod (operatorSpecN N)
            (fun t => (((1 : ℂ) + (c : ℂ) * Complex.I) - s - criticalLinePoint t)) := by
              simp [XiFinite, operatorCenterInvolutionArg, c]
    _ = Finset.prod (operatorSpecN N)
          (fun t => ((-1 : ℂ) * (s - criticalLinePoint (c - t)))) := by
          refine Finset.prod_congr rfl ?_
          intro t ht
          simp [criticalLinePoint, c]
          ring
    _ = (Finset.prod (operatorSpecN N) (fun _t => (-1 : ℂ))) *
          Finset.prod (operatorSpecN N) (fun t => (s - criticalLinePoint (c - t))) := by
          simpa using
            (Finset.prod_mul_distrib
              (s := operatorSpecN N)
              (f := fun _ : ℝ => (-1 : ℂ))
              (g := fun t : ℝ => s - criticalLinePoint (c - t)))
    _ = (-1 : ℂ) ^ (operatorSpecN N).card *
          Finset.prod (operatorSpecN N) (fun t => (s - criticalLinePoint (c - t))) := by
          simp
    _ = (-1 : ℂ) ^ (operatorSpecN N).card * XiFinite (operatorSpecN N) s := by
          have href :
              XiFinite (operatorSpecN N) s =
                Finset.prod (operatorSpecN N) (fun t => (s - criticalLinePoint (c - t))) := by
            simpa [c] using (XiFinite_operatorSpecN_reflectedFactors N s)
          rw [← href]

/-- Midpoint of the centered involution in the `s`-plane at finite level `N`. -/
def operatorCenterMidpoint (N : ℕ) : ℂ :=
  ((1 : ℂ) + (structuralCenterQ (N + 1) : ℂ) * Complex.I) / 2

/-- `(-1)^k` depends only on parity (`k mod 2`). -/
theorem negOne_pow_eq_negOne_pow_mod2 (k : ℕ) :
    (-1 : ℂ) ^ k = (-1 : ℂ) ^ (k % 2) := by
  have hkdecomp : k = k % 2 + 2 * (k / 2) := by
    simpa [Nat.add_comm, Nat.mul_comm] using (Nat.mod_add_div k 2).symm
  calc
    (-1 : ℂ) ^ k = (-1 : ℂ) ^ (k % 2 + 2 * (k / 2)) := by
      exact congrArg (fun n : ℕ => (-1 : ℂ) ^ n) hkdecomp
    _ = (-1 : ℂ) ^ (2 * (k / 2)) * (-1 : ℂ) ^ (k % 2) := by
      rw [Nat.add_comm, pow_add]
    _ = (-1 : ℂ) ^ (k % 2) := by
      simp [pow_mul]

/-- Mod-2 parity square for `-1` is always `1`. -/
theorem negOne_pow_mod2_mul_self (k : ℕ) :
    (-1 : ℂ) ^ (k % 2) * (-1 : ℂ) ^ (k % 2) = 1 := by
  rcases Nat.mod_two_eq_zero_or_one k with h0 | h1
  · simp [h0]
  · simp [h1]

/-- Parity cancellation identity:
the mod-2 correction exactly cancels the finite-level parity phase. -/
theorem negOne_pow_mod2_mul_pow_card_cancel (k : ℕ) :
    (-1 : ℂ) ^ (k % 2) * (-1 : ℂ) ^ k = 1 := by
  have hk : (-1 : ℂ) ^ k = (-1 : ℂ) ^ (k % 2) :=
    negOne_pow_eq_negOne_pow_mod2 k
  calc
    (-1 : ℂ) ^ (k % 2) * (-1 : ℂ) ^ k
        = (-1 : ℂ) ^ (k % 2) * (-1 : ℂ) ^ (k % 2) := by rw [hk]
    _ = 1 := negOne_pow_mod2_mul_self k

/-- Operator-native finite ladder witness:
each finite level list is exactly the real ordinates detected as operator
eigenvalues for that level. -/
structure OperatorSpecLadder where
  specN : ℕ → Finset ℝ
  eigenvalue_exact :
    ∀ N : ℕ, ∀ t : ℝ, t ∈ specN N ↔ ordinateIsEigenvalue (N + 1) t

/-- Canonical operator ladder witness coming directly from `operatorSpecN`. -/
def canonicalOperatorSpecLadder : OperatorSpecLadder where
  specN := operatorSpecN
  eigenvalue_exact := mem_operatorSpecN_iff_ordinateIsEigenvalue

/-- Obligation-1 in operator-native form. -/
def OperatorNontrivialCapture (hOp : OperatorSpecLadder) : Prop :=
  RiemannNontrivialLadderZeroCapture hOp.specN

/-- Direct operator-set capture surface (no external reference data). -/
def OperatorSetNontrivialCapture : Prop :=
  ∀ s : ℂ, riemannZeta s = 0 →
    (¬ ∃ n : ℕ, s = -2 * (n + 1)) →
    s ≠ 1 →
    ∃ N : ℕ, ∃ t : ℝ, t ∈ operatorSpecSet N ∧ s = criticalLinePoint t

/-- Concrete finite-ladder capture obligation specialized to the
operator-constructed `operatorSpecN`. -/
def OperatorEnumeratedNontrivialCapture : Prop :=
  RiemannNontrivialLadderZeroCapture operatorSpecN

/-- Operator-native equivalence of endgame capture surfaces. -/
theorem operator_nontrivial_capture_iff_xiTarget_capture
    (hOp : OperatorSpecLadder) :
    OperatorNontrivialCapture hOp ↔ XiTargetLadderZeroCapture hOp.specN := by
  exact nontrivial_capture_iff_xiTarget_capture hOp.specN

/-- Operator-native RH closure:
once nontrivial-`ζ` capture is established for the operator ladder, RH follows. -/
theorem mathlibRH_of_operator_nontrivial_capture
    (hOp : OperatorSpecLadder)
    (hCap : OperatorNontrivialCapture hOp) :
    RiemannHypothesis := by
  exact mathlibRH_of_nontrivial_capture hOp.specN hCap

/-- RH closure directly from nontrivial capture over the concrete operator
eigenvalue ladder `operatorSpecN`. -/
theorem mathlibRH_of_operator_enumerated_nontrivial_capture
    (hCap : OperatorEnumeratedNontrivialCapture) :
    RiemannHypothesis := by
  exact mathlibRH_of_nontrivial_capture operatorSpecN hCap

/-- Concrete operator finite-product ladder used in the convergence lane. -/
def operatorXiFiniteLadder : ℕ → (ℂ → ℂ) :=
  XiFiniteLadder operatorSpecN

/-- Concrete operator Hadamard-normalized finite-product ladder. -/
def operatorXiFiniteHadamardLadder : ℕ → (ℂ → ℂ) :=
  fun N => XiFiniteHadamard (operatorSpecN N)

/-- One-step Hadamard increment bound specialized to the concrete operator lane
when adding a fresh ordinate to the level-`N` spectrum. -/
theorem norm_operatorXiFiniteHadamard_insert_sub_le
    (N : ℕ) {t : ℝ} (ht : t ∉ operatorSpecN N) (s : ℂ) :
    ‖XiFiniteHadamard (insert t (operatorSpecN N)) s
      - operatorXiFiniteHadamardLadder N s‖ ≤
      ‖hadamardFactor t s - 1‖ * ‖operatorXiFiniteHadamardLadder N s‖ := by
  simpa [operatorXiFiniteHadamardLadder] using
    (norm_XiFiniteHadamard_insert_sub_le (operatorSpecN N) ht s)

/-- One-step Hadamard increment bound with explicit second-order control,
specialized to the concrete operator lane. -/
theorem norm_operatorXiFiniteHadamard_insert_sub_le_three_mul_sq
    (N : ℕ) {t : ℝ} (ht : t ∉ operatorSpecN N) (s : ℂ)
    (hz : ‖s / criticalLinePoint t‖ ≤ 1) :
    ‖XiFiniteHadamard (insert t (operatorSpecN N)) s
      - operatorXiFiniteHadamardLadder N s‖ ≤
      (3 * ‖s / criticalLinePoint t‖ ^ 2) * ‖operatorXiFiniteHadamardLadder N s‖ := by
  simpa [operatorXiFiniteHadamardLadder] using
    (norm_XiFiniteHadamard_insert_sub_le_three_mul_sq (operatorSpecN N) ht s hz)

/-- Canonical quadratic majorant profile for Hadamard-step increments. -/
def hadamardQuadraticProfile (R M c : ℝ) : ℕ → ℝ :=
  fun n => (3 * (R / (c * (n + 1 : ℝ))) ^ (2 : ℕ)) * M

/-- Summability of the canonical quadratic Hadamard profile. -/
theorem summable_hadamardQuadraticProfile
    (R M c : ℝ) (hc : c ≠ 0) :
    Summable (hadamardQuadraticProfile R M c) := by
  unfold hadamardQuadraticProfile
  have hbase : Summable (fun n : ℕ => (1 : ℝ) / ((n : ℝ) ^ (2 : ℕ))) :=
    (Real.summable_one_div_nat_pow).2 (by norm_num)
  have hshift : Summable (fun n : ℕ => (1 : ℝ) / (((n + 1 : ℕ) : ℝ) ^ (2 : ℕ))) := by
    simpa [Nat.cast_add, add_assoc, add_comm, add_left_comm] using
      (summable_nat_add_iff 1).2 hbase
  have hcongr :
      (fun n : ℕ => (3 * (R / (c * (n + 1 : ℝ))) ^ (2 : ℕ)) * M)
        = (fun n : ℕ =>
            ((3 * (R / c) ^ (2 : ℕ)) * M) *
              ((1 : ℝ) / (((n + 1 : ℕ) : ℝ) ^ (2 : ℕ)))) := by
    funext n
    have hn1 : ((n + 1 : ℕ) : ℝ) ≠ 0 := by
      exact_mod_cast Nat.succ_ne_zero n
    have hnc : c * (((n + 1 : ℕ) : ℝ)) ≠ 0 := mul_ne_zero hc hn1
    field_simp [hc, hn1, hnc]
    norm_num [Nat.cast_add] at *
    ring_nf
  rw [hcongr]
  exact hshift.mul_left ((3 * (R / c) ^ (2 : ℕ)) * M)

/-- API-level centered involution identity for the concrete operator Xi ladder. -/
theorem operatorXiFiniteLadder_centered_involution
    (N : ℕ) (s : ℂ) :
    operatorXiFiniteLadder N (operatorCenterInvolutionArg N s) =
      (-1 : ℂ) ^ (operatorSpecN N).card * operatorXiFiniteLadder N s := by
  simpa [operatorXiFiniteLadder, XiFiniteLadder] using
    (XiFinite_operatorSpecN_centered_involution N s)

/-- Centered involution preserves the norm of finite operator Xi products. -/
theorem norm_operatorXiFiniteLadder_centered_involution
    (N : ℕ) (s : ℂ) :
    ‖operatorXiFiniteLadder N (operatorCenterInvolutionArg N s)‖ =
      ‖operatorXiFiniteLadder N s‖ := by
  calc
    ‖operatorXiFiniteLadder N (operatorCenterInvolutionArg N s)‖
        = ‖(-1 : ℂ) ^ (operatorSpecN N).card * operatorXiFiniteLadder N s‖ := by
            rw [operatorXiFiniteLadder_centered_involution]
    _ = ‖(-1 : ℂ) ^ (operatorSpecN N).card‖ * ‖operatorXiFiniteLadder N s‖ := norm_mul _ _
    _ = ‖operatorXiFiniteLadder N s‖ := by simp

/-- Phase-normalized finite operator Xi ladder:
the parity correction removes the `(-1)^|specN|` phase from centered involution. -/
def operatorPhaseNormalizedXiLadder (N : ℕ) (s : ℂ) : ℂ :=
  (s - operatorCenterMidpoint N) ^ ((operatorSpecN N).card % 2) * operatorXiFiniteLadder N s

/-- Strict centered involution after phase normalization. -/
theorem operatorPhaseNormalizedXiLadder_centered_involution
    (N : ℕ) (s : ℂ) :
    operatorPhaseNormalizedXiLadder N (operatorCenterInvolutionArg N s) =
      operatorPhaseNormalizedXiLadder N s := by
  let k : ℕ := (operatorSpecN N).card
  let p : ℕ := k % 2
  let m : ℂ := operatorCenterMidpoint N
  have haff :
      operatorCenterInvolutionArg N s - m = (-1 : ℂ) * (s - m) := by
    simp [operatorCenterInvolutionArg, operatorCenterMidpoint, m]
    ring
  have hpow :
      (operatorCenterInvolutionArg N s - m) ^ p =
        (-1 : ℂ) ^ p * (s - m) ^ p := by
    rw [haff, mul_pow]
  have hXi :
      operatorXiFiniteLadder N (operatorCenterInvolutionArg N s) =
        (-1 : ℂ) ^ k * operatorXiFiniteLadder N s := by
    simpa [k] using (operatorXiFiniteLadder_centered_involution N s)
  calc
    operatorPhaseNormalizedXiLadder N (operatorCenterInvolutionArg N s)
        = (operatorCenterInvolutionArg N s - m) ^ p *
            operatorXiFiniteLadder N (operatorCenterInvolutionArg N s) := by
              simp [operatorPhaseNormalizedXiLadder, k, p, m]
    _ = ((-1 : ℂ) ^ p * (s - m) ^ p) *
          ((-1 : ℂ) ^ k * operatorXiFiniteLadder N s) := by rw [hpow, hXi]
    _ = (((-1 : ℂ) ^ p * (-1 : ℂ) ^ k) * ((s - m) ^ p * operatorXiFiniteLadder N s)) := by
          ring
    _ = (s - m) ^ p * operatorXiFiniteLadder N s := by
          rw [negOne_pow_mod2_mul_pow_card_cancel k]
          simp
    _ = operatorPhaseNormalizedXiLadder N s := by
          simp [operatorPhaseNormalizedXiLadder, k, p, m]

/-- Recentered phase-normalized ladder around the finite midpoint.
In this coordinate, the involution becomes `z ↦ -z`. -/
def operatorCenteredPhaseNormalizedXiLadder (N : ℕ) (z : ℂ) : ℂ :=
  operatorPhaseNormalizedXiLadder N (z + operatorCenterMidpoint N)

/-- Finite-level evenness of the recentered normalized ladder. -/
theorem operatorCenteredPhaseNormalizedXiLadder_even
    (N : ℕ) (z : ℂ) :
    operatorCenteredPhaseNormalizedXiLadder N (-z) =
      operatorCenteredPhaseNormalizedXiLadder N z := by
  let m : ℂ := operatorCenterMidpoint N
  have hmap : (-z) + m = operatorCenterInvolutionArg N (z + m) := by
    simp [operatorCenterInvolutionArg, operatorCenterMidpoint, m]
    ring
  calc
    operatorCenteredPhaseNormalizedXiLadder N (-z)
        = operatorPhaseNormalizedXiLadder N ((-z) + m) := by
            simp [operatorCenteredPhaseNormalizedXiLadder, m]
    _ = operatorPhaseNormalizedXiLadder N (operatorCenterInvolutionArg N (z + m)) := by
          rw [hmap]
    _ = operatorPhaseNormalizedXiLadder N (z + m) := by
          simpa using (operatorPhaseNormalizedXiLadder_centered_involution N (z + m))
    _ = operatorCenteredPhaseNormalizedXiLadder N z := by
          simp [operatorCenteredPhaseNormalizedXiLadder, m]

/-- Local-uniform limits of the recentered normalized ladder inherit exact
even symmetry. -/
theorem even_limit_of_locallyUniform_operatorCenteredPhaseNormalized
    {F : ℂ → ℂ}
    (hconv : TendstoLocallyUniformly operatorCenteredPhaseNormalizedXiLadder F
      (Filter.atTop : Filter ℕ)) :
    ∀ z : ℂ, F (-z) = F z := by
  intro z
  have hconvOn : TendstoLocallyUniformlyOn operatorCenteredPhaseNormalizedXiLadder F
      (Filter.atTop : Filter ℕ) Set.univ := by
    simpa [tendstoLocallyUniformlyOn_univ] using hconv
  have hz : Filter.Tendsto (fun N : ℕ => operatorCenteredPhaseNormalizedXiLadder N z)
      (Filter.atTop : Filter ℕ) (𝓝 (F z)) :=
    hconvOn.tendsto_at (by simp)
  have hneg : Filter.Tendsto (fun N : ℕ => operatorCenteredPhaseNormalizedXiLadder N (-z))
      (Filter.atTop : Filter ℕ) (𝓝 (F (-z))) :=
    hconvOn.tendsto_at (by simp)
  have hseqEq :
      (fun N : ℕ => operatorCenteredPhaseNormalizedXiLadder N (-z)) =
        (fun N : ℕ => operatorCenteredPhaseNormalizedXiLadder N z) := by
    funext N
    exact operatorCenteredPhaseNormalizedXiLadder_even N z
  have hnegOnZ :
      Filter.Tendsto (fun N : ℕ => operatorCenteredPhaseNormalizedXiLadder N z)
        (Filter.atTop : Filter ℕ) (𝓝 (F (-z))) := by
    simpa [hseqEq] using hneg
  exact (tendsto_nhds_unique hz hnegOnZ).symm

/-- The phase normalization contributes at most a linear factor in the centered
distance. This keeps growth order unchanged. -/
theorem norm_operatorPhaseNormalizedXiLadder_le_linear
    (N : ℕ) (s : ℂ) :
    ‖operatorPhaseNormalizedXiLadder N s‖ ≤
      (1 + ‖s - operatorCenterMidpoint N‖) * ‖operatorXiFiniteLadder N s‖ := by
  let p : ℕ := (operatorSpecN N).card % 2
  let m : ℂ := operatorCenterMidpoint N
  have hpow : ‖s - m‖ ^ p ≤ 1 + ‖s - m‖ := by
    rcases Nat.mod_two_eq_zero_or_one ((operatorSpecN N).card) with hp0 | hp1
    · have hpz : p = 0 := by simpa [p] using hp0
      simp [hpz]
    · have hpo : p = 1 := by simpa [p] using hp1
      simp [hpo]
  calc
    ‖operatorPhaseNormalizedXiLadder N s‖
        = ‖(s - m) ^ p * operatorXiFiniteLadder N s‖ := by
            simp [operatorPhaseNormalizedXiLadder, p, m]
    _ = ‖(s - m) ^ p‖ * ‖operatorXiFiniteLadder N s‖ := norm_mul _ _
    _ = ‖s - m‖ ^ p * ‖operatorXiFiniteLadder N s‖ := by
          simp [norm_pow]
    _ ≤ (1 + ‖s - m‖) * ‖operatorXiFiniteLadder N s‖ := by
          exact mul_le_mul_of_nonneg_right hpow (norm_nonneg _)
    _ = (1 + ‖s - operatorCenterMidpoint N‖) * ‖operatorXiFiniteLadder N s‖ := by
          simp [m]

/-- Recentered normalized ladder bound: at level `N`, only a linear factor in
`‖z‖` is added on top of the base finite-product magnitude. -/
theorem norm_operatorCenteredPhaseNormalizedXiLadder_le_linear
    (N : ℕ) (z : ℂ) :
    ‖operatorCenteredPhaseNormalizedXiLadder N z‖ ≤
      (1 + ‖z‖) *
        ‖operatorXiFiniteLadder N (z + operatorCenterMidpoint N)‖ := by
  simpa [operatorCenteredPhaseNormalizedXiLadder, sub_eq_add_neg, add_assoc, add_left_comm,
    add_comm] using
    (norm_operatorPhaseNormalizedXiLadder_le_linear N (z + operatorCenterMidpoint N))

/-- Telescoping norm control for complex sequences over a finite step range. -/
theorem norm_sub_le_sum_steps
    (f : ℕ → ℂ) (n m : ℕ) :
    ‖f (n + m) - f n‖ ≤
      Finset.sum (Finset.range m) (fun k => ‖f (n + (k + 1)) - f (n + k)‖) := by
  induction m with
  | zero =>
      simp
  | succ m ih =>
      have hsplit :
          f (n + (m + 1)) - f n =
            (f (n + (m + 1)) - f (n + m)) + (f (n + m) - f n) := by
        ring
      calc
        ‖f (n + (m + 1)) - f n‖
            = ‖(f (n + (m + 1)) - f (n + m)) + (f (n + m) - f n)‖ := by
                rw [hsplit]
        _ ≤ ‖f (n + (m + 1)) - f (n + m)‖ + ‖f (n + m) - f n‖ := by
              exact norm_add_le _ _
        _ ≤ ‖f (n + (m + 1)) - f (n + m)‖ +
              Finset.sum (Finset.range m) (fun k => ‖f (n + (k + 1)) - f (n + k)‖) := by
              simpa [add_comm, add_left_comm, add_assoc] using
                (add_le_add_right ih ‖f (n + (m + 1)) - f (n + m)‖)
        _ = Finset.sum (Finset.range (m + 1))
              (fun k => ‖f (n + (k + 1)) - f (n + k)‖) := by
              simp [Finset.sum_range_succ, add_comm]

/-- Telescoping norm control from explicit per-step bounds. -/
theorem norm_sub_le_sum_stepBounds
    (f : ℕ → ℂ) (a : ℕ → ℝ)
    (hstep : ∀ j : ℕ, ‖f (j + 1) - f j‖ ≤ a j)
    (n m : ℕ) :
    ‖f (n + m) - f n‖ ≤
      Finset.sum (Finset.range m) (fun k => a (n + k)) := by
  have htele := norm_sub_le_sum_steps f n m
  have hsum :
      (Finset.sum (Finset.range m) (fun k => ‖f (n + (k + 1)) - f (n + k)‖))
        ≤ Finset.sum (Finset.range m) (fun k => a (n + k)) := by
    exact Finset.sum_le_sum (fun k hk => by
      simpa [Nat.add_assoc] using hstep (n + k))
  exact le_trans htele hsum

/-- Telescoping control specialized to the recentered phase-normalized ladder. -/
theorem norm_operatorCenteredPhaseNormalizedXiLadder_sub_le_sum_stepBounds
    (z : ℂ) (a : ℕ → ℝ)
    (hstep : ∀ j : ℕ,
      ‖operatorCenteredPhaseNormalizedXiLadder (j + 1) z
        - operatorCenteredPhaseNormalizedXiLadder j z‖ ≤ a j)
    (n m : ℕ) :
    ‖operatorCenteredPhaseNormalizedXiLadder (n + m) z
      - operatorCenteredPhaseNormalizedXiLadder n z‖ ≤
      Finset.sum (Finset.range m) (fun k => a (n + k)) := by
  simpa using
    (norm_sub_le_sum_stepBounds
      (f := fun j => operatorCenteredPhaseNormalizedXiLadder j z)
      (a := a) hstep n m)

/-- Uniform Cauchy-on-closed-ball criterion for the recentered phase-normalized
ladder from global step-tail bounds. -/
theorem uniform_cauchyOnClosedBall_of_stepTail_operatorCenteredPhaseNormalized
    (R : ℝ) (a : ℕ → ℝ)
    (hstep : ∀ j : ℕ, ∀ z : ℂ, z ∈ Metric.closedBall (0 : ℂ) R →
      ‖operatorCenteredPhaseNormalizedXiLadder (j + 1) z
        - operatorCenteredPhaseNormalizedXiLadder j z‖ ≤ a j)
    (htail : ∀ ε : ℝ, 0 < ε → ∃ n0 : ℕ,
      ∀ n : ℕ, n0 ≤ n → ∀ m : ℕ,
        Finset.sum (Finset.range m) (fun k => a (n + k)) < ε) :
    ∀ ε : ℝ, 0 < ε → ∃ n0 : ℕ,
      ∀ n : ℕ, n0 ≤ n → ∀ m : ℕ, ∀ z : ℂ,
        z ∈ Metric.closedBall (0 : ℂ) R →
          ‖operatorCenteredPhaseNormalizedXiLadder (n + m) z
            - operatorCenteredPhaseNormalizedXiLadder n z‖ < ε := by
  intro ε hε
  rcases htail ε hε with ⟨n0, hn0⟩
  refine ⟨n0, ?_⟩
  intro n hn m z hz
  have hsum_lt : Finset.sum (Finset.range m) (fun k => a (n + k)) < ε :=
    hn0 n hn m
  have hbound :
      ‖operatorCenteredPhaseNormalizedXiLadder (n + m) z
        - operatorCenteredPhaseNormalizedXiLadder n z‖ ≤
      Finset.sum (Finset.range m) (fun k => a (n + k)) :=
    norm_operatorCenteredPhaseNormalizedXiLadder_sub_le_sum_stepBounds z a
      (fun j => hstep j z hz) n m
  exact lt_of_le_of_lt hbound hsum_lt

/-- Summable nonnegative step profiles give vanishing shifted finite tails. -/
theorem stepTail_of_summable_nonneg
    (a : ℕ → ℝ)
    (ha_nonneg : ∀ j : ℕ, 0 ≤ a j)
    (ha_sum : Summable a) :
    ∀ ε : ℝ, 0 < ε → ∃ n0 : ℕ,
      ∀ n : ℕ, n0 ≤ n → ∀ m : ℕ,
        Finset.sum (Finset.range m) (fun k => a (n + k)) < ε := by
  intro ε hε
  have htend : Filter.Tendsto (fun n : ℕ => ∑' k : ℕ, a (n + k))
      (Filter.atTop : Filter ℕ) (𝓝 (0 : ℝ)) := by
    simpa [Nat.add_comm] using (_root_.tendsto_sum_nat_add a)
  have hEventually :
      ∀ᶠ n : ℕ in (Filter.atTop : Filter ℕ),
        dist (∑' k : ℕ, a (n + k)) 0 < ε := by
    simpa [Metric.ball, dist_eq_norm] using
      (htend (Metric.ball_mem_nhds (0 : ℝ) hε))
  rcases hEventually.exists_forall_of_atTop with ⟨n0, hn0⟩
  refine ⟨n0, ?_⟩
  intro n hn m
  have htsum_lt_norm : ‖∑' k : ℕ, a (n + k)‖ < ε := by
    simpa [dist_eq_norm] using (hn0 n hn)
  have hsum_nonneg : 0 ≤ Finset.sum (Finset.range m) (fun k => a (n + k)) :=
    Finset.sum_nonneg (fun k _hk => ha_nonneg (n + k))
  have hsumNat : Summable (fun k : ℕ => a (n + k)) := by
    simpa [Nat.add_comm] using ((_root_.summable_nat_add_iff n).2 ha_sum)
  have hsum_le_tsum :
      Finset.sum (Finset.range m) (fun k => a (n + k))
        ≤ ∑' k : ℕ, a (n + k) := by
    exact hsumNat.sum_le_tsum (Finset.range m) (fun k _hk => ha_nonneg (n + k))
  have htsum_nonneg : 0 ≤ ∑' k : ℕ, a (n + k) := by
    exact tsum_nonneg (fun k => ha_nonneg (n + k))
  have htsum_lt : (∑' k : ℕ, a (n + k)) < ε := by
    simpa [Real.norm_of_nonneg htsum_nonneg] using htsum_lt_norm
  exact lt_of_le_of_lt hsum_le_tsum htsum_lt

/-- Closed-ball uniform convergence for the recentered normalized ladder from
global step-tail control plus pointwise convergence on that closed ball. -/
theorem tendstoUniformlyOn_closedBall_operatorCenteredPhaseNormalized_of_stepTail_and_pointwise
    (R : ℝ) (a : ℕ → ℝ)
    (hstep : ∀ j : ℕ, ∀ z : ℂ, z ∈ Metric.closedBall (0 : ℂ) R →
      ‖operatorCenteredPhaseNormalizedXiLadder (j + 1) z
        - operatorCenteredPhaseNormalizedXiLadder j z‖ ≤ a j)
    (htail : ∀ ε : ℝ, 0 < ε → ∃ n0 : ℕ,
      ∀ n : ℕ, n0 ≤ n → ∀ m : ℕ,
        Finset.sum (Finset.range m) (fun k => a (n + k)) < ε)
    {F : ℂ → ℂ}
    (hpt : ∀ z : ℂ, z ∈ Metric.closedBall (0 : ℂ) R →
      Filter.Tendsto (fun N : ℕ => operatorCenteredPhaseNormalizedXiLadder N z)
        (Filter.atTop : Filter ℕ) (𝓝 (F z))) :
    TendstoUniformlyOn operatorCenteredPhaseNormalizedXiLadder F
      (Filter.atTop : Filter ℕ) (Metric.closedBall (0 : ℂ) R) := by
  have htailCauchy :=
    uniform_cauchyOnClosedBall_of_stepTail_operatorCenteredPhaseNormalized
      R a hstep htail
  have hUC :
      UniformCauchySeqOn operatorCenteredPhaseNormalizedXiLadder
        (Filter.atTop : Filter ℕ) (Metric.closedBall (0 : ℂ) R) := by
    intro u hu
    rcases (Metric.mem_uniformity_dist).1 hu with ⟨ε, hε, hεsub⟩
    rcases htailCauchy ε hε with ⟨n0, hn0⟩
    rw [Filter.eventually_prod_iff]
    refine ⟨{i : ℕ | n0 ≤ i}, Filter.mem_atTop_sets.mpr ⟨n0, by intro i hi; exact hi⟩,
      {j : ℕ | n0 ≤ j}, Filter.mem_atTop_sets.mpr ⟨n0, by intro j hj; exact hj⟩, ?_⟩
    intro i hi j hj z hz
    have hdist : dist (operatorCenteredPhaseNormalizedXiLadder i z)
        (operatorCenteredPhaseNormalizedXiLadder j z) < ε := by
      by_cases hij : i ≤ j
      · let m : ℕ := j - i
        have hjEq : j = i + m := by
          dsimp [m]
          exact (Nat.add_sub_of_le hij).symm
        have hlt : ‖operatorCenteredPhaseNormalizedXiLadder j z
            - operatorCenteredPhaseNormalizedXiLadder i z‖ < ε := by
          simpa [hjEq] using (hn0 i hi m z hz)
        simpa [dist_eq_norm, norm_sub_rev] using hlt
      · have hji : j ≤ i := le_of_not_ge hij
        let m : ℕ := i - j
        have hiEq : i = j + m := by
          dsimp [m]
          exact (Nat.add_sub_of_le hji).symm
        have hlt : ‖operatorCenteredPhaseNormalizedXiLadder i z
            - operatorCenteredPhaseNormalizedXiLadder j z‖ < ε := by
          simpa [hiEq] using (hn0 j hj m z hz)
        simpa [dist_eq_norm] using hlt
    exact hεsub (by simpa using hdist)
  exact hUC.tendstoUniformlyOn_of_tendsto (fun z hz => hpt z hz)

/-- Local-uniform convergence of the recentered normalized ladder from
closed-ball step-tail control and pointwise convergence. -/
theorem tendstoLocallyUniformly_operatorCenteredPhaseNormalized_of_stepTail_and_pointwise
    (a : ℕ → ℝ)
    (hstep : ∀ R : ℝ, ∀ j : ℕ, ∀ z : ℂ, z ∈ Metric.closedBall (0 : ℂ) R →
      ‖operatorCenteredPhaseNormalizedXiLadder (j + 1) z
        - operatorCenteredPhaseNormalizedXiLadder j z‖ ≤ a j)
    (htail : ∀ R : ℝ, ∀ ε : ℝ, 0 < ε → ∃ n0 : ℕ,
      ∀ n : ℕ, n0 ≤ n → ∀ m : ℕ,
        Finset.sum (Finset.range m) (fun k => a (n + k)) < ε)
    {F : ℂ → ℂ}
    (hpt : ∀ z : ℂ,
      Filter.Tendsto (fun N : ℕ => operatorCenteredPhaseNormalizedXiLadder N z)
        (Filter.atTop : Filter ℕ) (𝓝 (F z))) :
    TendstoLocallyUniformly operatorCenteredPhaseNormalizedXiLadder F
      (Filter.atTop : Filter ℕ) := by
  rw [Metric.tendstoLocallyUniformly_iff]
  intro ε hε x
  let R : ℝ := ‖x‖ + 1
  have hUnif :
      TendstoUniformlyOn operatorCenteredPhaseNormalizedXiLadder F
        (Filter.atTop : Filter ℕ) (Metric.closedBall (0 : ℂ) R) :=
    tendstoUniformlyOn_closedBall_operatorCenteredPhaseNormalized_of_stepTail_and_pointwise
      R a
      (hstep R)
      (htail R)
      (fun z _hz => hpt z)
  have hBall :
      Metric.ball x 1 ⊆ Metric.closedBall (0 : ℂ) R := by
    intro y hy
    have hxy : ‖y - x‖ < 1 := by
      simpa [Metric.mem_ball, dist_eq_norm] using hy
    have hyNorm : ‖y‖ ≤ R := by
      have hyLt : ‖y‖ < R := by
        calc
          ‖y‖ = ‖(y - x) + x‖ := by ring
          _ ≤ ‖y - x‖ + ‖x‖ := norm_add_le _ _
          _ < 1 + ‖x‖ := by linarith
          _ = R := by simp [R, add_comm]
      exact le_of_lt hyLt
    simpa [Metric.mem_closedBall, dist_eq_norm, R] using hyNorm
  refine ⟨Metric.ball x 1, Metric.ball_mem_nhds x zero_lt_one, ?_⟩
  have hUnifε :
      ∀ᶠ n : ℕ in (Filter.atTop : Filter ℕ),
        ∀ y : ℂ, y ∈ Metric.closedBall (0 : ℂ) R →
          dist (F y) (operatorCenteredPhaseNormalizedXiLadder n y) < ε :=
    (Metric.tendstoUniformlyOn_iff.1 hUnif) ε hε
  exact hUnifε.mono (fun n hn y hy => hn y (hBall hy))

/-- Local-uniform convergence of the recentered normalized ladder from
pointwise convergence plus a summable nonnegative step profile. -/
theorem tendstoLocallyUniformly_operatorCenteredPhaseNormalized_of_stepSummable_and_pointwise
    (a : ℕ → ℝ)
    (ha_nonneg : ∀ j : ℕ, 0 ≤ a j)
    (ha_sum : Summable a)
    (hstep : ∀ R : ℝ, ∀ j : ℕ, ∀ z : ℂ, z ∈ Metric.closedBall (0 : ℂ) R →
      ‖operatorCenteredPhaseNormalizedXiLadder (j + 1) z
        - operatorCenteredPhaseNormalizedXiLadder j z‖ ≤ a j)
    {F : ℂ → ℂ}
    (hpt : ∀ z : ℂ,
      Filter.Tendsto (fun N : ℕ => operatorCenteredPhaseNormalizedXiLadder N z)
        (Filter.atTop : Filter ℕ) (𝓝 (F z))) :
    TendstoLocallyUniformly operatorCenteredPhaseNormalizedXiLadder F
      (Filter.atTop : Filter ℕ) := by
  exact
    tendstoLocallyUniformly_operatorCenteredPhaseNormalized_of_stepTail_and_pointwise
      a hstep
      (fun _R ε hε => stepTail_of_summable_nonneg a ha_nonneg ha_sum ε hε)
      hpt

/-- Closed-ball uniform convergence for the concrete operator Xi ladder from
global step-tail control plus pointwise convergence on that closed ball. -/
theorem tendstoUniformlyOn_closedBall_operatorXiFiniteLadder_of_stepTail_and_pointwise
    (R : ℝ) (a : ℕ → ℝ)
    (hstep : ∀ j : ℕ, ∀ z : ℂ, z ∈ Metric.closedBall (0 : ℂ) R →
      ‖operatorXiFiniteLadder (j + 1) z - operatorXiFiniteLadder j z‖ ≤ a j)
    (htail : ∀ ε : ℝ, 0 < ε → ∃ n0 : ℕ,
      ∀ n : ℕ, n0 ≤ n → ∀ m : ℕ,
        Finset.sum (Finset.range m) (fun k => a (n + k)) < ε)
    {F : ℂ → ℂ}
    (hpt : ∀ z : ℂ, z ∈ Metric.closedBall (0 : ℂ) R →
      Filter.Tendsto (fun N : ℕ => operatorXiFiniteLadder N z)
        (Filter.atTop : Filter ℕ) (𝓝 (F z))) :
    TendstoUniformlyOn operatorXiFiniteLadder F
      (Filter.atTop : Filter ℕ) (Metric.closedBall (0 : ℂ) R) := by
  have hUC :
      UniformCauchySeqOn operatorXiFiniteLadder
        (Filter.atTop : Filter ℕ) (Metric.closedBall (0 : ℂ) R) := by
    intro u hu
    rcases (Metric.mem_uniformity_dist).1 hu with ⟨ε, hε, hεsub⟩
    rcases htail ε hε with ⟨n0, hn0⟩
    rw [Filter.eventually_prod_iff]
    refine ⟨{i : ℕ | n0 ≤ i}, Filter.mem_atTop_sets.mpr ⟨n0, by intro i hi; exact hi⟩,
      {j : ℕ | n0 ≤ j}, Filter.mem_atTop_sets.mpr ⟨n0, by intro j hj; exact hj⟩, ?_⟩
    intro i hi j hj z hz
    have hdist : dist (operatorXiFiniteLadder i z) (operatorXiFiniteLadder j z) < ε := by
      by_cases hij : i ≤ j
      · let m : ℕ := j - i
        have hjEq : j = i + m := by
          dsimp [m]
          exact (Nat.add_sub_of_le hij).symm
        have hbound :
            ‖operatorXiFiniteLadder j z - operatorXiFiniteLadder i z‖ ≤
              Finset.sum (Finset.range m) (fun k => a (i + k)) := by
          have hle := norm_sub_le_sum_stepBounds
            (fun j => operatorXiFiniteLadder j z) a (fun j => hstep j z hz) i m
          exact by simpa [hjEq] using hle
        have hlt : ‖operatorXiFiniteLadder j z - operatorXiFiniteLadder i z‖ < ε :=
          lt_of_le_of_lt hbound (hn0 i hi m)
        simpa [dist_eq_norm, norm_sub_rev] using hlt
      · have hji : j ≤ i := le_of_not_ge hij
        let m : ℕ := i - j
        have hiEq : i = j + m := by
          dsimp [m]
          exact (Nat.add_sub_of_le hji).symm
        have hbound :
            ‖operatorXiFiniteLadder i z - operatorXiFiniteLadder j z‖ ≤
              Finset.sum (Finset.range m) (fun k => a (j + k)) := by
          have hle := norm_sub_le_sum_stepBounds
            (fun j => operatorXiFiniteLadder j z) a (fun j => hstep j z hz) j m
          exact by simpa [hiEq] using hle
        have hlt : ‖operatorXiFiniteLadder i z - operatorXiFiniteLadder j z‖ < ε :=
          lt_of_le_of_lt hbound (hn0 j hj m)
        simpa [dist_eq_norm] using hlt
    exact hεsub (by simpa using hdist)
  exact hUC.tendstoUniformlyOn_of_tendsto (fun z hz => hpt z hz)

/-- Local-uniform convergence of the concrete operator Xi ladder from
closed-ball step-tail control and pointwise convergence. -/
theorem tendstoLocallyUniformly_operatorXiFiniteLadder_of_stepTail_and_pointwise
    (a : ℕ → ℝ)
    (hstep : ∀ R : ℝ, ∀ j : ℕ, ∀ z : ℂ, z ∈ Metric.closedBall (0 : ℂ) R →
      ‖operatorXiFiniteLadder (j + 1) z - operatorXiFiniteLadder j z‖ ≤ a j)
    (htail : ∀ R : ℝ, ∀ ε : ℝ, 0 < ε → ∃ n0 : ℕ,
      ∀ n : ℕ, n0 ≤ n → ∀ m : ℕ,
        Finset.sum (Finset.range m) (fun k => a (n + k)) < ε)
    {F : ℂ → ℂ}
    (hpt : ∀ z : ℂ,
      Filter.Tendsto (fun N : ℕ => operatorXiFiniteLadder N z)
        (Filter.atTop : Filter ℕ) (𝓝 (F z))) :
    TendstoLocallyUniformly operatorXiFiniteLadder F
      (Filter.atTop : Filter ℕ) := by
  rw [Metric.tendstoLocallyUniformly_iff]
  intro ε hε x
  let R : ℝ := ‖x‖ + 1
  have hUnif :
      TendstoUniformlyOn operatorXiFiniteLadder F
        (Filter.atTop : Filter ℕ) (Metric.closedBall (0 : ℂ) R) :=
    tendstoUniformlyOn_closedBall_operatorXiFiniteLadder_of_stepTail_and_pointwise
      R a (hstep R) (htail R) (fun z _hz => hpt z)
  have hBall : Metric.ball x 1 ⊆ Metric.closedBall (0 : ℂ) R := by
    intro y hy
    have hxy : ‖y - x‖ < 1 := by
      simpa [Metric.mem_ball, dist_eq_norm] using hy
    have hyNorm : ‖y‖ ≤ R := by
      have hyLt : ‖y‖ < R := by
        calc
          ‖y‖ = ‖(y - x) + x‖ := by ring
          _ ≤ ‖y - x‖ + ‖x‖ := norm_add_le _ _
          _ < 1 + ‖x‖ := by linarith
          _ = R := by simp [R, add_comm]
      exact le_of_lt hyLt
    simpa [Metric.mem_closedBall, dist_eq_norm, R] using hyNorm
  refine ⟨Metric.ball x 1, Metric.ball_mem_nhds x zero_lt_one, ?_⟩
  have hUnifε :
      ∀ᶠ n : ℕ in (Filter.atTop : Filter ℕ),
        ∀ y : ℂ, y ∈ Metric.closedBall (0 : ℂ) R →
          dist (F y) (operatorXiFiniteLadder n y) < ε :=
    (Metric.tendstoUniformlyOn_iff.1 hUnif) ε hε
  exact hUnifε.mono (fun n hn y hy => hn y (hBall hy))

/-- Local-uniform convergence of the concrete operator Xi ladder from
pointwise convergence plus a summable nonnegative step profile. -/
theorem tendstoLocallyUniformly_operatorXiFiniteLadder_of_stepSummable_and_pointwise
    (a : ℕ → ℝ)
    (ha_nonneg : ∀ j : ℕ, 0 ≤ a j)
    (ha_sum : Summable a)
    (hstep : ∀ R : ℝ, ∀ j : ℕ, ∀ z : ℂ, z ∈ Metric.closedBall (0 : ℂ) R →
      ‖operatorXiFiniteLadder (j + 1) z - operatorXiFiniteLadder j z‖ ≤ a j)
    {F : ℂ → ℂ}
    (hpt : ∀ z : ℂ,
      Filter.Tendsto (fun N : ℕ => operatorXiFiniteLadder N z)
        (Filter.atTop : Filter ℕ) (𝓝 (F z))) :
    TendstoLocallyUniformly operatorXiFiniteLadder F
      (Filter.atTop : Filter ℕ) := by
  exact
    tendstoLocallyUniformly_operatorXiFiniteLadder_of_stepTail_and_pointwise
      a hstep
      (fun _R ε hε => stepTail_of_summable_nonneg a ha_nonneg ha_sum ε hε)
      hpt

/-- Critical-line points are never zero. -/
theorem criticalLinePoint_ne_zero (t : ℝ) :
    criticalLinePoint t ≠ 0 := by
  intro h
  have hre : (criticalLinePoint t).re = 0 := by
    simpa [h] using congrArg Complex.re h
  have hhalf : (1 / 2 : ℝ) = 0 := by
    simpa [criticalLinePoint_re] using hre
  norm_num at hhalf

/-- The finite model `XiFinite spec` is complex-differentiable everywhere. -/
theorem differentiable_XiFinite
    (spec : Finset ℝ) :
    Differentiable ℂ (XiFinite spec) := by
  intro z
  unfold XiFinite
  let f : ℝ → ℂ → ℂ := fun t s => s - criticalLinePoint t
  have hf : ∀ t ∈ spec, DifferentiableAt ℂ (f t) z := by
    intro t ht
    simpa [f] using (differentiableAt_id.sub_const (criticalLinePoint t))
  simpa [f] using (DifferentiableAt.fun_finset_prod (u := spec) (f := f) hf)

/-- Each finite level of the concrete operator Xi ladder is differentiable. -/
theorem differentiable_operatorXiFiniteLadder (N : ℕ) :
    Differentiable ℂ (operatorXiFiniteLadder N) := by
  simpa [operatorXiFiniteLadder, XiFiniteLadder] using
    (differentiable_XiFinite (operatorSpecN N))

/-- The finite operator Xi products are nontrivial: they never vanish at `0`. -/
theorem xiFinite_operatorSpecN_zero_ne (N : ℕ) :
    XiFinite (operatorSpecN N) 0 ≠ 0 := by
  unfold XiFinite
  refine Finset.prod_ne_zero_iff.mpr ?_
  intro t ht
  simp [criticalLinePoint_ne_zero t]

/-- Concrete operator-level envelope bound for finite Xi factors. -/
theorem norm_operatorXiFiniteLadder_le_factorized_envelope
    (N : ℕ) (s : ℂ) :
    ‖operatorXiFiniteLadder N s‖ ≤
      Finset.prod (operatorSpecN N) (fun t => (‖s‖ + ‖criticalLinePoint t‖)) := by
  simpa [operatorXiFiniteLadder, XiFiniteLadder] using
    (norm_XiFinite_le_factorized_envelope (operatorSpecN N) s)

/-- Concrete operator-level cardinality growth bound under a uniform
ordinate envelope at level `N`. -/
theorem norm_operatorXiFiniteLadder_le_pow_card_of_ordinate_bound
    (N : ℕ) (s : ℂ) (B : ℝ)
    (hB : ∀ t ∈ operatorSpecN N, ‖criticalLinePoint t‖ ≤ B) :
    ‖operatorXiFiniteLadder N s‖ ≤ (‖s‖ + B) ^ (operatorSpecN N).card := by
  simpa [operatorXiFiniteLadder, XiFiniteLadder] using
    (norm_XiFinite_le_pow_card_of_ordinate_bound
      (operatorSpecN N) s B hB)

/-- One-step increment bound specialized to the concrete operator finite product
at level `N` when adding a fresh ordinate. -/
theorem norm_operatorXiFinite_insert_sub_le
    (N : ℕ) {t : ℝ} (ht : t ∉ operatorSpecN N) (s : ℂ) :
    ‖XiFinite (insert t (operatorSpecN N)) s - operatorXiFiniteLadder N s‖ ≤
      ‖s - criticalLinePoint t - 1‖ *
        Finset.prod (operatorSpecN N) (fun u => (‖s‖ + ‖criticalLinePoint u‖)) := by
  simpa [operatorXiFiniteLadder, XiFiniteLadder] using
    (norm_XiFinite_insert_sub_le (operatorSpecN N) ht s)

/-- Concrete Hadamard step bound from an insert-step model plus linear lower growth
of inserted ordinates and a uniform previous-level norm bound on the closed ball. -/
theorem norm_operatorXiFiniteHadamard_step_le_profile_of_insert_linear_growth
    (R M c : ℝ)
    (hR0 : 0 ≤ R) (hM0 : 0 ≤ M) (hc0 : 0 < c) (hRc : R ≤ c)
    (hinsert : ∀ n : ℕ, ∃ t : ℝ,
      t ∉ operatorSpecN n ∧
      operatorSpecN (n + 1) = insert t (operatorSpecN n) ∧
      c * (n + 1 : ℝ) ≤ ‖criticalLinePoint t‖)
    (hprod : ∀ n : ℕ, ∀ z : ℂ, z ∈ Metric.closedBall (0 : ℂ) R →
      ‖operatorXiFiniteHadamardLadder n z‖ ≤ M) :
    ∀ n : ℕ, ∀ z : ℂ, z ∈ Metric.closedBall (0 : ℂ) R →
      ‖operatorXiFiniteHadamardLadder (n + 1) z - operatorXiFiniteHadamardLadder n z‖
        ≤ hadamardQuadraticProfile R M c n := by
  intro n z hz
  rcases hinsert n with ⟨t, htNot, hstepEq, hgrow⟩
  have hsR : ‖z‖ ≤ R := by
    simpa [Metric.mem_closedBall, dist_eq_norm] using hz
  have hn1pos : (0 : ℝ) < (n + 1 : ℝ) := by
    exact_mod_cast Nat.succ_pos n
  have hcmul_pos : 0 < c * (n + 1 : ℝ) := mul_pos hc0 hn1pos
  have hden_pos : 0 < ‖criticalLinePoint t‖ := lt_of_lt_of_le hcmul_pos hgrow
  have hc_le_cmul : c ≤ c * (n + 1 : ℝ) := by
    have h1le : (1 : ℝ) ≤ (n + 1 : ℝ) := by linarith [hn1pos]
    nlinarith [hc0.le, h1le]
  have hR_le_den : R ≤ ‖criticalLinePoint t‖ := by
    exact le_trans hRc (le_trans hc_le_cmul hgrow)
  have hratio_le_one : ‖z / criticalLinePoint t‖ ≤ 1 := by
    rw [norm_div]
    exact div_le_one_of_le₀ (by linarith [hsR, hR_le_den]) (norm_nonneg _)
  have hratio_le_linear : ‖z / criticalLinePoint t‖ ≤ R / (c * (n + 1 : ℝ)) := by
    rw [norm_div]
    have hz_over_den : ‖z‖ / ‖criticalLinePoint t‖ ≤ R / ‖criticalLinePoint t‖ :=
      div_le_div_of_nonneg_right hsR (norm_nonneg _)
    have hR_div :
        R / ‖criticalLinePoint t‖ ≤ R / (c * (n + 1 : ℝ)) := by
      have hInv :
          (1 / ‖criticalLinePoint t‖) ≤ (1 / (c * (n + 1 : ℝ))) := by
        exact one_div_le_one_div_of_le hcmul_pos hgrow
      simpa [div_eq_mul_inv, mul_assoc, mul_left_comm, mul_comm] using
        (mul_le_mul_of_nonneg_left hInv hR0)
    exact le_trans hz_over_den hR_div
  have hratio_sq :
      ‖z / criticalLinePoint t‖ ^ (2 : ℕ) ≤
        (R / (c * (n + 1 : ℝ))) ^ (2 : ℕ) := by
    nlinarith [norm_nonneg (z / criticalLinePoint t), hratio_le_linear]
  have hstepInsert :
      ‖XiFiniteHadamard (insert t (operatorSpecN n)) z - operatorXiFiniteHadamardLadder n z‖
        ≤ (3 * ‖z / criticalLinePoint t‖ ^ (2 : ℕ)) * ‖operatorXiFiniteHadamardLadder n z‖ :=
    norm_operatorXiFiniteHadamard_insert_sub_le_three_mul_sq n htNot z hratio_le_one
  have hstepLadder :
      ‖operatorXiFiniteHadamardLadder (n + 1) z - operatorXiFiniteHadamardLadder n z‖
        ≤ (3 * ‖z / criticalLinePoint t‖ ^ (2 : ℕ)) * ‖operatorXiFiniteHadamardLadder n z‖ := by
    simpa [operatorXiFiniteHadamardLadder, hstepEq] using hstepInsert
  have hprodN : ‖operatorXiFiniteHadamardLadder n z‖ ≤ M := hprod n z hz
  have hfac_nonneg : 0 ≤ 3 * ‖z / criticalLinePoint t‖ ^ (2 : ℕ) := by positivity
  have hfac_le :
      3 * ‖z / criticalLinePoint t‖ ^ (2 : ℕ) ≤
        3 * (R / (c * (n + 1 : ℝ))) ^ (2 : ℕ) := by
    gcongr
  have hfacR_nonneg : 0 ≤ 3 * (R / (c * (n + 1 : ℝ))) ^ (2 : ℕ) := by positivity
  have hmul_le :
      (3 * ‖z / criticalLinePoint t‖ ^ (2 : ℕ)) * ‖operatorXiFiniteHadamardLadder n z‖ ≤
        (3 * (R / (c * (n + 1 : ℝ))) ^ (2 : ℕ)) * M := by
    exact mul_le_mul hfac_le hprodN (norm_nonneg _) hfacR_nonneg
  exact le_trans hstepLadder (by simpa [hadamardQuadraticProfile] using hmul_le)

/-- The concrete operator spectral ladder at level `N` is nonempty. -/
theorem operatorSpecN_nonempty (N : ℕ) :
    (operatorSpecN N).Nonempty := by
  classical
  have huniv : (Finset.univ : Finset (Fin (N + 1))).Nonempty := Finset.univ_nonempty
  rcases huniv with ⟨i, hi⟩
  refine ⟨operatorEigenvalues N i, ?_⟩
  exact Finset.mem_image.mpr ⟨i, hi, rfl⟩

/-- Existence of a concrete operator ordinate at level `N` that lies on or above
the structural center midpoint. -/
theorem exists_mem_operatorSpecN_ge_center_midpoint (N : ℕ) :
    ∃ t : ℝ, t ∈ operatorSpecN N ∧
      ((structuralCenterQ (N + 1) : ℝ) / 2 ≤ t) := by
  rcases operatorSpecN_nonempty N with ⟨t0, ht0⟩
  let c : ℝ := (structuralCenterQ (N + 1) : ℝ)
  let t1 : ℝ := c - t0
  have ht1 : t1 ∈ operatorSpecN N := by
    simpa [t1, c] using
      (mem_operatorSpecN_reflect_structuralCenter (N := N) (t := t0) ht0)
  by_cases h0 : c / 2 ≤ t0
  · exact ⟨t0, ht0, h0⟩
  · have h0lt : t0 < c / 2 := lt_of_not_ge h0
    have h1ge : c / 2 ≤ t1 := by
      dsimp [t1, c]
      linarith
    exact ⟨t1, ht1, h1ge⟩

/-- Existence of a concrete operator ordinate at level `N` with explicit linear
lower growth `(N+1)/2`. -/
theorem exists_mem_operatorSpecN_ge_linear_half (N : ℕ) :
    ∃ t : ℝ, t ∈ operatorSpecN N ∧ (((N + 1 : ℕ) : ℝ) / 2 ≤ t) := by
  rcases exists_mem_operatorSpecN_ge_center_midpoint N with ⟨t, ht, hmid⟩
  refine ⟨t, ht, ?_⟩
  have hmid_ge :
      (((N + 1 : ℕ) : ℝ) / 2) ≤ ((structuralCenterQ (N + 1) : ℝ) / 2) := by
    norm_num [structuralCenterQ, timelikeOffsetQ]
    nlinarith
  exact le_trans hmid_ge hmid

/-- Canonical concrete large-ordinate selector from the operator ladder at each level. -/
noncomputable def operatorLargeOrdinate (N : ℕ) : ℝ :=
  Classical.choose (exists_mem_operatorSpecN_ge_linear_half N)

/-- Membership of the canonical large-ordinate selector in the concrete operator ladder. -/
theorem operatorLargeOrdinate_mem (N : ℕ) :
    operatorLargeOrdinate N ∈ operatorSpecN N := by
  exact (Classical.choose_spec (exists_mem_operatorSpecN_ge_linear_half N)).1

/-- Explicit linear lower bound for the canonical large-ordinate selector. -/
theorem operatorLargeOrdinate_lower (N : ℕ) :
    (((N + 1 : ℕ) : ℝ) / 2) ≤ operatorLargeOrdinate N := by
  exact (Classical.choose_spec (exists_mem_operatorSpecN_ge_linear_half N)).2

/-- The complex critical-line norm at the canonical large ordinate is at least
the same explicit linear half-growth lower bound. -/
theorem operatorLargeOrdinate_criticalLine_norm_lower (N : ℕ) :
    (((N + 1 : ℕ) : ℝ) / 2) ≤ ‖criticalLinePoint (operatorLargeOrdinate N)‖ := by
  have hlin : (((N + 1 : ℕ) : ℝ) / 2) ≤ operatorLargeOrdinate N :=
    operatorLargeOrdinate_lower N
  have hnonneg_t : 0 ≤ operatorLargeOrdinate N := by
    have hhalf_nonneg : 0 ≤ (((N + 1 : ℕ) : ℝ) / 2) := by positivity
    exact le_trans hhalf_nonneg hlin
  have him_to_norm :
      operatorLargeOrdinate N ≤ ‖criticalLinePoint (operatorLargeOrdinate N)‖ := by
    have habs_im : |(criticalLinePoint (operatorLargeOrdinate N)).im|
        ≤ ‖criticalLinePoint (operatorLargeOrdinate N)‖ :=
      Complex.abs_im_le_norm (criticalLinePoint (operatorLargeOrdinate N))
    have habs_im_eq :
        |(criticalLinePoint (operatorLargeOrdinate N)).im| = operatorLargeOrdinate N := by
      simpa [criticalLinePoint_im, abs_of_nonneg hnonneg_t]
    simpa [habs_im_eq] using habs_im
  exact le_trans hlin him_to_norm

/-- The inverse-square series over the canonical large-ordinate selector is summable. -/
theorem summable_one_div_abs_sq_operatorLargeOrdinate :
    Summable (fun N : ℕ => (1 : ℝ) / (|operatorLargeOrdinate N| ^ (2 : ℕ))) := by
  have hbase : Summable (fun n : ℕ => (1 : ℝ) / ((n : ℝ) ^ (2 : ℕ))) :=
    (Real.summable_one_div_nat_pow).2 (by norm_num)
  have hshift : Summable (fun n : ℕ => (1 : ℝ) / (((n + 1 : ℕ) : ℝ) ^ (2 : ℕ))) := by
    simpa [Nat.cast_add, add_assoc, add_comm, add_left_comm] using
      (summable_nat_add_iff 1).2 hbase
  have hmajor : Summable (fun n : ℕ =>
      (4 : ℝ) * ((1 : ℝ) / (((n + 1 : ℕ) : ℝ) ^ (2 : ℕ)))) := by
    exact hshift.mul_left (4 : ℝ)
  have hnonneg :
      ∀ n : ℕ, 0 ≤ (1 : ℝ) / (|operatorLargeOrdinate n| ^ (2 : ℕ)) := by
    intro n
    positivity
  have hle :
      ∀ n : ℕ,
        (1 : ℝ) / (|operatorLargeOrdinate n| ^ (2 : ℕ))
          ≤ (4 : ℝ) * ((1 : ℝ) / (((n + 1 : ℕ) : ℝ) ^ (2 : ℕ))) := by
    intro n
    have hlin : (((n + 1 : ℕ) : ℝ) / 2) ≤ operatorLargeOrdinate n :=
      operatorLargeOrdinate_lower n
    have hnonneg_t : 0 ≤ operatorLargeOrdinate n := by
      have hhalf_nonneg : 0 ≤ (((n + 1 : ℕ) : ℝ) / 2) := by positivity
      exact le_trans hhalf_nonneg hlin
    have habs_lower : (((n + 1 : ℕ) : ℝ) / 2) ≤ |operatorLargeOrdinate n| := by
      simpa [abs_of_nonneg hnonneg_t] using hlin
    have ha_nonneg : 0 ≤ (((n + 1 : ℕ) : ℝ) / 2) := by
      positivity
    have hb_nonneg : 0 ≤ |operatorLargeOrdinate n| := abs_nonneg _
    have hpow_le :
        ((((n + 1 : ℕ) : ℝ) / 2) ^ (2 : ℕ)) ≤ (|operatorLargeOrdinate n| ^ (2 : ℕ)) := by
      nlinarith [habs_lower, ha_nonneg, hb_nonneg]
    have hhalf_sq_pos : 0 < ((((n + 1 : ℕ) : ℝ) / 2) ^ (2 : ℕ)) := by
      positivity
    have hrecip :
        (1 : ℝ) / (|operatorLargeOrdinate n| ^ (2 : ℕ))
          ≤ (1 : ℝ) / ((((n + 1 : ℕ) : ℝ) / 2) ^ (2 : ℕ)) := by
      exact one_div_le_one_div_of_le hhalf_sq_pos hpow_le
    have hn1 : (((n + 1 : ℕ) : ℝ)) ≠ 0 := by
      exact_mod_cast Nat.succ_ne_zero n
    have hhalf_rewrite :
        (1 : ℝ) / ((((n + 1 : ℕ) : ℝ) / 2) ^ (2 : ℕ))
          = (4 : ℝ) * ((1 : ℝ) / (((n + 1 : ℕ) : ℝ) ^ (2 : ℕ))) := by
      field_simp [hn1]
      ring
    calc
      (1 : ℝ) / (|operatorLargeOrdinate n| ^ (2 : ℕ))
          ≤ (1 : ℝ) / ((((n + 1 : ℕ) : ℝ) / 2) ^ (2 : ℕ)) := hrecip
      _ = (4 : ℝ) * ((1 : ℝ) / (((n + 1 : ℕ) : ℝ) ^ (2 : ℕ))) := hhalf_rewrite
  exact Summable.of_nonneg_of_le hnonneg hle hmajor

/-- Uniform inverse-square summability through critical-line norms at the
canonical large-ordinate selector (Path-B ready majorant). -/
theorem summable_one_div_normSq_criticalLine_operatorLargeOrdinate :
    Summable (fun N : ℕ =>
      (1 : ℝ) / (‖criticalLinePoint (operatorLargeOrdinate N)‖ ^ (2 : ℕ))) := by
  have hnonneg :
      ∀ N : ℕ,
        0 ≤ (1 : ℝ) / (‖criticalLinePoint (operatorLargeOrdinate N)‖ ^ (2 : ℕ)) := by
    intro N
    positivity
  have hle :
      ∀ N : ℕ,
        (1 : ℝ) / (‖criticalLinePoint (operatorLargeOrdinate N)‖ ^ (2 : ℕ))
          ≤ (1 : ℝ) / (|operatorLargeOrdinate N| ^ (2 : ℕ)) := by
    intro N
    let t : ℝ := operatorLargeOrdinate N
    have hnorm_ge_abs : |t| ≤ ‖criticalLinePoint t‖ := by
      have habs_im : |(criticalLinePoint t).im| ≤ ‖criticalLinePoint t‖ :=
        Complex.abs_im_le_norm (criticalLinePoint t)
      simpa [criticalLinePoint_im] using habs_im
    have hpow : |t| ^ (2 : ℕ) ≤ ‖criticalLinePoint t‖ ^ (2 : ℕ) := by
      nlinarith [hnorm_ge_abs, abs_nonneg t, norm_nonneg (criticalLinePoint t)]
    have hlin : (((N + 1 : ℕ) : ℝ) / 2) ≤ t := by
      simpa [t] using operatorLargeOrdinate_lower N
    have hhalf_pos : 0 < (((N + 1 : ℕ) : ℝ) / 2) := by positivity
    have ht_pos : 0 < |t| := by
      have ht_nonzero : t ≠ 0 := by
        linarith [hlin, hhalf_pos]
      exact abs_pos.mpr ht_nonzero
    have hpow_pos : 0 < |t| ^ (2 : ℕ) := by positivity
    have hrecip :
        (1 : ℝ) / (‖criticalLinePoint t‖ ^ (2 : ℕ))
          ≤ (1 : ℝ) / (|t| ^ (2 : ℕ)) := by
      exact one_div_le_one_div_of_le hpow_pos hpow
    simpa [t] using hrecip
  exact Summable.of_nonneg_of_le hnonneg hle
    summable_one_div_abs_sq_operatorLargeOrdinate

/-- Path-B finite-level full-spectrum bound under indexed linear growth:
if each indexed operator eigenvalue at level `M+1` has linear lower growth,
then the complete inverse-square sum over `operatorSpecN M` is uniformly bounded in `M`. -/
theorem exists_uniform_bound_sum_one_div_normSq_operatorSpecN_of_indexed_linear_growth
    (hlin : ∀ M : ℕ, ∀ i : Fin (M + 1),
      (((i.1 + 1 : ℕ) : ℝ) / 2) ≤ operatorEigenvalues M i) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ M : ℕ,
      Finset.sum (operatorSpecN M) (fun t => (1 : ℝ) / (‖criticalLinePoint t‖ ^ (2 : ℕ))) ≤ C := by
  let a : ℕ → ℝ := fun n => (4 : ℝ) * ((1 : ℝ) / (((n + 1 : ℕ) : ℝ) ^ (2 : ℕ)))
  have ha_nonneg : ∀ n : ℕ, 0 ≤ a n := by
    intro n
    dsimp [a]
    positivity
  have hsum_a : Summable a := by
    have hbase0 : Summable (fun n : ℕ => (1 : ℝ) / ((n : ℝ) ^ (2 : ℕ))) :=
      (Real.summable_one_div_nat_pow).2 (by norm_num)
    have hshift : Summable (fun n : ℕ => (1 : ℝ) / (((n + 1 : ℕ) : ℝ) ^ (2 : ℕ))) := by
      simpa [Nat.cast_add, add_assoc, add_comm, add_left_comm] using
        (summable_nat_add_iff 1).2 hbase0
    simpa [a] using hshift.mul_left (4 : ℝ)
  refine ⟨∑' n : ℕ, a n, tsum_nonneg ha_nonneg, ?_⟩
  intro M
  classical
  have himage_le :
      Finset.sum (operatorSpecN M) (fun t => (1 : ℝ) / (‖criticalLinePoint t‖ ^ (2 : ℕ)))
        ≤ Finset.sum Finset.univ
            (fun i : Fin (M + 1) =>
              (1 : ℝ) / (‖criticalLinePoint (operatorEigenvalues M i)‖ ^ (2 : ℕ))) := by
    unfold operatorSpecN
    simpa using
      (Finset.sum_image_le_of_nonneg
        (s := (Finset.univ : Finset (Fin (M + 1))))
        (g := fun i : Fin (M + 1) => operatorEigenvalues M i)
        (f := fun t : ℝ => (1 : ℝ) / (‖criticalLinePoint t‖ ^ (2 : ℕ)))
        (by intro u hu; positivity))
  have hindex_le :
      Finset.sum Finset.univ
          (fun i : Fin (M + 1) =>
            (1 : ℝ) / (‖criticalLinePoint (operatorEigenvalues M i)‖ ^ (2 : ℕ)))
        ≤ Finset.sum Finset.univ (fun i : Fin (M + 1) => a i.1) := by
    refine Finset.sum_le_sum (fun i _ => ?_)
    let t : ℝ := operatorEigenvalues M i
    have hlin_i : (((i.1 + 1 : ℕ) : ℝ) / 2) ≤ t := hlin M i
    have hnonneg_half : 0 ≤ (((i.1 + 1 : ℕ) : ℝ) / 2) := by positivity
    have hnonneg_t : 0 ≤ t := le_trans hnonneg_half hlin_i
    have hnorm_ge_t : t ≤ ‖criticalLinePoint t‖ := by
      have him_abs : |(criticalLinePoint t).im| ≤ ‖criticalLinePoint t‖ :=
        Complex.abs_im_le_norm (criticalLinePoint t)
      have him_eq : |(criticalLinePoint t).im| = t := by
        simpa [criticalLinePoint_im, abs_of_nonneg hnonneg_t]
      simpa [him_eq] using him_abs
    have hnorm_lower : (((i.1 + 1 : ℕ) : ℝ) / 2) ≤ ‖criticalLinePoint t‖ :=
      le_trans hlin_i hnorm_ge_t
    have hpow_lower :
        ((((i.1 + 1 : ℕ) : ℝ) / 2) ^ (2 : ℕ))
          ≤ ‖criticalLinePoint t‖ ^ (2 : ℕ) := by
      exact pow_le_pow_left₀ hnonneg_half hnorm_lower (2 : ℕ)
    have hpow_half_pos : 0 < ((((i.1 + 1 : ℕ) : ℝ) / 2) ^ (2 : ℕ)) := by
      positivity
    have hrecip :
        (1 : ℝ) / (‖criticalLinePoint t‖ ^ (2 : ℕ))
          ≤ (1 : ℝ) / ((((i.1 + 1 : ℕ) : ℝ) / 2) ^ (2 : ℕ)) := by
      exact one_div_le_one_div_of_le hpow_half_pos hpow_lower
    have hi1_ne : (((i.1 + 1 : ℕ) : ℝ)) ≠ 0 := by
      exact_mod_cast Nat.succ_ne_zero i.1
    have hrewrite :
        (1 : ℝ) / ((((i.1 + 1 : ℕ) : ℝ) / 2) ^ (2 : ℕ)) = a i.1 := by
      dsimp [a]
      field_simp [hi1_ne]
      ring
    have hmain :
        (1 : ℝ) / (‖criticalLinePoint t‖ ^ (2 : ℕ)) ≤ a i.1 := by
      rw [← hrewrite]
      exact hrecip
    simpa [t] using hmain
  have hsum_range :
      Finset.sum Finset.univ (fun i : Fin (M + 1) => a i.1)
        = Finset.sum (Finset.range (M + 1)) a := by
    simpa using (Fin.sum_univ_eq_sum_range (f := fun n : ℕ => a n) (n := M + 1))
  have hrange_le_tsum :
      Finset.sum (Finset.range (M + 1)) a ≤ ∑' n : ℕ, a n := by
    exact Summable.sum_le_tsum (s := Finset.range (M + 1)) (f := a)
      (hs := by intro n hn; exact ha_nonneg n) hsum_a
  calc
    Finset.sum (operatorSpecN M) (fun t => (1 : ℝ) / (‖criticalLinePoint t‖ ^ (2 : ℕ)))
        ≤ Finset.sum Finset.univ
            (fun i : Fin (M + 1) =>
              (1 : ℝ) / (‖criticalLinePoint (operatorEigenvalues M i)‖ ^ (2 : ℕ))) := himage_le
    _ ≤ Finset.sum Finset.univ (fun i : Fin (M + 1) => a i.1) := hindex_le
    _ = Finset.sum (Finset.range (M + 1)) a := hsum_range
    _ ≤ ∑' n : ℕ, a n := hrange_le_tsum

/-- Path-B finite-level full-spectrum bound under permutation-invariant center-gap:
if at each level eigenvalues can be bijectively paired to diagonal centers
`k + 29/16` with uniform gap `≤ 12/11`, then the complete inverse-square sum over
`operatorSpecN M` is uniformly bounded in `M`. -/
theorem exists_uniform_bound_sum_one_div_normSq_operatorSpecN_of_centerGapPermutationInvariant
    (hP : OperatorCenterGapPermutationInvariant) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ M : ℕ,
      Finset.sum (operatorSpecN M) (fun t => (1 : ℝ) / (‖criticalLinePoint t‖ ^ (2 : ℕ))) ≤ C := by
  let a : ℕ → ℝ := fun n => (4 : ℝ) * ((1 : ℝ) / (((n + 1 : ℕ) : ℝ) ^ (2 : ℕ)))
  have ha_nonneg : ∀ n : ℕ, 0 ≤ a n := by
    intro n
    dsimp [a]
    positivity
  have hsum_a : Summable a := by
    have hbase0 : Summable (fun n : ℕ => (1 : ℝ) / ((n : ℝ) ^ (2 : ℕ))) :=
      (Real.summable_one_div_nat_pow).2 (by norm_num)
    have hshift : Summable (fun n : ℕ => (1 : ℝ) / (((n + 1 : ℕ) : ℝ) ^ (2 : ℕ))) := by
      simpa [Nat.cast_add, add_assoc, add_comm, add_left_comm] using
        (summable_nat_add_iff 1).2 hbase0
    simpa [a] using hshift.mul_left (4 : ℝ)
  refine ⟨∑' n : ℕ, a n, tsum_nonneg ha_nonneg, ?_⟩
  intro M
  classical
  rcases hP M with ⟨σ, hσ⟩
  have himage_le :
      Finset.sum (operatorSpecN M) (fun t => (1 : ℝ) / (‖criticalLinePoint t‖ ^ (2 : ℕ)))
        ≤ Finset.sum Finset.univ
            (fun i : Fin (M + 1) =>
              (1 : ℝ) / (‖criticalLinePoint (operatorEigenvalues M i)‖ ^ (2 : ℕ))) := by
    unfold operatorSpecN
    simpa using
      (Finset.sum_image_le_of_nonneg
        (s := (Finset.univ : Finset (Fin (M + 1))))
        (g := fun i : Fin (M + 1) => operatorEigenvalues M i)
        (f := fun t : ℝ => (1 : ℝ) / (‖criticalLinePoint t‖ ^ (2 : ℕ)))
        (by intro u hu; positivity))
  have hindex_le :
      Finset.sum Finset.univ
          (fun i : Fin (M + 1) =>
            (1 : ℝ) / (‖criticalLinePoint (operatorEigenvalues M i)‖ ^ (2 : ℕ)))
        ≤ Finset.sum Finset.univ (fun i : Fin (M + 1) => a ((σ i).1)) := by
    refine Finset.sum_le_sum (fun i _ => ?_)
    let t : ℝ := operatorEigenvalues M i
    let k : ℕ := (σ i).1
    have hgap_i : |t - ((k : ℝ) + (29 : ℝ) / 16)| ≤ (12 : ℝ) / 11 := by
      simpa [t, k] using hσ i
    have hpair := abs_le.mp hgap_i
    have hlower :
        (k : ℝ) + (127 : ℝ) / 176 ≤ t := by
      have hconst : (k : ℝ) + (29 : ℝ) / 16 - (12 : ℝ) / 11
          = (k : ℝ) + (127 : ℝ) / 176 := by
        ring_nf
      linarith [hpair.1, hconst]
    have hhalf_le : (((k + 1 : ℕ) : ℝ) / 2) ≤ t := by
      have hle_half : (((k + 1 : ℕ) : ℝ) / 2) ≤ (k : ℝ) + (1 : ℝ) / 2 := by
        calc
          (((k + 1 : ℕ) : ℝ) / 2)
              = ((k : ℝ) + 1) / 2 := by norm_num [Nat.cast_add]
          _ = (k : ℝ) / 2 + (1 : ℝ) / 2 := by ring
          _ ≤ (k : ℝ) + (1 : ℝ) / 2 := by nlinarith
      have hconst : (1 : ℝ) / 2 ≤ (127 : ℝ) / 176 := by norm_num
      have hle_const : (k : ℝ) + (1 : ℝ) / 2 ≤ (k : ℝ) + (127 : ℝ) / 176 := by
        simpa [add_assoc, add_comm, add_left_comm] using add_le_add_left hconst (k : ℝ)
      exact le_trans (le_trans hle_half hle_const) hlower
    have hnonneg_half : 0 ≤ (((k + 1 : ℕ) : ℝ) / 2) := by positivity
    have hnonneg_t : 0 ≤ t := le_trans hnonneg_half hhalf_le
    have hnorm_ge_t : t ≤ ‖criticalLinePoint t‖ := by
      have him_abs : |(criticalLinePoint t).im| ≤ ‖criticalLinePoint t‖ :=
        Complex.abs_im_le_norm (criticalLinePoint t)
      have him_eq : |(criticalLinePoint t).im| = t := by
        simpa [criticalLinePoint_im, abs_of_nonneg hnonneg_t]
      simpa [him_eq] using him_abs
    have hnorm_lower : (((k + 1 : ℕ) : ℝ) / 2) ≤ ‖criticalLinePoint t‖ :=
      le_trans hhalf_le hnorm_ge_t
    have hpow_lower :
        ((((k + 1 : ℕ) : ℝ) / 2) ^ (2 : ℕ))
          ≤ ‖criticalLinePoint t‖ ^ (2 : ℕ) := by
      exact pow_le_pow_left₀ hnonneg_half hnorm_lower (2 : ℕ)
    have hpow_half_pos : 0 < ((((k + 1 : ℕ) : ℝ) / 2) ^ (2 : ℕ)) := by
      positivity
    have hrecip :
        (1 : ℝ) / (‖criticalLinePoint t‖ ^ (2 : ℕ))
          ≤ (1 : ℝ) / ((((k + 1 : ℕ) : ℝ) / 2) ^ (2 : ℕ)) := by
      exact one_div_le_one_div_of_le hpow_half_pos hpow_lower
    have hk1_ne : (((k + 1 : ℕ) : ℝ)) ≠ 0 := by
      exact_mod_cast Nat.succ_ne_zero k
    have hrewrite :
        (1 : ℝ) / ((((k + 1 : ℕ) : ℝ) / 2) ^ (2 : ℕ)) = a k := by
      dsimp [a]
      field_simp [hk1_ne]
      ring
    have hmain :
        (1 : ℝ) / (‖criticalLinePoint t‖ ^ (2 : ℕ)) ≤ a k := by
      rw [← hrewrite]
      exact hrecip
    simpa [t, k] using hmain
  have hsum_perm :
      Finset.sum Finset.univ (fun i : Fin (M + 1) => a ((σ i).1))
        = Finset.sum Finset.univ (fun i : Fin (M + 1) => a i.1) := by
    simpa using (Equiv.sum_comp σ (fun i : Fin (M + 1) => a i.1))
  have hsum_range :
      Finset.sum Finset.univ (fun i : Fin (M + 1) => a i.1)
        = Finset.sum (Finset.range (M + 1)) a := by
    simpa using (Fin.sum_univ_eq_sum_range (f := fun n : ℕ => a n) (n := M + 1))
  have hrange_le_tsum :
      Finset.sum (Finset.range (M + 1)) a ≤ ∑' n : ℕ, a n := by
    exact Summable.sum_le_tsum (s := Finset.range (M + 1)) (f := a)
      (hs := by intro n hn; exact ha_nonneg n) hsum_a
  calc
    Finset.sum (operatorSpecN M) (fun t => (1 : ℝ) / (‖criticalLinePoint t‖ ^ (2 : ℕ)))
        ≤ Finset.sum Finset.univ
            (fun i : Fin (M + 1) =>
              (1 : ℝ) / (‖criticalLinePoint (operatorEigenvalues M i)‖ ^ (2 : ℕ))) := himage_le
    _ ≤ Finset.sum Finset.univ (fun i : Fin (M + 1) => a ((σ i).1)) := hindex_le
    _ = Finset.sum Finset.univ (fun i : Fin (M + 1) => a i.1) := hsum_perm
    _ = Finset.sum (Finset.range (M + 1)) a := hsum_range
    _ ≤ ∑' n : ℕ, a n := hrange_le_tsum

/-- Static (non-recurrence) summability route:
the finite-level uniform inverse-square bound follows directly from
permutation-invariant center-gap geometry, with no Sturm step-exclusion
induction required. -/
theorem exists_uniform_bound_sum_one_div_normSq_operatorSpecN_of_static_center_gap_geometry
    (hP : OperatorCenterGapPermutationInvariant) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ M : ℕ,
      Finset.sum (operatorSpecN M) (fun t => (1 : ℝ) / (‖criticalLinePoint t‖ ^ (2 : ℕ))) ≤ C := by
  exact exists_uniform_bound_sum_one_div_normSq_operatorSpecN_of_centerGapPermutationInvariant hP

/-- One-assumption Weyl reduction:
the Weyl center-gap contract directly yields the finite-level uniform inverse-square bound. -/
theorem exists_uniform_bound_sum_one_div_normSq_operatorSpecN_of_weylCenterGap
    (hW : OperatorWeylCenterGap) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ M : ℕ,
      Finset.sum (operatorSpecN M) (fun t => (1 : ℝ) / (‖criticalLinePoint t‖ ^ (2 : ℕ))) ≤ C := by
  exact exists_uniform_bound_sum_one_div_normSq_operatorSpecN_of_centerGapPermutationInvariant
    (operatorCenterGapPermutationInvariant_of_weylCenterGap hW)

/-- Sturm-route finite-level uniform inverse-square bound:
if a Sturm contract is available together with a bridge to permutation-invariant
center-gap pairing, we obtain the same Path-B uniform summability bound. -/
theorem exists_uniform_bound_sum_one_div_normSq_operatorSpecN_of_sturm
    (hS : OperatorSturmCountContract)
    (hBridge : OperatorSturmCountContract → OperatorCenterGapPermutationInvariant) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ M : ℕ,
      Finset.sum (operatorSpecN M) (fun t => (1 : ℝ) / (‖criticalLinePoint t‖ ^ (2 : ℕ))) ≤ C := by
  exact exists_uniform_bound_sum_one_div_normSq_operatorSpecN_of_centerGapPermutationInvariant
    (hBridge hS)

/-- Direct Sturm+Weyl reduction:
if the corrected Sturm-route counting contract holds and Weyl center-gap is
available, then the finite-level uniform inverse-square bound follows without an
explicit bridge parameter. -/
theorem exists_uniform_bound_sum_one_div_normSq_operatorSpecN_of_sturm_and_weylCenterGap
    (hS : OperatorSturmCountContract)
    (hW : OperatorWeylCenterGap) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ M : ℕ,
      Finset.sum (operatorSpecN M) (fun t => (1 : ℝ) / (‖criticalLinePoint t‖ ^ (2 : ℕ))) ≤ C := by
  exact exists_uniform_bound_sum_one_div_normSq_operatorSpecN_of_sturm hS
    (fun _ => operatorCenterGapPermutationInvariant_of_weylCenterGap hW)

/-- Sturm-to-summability route (minimal step-compatibility form):
if eigenvalue counting matches Sturm sign-variation counting and the recurrence-step
compatibility holds, then any bridge from Sturm contract to permutation-invariant
center-gap pairing yields the Path-B uniform inverse-square bound. -/
theorem exists_uniform_bound_sum_one_div_normSq_operatorSpecN_of_sturm_stepCompatibility
    (hEigSturm :
      ∀ M : ℕ, ∀ x : ℝ,
        operatorEigenvalueCountLE M x = operatorSturmSignVariationCount M x)
    (hCompat : OperatorSturmStepCompatibility)
    (hBridge : OperatorSturmCountContract → OperatorCenterGapPermutationInvariant) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ M : ℕ,
      Finset.sum (operatorSpecN M) (fun t => (1 : ℝ) / (‖criticalLinePoint t‖ ^ (2 : ℕ))) ≤ C := by
  exact exists_uniform_bound_sum_one_div_normSq_operatorSpecN_of_sturm
    (operatorSturmCountContract_of_signVariationBridge_and_stepCompatibility hEigSturm hCompat)
    hBridge

/-- Step-compatibility route with Weyl center-gap:
the explicit bridge argument is eliminated by the canonical Weyl→permutation
map. -/
theorem exists_uniform_bound_sum_one_div_normSq_operatorSpecN_of_sturm_stepCompatibility_and_weylCenterGap
    (hEigSturm :
      ∀ M : ℕ, ∀ x : ℝ,
        operatorEigenvalueCountLE M x = operatorSturmSignVariationCount M x)
    (hCompat : OperatorSturmStepCompatibility)
    (hW : OperatorWeylCenterGap) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ M : ℕ,
      Finset.sum (operatorSpecN M) (fun t => (1 : ℝ) / (‖criticalLinePoint t‖ ^ (2 : ℕ))) ≤ C := by
  exact exists_uniform_bound_sum_one_div_normSq_operatorSpecN_of_sturm_stepCompatibility
    hEigSturm hCompat
    (fun _ => operatorCenterGapPermutationInvariant_of_weylCenterGap hW)

/-- Sturm-to-summability route (edge-lock form):
if eigenvalue counting matches Sturm sign-variation counting and edge-lock holds
for the recurrence increments, then any bridge from Sturm contract to
permutation-invariant center-gap pairing yields the Path-B uniform inverse-square
bound. -/
theorem exists_uniform_bound_sum_one_div_normSq_operatorSpecN_of_sturm_edgeLock
    (hEigSturm :
      ∀ M : ℕ, ∀ x : ℝ,
        operatorEigenvalueCountLE M x = operatorSturmSignVariationCount M x)
    (hEdgeLock :
      ∀ M : ℕ, ∀ x : ℝ,
        (operatorSturmSignVariationCount M x = operatorCenterCountLE M x + 1 →
          operatorCenterAt (M + 1) ≤ x) ∧
        (operatorCenterCountLE M x = operatorSturmSignVariationCount M x + 1 →
          operatorSturmSign (operatorSturmP (M + 1) x) ≠
            operatorSturmSign (operatorSturmP (M + 2) x)))
    (hBridge : OperatorSturmCountContract → OperatorCenterGapPermutationInvariant) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ M : ℕ,
      Finset.sum (operatorSpecN M) (fun t => (1 : ℝ) / (‖criticalLinePoint t‖ ^ (2 : ℕ))) ≤ C := by
  exact exists_uniform_bound_sum_one_div_normSq_operatorSpecN_of_sturm
    (operatorSturmCountContract_of_signVariationBridge_and_edgeLock hEigSturm hEdgeLock)
    hBridge

/-- Edge-lock route with Weyl center-gap:
the explicit bridge argument is eliminated by the canonical Weyl→permutation
map. -/
theorem exists_uniform_bound_sum_one_div_normSq_operatorSpecN_of_sturm_edgeLock_and_weylCenterGap
    (hEigSturm :
      ∀ M : ℕ, ∀ x : ℝ,
        operatorEigenvalueCountLE M x = operatorSturmSignVariationCount M x)
    (hEdgeLock :
      ∀ M : ℕ, ∀ x : ℝ,
        (operatorSturmSignVariationCount M x = operatorCenterCountLE M x + 1 →
          operatorCenterAt (M + 1) ≤ x) ∧
        (operatorCenterCountLE M x = operatorSturmSignVariationCount M x + 1 →
          operatorSturmSign (operatorSturmP (M + 1) x) ≠
            operatorSturmSign (operatorSturmP (M + 2) x)))
    (hW : OperatorWeylCenterGap) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ M : ℕ,
      Finset.sum (operatorSpecN M) (fun t => (1 : ℝ) / (‖criticalLinePoint t‖ ^ (2 : ℕ))) ≤ C := by
  exact exists_uniform_bound_sum_one_div_normSq_operatorSpecN_of_sturm_edgeLock
    hEigSturm hEdgeLock
    (fun _ => operatorCenterGapPermutationInvariant_of_weylCenterGap hW)

/-- Cardinality growth control for the concrete operator spectral ladder. -/
theorem card_operatorSpecN_le (N : ℕ) :
    (operatorSpecN N).card ≤ N + 1 := by
  classical
  unfold operatorSpecN
  simpa using (Finset.card_image_le (f := operatorEigenvalues N) (s := Finset.univ))

/-- There is always a finite ordinate envelope at each level `N`. -/
theorem exists_operator_ordinate_envelope (N : ℕ) :
    ∃ B : ℝ, 0 ≤ B ∧ ∀ t ∈ operatorSpecN N, ‖criticalLinePoint t‖ ≤ B := by
  classical
  let B :=
    (operatorSpecN N).sup' (operatorSpecN_nonempty N) (fun t => ‖criticalLinePoint t‖)
  refine ⟨B, ?_, ?_⟩
  · rcases operatorSpecN_nonempty N with ⟨t0, ht0⟩
    have ht0le : ‖criticalLinePoint t0‖ ≤ B := by
      exact Finset.le_sup' (s := operatorSpecN N) (f := fun t => ‖criticalLinePoint t‖) ht0
    exact le_trans (norm_nonneg _) ht0le
  · intro t ht
    exact Finset.le_sup' (s := operatorSpecN N) (f := fun t => ‖criticalLinePoint t‖) ht

/-- Uniform closed-ball norm control at level `N` from an ordinate envelope. -/
theorem norm_operatorXiFiniteLadder_le_on_closedBall
    (N : ℕ) (R B : ℝ)
    (hB0 : 0 ≤ B)
    (hB : ∀ t ∈ operatorSpecN N, ‖criticalLinePoint t‖ ≤ B)
    {s : ℂ} (hs : s ∈ Metric.closedBall (0 : ℂ) R) :
    ‖operatorXiFiniteLadder N s‖ ≤ (R + B) ^ (operatorSpecN N).card := by
  have hsR : ‖s‖ ≤ R := by
    simpa [Metric.mem_closedBall, dist_eq_norm] using hs
  have hbase_nonneg : 0 ≤ ‖s‖ + B := by
    exact add_nonneg (norm_nonneg _) hB0
  have hbase_le : ‖s‖ + B ≤ R + B := by
    linarith [hsR]
  have hpow_le :
      (‖s‖ + B) ^ (operatorSpecN N).card ≤ (R + B) ^ (operatorSpecN N).card := by
    exact pow_le_pow_left₀ hbase_nonneg hbase_le _
  have hN :=
    norm_operatorXiFiniteLadder_le_pow_card_of_ordinate_bound N s B hB
  exact le_trans hN hpow_le

/-- Existence form of uniform closed-ball operator control at level `N`. -/
theorem exists_closedBall_uniform_operator_bound
    (N : ℕ) (R : ℝ) :
    ∃ B : ℝ, 0 ≤ B ∧
      ∀ s ∈ Metric.closedBall (0 : ℂ) R,
        ‖operatorXiFiniteLadder N s‖ ≤ (R + B) ^ (operatorSpecN N).card := by
  rcases exists_operator_ordinate_envelope N with ⟨B, hB0, hB⟩
  refine ⟨B, hB0, ?_⟩
  intro s hs
  exact norm_operatorXiFiniteLadder_le_on_closedBall N R B hB0 hB hs

/-- Operator-ladder telescoping norm control over `m` successive levels. -/
theorem norm_operatorXiFiniteLadder_sub_le_sum_steps
    (s : ℂ) (n m : ℕ) :
    ‖operatorXiFiniteLadder (n + m) s - operatorXiFiniteLadder n s‖ ≤
      Finset.sum (Finset.range m)
        (fun k => ‖operatorXiFiniteLadder (n + (k + 1)) s - operatorXiFiniteLadder (n + k) s‖) := by
  simpa using norm_sub_le_sum_steps (fun j => operatorXiFiniteLadder j s) n m

/-- Operator-ladder telescoping control from explicit per-step bounds. -/
theorem norm_operatorXiFiniteLadder_sub_le_sum_stepBounds
    (s : ℂ) (a : ℕ → ℝ)
    (hstep : ∀ j : ℕ,
      ‖operatorXiFiniteLadder (j + 1) s - operatorXiFiniteLadder j s‖ ≤ a j)
    (n m : ℕ) :
    ‖operatorXiFiniteLadder (n + m) s - operatorXiFiniteLadder n s‖ ≤
      Finset.sum (Finset.range m) (fun k => a (n + k)) := by
  simpa using norm_sub_le_sum_stepBounds (fun j => operatorXiFiniteLadder j s) a hstep n m

/-- Finite zero witnesses for the operator Xi ladder:
every finite-level zero is exactly a critical-line point listed in `operatorSpecN N`. -/
theorem finiteZeroWitness_operatorXiFiniteLadder :
    FiniteZeroWitness operatorXiFiniteLadder operatorSpecN := by
  intro N z hz0
  exact (XiFinite_zero_iff_exists (operatorSpecN N) z).1 hz0

/-- Constant tolerance profile for the concrete operator ladder. -/
def operatorTolConst (τ : ℝ) : ℕ → ℝ := fun _ => τ

/-- Local-uniform convergence implies nonzero-tolerance approximate zero capture
at every target zero. -/
theorem operatorApprox_of_locallyUniform_constTol
    (hconv : TendstoLocallyUniformly operatorXiFiniteLadder XiTarget
      (Filter.atTop : Filter ℕ))
    {τ : ℝ} (hτ : 0 < τ) :
    ApproxZeroConvergence XiTarget operatorXiFiniteLadder (operatorTolConst τ) := by
  intro s hsXi
  have hconvOn : TendstoLocallyUniformlyOn operatorXiFiniteLadder XiTarget
      (Filter.atTop : Filter ℕ) (Set.univ) := by
    simpa [tendstoLocallyUniformlyOn_univ] using hconv
  have hpt : Filter.Tendsto (fun N : ℕ => operatorXiFiniteLadder N s)
      (Filter.atTop : Filter ℕ) (𝓝 (XiTarget s)) :=
    hconvOn.tendsto_at (by simp)
  have hball : ∀ᶠ N : ℕ in (Filter.atTop : Filter ℕ),
      ‖operatorXiFiniteLadder N s - XiTarget s‖ < τ := by
    simpa [Metric.ball, dist_eq_norm] using
      (hpt (Metric.ball_mem_nhds (XiTarget s) hτ))
  rcases hball.exists_forall_of_atTop with ⟨N0, hN0⟩
  refine ⟨N0, ?_⟩
  have hlt : ‖operatorXiFiniteLadder N0 s - XiTarget s‖ < τ := hN0 N0 (le_rfl)
  have hlt0 : ‖operatorXiFiniteLadder N0 s‖ < τ := by
    simpa [hsXi] using hlt
  exact le_of_lt hlt0

/-- If local-uniform convergence is available and the constant tolerance level
is rigid, one gets zero-forward transfer for the operator Xi ladder. -/
theorem operatorZeroForward_of_locallyUniform_and_constRigidity
    (hconv : TendstoLocallyUniformly operatorXiFiniteLadder XiTarget
      (Filter.atTop : Filter ℕ))
    {τ : ℝ} (hτ : 0 < τ)
    (hRigid : SpectralRigidity operatorXiFiniteLadder (operatorTolConst τ)) :
    ZeroForwardTransfer XiTarget operatorXiFiniteLadder := by
  exact zeroForward_of_convergence_and_rigidity
    XiTarget
    operatorXiFiniteLadder
    (operatorTolConst τ)
    (operatorApprox_of_locallyUniform_constTol hconv hτ)
    hRigid

/-- Operator-ladder convergence obligation at zero tolerance. -/
def OperatorApproxZeroConvergence : Prop :=
  ApproxZeroConvergence
    XiTarget
    operatorXiFiniteLadder
    Gutoe.RiemannTargetFiniteLadder.tolZero

/-- For the concrete operator ladder, zero-tolerance approximation is exactly
finite-ladder `XiTarget` zero-capture. -/
theorem operatorApproxZero_iff_xiTarget_capture :
    OperatorApproxZeroConvergence ↔ XiTargetLadderZeroCapture operatorSpecN := by
  simpa [OperatorApproxZeroConvergence, operatorXiFiniteLadder] using
    (approxZero_tolZero_iff_zeroCapture operatorSpecN)

/-- Convergence closure in the operator lane:
zero-tolerance approximation implies nontrivial-`ζ` capture for `operatorSpecN`. -/
theorem operator_enumerated_capture_of_approxZero
    (hApprox : OperatorApproxZeroConvergence) :
    OperatorEnumeratedNontrivialCapture := by
  exact (nontrivial_capture_iff_xiTarget_capture operatorSpecN).2
    ((operatorApproxZero_iff_xiTarget_capture).1 hApprox)

/-- RH closure from operator-ladder convergence obligation. -/
theorem mathlibRH_of_operator_approxZero
    (hApprox : OperatorApproxZeroConvergence) :
    RiemannHypothesis := by
  exact mathlibRH_of_operator_enumerated_nontrivial_capture
    (operator_enumerated_capture_of_approxZero hApprox)

/-- Positive-threshold rigidity surface: there exists at least one strictly
positive tolerance level at which small values force exact finite-level zeros. -/
def OperatorPositiveRigidity : Prop :=
  ∃ τ : ℝ, 0 < τ ∧ SpectralRigidity operatorXiFiniteLadder (operatorTolConst τ)

/-- Local-uniform convergence plus any positive-threshold rigidity already
forces the zero-tolerance operator approximation obligation. -/
theorem operatorApproxZero_of_locallyUniform_and_positiveRigidity
    (hconv : TendstoLocallyUniformly operatorXiFiniteLadder XiTarget
      (Filter.atTop : Filter ℕ))
    (hPosRigid : OperatorPositiveRigidity) :
    OperatorApproxZeroConvergence := by
  rcases hPosRigid with ⟨τ, hτ, hRigid⟩
  intro s hsXi
  have hApproxτ :
      ApproxZeroConvergence XiTarget operatorXiFiniteLadder (operatorTolConst τ) :=
    operatorApprox_of_locallyUniform_constTol hconv hτ
  rcases hApproxτ s hsXi with ⟨N, hNτ⟩
  have hZero : operatorXiFiniteLadder N s = 0 := hRigid N s hNτ
  refine ⟨N, ?_⟩
  simpa [Gutoe.RiemannTargetFiniteLadder.tolZero, hZero]

/-- RH closure from the analytic surface
`local-uniform convergence + positive-threshold rigidity`. -/
theorem mathlibRH_of_locallyUniform_and_positiveRigidity
    (hconv : TendstoLocallyUniformly operatorXiFiniteLadder XiTarget
      (Filter.atTop : Filter ℕ))
    (hPosRigid : OperatorPositiveRigidity) :
    RiemannHypothesis := by
  exact mathlibRH_of_operator_approxZero
    (operatorApproxZero_of_locallyUniform_and_positiveRigidity hconv hPosRigid)

/-- The two operator capture surfaces are equivalent:
finite ladder capture (`operatorSpecN`) and direct spectrum-set capture (`operatorSpecSet`). -/
theorem operator_enumerated_capture_iff_operator_set_capture :
    OperatorEnumeratedNontrivialCapture ↔ OperatorSetNontrivialCapture := by
  constructor
  · intro hCap
    intro s hs htriv h1
    rcases hCap s hs htriv h1 with ⟨N, t, ht, hsEq⟩
    exact ⟨N, t, (mem_operatorSpecN_iff_mem_operatorSpecSet N t).1 ht, hsEq⟩
  · intro hCap
    intro s hs htriv h1
    rcases hCap s hs htriv h1 with ⟨N, t, ht, hsEq⟩
    exact ⟨N, t, (mem_operatorSpecN_iff_mem_operatorSpecSet N t).2 ht, hsEq⟩

/-- Endgame equivalence in operator-native form:
the convergence obligation is equivalent to direct operator-spectrum capture. -/
theorem operatorApproxZero_iff_operator_set_capture :
    OperatorApproxZeroConvergence ↔ OperatorSetNontrivialCapture := by
  exact (operatorApproxZero_iff_xiTarget_capture.trans
    (nontrivial_capture_iff_xiTarget_capture operatorSpecN).symm).trans
      operator_enumerated_capture_iff_operator_set_capture

/-- A finite operator-ladder witness immediately yields set-capture on the
explicit operator spectrum lane. -/
theorem operator_set_capture_of_operator_ladder
    (hOp : OperatorSpecLadder)
    (hCap : OperatorNontrivialCapture hOp) :
    OperatorSetNontrivialCapture := by
  intro s hs htriv h1
  rcases hCap s hs htriv h1 with ⟨N, t, ht, hsEq⟩
  refine ⟨N, t, ?_, hsEq⟩
  exact (ordinateIsEigenvalue_iff_mem_operatorSpecSet N t).1
    ((hOp.eigenvalue_exact N t).1 ht)

/-- RH closure from the direct operator spectrum set-capture surface. -/
theorem mathlibRH_of_operator_set_nontrivial_capture
    (hCap : OperatorSetNontrivialCapture) :
    RiemannHypothesis := by
  intro s hs htriv h1
  rcases hCap s hs htriv h1 with ⟨_N, t, _ht, hsEq⟩
  simpa [hsEq, onCriticalLine, criticalLinePoint_re] using (criticalLinePoint_re t)

/-- Approximate operator-spectrum capture:
every `XiTarget` zero can be approximated arbitrarily well by critical-line
points coming from finite operator spectra. -/
def OperatorApproximateCapture : Prop :=
  ∀ s : ℂ, XiTarget s = 0 → ∀ ε : ℝ, 0 < ε →
    ∃ N : ℕ, ∃ t : ℝ, t ∈ operatorSpecSet N ∧ ‖s - criticalLinePoint t‖ < ε

/-- Hurwitz-output surface specialized to the concrete operator lane:
every `XiTarget` zero can be approximated by an actual finite-level zero of the
operator finite product. -/
def OperatorHurwitzZeroApproxTransfer : Prop :=
  ∀ s : ℂ, XiTarget s = 0 → ∀ ε : ℝ, 0 < ε →
    ∃ N : ℕ, ∃ z : ℂ,
      XiFinite (operatorSpecN N) z = 0 ∧ ‖s - z‖ < ε

/-- Concrete operator-lane Hurwitz kernel obligation. -/
def OperatorHurwitzKernel : Prop :=
  HurwitzZeroApproxKernel XiTarget operatorXiFiniteLadder

/-- Instantiation of the abstract Hurwitz kernel on the concrete operator lane:
all regularity/nonvanishing side conditions are discharged from existing lemmas. -/
theorem operatorHurwitzTransfer_of_kernel
    (hKernel : OperatorHurwitzKernel)
    (hconv : TendstoLocallyUniformly operatorXiFiniteLadder XiTarget
      (Filter.atTop : Filter ℕ)) :
    OperatorHurwitzZeroApproxTransfer := by
  exact hurwitzTransfer_of_kernel
    (Xi := XiTarget)
    (XiN := operatorXiFiniteLadder)
    hKernel
    hconv
    differentiable_operatorXiFiniteLadder
    (fun s hs => differentiableAt_XiTarget_of_zero hs)
    xiFinite_operatorSpecN_zero_ne

/-- The operator Hurwitz-output surface implies the internal approximate-capture
surface used by the RH closure theorem. -/
theorem operatorApproximateCapture_of_hurwitzTransfer
    (hHurwitz : OperatorHurwitzZeroApproxTransfer) :
    OperatorApproximateCapture := by
  have hApproxN : ApproximateCriticalCapture XiTarget operatorSpecN :=
    approximateCriticalCapture_of_hurwitzTransfer_and_witness
      (Xi := XiTarget)
      (XiN := operatorXiFiniteLadder)
      (specN := operatorSpecN)
      hHurwitz
      finiteZeroWitness_operatorXiFiniteLadder
  intro s hsXi ε hε
  rcases hApproxN s hsXi ε hε with ⟨N, t, htN, hdist⟩
  exact ⟨N, t, (mem_operatorSpecN_iff_mem_operatorSpecSet N t).1 htN, hdist⟩

/-- The approximate-capture surface also yields the Hurwitz-output surface:
finite-level critical-line witnesses are concrete finite-product zeros. -/
theorem operatorHurwitzTransfer_of_operatorApproximateCapture
    (hApprox : OperatorApproximateCapture) :
    OperatorHurwitzZeroApproxTransfer := by
  intro s hsXi ε hε
  rcases hApprox s hsXi ε hε with ⟨N, t, htSet, hdist⟩
  have htN : t ∈ operatorSpecN N :=
    (mem_operatorSpecN_iff_mem_operatorSpecSet N t).2 htSet
  refine ⟨N, criticalLinePoint t, XiFinite_zero_of_mem (operatorSpecN N) htN, ?_⟩
  simpa using hdist

/-- Exact equivalence between the operator Hurwitz-output and approximate-capture
surfaces. This isolates the remaining analytic gap to either form. -/
theorem operatorHurwitz_iff_operatorApproximateCapture :
    OperatorHurwitzZeroApproxTransfer ↔ OperatorApproximateCapture := by
  constructor
  · exact operatorApproximateCapture_of_hurwitzTransfer
  · exact operatorHurwitzTransfer_of_operatorApproximateCapture

/-- RH closure from approximate operator-spectrum capture.
This is the exact interface needed for a future Hurwitz-style transfer lemma. -/
theorem mathlibRH_of_operator_approximate_capture
    (hApprox : OperatorApproximateCapture) :
    RiemannHypothesis := by
  intro s hs htriv h1
  have hsXi : XiTarget s = 0 := nontrivialZeroTransferToXiTarget s hs htriv h1
  have hsmall : ∀ ε : ℝ, 0 < ε → |s.re - (1 / 2 : ℝ)| < ε := by
    intro ε hε
    rcases hApprox s hsXi ε hε with ⟨N, t, _ht, hdist⟩
    have hreLe : |(s - criticalLinePoint t).re| ≤ ‖s - criticalLinePoint t‖ := by
      simpa using (Complex.abs_re_le_norm (s - criticalLinePoint t))
    have hreLt : |(s - criticalLinePoint t).re| < ε := lt_of_le_of_lt hreLe hdist
    simpa [criticalLinePoint_re, Complex.sub_re] using hreLt
  have hEq : s.re = (1 / 2 : ℝ) := by
    by_contra hne
    let δ : ℝ := |s.re - (1 / 2 : ℝ)| / 2
    have hδpos : 0 < δ := by
      unfold δ
      exact half_pos (abs_pos.mpr (sub_ne_zero.mpr hne))
    have hlt : |s.re - (1 / 2 : ℝ)| < δ := hsmall δ hδpos
    have hge : δ ≤ |s.re - (1 / 2 : ℝ)| := by
      unfold δ
      nlinarith [abs_nonneg (s.re - (1 / 2 : ℝ))]
    exact (not_lt_of_ge hge) hlt
  exact hEq

/-- RH closure from the operator Hurwitz-output obligation. -/
theorem mathlibRH_of_operator_hurwitz_zero_approx
    (hHurwitz : OperatorHurwitzZeroApproxTransfer) :
    RiemannHypothesis := by
  exact mathlibRH_of_operator_approximate_capture
    (operatorApproximateCapture_of_hurwitzTransfer hHurwitz)

/-- RH closure from the instantiated concrete operator Hurwitz kernel plus
local-uniform convergence. -/
theorem mathlibRH_of_operator_hurwitzKernel_and_locallyUniform
    (hKernel : OperatorHurwitzKernel)
    (hconv : TendstoLocallyUniformly operatorXiFiniteLadder XiTarget
      (Filter.atTop : Filter ℕ)) :
    RiemannHypothesis := by
  exact mathlibRH_of_operator_hurwitz_zero_approx
    (operatorHurwitzTransfer_of_kernel hKernel hconv)

/-- Named local-uniform convergence obligation on the concrete operator lane. -/
def OperatorXiFiniteLocallyUniformConvergence : Prop :=
  TendstoLocallyUniformly operatorXiFiniteLadder XiTarget
    (Filter.atTop : Filter ℕ)

/-- Existence of a locally-uniform limit for the concrete operator finite ladder
along `atTop`. -/
def OperatorXiFiniteLocallyUniformLimitExists : Prop :=
  ∃ F : ℂ → ℂ,
    TendstoLocallyUniformly operatorXiFiniteLadder F (Filter.atTop : Filter ℕ)

/-- Overconstrained uniqueness surface:
any locally-uniform limit of the concrete operator finite ladder is forced to be
`XiTarget`. This encodes the "single occupant in an overconstrained function
class" principle. -/
def OperatorXiFiniteLocallyUniformLimitUniqueness : Prop :=
  ∀ F : ℂ → ℂ,
    TendstoLocallyUniformly operatorXiFiniteLadder F (Filter.atTop : Filter ℕ) →
      F = XiTarget

/-- If a locally-uniform limit exists and the locally-uniform limit is unique in
the overconstrained class (hence equals `XiTarget`), then the concrete
local-uniform convergence obligation follows automatically. -/
theorem operatorXiFiniteLocallyUniformConvergence_of_limitExists_and_uniqueness
    (hExist : OperatorXiFiniteLocallyUniformLimitExists)
    (hUnique : OperatorXiFiniteLocallyUniformLimitUniqueness) :
    OperatorXiFiniteLocallyUniformConvergence := by
  rcases hExist with ⟨F, hF⟩
  have hEq : F = XiTarget := hUnique F hF
  simpa [OperatorXiFiniteLocallyUniformConvergence, hEq] using hF

/-- Named three-obligation RH surface on the concrete operator lane.
1) permutation-invariant center-gap (summability geometry),
2) Hurwitz kernel, and
3) local-uniform convergence of the operator finite ladder to `XiTarget`. -/
def OperatorRHThreeConditionalObligations : Prop :=
  OperatorCenterGapPermutationInvariant ∧
  OperatorHurwitzKernel ∧
  OperatorXiFiniteLocallyUniformConvergence

/-- Minimal RH surface currently used by the compiled closure:
Hurwitz kernel + local-uniform convergence on the concrete operator lane. -/
def OperatorRHTwoConditionalObligations : Prop :=
  OperatorHurwitzKernel ∧
  OperatorXiFiniteLocallyUniformConvergence

/-- RH closure from the minimal two-obligation surface. -/
theorem mathlibRH_of_operator_two_conditional_obligations
    (h2 : OperatorRHTwoConditionalObligations) :
    RiemannHypothesis := by
  rcases h2 with ⟨hKernel, hConv⟩
  exact mathlibRH_of_operator_hurwitzKernel_and_locallyUniform hKernel hConv

/-- The historical three-obligation surface implies the minimal two-obligation
surface by forgetting the extra geometry slot. -/
theorem operator_two_conditional_obligations_of_three
    (h3 : OperatorRHThreeConditionalObligations) :
    OperatorRHTwoConditionalObligations := by
  rcases h3 with ⟨_hGap, hKernel, hConv⟩
  exact ⟨hKernel, hConv⟩

/-- RH closure from the named three-obligation surface. -/
theorem mathlibRH_of_operator_three_conditional_obligations
    (h3 : OperatorRHThreeConditionalObligations) :
    RiemannHypothesis := by
  exact mathlibRH_of_operator_two_conditional_obligations
    (operator_two_conditional_obligations_of_three h3)

/-- RH closure adapter from:
1) explicit same-min endgame obligations, and
2) the minimal operator two-obligation surface. -/
theorem mathlibRH_of_sameMinEndgame_and_two_conditional
    (hMaxAboveInSameMin : OperatorSameMinMaxAboveObligation)
    (hNoPrevPlusOneInSameMin : OperatorSameMinPlusOneNoPrevObligation)
    (h2 : OperatorRHTwoConditionalObligations) :
    RiemannHypothesis := by
  rcases h2 with ⟨hKernel, hConv⟩
  have hGap : OperatorCenterGapPermutationInvariant :=
    operatorCenterGapPermutationInvariant_of_sameMinEndgameObligations
      hMaxAboveInSameMin hNoPrevPlusOneInSameMin
  exact mathlibRH_of_operator_three_conditional_obligations ⟨hGap, hKernel, hConv⟩

/-- Full Clifford-operator endgame bundle:
1) Clifford same-min symmetry contract (geometry lane),
2) concrete operator Hurwitz kernel, and
3) local-uniform convergence of finite operator products to `XiTarget`. -/
def OperatorCliffordRHEndgameObligation : Prop :=
  OperatorCliffordSameMinSymmetryContract ∧
  OperatorHurwitzKernel ∧
  OperatorXiFiniteLocallyUniformConvergence

/-- RH closure from the Clifford-operator endgame bundle. -/
theorem mathlibRH_of_clifford_operator_endgame
    (hEnd : OperatorCliffordRHEndgameObligation) :
    RiemannHypothesis := by
  rcases hEnd with ⟨hCliff, hKernel, hConv⟩
  have hGap : OperatorCenterGapPermutationInvariant :=
    operatorCenterGapPermutationInvariant_of_cliffordSameMinSymmetryContract hCliff
  exact mathlibRH_of_operator_three_conditional_obligations ⟨hGap, hKernel, hConv⟩

/-- Obligation instantiator: from Clifford same-min symmetry + the minimal
operator two-obligation surface, recover the Clifford RH endgame bundle. -/
theorem operatorCliffordRHEndgameObligation_of_clifford_and_two_conditional
    (hCliff : OperatorCliffordSameMinSymmetryContract)
    (h2 : OperatorRHTwoConditionalObligations) :
    OperatorCliffordRHEndgameObligation := by
  rcases h2 with ⟨hKernel, hConv⟩
  exact ⟨hCliff, hKernel, hConv⟩

/-- RH closure adapter: Clifford same-min symmetry + the minimal operator
two-obligation surface. -/
theorem mathlibRH_of_clifford_and_two_conditional
    (hCliff : OperatorCliffordSameMinSymmetryContract)
    (h2 : OperatorRHTwoConditionalObligations) :
    RiemannHypothesis := by
  exact mathlibRH_of_clifford_operator_endgame
    (operatorCliffordRHEndgameObligation_of_clifford_and_two_conditional hCliff h2)

/-- RH closure from the concrete operator Hurwitz kernel plus a summable
nonnegative step profile and pointwise convergence to `XiTarget`. This packages
the local-uniform promotion in one theorem so the remaining analytic inputs are
exactly `(ha_sum, hpt)`. -/
theorem mathlibRH_of_operator_hurwitzKernel_and_stepSummable_pointwise
    (hKernel : OperatorHurwitzKernel)
    (a : ℕ → ℝ)
    (ha_nonneg : ∀ j : ℕ, 0 ≤ a j)
    (ha_sum : Summable a)
    (hstep : ∀ R : ℝ, ∀ j : ℕ, ∀ z : ℂ, z ∈ Metric.closedBall (0 : ℂ) R →
      ‖operatorXiFiniteLadder (j + 1) z - operatorXiFiniteLadder j z‖ ≤ a j)
    (hpt : ∀ z : ℂ,
      Filter.Tendsto (fun N : ℕ => operatorXiFiniteLadder N z)
        (Filter.atTop : Filter ℕ) (𝓝 (XiTarget z))) :
    RiemannHypothesis := by
  exact mathlibRH_of_operator_hurwitzKernel_and_locallyUniform hKernel
    (tendstoLocallyUniformly_operatorXiFiniteLadder_of_stepSummable_and_pointwise
      a ha_nonneg ha_sum hstep hpt)

/-- RH closure from:
1) concrete operator Hurwitz kernel,
2) existence of a locally-uniform limit of the finite operator ladder, and
3) uniqueness identifying any such limit with `XiTarget`.

This removes "pointwise-to-`XiTarget`" as a primitive obligation and replaces it
with an overconstrained uniqueness obligation. -/
theorem mathlibRH_of_operator_hurwitzKernel_and_limitExists_uniqueness
    (hKernel : OperatorHurwitzKernel)
    (hExist : OperatorXiFiniteLocallyUniformLimitExists)
    (hUnique : OperatorXiFiniteLocallyUniformLimitUniqueness) :
    RiemannHypothesis := by
  exact mathlibRH_of_operator_hurwitzKernel_and_locallyUniform hKernel
    (operatorXiFiniteLocallyUniformConvergence_of_limitExists_and_uniqueness
      hExist hUnique)

/-- Endgame surface (overconstrained-limit route):
1) concrete operator Hurwitz kernel, and
2) existence + uniqueness of a locally-uniform limit of the operator finite ladder.

This packages the final analytic burden in one named obligation bundle. -/
def OperatorOverconstrainedEndgameObligation : Prop :=
  OperatorHurwitzKernel ∧
  OperatorXiFiniteLocallyUniformLimitExists ∧
  OperatorXiFiniteLocallyUniformLimitUniqueness

/-- RH closure from the overconstrained endgame surface. -/
theorem mathlibRH_of_operator_overconstrained_endgame
    (hEnd : OperatorOverconstrainedEndgameObligation) :
    RiemannHypothesis := by
  rcases hEnd with ⟨hKernel, hExist, hUnique⟩
  exact mathlibRH_of_operator_hurwitzKernel_and_limitExists_uniqueness
    hKernel hExist hUnique

/-- RH closure adapter from:
1) explicit same-min endgame obligations, and
2) operator overconstrained endgame obligations. -/
theorem mathlibRH_of_sameMinEndgame_and_overconstrained
    (hMaxAboveInSameMin : OperatorSameMinMaxAboveObligation)
    (hNoPrevPlusOneInSameMin : OperatorSameMinPlusOneNoPrevObligation)
    (hEnd : OperatorOverconstrainedEndgameObligation) :
    RiemannHypothesis := by
  rcases hEnd with ⟨hKernel, hExist, hUnique⟩
  have hConv : OperatorXiFiniteLocallyUniformConvergence :=
    operatorXiFiniteLocallyUniformConvergence_of_limitExists_and_uniqueness
      hExist hUnique
  exact mathlibRH_of_sameMinEndgame_and_two_conditional
    hMaxAboveInSameMin hNoPrevPlusOneInSameMin ⟨hKernel, hConv⟩

/-- Obligation instantiator: Clifford same-min symmetry plus the operator
overconstrained endgame obligations instantiate the Clifford RH endgame bundle. -/
theorem operatorCliffordRHEndgameObligation_of_clifford_and_overconstrained
    (hCliff : OperatorCliffordSameMinSymmetryContract)
    (hEnd : OperatorOverconstrainedEndgameObligation) :
    OperatorCliffordRHEndgameObligation := by
  rcases hEnd with ⟨hKernel, hExist, hUnique⟩
  have hConv : OperatorXiFiniteLocallyUniformConvergence :=
    operatorXiFiniteLocallyUniformConvergence_of_limitExists_and_uniqueness
      hExist hUnique
  exact ⟨hCliff, hKernel, hConv⟩

/-- RH closure adapter: Clifford same-min symmetry plus operator overconstrained
endgame obligations. -/
theorem mathlibRH_of_clifford_and_overconstrained_endgame
    (hCliff : OperatorCliffordSameMinSymmetryContract)
    (hEnd : OperatorOverconstrainedEndgameObligation) :
    RiemannHypothesis := by
  exact mathlibRH_of_clifford_operator_endgame
    (operatorCliffordRHEndgameObligation_of_clifford_and_overconstrained
      hCliff hEnd)

/-- RH closure adapter from the Hodge parity contract plus the operator
overconstrained endgame obligations. -/
theorem mathlibRH_of_hodgeParity_and_overconstrained_endgame
    (hHP : OperatorHodgeParityContract)
    (hEnd : OperatorOverconstrainedEndgameObligation) :
    RiemannHypothesis := by
  exact mathlibRH_of_clifford_and_overconstrained_endgame
    (operatorCliffordSameMinSymmetryContract_of_hodgeParity hHP) hEnd

/-- Single hard obligation once Hurwitz kernel is fixed:
existence + uniqueness of the locally-uniform operator-ladder limit. -/
def OperatorSingleHardObligation : Prop :=
  OperatorXiFiniteLocallyUniformLimitExists ∧
  OperatorXiFiniteLocallyUniformLimitUniqueness

/-- RH closure from a fixed Hurwitz kernel plus the single hard obligation. -/
theorem mathlibRH_of_operator_hurwitzKernel_and_single_hard_obligation
    (hKernel : OperatorHurwitzKernel)
    (hHard : OperatorSingleHardObligation) :
    RiemannHypothesis := by
  rcases hHard with ⟨hExist, hUnique⟩
  exact mathlibRH_of_operator_hurwitzKernel_and_limitExists_uniqueness
    hKernel hExist hUnique

/-- Zero-tolerance operator approximation already implies the Hurwitz-output
surface (with exact finite-level zeros at the target point itself). -/
theorem operatorHurwitz_of_operatorApproxZero
    (hApprox0 : OperatorApproxZeroConvergence) :
    OperatorHurwitzZeroApproxTransfer := by
  have hCapXi : XiTargetLadderZeroCapture operatorSpecN :=
    (operatorApproxZero_iff_xiTarget_capture).1 hApprox0
  intro s hsXi ε hε
  rcases hCapXi s hsXi with ⟨N, t, htN, hsEq⟩
  refine ⟨N, criticalLinePoint t, XiFinite_zero_of_mem (operatorSpecN N) htN, ?_⟩
  subst hsEq
  simpa using hε

/-- RH closure from the zero-tolerance operator approximation routed through
the Hurwitz-output surface. -/
theorem mathlibRH_of_operator_approxZero_via_hurwitz
    (hApprox0 : OperatorApproxZeroConvergence) :
    RiemannHypothesis := by
  exact mathlibRH_of_operator_hurwitz_zero_approx
    (operatorHurwitz_of_operatorApproxZero hApprox0)

end

end Gutoe.RiemannOperatorLadder
