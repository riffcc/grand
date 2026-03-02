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

/-- Ordered Weyl center-gap implies the permutation-invariant center-gap
via the identity permutation. -/
theorem operatorCenterGapPermutationInvariant_of_weylCenterGap
    (hW : OperatorWeylCenterGap) :
    OperatorCenterGapPermutationInvariant := by
  intro M
  refine ⟨Equiv.refl _, ?_⟩
  intro i
  simpa using hW M i

/-- Structural center on the real axis for the `k`-th Gershgorin lane. -/
def operatorCenterAt (k : ℕ) : ℝ := (k : ℝ) + (29 : ℝ) / 16

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

/-- Sturm-route contract:
eigenvalue and center counting functions agree at every threshold. This is the
counting-surface targeted to discharge permutation-invariant center-gap without
index-order assumptions. -/
def OperatorSturmCountContract : Prop :=
  ∀ M : ℕ, ∀ x : ℝ, operatorEigenvalueCountLE M x = operatorCenterCountLE M x

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

/-- Sturm-route counting consequence:
if eigenvalue and center counting functions agree at every threshold, then the
eigenvalue counting function also increases by at most `3` across the structural
`±12/11` window around each center. -/
theorem operatorEigenvalueCountLE_window_sub_le_three_of_sturm
    (hS : OperatorSturmCountContract) (M : ℕ) (k : Fin (M + 1)) :
    operatorEigenvalueCountLE M (operatorCenterAt k.1 + (12 : ℝ) / 11) -
      operatorEigenvalueCountLE M (operatorCenterAt k.1 - (12 : ℝ) / 11) ≤ 3 := by
  have hUpper :
      operatorEigenvalueCountLE M (operatorCenterAt k.1 + (12 : ℝ) / 11) =
        operatorCenterCountLE M (operatorCenterAt k.1 + (12 : ℝ) / 11) := by
    exact hS M (operatorCenterAt k.1 + (12 : ℝ) / 11)
  have hLower :
      operatorEigenvalueCountLE M (operatorCenterAt k.1 - (12 : ℝ) / 11) =
        operatorCenterCountLE M (operatorCenterAt k.1 - (12 : ℝ) / 11) := by
    exact hS M (operatorCenterAt k.1 - (12 : ℝ) / 11)
  rw [hUpper, hLower]
  exact operatorCenterCountLE_window_sub_le_three M k

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

/-- One-assumption Weyl reduction:
the Weyl center-gap contract directly yields the finite-level uniform inverse-square bound. -/
theorem exists_uniform_bound_sum_one_div_normSq_operatorSpecN_of_weylCenterGap
    (hW : OperatorWeylCenterGap) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ M : ℕ,
      Finset.sum (operatorSpecN M) (fun t => (1 : ℝ) / (‖criticalLinePoint t‖ ^ (2 : ℕ))) ≤ C := by
  exact exists_uniform_bound_sum_one_div_normSq_operatorSpecN_of_centerGapPermutationInvariant
    (operatorCenterGapPermutationInvariant_of_weylCenterGap hW)

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

/-- Named three-obligation RH surface on the concrete operator lane.
1) permutation-invariant center-gap (summability geometry),
2) Hurwitz kernel, and
3) local-uniform convergence of the operator finite ladder to `XiTarget`. -/
def OperatorRHThreeConditionalObligations : Prop :=
  OperatorCenterGapPermutationInvariant ∧
  OperatorHurwitzKernel ∧
  OperatorXiFiniteLocallyUniformConvergence

/-- RH closure from the named three-obligation surface. -/
theorem mathlibRH_of_operator_three_conditional_obligations
    (h3 : OperatorRHThreeConditionalObligations) :
    RiemannHypothesis := by
  rcases h3 with ⟨_hGap, hKernel, hConv⟩
  exact mathlibRH_of_operator_hurwitzKernel_and_locallyUniform hKernel hConv

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
