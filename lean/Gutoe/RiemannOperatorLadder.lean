import Mathlib
import Gutoe.RiemannCore
import Gutoe.RiemannSelfAdjoint
import Gutoe.RiemannTargetFiniteLadder
import Gutoe.RiemannLimitBridge
import Gutoe.RiemannConvergenceTransfer
import Gutoe.RiemannFinalTarget
import Gutoe.RiemannFiniteXiModel
import Gutoe.RiemannHurwitzKernel

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

/-- The concrete operator spectral ladder at level `N` is nonempty. -/
theorem operatorSpecN_nonempty (N : ℕ) :
    (operatorSpecN N).Nonempty := by
  classical
  have huniv : (Finset.univ : Finset (Fin (N + 1))).Nonempty := Finset.univ_nonempty
  rcases huniv with ⟨i, hi⟩
  refine ⟨operatorEigenvalues N i, ?_⟩
  exact Finset.mem_image.mpr ⟨i, hi, rfl⟩

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
