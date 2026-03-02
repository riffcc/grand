import Mathlib
import Gutoe.RiemannCore
import Gutoe.RiemannSelfAdjoint
import Gutoe.RiemannTargetFiniteLadder
import Gutoe.RiemannConvergenceTransfer
import Gutoe.RiemannFinalTarget

namespace Gutoe.RiemannOperatorLadder

open Gutoe.RiemannCore
open Gutoe.RiemannSelfAdjoint
open Gutoe.RiemannTargetFiniteLadder
open Gutoe.RiemannConvergenceTransfer
open Gutoe.RiemannFinalTarget

noncomputable section

/-- Complex-lifted structural matrix lane for spectral statements. -/
def structuralRiemannMatrixC (n : ℕ) : Matrix (Fin n) (Fin n) ℂ :=
  fun i j => (structuralRiemannMatrix n i j : ℂ)

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

end

end Gutoe.RiemannOperatorLadder
