import Mathlib
import Gutoe.RiemannCore
import Gutoe.RiemannSelfAdjoint
import Gutoe.RiemannTargetFiniteLadder

namespace Gutoe.RiemannOperatorLadder

open Gutoe.RiemannCore
open Gutoe.RiemannSelfAdjoint
open Gutoe.RiemannTargetFiniteLadder

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
  Module.End.HasEigenvalue (structuralRiemannEnd n) (t : ℂ)

/-- Direct operator-defined spectrum lane at level `N`:
real ordinates whose complex lift lies in the matrix spectrum. -/
def operatorSpecSet (N : ℕ) : Set ℝ :=
  fun t => (t : ℂ) ∈ spectrum ℂ (structuralRiemannEnd (N + 1))

theorem ordinateIsEigenvalue_iff_mem_operatorSpecSet
    (N : ℕ) (t : ℝ) :
    ordinateIsEigenvalue (N + 1) t ↔ t ∈ operatorSpecSet N := by
  unfold ordinateIsEigenvalue operatorSpecSet structuralRiemannEnd
  simpa using
    (Module.End.hasEigenvalue_iff_mem_spectrum
      (f := Matrix.mulVecLin (structuralRiemannMatrixC (N + 1)))
      (μ := (t : ℂ)))

/-- Operator-native finite ladder witness:
each finite level list is exactly the real ordinates detected as operator
eigenvalues for that level. -/
structure OperatorSpecLadder where
  specN : ℕ → Finset ℝ
  eigenvalue_exact :
    ∀ N : ℕ, ∀ t : ℝ, t ∈ specN N ↔ ordinateIsEigenvalue (N + 1) t

/-- Obligation-1 in operator-native form. -/
def OperatorNontrivialCapture (hOp : OperatorSpecLadder) : Prop :=
  RiemannNontrivialLadderZeroCapture hOp.specN

/-- Direct operator-set capture surface (no external reference data). -/
def OperatorSetNontrivialCapture : Prop :=
  ∀ s : ℂ, riemannZeta s = 0 →
    (¬ ∃ n : ℕ, s = -2 * (n + 1)) →
    s ≠ 1 →
    ∃ N : ℕ, ∃ t : ℝ, t ∈ operatorSpecSet N ∧ s = criticalLinePoint t

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
