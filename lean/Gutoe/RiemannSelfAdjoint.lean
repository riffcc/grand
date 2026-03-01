import Mathlib
import Gutoe.FineStructure
import Gutoe.RiemannCore

namespace Gutoe.RiemannSelfAdjoint

open Gutoe.FineStructure

/-- Structural diagonal offset `13/16` (non-void / full Clifford ratio). -/
def timelikeOffsetQ : ℚ := (13 : ℚ) / 16

/-- Structural nearest-neighbor hop `6/11` from grade-2 and `α⁻¹(d=2)=11`. -/
def structuralHopQ : ℚ := (Nat.choose 4 2 : ℚ) / (alphaInverse 2 : ℚ)

theorem timelikeOffsetQ_pos : 0 < timelikeOffsetQ := by
  norm_num [timelikeOffsetQ]

theorem structuralHopQ_eq_six_over_eleven :
    structuralHopQ = (6 : ℚ) / 11 := by
  have hchooseNat : Nat.choose 4 2 = 6 := by
    native_decide
  have hchooseQ : (Nat.choose 4 2 : ℚ) = 6 := by
    exact_mod_cast hchooseNat
  simp [structuralHopQ, alpha_inverse_d2, hchooseQ]

theorem structuralHopQ_pos : 0 < structuralHopQ := by
  rw [structuralHopQ_eq_six_over_eleven]
  norm_num

/-- Structural finite-dimensional RH candidate operator (tridiagonal and symmetric). -/
def structuralRiemannMatrix (n : ℕ) : Matrix (Fin n) (Fin n) ℚ :=
  fun i j =>
    if i = j then
      ((i.1 + 1 : ℚ) + timelikeOffsetQ)
    else if i.1 + 1 = j.1 ∨ j.1 + 1 = i.1 then
      structuralHopQ
    else
      0

/-- Finite self-adjointness proxy in this rational matrix lane. -/
def finiteSelfAdjoint {n : ℕ} (A : Matrix (Fin n) (Fin n) ℚ) : Prop := A.IsSymm

theorem structuralRiemannMatrix_diag (n : ℕ) (i : Fin n) :
    structuralRiemannMatrix n i i = ((i.1 + 1 : ℚ) + timelikeOffsetQ) := by
  simp [structuralRiemannMatrix]

theorem structuralRiemannMatrix_isSymm (n : ℕ) :
    (structuralRiemannMatrix n).IsSymm := by
  refine Matrix.IsSymm.ext ?_
  intro i j
  by_cases h : j = i
  · simp [structuralRiemannMatrix, h]
  · have h' : i ≠ j := by simpa [eq_comm] using h
    simp [structuralRiemannMatrix, h, h', or_comm, or_left_comm, or_assoc]

theorem structuralRiemannMatrix_finiteSelfAdjoint (n : ℕ) :
    finiteSelfAdjoint (structuralRiemannMatrix n) := by
  simpa [finiteSelfAdjoint] using structuralRiemannMatrix_isSymm n

/-- Off-band entries vanish in the structural tridiagonal model. -/
theorem structuralRiemannMatrix_offband_zero
    (n : ℕ) (i j : Fin n)
    (hdiag : i ≠ j)
    (hband : ¬ (i.1 + 1 = j.1 ∨ j.1 + 1 = i.1)) :
    structuralRiemannMatrix n i j = 0 := by
  simp [structuralRiemannMatrix, hdiag, hband]

end Gutoe.RiemannSelfAdjoint
