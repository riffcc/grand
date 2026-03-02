import Mathlib
import Gutoe.FineStructure
import Gutoe.RiemannCore

namespace Gutoe.RiemannSelfAdjoint

open Gutoe.FineStructure

/-- Structural diagonal offset `13/16` (non-void / full Clifford ratio). -/
def timelikeOffsetQ : ℚ := (13 : ℚ) / 16

/-- Structural nearest-neighbor hop `6/11` from grade-2 and `α⁻¹(d=2)=11`. -/
def structuralHopQ : ℚ := (Nat.choose 4 2 : ℚ) / (alphaInverse 2 : ℚ)

/-- Finite-level spectral center candidate coming from reversal symmetry. -/
def structuralCenterQ (n : ℕ) : ℚ := (n + 1 : ℚ) + 2 * timelikeOffsetQ

/-- Parity sign used for alternating-sign conjugation. -/
def paritySignQ (k : ℕ) : ℚ := if Even k then 1 else -1

/-- Parity sign lifted to finite indices. -/
def finParitySignQ {n : ℕ} (i : Fin n) : ℚ := paritySignQ i.1

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

theorem paritySignQ_sq_one (k : ℕ) :
    paritySignQ k * paritySignQ k = 1 := by
  by_cases hk : Even k
  · simp [paritySignQ, hk]
  · simp [paritySignQ, hk]

theorem paritySignQ_succ_mul (k : ℕ) :
    paritySignQ k * paritySignQ (k + 1) = -1 := by
  by_cases hk : Even k
  · have hk1 : ¬ Even (k + 1) := by
      simpa [Nat.even_add_one] using hk
    simp [paritySignQ, hk, hk1]
  · have hk1 : Even (k + 1) := by
      simpa [Nat.even_add_one] using hk
    simp [paritySignQ, hk, hk1]

theorem finParitySignQ_sq_one {n : ℕ} (i : Fin n) :
    finParitySignQ i * finParitySignQ i = 1 := by
  simpa [finParitySignQ] using paritySignQ_sq_one i.1

theorem finParitySignQ_mul_of_succ_eq {n : ℕ} {i j : Fin n}
    (h : i.1 + 1 = j.1) :
    finParitySignQ i * finParitySignQ j = -1 := by
  simpa [finParitySignQ, h] using paritySignQ_succ_mul i.1

theorem finParitySignQ_mul_of_succ_eq_rev {n : ℕ} {i j : Fin n}
    (h : j.1 + 1 = i.1) :
    finParitySignQ i * finParitySignQ j = -1 := by
  simpa [mul_comm, finParitySignQ, h] using paritySignQ_succ_mul j.1

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

theorem rev_band_forward_iff {n : ℕ} (i j : Fin n) :
    ((Fin.rev i).1 + 1 = (Fin.rev j).1) ↔ (j.1 + 1 = i.1) := by
  have hri : (Fin.rev i).1 = n - (i.1 + 1) := by
    simp [Fin.rev]
  have hrj : (Fin.rev j).1 = n - (j.1 + 1) := by
    simp [Fin.rev]
  constructor <;> intro h <;> omega

theorem rev_band_backward_iff {n : ℕ} (i j : Fin n) :
    ((Fin.rev j).1 + 1 = (Fin.rev i).1) ↔ (i.1 + 1 = j.1) := by
  have hri : (Fin.rev i).1 = n - (i.1 + 1) := by
    simp [Fin.rev]
  have hrj : (Fin.rev j).1 = n - (j.1 + 1) := by
    simp [Fin.rev]
  constructor <;> intro h <;> omega

theorem rev_band_iff {n : ℕ} (i j : Fin n) :
    ((Fin.rev i).1 + 1 = (Fin.rev j).1 ∨ (Fin.rev j).1 + 1 = (Fin.rev i).1)
      ↔ (i.1 + 1 = j.1 ∨ j.1 + 1 = i.1) := by
  constructor
  · intro h
    rcases h with h | h
    · exact Or.inr ((rev_band_forward_iff i j).1 h)
    · exact Or.inl ((rev_band_backward_iff i j).1 h)
  · intro h
    rcases h with h | h
    · exact Or.inr ((rev_band_backward_iff i j).2 h)
    · exact Or.inl ((rev_band_forward_iff i j).2 h)

theorem structuralRiemannMatrix_rev_diag_sum (n : ℕ) (i : Fin n) :
    structuralRiemannMatrix n i i +
      structuralRiemannMatrix n (Fin.rev i) (Fin.rev i)
      = structuralCenterQ n := by
  have hri1 : (Fin.rev i).1 + 1 = n - i.1 := by
    have hri : (Fin.rev i).1 = n - (i.1 + 1) := by
      simp [Fin.rev]
    omega
  have hri1Q : (((Fin.rev i).1 + 1 : ℕ) : ℚ) = ((n - i.1 : ℕ) : ℚ) := by
    exact_mod_cast hri1
  have hsumNat : (n - i.1) + i.1 = n := Nat.sub_add_cancel (Nat.le_of_lt i.2)
  have hsumQ : ((n - i.1 : ℕ) : ℚ) + (i.1 : ℚ) = (n : ℚ) := by
    exact_mod_cast hsumNat
  calc
    structuralRiemannMatrix n i i +
      structuralRiemannMatrix n (Fin.rev i) (Fin.rev i)
        = (((i.1 + 1 : ℚ) + timelikeOffsetQ) +
            (((Fin.rev i).1 + 1 : ℚ) + timelikeOffsetQ)) := by
              simp [structuralRiemannMatrix]
    _ = (((i.1 + 1 : ℚ) + timelikeOffsetQ) + (((n - i.1 : ℕ) : ℚ) + timelikeOffsetQ)) := by
          simpa using congrArg (fun x : ℚ => ((i.1 + 1 : ℚ) + timelikeOffsetQ) + (x + timelikeOffsetQ)) hri1Q
    _ = (((n : ℚ) + 1) + 2 * timelikeOffsetQ) := by
          have hsumQ' :
              ((i.1 : ℚ) + ((n - i.1 : ℕ) : ℚ)) = (n : ℚ) := by
            simpa [add_comm] using hsumQ
          nlinarith
    _ = structuralCenterQ n := by
          simp [structuralCenterQ]

theorem structuralRiemannMatrix_isSymm (n : ℕ) :
    (structuralRiemannMatrix n).IsSymm := by
  refine Matrix.IsSymm.ext ?_
  intro i j
  by_cases h : j = i
  · simp [structuralRiemannMatrix, h]
  · have h' : i ≠ j := by simpa [eq_comm] using h
    simp [structuralRiemannMatrix, h, h', or_comm]

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

/-- Entrywise reversal/parity balance identity for the structural operator:
diagonal entries center at `structuralCenterQ n`, while off-diagonal entries
cancel under reversal with alternating-sign conjugation. -/
theorem structuralRiemannMatrix_rev_parity_balance
    (n : ℕ) (i j : Fin n) :
    structuralRiemannMatrix n i j +
      finParitySignQ i * structuralRiemannMatrix n (Fin.rev i) (Fin.rev j) * finParitySignQ j
      = if i = j then structuralCenterQ n else 0 := by
  by_cases hij : i = j
  · subst hij
    calc
      structuralRiemannMatrix n i i +
          finParitySignQ i * structuralRiemannMatrix n (Fin.rev i) (Fin.rev i) * finParitySignQ i
          = structuralRiemannMatrix n i i +
              structuralRiemannMatrix n (Fin.rev i) (Fin.rev i) * (finParitySignQ i * finParitySignQ i) := by
                ring
      _ = structuralRiemannMatrix n i i +
            structuralRiemannMatrix n (Fin.rev i) (Fin.rev i) := by
              simp [finParitySignQ_sq_one]
      _ = structuralCenterQ n := structuralRiemannMatrix_rev_diag_sum n i
      _ = if i = i then structuralCenterQ n else 0 := by simp
  · by_cases hband : (i.1 + 1 = j.1 ∨ j.1 + 1 = i.1)
    · have hAij : structuralRiemannMatrix n i j = structuralHopQ := by
        rw [structuralRiemannMatrix, if_neg hij, if_pos hband]
      have hrevij : Fin.rev i ≠ Fin.rev j := by
        intro h
        exact hij (Fin.rev_injective h)
      have hrevband :
          ((Fin.rev i).1 + 1 = (Fin.rev j).1 ∨ (Fin.rev j).1 + 1 = (Fin.rev i).1) :=
        (rev_band_iff i j).2 hband
      have hArev : structuralRiemannMatrix n (Fin.rev i) (Fin.rev j) = structuralHopQ := by
        rw [structuralRiemannMatrix, if_neg hrevij, if_pos hrevband]
      have hsign : finParitySignQ i * finParitySignQ j = -1 := by
        rcases hband with hsucc | hpred
        · exact finParitySignQ_mul_of_succ_eq hsucc
        · exact finParitySignQ_mul_of_succ_eq_rev hpred
      calc
        structuralRiemannMatrix n i j +
            finParitySignQ i * structuralRiemannMatrix n (Fin.rev i) (Fin.rev j) * finParitySignQ j
            = structuralHopQ + finParitySignQ i * structuralHopQ * finParitySignQ j := by
                simp [hAij, hArev]
        _ = structuralHopQ + structuralHopQ * (finParitySignQ i * finParitySignQ j) := by
              ring
        _ = structuralHopQ + structuralHopQ * (-1) := by simp [hsign]
        _ = if i = j then structuralCenterQ n else 0 := by simp [hij]
    · have hAij0 : structuralRiemannMatrix n i j = 0 :=
        structuralRiemannMatrix_offband_zero n i j hij hband
      have hrevij : Fin.rev i ≠ Fin.rev j := by
        intro h
        exact hij (Fin.rev_injective h)
      have hrevband :
          ¬ ((Fin.rev i).1 + 1 = (Fin.rev j).1 ∨ (Fin.rev j).1 + 1 = (Fin.rev i).1) := by
        intro h
        exact hband ((rev_band_iff i j).1 h)
      have hArev0 : structuralRiemannMatrix n (Fin.rev i) (Fin.rev j) = 0 :=
        structuralRiemannMatrix_offband_zero n (Fin.rev i) (Fin.rev j) hrevij hrevband
      calc
        structuralRiemannMatrix n i j +
            finParitySignQ i * structuralRiemannMatrix n (Fin.rev i) (Fin.rev j) * finParitySignQ j
            = 0 := by simp [hAij0, hArev0]
        _ = if i = j then structuralCenterQ n else 0 := by simp [hij]

end Gutoe.RiemannSelfAdjoint
