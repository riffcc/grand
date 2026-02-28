import Mathlib
import Gutoe.TriangulatedConstants
import Gutoe.TriangulatedTermProvenance

namespace Gutoe.TriangulatedClosureUniqueness

open Gutoe.TriangulatedConstants
open Gutoe.TriangulatedTermProvenance

/-!
Constrained-grammar uniqueness closure for triangulated constants.

We enumerate small structural coefficient families and prove that only one
coefficient choice lands inside the frozen tolerance windows.
-/

/-- Three-sign coefficient family `{-1, 0, 1}`. -/
def sign3 : Finset ℤ := {-1, 0, 1}

/-- Binary coefficient family `{0, 1}`. -/
def bit2 : Finset ℕ := {0, 1}

/-- `p` family around `137/10` with lattice-gauge correction step `1/(7*12)`. -/
def pFamilyQ (s : ℤ) : ℚ :=
  (137 / 10 : ℚ) + (s : ℚ) * (1 / (7 * 12 : ℚ))

/-- Good `p` coefficients: frozen-window closure condition. -/
def pGoodSigns : Finset ℤ :=
  sign3.filter (fun s => |pFamilyQ s - pFrozenQ| < (1 / 50000 : ℚ))

theorem p_good_signs_unique :
    pGoodSigns = ({-1} : Finset ℤ) := by
  native_decide

theorem p_family_selected_eq_candidate :
    pFamilyQ (-1) = pCandidateQ := by
  rw [p_candidate_closed_form]
  unfold pFamilyQ
  ring

/-- `kappa` family from three provenance terms with binary include/exclude knobs. -/
def kappaFamilyQ (a b c : ℕ) : ℚ :=
  (60 / 11 : ℚ) *
    ((a : ℚ) * (19 / 3 : ℚ)
      + (b : ℚ) * (1 / 36 : ℚ)
      + (c : ℚ) * (1 / (7 * 13 * 136 : ℚ)))

/-- Candidate tuples `(a,b,c)` with `a,b,c ∈ {0,1}`. -/
def kappaCoeffTuples : Finset (ℕ × (ℕ × ℕ)) :=
  bit2.product (bit2.product bit2)

/-- Good `kappa` tuples: frozen-window closure condition. -/
def kappaGoodTuples : Finset (ℕ × (ℕ × ℕ)) :=
  kappaCoeffTuples.filter (fun t =>
    let a := t.1
    let b := t.2.1
    let c := t.2.2
    |kappaFamilyQ a b c - kappaFrozenQ| < (1 / 50000 : ℚ))

theorem kappa_good_tuples_unique :
    kappaGoodTuples = ({(1, (1, 1))} : Finset (ℕ × (ℕ × ℕ))) := by
  native_decide

theorem kappa_family_selected_eq_candidate :
    kappaFamilyQ 1 1 1 = kappaCandidateQ := by
  rw [kappa_candidate_closed_form]
  unfold kappaFamilyQ
  ring

/-- EW coefficient family:
    baseline `8`, optional `6/13`, and signed finite correction `±1/(7*136)`. -/
def ewFamilyQ (b : ℕ) (s : ℤ) : ℚ :=
  (8 : ℚ) + (b : ℚ) * (6 / 13 : ℚ) + (s : ℚ) * (1 / (7 * 136 : ℚ))

/-- EW tuples `(b,s)` with `b ∈ {0,1}` and `s ∈ {-1,0,1}`. -/
def ewCoeffTuples : Finset (ℕ × ℤ) :=
  bit2.product sign3

/-- Good EW tuples: frozen-window closure condition. -/
def ewGoodTuples : Finset (ℕ × ℤ) :=
  ewCoeffTuples.filter (fun t =>
    let b := t.1
    let s := t.2
    |ewFamilyQ b s - ewCoeffFrozenQ| < (1 / 1000000 : ℚ))

theorem ew_good_tuples_unique :
    ewGoodTuples = ({(1, -1)} : Finset (ℕ × ℤ)) := by
  native_decide

theorem ew_family_selected_eq_candidate :
    ewFamilyQ 1 (-1) = ewCoeffCandidateQ := by
  rw [ew_coeff_candidate_closed_form]
  unfold ewFamilyQ
  ring

/-- Combined uniqueness closure for the constrained operator grammar. -/
theorem constrained_grammar_uniqueness_closure :
    pGoodSigns.card = 1 ∧
    kappaGoodTuples.card = 1 ∧
    ewGoodTuples.card = 1 := by
  rw [p_good_signs_unique, kappa_good_tuples_unique, ew_good_tuples_unique]
  native_decide

/-- Selected coefficients recovered by closure:
    `s_p = -1`, `(a,b,c) = (1,1,1)`, `(b_ew,s_ew) = (1,-1)`. -/
theorem constrained_grammar_selected_coefficients :
    (-1 : ℤ) ∈ pGoodSigns ∧
    (1, (1, 1)) ∈ kappaGoodTuples ∧
    (1, (-1 : ℤ)) ∈ ewGoodTuples := by
  rw [p_good_signs_unique, kappa_good_tuples_unique, ew_good_tuples_unique]
  native_decide

end Gutoe.TriangulatedClosureUniqueness
