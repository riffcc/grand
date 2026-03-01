import Mathlib
import Gutoe.CliffordStructure
import Gutoe.FineStructure

/-!
GUTOE — Even-Subalgebra Suppression for Origin Fan-In

This module formalizes the structural `1/2` suppression candidate from the
Cl(1,3) grade split:

- even grades: 0,2,4 with counts `1 + 6 + 1 = 8`
- odd grades: 1,3 with counts `4 + 4 = 8`
- total: `16`

So the even-grade filter is exactly `8/16 = 1/2`.
-/

namespace Gutoe.EvenSubalgebraSuppression

open Gutoe.CliffordStructure
open Gutoe.FineStructure

/-- Grade-`k` basis count in Cl(1,3) as `C(4,k)`. -/
def gradeCount (k : ℕ) : ℕ := Nat.choose 4 k

/-- Even-grade basis count: grades `0,2,4`. -/
def evenGradeCount : ℕ := gradeCount 0 + gradeCount 2 + gradeCount 4

/-- Odd-grade basis count: grades `1,3`. -/
def oddGradeCount : ℕ := gradeCount 1 + gradeCount 3

/-- Total basis count from even+odd split. -/
def totalGradeCount : ℕ := evenGradeCount + oddGradeCount

/-- Even-grade count is exactly `8` for Cl(1,3). -/
theorem even_grade_count_eq_8 : evenGradeCount = 8 := by
  unfold evenGradeCount gradeCount
  native_decide

/-- Odd-grade count is exactly `8` for Cl(1,3). -/
theorem odd_grade_count_eq_8 : oddGradeCount = 8 := by
  unfold oddGradeCount gradeCount
  native_decide

/-- Even+odd grade counts recover the full `16`-state basis. -/
theorem total_grade_count_eq_16 : totalGradeCount = 16 := by
  unfold totalGradeCount
  rw [even_grade_count_eq_8, odd_grade_count_eq_8]

/-- The even-grade suppression factor from Cl(1,3) grading. -/
def evenSuppressionQ : ℚ := (evenGradeCount : ℚ) / (totalGradeCount : ℚ)

/-- Exact structural closure: `evenSuppressionQ = 1/2`. -/
theorem even_suppression_eq_one_half : evenSuppressionQ = (1 : ℚ) / 2 := by
  unfold evenSuppressionQ
  rw [even_grade_count_eq_8, total_grade_count_eq_16]
  norm_num

/-- Canonical fan-in gain with `branching=2`, `merge=1`, and even filter. -/
def geffCanonicalQ : ℚ := (2 : ℚ) * (1 : ℚ) * evenSuppressionQ

/-- Under the canonical lane, even filtering pins gain exactly to `1`. -/
theorem geff_canonical_eq_one : geffCanonicalQ = 1 := by
  unfold geffCanonicalQ
  rw [even_suppression_eq_one_half]
  norm_num

/-- Z3 branching factor for the origin fan-in lane. -/
def branchingZ3Q : ℚ := 3

/-- Void merge fraction in the origin fan-in lane. -/
def mergeVoidQ : ℚ := (3 : ℚ) / 16

/-- Structural split candidate: visible vectors over bivectors = `4/6 = 2/3`. -/
def etaGrade1OverGrade2Q : ℚ := (2 : ℚ) / 3

/-- Structural split candidate: full basis over bivectors = `16/6 = 8/3`. -/
def infraBasisOverGrade2Q : ℚ := (8 : ℚ) / 3

/-- Product of the two split factors is exactly `16/9`. -/
theorem eta_infra_split_product_eq_16_over_9 :
    etaGrade1OverGrade2Q * infraBasisOverGrade2Q = (16 : ℚ) / 9 := by
  unfold etaGrade1OverGrade2Q infraBasisOverGrade2Q
  norm_num

/-- Full structural gain from `3 * (3/16) * (2/3) * (8/3)`. -/
def geffZ3VoidSplitQ : ℚ :=
  branchingZ3Q * mergeVoidQ * etaGrade1OverGrade2Q * infraBasisOverGrade2Q

/-- The structural split closes exactly at the knife-edge (`G_eff = 1`). -/
theorem geff_z3_void_split_eq_one : geffZ3VoidSplitQ = 1 := by
  unfold geffZ3VoidSplitQ branchingZ3Q mergeVoidQ
  calc
    (3 : ℚ) * ((3 : ℚ) / 16) * etaGrade1OverGrade2Q * infraBasisOverGrade2Q
        = ((9 : ℚ) / 16) * (etaGrade1OverGrade2Q * infraBasisOverGrade2Q) := by ring
    _ = ((9 : ℚ) / 16) * ((16 : ℚ) / 9) := by rw [eta_infra_split_product_eq_16_over_9]
    _ = 1 := by norm_num

/-- Rational form of the uncapped measured gain `1.9992`. -/
def uncappedGainQ : ℚ := (2499 : ℚ) / 1250

/-- Suppression needed to reduce `1.9992` to unit gain. -/
def suppressionForUnitGainQ : ℚ := (1 : ℚ) / uncappedGainQ

/-- Exact required suppression for `1.9992 → 1`: `1250/2499`. -/
theorem suppression_for_unit_gain_exact :
    suppressionForUnitGainQ = (1250 : ℚ) / 2499 := by
  unfold suppressionForUnitGainQ uncappedGainQ
  norm_num

/-- The required suppression sits just above `1/2` by `1/4998`. -/
theorem suppression_for_unit_gain_half_offset :
    suppressionForUnitGainQ - (1 : ℚ) / 2 = (1 : ℚ) / 4998 := by
  rw [suppression_for_unit_gain_exact]
  norm_num

end Gutoe.EvenSubalgebraSuppression
