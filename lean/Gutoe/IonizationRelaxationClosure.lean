import Mathlib
import Gutoe.DimensionalStructure
import Gutoe.Z3Uniqueness

namespace Gutoe.IonizationRelaxationClosure

open Gutoe.DimensionalStructure
open Gutoe.Z3Uniqueness

/-!
Lean closure for Koopmans-relaxation coefficients used in atomic IE correction.
-/

def dNat : ℕ := 2 ^ 4
def grade2Count : ℕ := Nat.choose 4 2
def su2Count : ℕ := magneticTriplet.card
def z3FixedCount : ℕ := (grade1_4d.filter (fun s => z3_4d s = s)).card

theorem d_nat_eq_16 : dNat = 16 := by native_decide
theorem grade2_count_eq_6 : grade2Count = 6 := by native_decide
theorem su2_count_eq_3 : su2Count = 3 := by native_decide
theorem z3_fixed_count_eq_1 : z3FixedCount = 1 := by
  simpa [z3FixedCount] using z3_grade1_fixed_count

/-- Nonmetal Koopmans relaxation gain. -/
def nonmetalRelaxGainQ : ℚ := (z3FixedCount : ℚ) / ((Nat.choose 4 1 : ℕ) : ℚ)

/-- Lanthanide spread gain. -/
def lanthanideSpreadGainQ : ℚ :=
  ((dNat + su2Count ^ 2 : ℕ) : ℚ) / (grade2Count : ℚ)

theorem ie_relaxation_coefficients_closed_form :
    nonmetalRelaxGainQ = (1 / 4 : ℚ) ∧
    lanthanideSpreadGainQ = (25 / 6 : ℚ) := by
  constructor
  · unfold nonmetalRelaxGainQ z3FixedCount
    native_decide
  · unfold lanthanideSpreadGainQ dNat su2Count grade2Count
    native_decide

/-- Finite numerator elimination for denominator-4 relaxation gain. -/
def numCandidates4 : Finset ℕ := Finset.range 5
def goodNums4 : Finset ℕ :=
  numCandidates4.filter (fun n => ((n : ℚ) / (4 : ℚ)) = nonmetalRelaxGainQ)

theorem nonmetal_relax_num_unique : goodNums4 = ({1} : Finset ℕ) := by
  unfold goodNums4 numCandidates4 nonmetalRelaxGainQ z3FixedCount
  native_decide

/-- Finite numerator elimination for denominator-6 lanthanide spread gain. -/
def numCandidates6 : Finset ℕ := Finset.range 31
def goodNums6 : Finset ℕ :=
  numCandidates6.filter (fun n => ((n : ℚ) / (6 : ℚ)) = lanthanideSpreadGainQ)

theorem lanthanide_spread_num_unique : goodNums6 = ({25} : Finset ℕ) := by
  unfold goodNums6 numCandidates6 lanthanideSpreadGainQ dNat su2Count grade2Count
  native_decide

end Gutoe.IonizationRelaxationClosure
