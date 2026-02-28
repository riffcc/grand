import Mathlib
import Gutoe.DimensionalStructure
import Gutoe.Z3Uniqueness

namespace Gutoe.ThermalEntropyClosure

open Gutoe.DimensionalStructure
open Gutoe.Z3Uniqueness

/-!
Lean closure for thermal entropy transduction coefficients used in chemical thermodynamics.
-/

def grade1Count : ℕ := grade1_4d.card
def grade2Count : ℕ := Nat.choose 4 2
def su2Count : ℕ := magneticTriplet.card
def z3FixedCount : ℕ := (grade1_4d.filter (fun s => z3_4d s = s)).card

theorem grade1_count_eq_4 : grade1Count = 4 := by native_decide
theorem grade2_count_eq_6 : grade2Count = 6 := by native_decide
theorem su2_count_eq_3 : su2Count = 3 := by native_decide
theorem z3_fixed_count_eq_1 : z3FixedCount = 1 := by
  simpa [z3FixedCount] using z3_grade1_fixed_count

/-- d-band fusion entropy gain coefficient. -/
def dFusionGainQ : ℚ := ((grade2Count - z3FixedCount : ℕ) : ℚ) / (grade2Count : ℚ)

/-- d-band vapor entropy gain coefficient. -/
def dVaporGainQ : ℚ := (su2Count : ℚ) / ((grade1Count + z3FixedCount : ℕ) : ℚ)

/-- Covalent-network entropy gain coefficient. -/
def covalentGainQ : ℚ := (grade1Count : ℚ) / ((grade1Count + z3FixedCount : ℕ) : ℚ)

/-- Metalloid entropy penalty coefficient. -/
def metalloidPenaltyQ : ℚ := ((grade2Count + z3FixedCount : ℕ) : ℚ) / (grade1Count : ℚ)

/-- Molecular entropy penalty coefficient. -/
def molecularPenaltyQ : ℚ := ((grade1Count + z3FixedCount : ℕ) : ℚ) / (grade1Count : ℚ)

theorem thermal_entropy_coefficients_closed_form :
    dFusionGainQ = (5 / 6 : ℚ) ∧
    dVaporGainQ = (3 / 5 : ℚ) ∧
    covalentGainQ = (4 / 5 : ℚ) ∧
    metalloidPenaltyQ = (7 / 4 : ℚ) ∧
    molecularPenaltyQ = (5 / 4 : ℚ) := by
  constructor
  · unfold dFusionGainQ grade2Count z3FixedCount
    native_decide
  constructor
  · unfold dVaporGainQ su2Count grade1Count z3FixedCount
    native_decide
  constructor
  · unfold covalentGainQ grade1Count z3FixedCount
    native_decide
  constructor
  · unfold metalloidPenaltyQ grade1Count grade2Count z3FixedCount
    native_decide
  · unfold molecularPenaltyQ grade1Count z3FixedCount
    native_decide

/-- Finite numerator search for a fixed denominator and target rational. -/
def goodNums (den maxNum : ℕ) (target : ℚ) : Finset ℕ :=
  (Finset.range (maxNum + 1)).filter (fun n => ((n : ℚ) / (den : ℚ)) = target)

theorem d_fusion_num_unique : goodNums 6 12 dFusionGainQ = ({5} : Finset ℕ) := by
  unfold goodNums dFusionGainQ grade2Count z3FixedCount
  native_decide

theorem d_vapor_num_unique : goodNums 5 12 dVaporGainQ = ({3} : Finset ℕ) := by
  unfold goodNums dVaporGainQ su2Count grade1Count z3FixedCount
  native_decide

theorem covalent_num_unique : goodNums 5 12 covalentGainQ = ({4} : Finset ℕ) := by
  unfold goodNums covalentGainQ grade1Count z3FixedCount
  native_decide

theorem metalloid_num_unique : goodNums 4 12 metalloidPenaltyQ = ({7} : Finset ℕ) := by
  unfold goodNums metalloidPenaltyQ grade1Count grade2Count z3FixedCount
  native_decide

theorem molecular_num_unique : goodNums 4 12 molecularPenaltyQ = ({5} : Finset ℕ) := by
  unfold goodNums molecularPenaltyQ grade1Count z3FixedCount
  native_decide

end Gutoe.ThermalEntropyClosure
