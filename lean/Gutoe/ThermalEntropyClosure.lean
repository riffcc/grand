import Mathlib
import Gutoe.DimensionalStructure
import Gutoe.Z3Uniqueness
import Gutoe.GaugeConstants

namespace Gutoe.ThermalEntropyClosure

open Gutoe.DimensionalStructure
open Gutoe.Z3Uniqueness
open Gutoe.GaugeConstants

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

/-- Total SM gauge-generator count `8 + 3 + 1 = 12`. -/
def gaugeGeneratorCount : ℕ := (3^2 - 1) + (2^2 - 1) + 1

/-- Spatial grade-1 count from grade-1 minus the unique fixed timelike slot. -/
def spatialVectorCount : ℕ := grade1Count - z3FixedCount

/-- Spatial grade-2 count from total grade-2 minus mixed timelike bivectors. -/
def spatialBivectorCount : ℕ := grade2Count - spatialVectorCount

/-- Odd-parity structural basis proxy: `3 vectors + 3 bivectors + 1 fixed slot = 7`. -/
def oddParityBasisCount : ℕ := spatialVectorCount + spatialBivectorCount + z3FixedCount

/-- Refractory suppression ratio from gauge/odd-parity structural balance. -/
def refractorySuppressionRatioQ : ℚ := (gaugeGeneratorCount : ℚ) / (oddParityBasisCount : ℚ)

theorem gauge_generator_count_eq_12 : gaugeGeneratorCount = 12 := by
  simpa [gaugeGeneratorCount] using total_gauge_bosons

theorem spatial_vector_count_eq_3 : spatialVectorCount = 3 := by
  unfold spatialVectorCount grade1Count z3FixedCount
  native_decide

theorem spatial_bivector_count_eq_3 : spatialBivectorCount = 3 := by
  unfold spatialBivectorCount grade2Count spatialVectorCount grade1Count z3FixedCount
  native_decide

theorem odd_parity_basis_count_eq_7 : oddParityBasisCount = 7 := by
  unfold oddParityBasisCount
  rw [spatial_vector_count_eq_3, spatial_bivector_count_eq_3, z3_fixed_count_eq_1]

theorem refractory_suppression_ratio_eq_12_over_7 :
    refractorySuppressionRatioQ = (12 / 7 : ℚ) := by
  unfold refractorySuppressionRatioQ
  rw [gauge_generator_count_eq_12, odd_parity_basis_count_eq_7]
  norm_num

end Gutoe.ThermalEntropyClosure
