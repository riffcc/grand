import Mathlib
import Gutoe.DimensionalStructure
import Gutoe.Z3Uniqueness

namespace Gutoe.CrystalCoefficientClosure

open Gutoe.DimensionalStructure
open Gutoe.Z3Uniqueness

/-!
Finite-elimination closure for crystal-channel coefficients.

All denominators are structural counts from Cl(1,3):
- d = 2^4 = 16
- |grade1| = 4
- |grade2| = C(4,2) = 6
- |SU(2)| = |magneticTriplet| = 3
- Z3 grade-1 fixed count = 1

Coefficients are represented as exact rationals and recovered by finite
numerator elimination over each structural denominator.
-/

def dNat : ℕ := 2 ^ 4
def grade1Count : ℕ := grade1_4d.card
def grade2Count : ℕ := Nat.choose 4 2
def su2Count : ℕ := magneticTriplet.card
def z3FixedCount : ℕ := (grade1_4d.filter (fun s => z3_4d s = s)).card

theorem d_nat_eq_16 : dNat = 16 := by native_decide
theorem grade1_count_eq_4 : grade1Count = 4 := by native_decide
theorem grade2_count_eq_6 : grade2Count = 6 := by native_decide
theorem su2_count_eq_3 : su2Count = 3 := by native_decide
theorem z3_fixed_count_eq_1 : z3FixedCount = 1 := by
  simpa [z3FixedCount] using z3_grade1_fixed_count

/-- Transition corridor blend weights (d-corridor, valence-corridor). -/
def corridorDen : ℕ := dNat + grade1Count

def corridorDWeightQ : ℚ := ((dNat - su2Count : ℕ) : ℚ) / (corridorDen : ℚ)
def corridorVWeightQ : ℚ := ((grade2Count + z3FixedCount : ℕ) : ℚ) / (corridorDen : ℚ)

theorem corridor_weights_closed_form :
    corridorDWeightQ = (13 / 20 : ℚ) ∧ corridorVWeightQ = (7 / 20 : ℚ) := by
  constructor
  · unfold corridorDWeightQ corridorDen dNat grade1Count su2Count
    native_decide
  · unfold corridorVWeightQ corridorDen dNat grade1Count grade2Count z3FixedCount
    native_decide

/-- Structural denominators for pack/radius coefficients. -/
def transitionDen : ℕ := grade2Count
def postTransitionDen : ℕ := dNat - grade1Count
def lanthanidePackDen : ℕ := dNat - grade1Count - (z3FixedCount + z3FixedCount)
def actinideDen : ℕ := dNat + su2Count ^ 2
def lanthanideRadiusDen : ℕ := 2 * actinideDen

/-- Exact coefficient rationals used by the crystal channel. -/
def transitionPackGainQ : ℚ := (su2Count : ℚ) / (transitionDen : ℚ)
def postTransitionPackGainQ : ℚ := ((grade1Count + z3FixedCount : ℕ) : ℚ) / (postTransitionDen : ℚ)
def lanthanidePackGainQ : ℚ := (su2Count : ℚ) / (lanthanidePackDen : ℚ)
def actinidePackGainQ : ℚ := (su2Count : ℚ) / (actinideDen : ℚ)
def lanthanideRadiusGainQ : ℚ := ((su2Count ^ 2 : ℕ) : ℚ) / (lanthanideRadiusDen : ℚ)
def actinideRadiusGainQ : ℚ := ((z3FixedCount + 1 : ℕ) : ℚ) / (actinideDen : ℚ)

theorem crystal_coefficients_closed_form :
    transitionPackGainQ = (1 / 2 : ℚ) ∧
    postTransitionPackGainQ = (5 / 12 : ℚ) ∧
    lanthanidePackGainQ = (3 / 10 : ℚ) ∧
    actinidePackGainQ = (3 / 25 : ℚ) ∧
    lanthanideRadiusGainQ = (9 / 50 : ℚ) ∧
    actinideRadiusGainQ = (2 / 25 : ℚ) := by
  constructor
  · unfold transitionPackGainQ transitionDen su2Count grade2Count
    native_decide
  constructor
  · unfold postTransitionPackGainQ postTransitionDen dNat grade1Count z3FixedCount
    native_decide
  constructor
  · unfold lanthanidePackGainQ lanthanidePackDen dNat grade1Count z3FixedCount su2Count
    native_decide
  constructor
  · unfold actinidePackGainQ actinideDen dNat su2Count
    native_decide
  constructor
  · unfold lanthanideRadiusGainQ lanthanideRadiusDen actinideDen dNat su2Count
    native_decide
  · unfold actinideRadiusGainQ actinideDen dNat z3FixedCount su2Count
    native_decide

/-- Numerator candidate family for a fixed denominator. -/
def numeratorCandidates (den : ℕ) : Finset ℕ := Finset.range (den + 1)

def goodNums (den : ℕ) (target : ℚ) : Finset ℕ :=
  (numeratorCandidates den).filter (fun n => ((n : ℚ) / (den : ℚ)) = target)

/-- Finite elimination: each structural target has exactly one numerator. -/
theorem transition_num_unique :
    goodNums transitionDen transitionPackGainQ = ({3} : Finset ℕ) := by
  unfold goodNums numeratorCandidates transitionDen transitionPackGainQ su2Count grade2Count
  native_decide

theorem post_transition_num_unique :
    goodNums postTransitionDen postTransitionPackGainQ = ({5} : Finset ℕ) := by
  unfold goodNums numeratorCandidates postTransitionDen postTransitionPackGainQ
  unfold dNat grade1Count z3FixedCount
  native_decide

theorem lanthanide_pack_num_unique :
    goodNums lanthanidePackDen lanthanidePackGainQ = ({3} : Finset ℕ) := by
  unfold goodNums numeratorCandidates lanthanidePackDen lanthanidePackGainQ
  unfold dNat grade1Count z3FixedCount su2Count
  native_decide

theorem actinide_pack_num_unique :
    goodNums actinideDen actinidePackGainQ = ({3} : Finset ℕ) := by
  unfold goodNums numeratorCandidates actinideDen actinidePackGainQ
  unfold dNat su2Count
  native_decide

theorem lanthanide_radius_num_unique :
    goodNums lanthanideRadiusDen lanthanideRadiusGainQ = ({9} : Finset ℕ) := by
  unfold goodNums numeratorCandidates lanthanideRadiusDen lanthanideRadiusGainQ
  unfold actinideDen dNat su2Count
  native_decide

theorem actinide_radius_num_unique :
    goodNums actinideDen actinideRadiusGainQ = ({2} : Finset ℕ) := by
  unfold goodNums numeratorCandidates actinideDen actinideRadiusGainQ
  unfold dNat su2Count z3FixedCount
  native_decide

end Gutoe.CrystalCoefficientClosure
