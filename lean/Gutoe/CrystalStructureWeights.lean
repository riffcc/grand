import Mathlib
import Gutoe.DimensionalStructure
import Gutoe.Z3Uniqueness

namespace Gutoe.CrystalStructureWeights

open Gutoe.DimensionalStructure
open Gutoe.Z3Uniqueness

/-!
Zero-parameter crystal-structure blend weights.

We force the family-vs-crystal blending numerators from shared Cl(1,3)
combinatorics:
- family lane numerator = `grade1_4d.card`
- crystal lane numerator = `magneticTriplet.card`

No fitted decimals: weights are exact rationals from finite cardinalities.
-/

/-- Family-lane numerator from grade-1 state count. -/
def familyCount : ℕ := grade1_4d.card

/-- Crystal-lane numerator from spatial bivector triplet count. -/
def crystalCount : ℕ := magneticTriplet.card

/-- Total numerator budget for convex blending. -/
def totalCount : ℕ := familyCount + crystalCount

/-- Family blend weight as a rational. -/
def familyWeightQ : ℚ := (familyCount : ℚ) / (totalCount : ℚ)

/-- Crystal blend weight as a rational. -/
def crystalWeightQ : ℚ := (crystalCount : ℚ) / (totalCount : ℚ)

/-- Shared primitive cardinalities are exactly 4 and 3. -/
theorem family_count_eq_four : familyCount = 4 := by native_decide

theorem crystal_count_eq_three : crystalCount = 3 := by native_decide

theorem total_count_eq_seven : totalCount = 7 := by native_decide

/-- Weights are exactly 4/7 and 3/7. -/
theorem blend_weights_exact :
    familyWeightQ = (4 / 7 : ℚ) ∧ crystalWeightQ = (3 / 7 : ℚ) := by
  constructor
  · unfold familyWeightQ totalCount familyCount crystalCount
    native_decide
  · unfold crystalWeightQ totalCount familyCount crystalCount
    native_decide

/-- Candidate numerator set `0..totalCount`. -/
def numeratorCandidates : Finset ℕ := Finset.range (totalCount + 1)

/-- Good numerator pairs satisfy both:
    (A) convex budget `nF + nC = totalCount`
    (B) Clifford ratio lock `nF * crystalCount = nC * familyCount`. -/
def goodNumeratorPairs : Finset (ℕ × ℕ) :=
  (numeratorCandidates.product numeratorCandidates).filter (fun t =>
    t.1 + t.2 = totalCount ∧ t.1 * crystalCount = t.2 * familyCount)

/-- Finite elimination: only `(4,3)` survives the dual constraints. -/
theorem good_numerator_pairs_unique :
    goodNumeratorPairs = ({(4, 3)} : Finset (ℕ × ℕ)) := by
  unfold goodNumeratorPairs numeratorCandidates totalCount familyCount crystalCount
  native_decide

/-- Uniqueness consequence in proposition form. -/
theorem blend_numerators_unique :
    ∀ nF nC : ℕ,
      (nF, nC) ∈ goodNumeratorPairs → nF = 4 ∧ nC = 3 := by
  intro nF nC hmem
  have hset : goodNumeratorPairs = ({(4, 3)} : Finset (ℕ × ℕ)) :=
    good_numerator_pairs_unique
  rw [hset] at hmem
  simpa using hmem

end Gutoe.CrystalStructureWeights
