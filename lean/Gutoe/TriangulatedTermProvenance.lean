import Mathlib
import Gutoe.FineStructure
import Gutoe.Z3Uniqueness
import Gutoe.DarkMatterSector
import Gutoe.TriangulatedConstants

namespace Gutoe.TriangulatedTermProvenance

open Gutoe.FineStructure
open Gutoe.Z3Uniqueness
open Gutoe.DarkMatterSector
open Gutoe.TriangulatedConstants

/-!
Term-provenance theorems for triangulated constants.

Goal: move from "familiar numbers" to explicit structural origins:

- `7` from `|grade₂| + 1 = C(4,2) + 1`
- `19/3` from `(d + |SU(2)|)/|SU(2)| = (16 + 3)/3`
- correction factors from finite Cl(1,3) closure terms
-/

/-- Structural lattice shift reused across PMNS/triangulation lanes:
    `|grade₂| + 1 = C(4,2) + 1`. -/
def latticeShiftQ : ℚ := ((Nat.choose 4 2 + 1 : ℕ) : ℚ)

theorem lattice_shift_eq_7 : latticeShiftQ = 7 := by
  native_decide

/-- Weak/color balance ratio from Clifford dimension and SU(2) orbit size:
    `(d + |SU(2)|)/|SU(2)| = (16+3)/3`. -/
def weakColorBalanceQ : ℚ :=
  (((2 ^ 4) + magneticTriplet.card : ℕ) : ℚ) / (magneticTriplet.card : ℚ)

theorem weak_color_balance_eq_19_over_3 :
    weakColorBalanceQ = 19 / 3 := by
  unfold weakColorBalanceQ
  rw [su2_dim]
  native_decide

/-- Grade-2 self-coupling suppression from the bivector count:
    `1/|grade₂|² = 1/36`. -/
def bivectorSelfCouplingQ : ℚ :=
  1 / (((Nat.choose 4 2 : ℕ) : ℚ) ^ 2)

theorem bivector_self_coupling_eq_1_over_36 :
    bivectorSelfCouplingQ = 1 / 36 := by
  unfold bivectorSelfCouplingQ
  native_decide

/-- Finite-closure coupling from the product
    `( |grade₂|+1 ) * (d-|SU(2)|) * T(16)` in the denominator. -/
def finiteClosureCouplingQ : ℚ :=
  1 / (latticeShiftQ
      * (((2 ^ 4) - magneticTriplet.card : ℕ) : ℚ)
      * (triangularNumber (2 ^ 4) : ℚ))

theorem finite_closure_coupling_eq_1_over_7_13_136 :
    finiteClosureCouplingQ = 1 / (7 * 13 * 136 : ℚ) := by
  unfold finiteClosureCouplingQ latticeShiftQ
  have hs : magneticTriplet.card = 3 := su2_dim
  have hpow : (2 ^ 4 : ℕ) = 16 := by native_decide
  have hT : triangularNumber (2 ^ 4) = 136 := by
    rw [hpow, T16_eq_136]
  rw [hs, hT]
  native_decide

/-- Candidate `kappa` decomposition into provenance terms. -/
def kappaProvenanceDecompQ : ℚ :=
  geometricDarkToVisibleRatio *
    (weakColorBalanceQ + bivectorSelfCouplingQ + finiteClosureCouplingQ)

theorem kappa_provenance_decomp_eq_candidate :
    kappaProvenanceDecompQ = kappaCandidateQ := by
  unfold kappaProvenanceDecompQ kappaCandidateQ weakColorBalanceQ bivectorSelfCouplingQ
  unfold finiteClosureCouplingQ latticeShiftQ
  ring

theorem kappa_provenance_closed_form :
    kappaProvenanceDecompQ = (60 / 11 : ℚ) * ((19 / 3 : ℚ) + 1 / 36 + 1 / (7 * 13 * 136 : ℚ)) := by
  rw [kappa_provenance_decomp_eq_candidate]
  exact kappa_candidate_closed_form

/-- EW uplift provenance term:
    `|grade₂|/(d-|SU(2)|) - 1/((|grade₂|+1)T(16)) = 6/13 - 1/(7*136)`. -/
def ewUpliftProvenanceQ : ℚ :=
  ((Nat.choose 4 2 : ℕ) : ℚ) / (((2 ^ 4) - magneticTriplet.card : ℕ) : ℚ)
    - 1 / (latticeShiftQ * (triangularNumber (2 ^ 4) : ℚ))

theorem ew_uplift_provenance_closed_form :
    ewUpliftProvenanceQ = 6 / 13 - 1 / (7 * 136 : ℚ) := by
  unfold ewUpliftProvenanceQ latticeShiftQ
  have hs : magneticTriplet.card = 3 := su2_dim
  have hpow : (2 ^ 4 : ℕ) = 16 := by native_decide
  have hT : triangularNumber (2 ^ 4) = 136 := by
    rw [hpow, T16_eq_136]
  rw [hs, hT]
  native_decide

/-- EW candidate coefficient is exactly `d/2 +` the provenance uplift. -/
theorem ew_candidate_from_provenance :
    ewCoeffCandidateQ = ((2 ^ 4 : ℚ) / 2) + ewUpliftProvenanceQ := by
  unfold ewCoeffCandidateQ ewUpliftProvenanceQ latticeShiftQ
  ring

/-- `p` candidate provenance: structural baseline `α⁻¹/10` plus a finite-lattice
    subtraction from `(|grade₂|+1)*N_gauge`. -/
theorem p_candidate_from_baseline_and_finite_subtraction :
    pCandidateQ =
      (alphaInverse 4 : ℚ) / ((Nat.choose 4 1 + Nat.choose 4 2 : ℕ) : ℚ)
      - 1 / (latticeShiftQ * totalGaugeGeneratorsQ) := by
  unfold pCandidateQ latticeShiftQ
  ring

end Gutoe.TriangulatedTermProvenance
