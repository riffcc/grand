import Mathlib
import Gutoe.FineStructure
import Gutoe.Z3Uniqueness
import Gutoe.DarkMatterSector

namespace Gutoe.TriangulatedConstants

open Gutoe.FineStructure
open Gutoe.Z3Uniqueness
open Gutoe.DarkMatterSector

/-!
Formal candidate derivations for the three triangulated constants.

These are Cl(1,3)-count formulas (no floating arithmetic in the statements),
followed by strict rational proximity checks to the frozen runtime anchors.
-/

/-- Shared total SM gauge-generator count from group dimensions: 8 + 3 + 1 = 12. -/
def totalGaugeGeneratorsQ : ℚ :=
  (((3 ^ 2 - 1) + (2 ^ 2 - 1) + 1 : ℕ) : ℚ)

theorem total_gauge_generators_eq :
    totalGaugeGeneratorsQ = 12 := by
  native_decide

/-- Frozen runtime anchor for `p_ratio` (triangulation lane). -/
def pFrozenQ : ℚ := 13688110433760 / 1000000000000

/-- Frozen runtime anchor for `kappa_geo` (triangulation lane). -/
def kappaFrozenQ : ℚ := 34697396055505 / 1000000000000

/-- Frozen runtime anchor for EW bridge coefficient (triangulation lane). -/
def ewCoeffFrozenQ : ℚ := 8460487692308 / 1000000000000

/-- Candidate exponent from shared Cl(1,3) counts:
    `p = α⁻¹/(|grade₁|+|grade₂|) - 1/((|grade₂|+1)*N_gauge)`. -/
def pCandidateQ : ℚ :=
  (alphaInverse 4 : ℚ) / ((Nat.choose 4 1 + Nat.choose 4 2 : ℕ) : ℚ)
    - 1 / (((Nat.choose 4 2 + 1 : ℕ) : ℚ) * totalGaugeGeneratorsQ)

theorem p_candidate_closed_form :
    pCandidateQ = (137 / 10 : ℚ) - 1 / (7 * 12 : ℚ) := by
  unfold pCandidateQ
  rw [alpha_inverse_d4, total_gauge_generators_eq]
  native_decide

theorem p_candidate_eq_5749_over_420 :
    pCandidateQ = 5749 / 420 := by
  rw [p_candidate_closed_form]
  native_decide

/-- Candidate neutrino normalization uplift from shared counts:
    `κ = (60/11) * (19/3 + 1/36 + 1/(7*13*136))` in closed form. -/
def kappaCandidateQ : ℚ :=
  geometricDarkToVisibleRatio *
    ((((2 ^ 4) + magneticTriplet.card : ℕ) : ℚ) / (magneticTriplet.card : ℚ)
      + 1 / (((Nat.choose 4 2 : ℕ) : ℚ) ^ 2)
      + 1 / (((Nat.choose 4 2 + 1 : ℕ) : ℚ)
            * (((2 ^ 4) - magneticTriplet.card : ℕ) : ℚ)
            * (triangularNumber (2 ^ 4) : ℚ)))

theorem kappa_candidate_closed_form :
    kappaCandidateQ = (60 / 11 : ℚ) * ((19 / 3 : ℚ) + 1 / 36 + 1 / (7 * 13 * 136 : ℚ)) := by
  unfold kappaCandidateQ
  rw [geometric_dark_to_visible_ratio_eq]
  have hs : magneticTriplet.card = 3 := su2_dim
  have hpow : (2 ^ 4 : ℕ) = 16 := by native_decide
  have hT : triangularNumber (2 ^ 4) = 136 := by
    rw [hpow, T16_eq_136]
  rw [hs, hT]
  native_decide

/-- Candidate EW bridge uplift coefficient from shared counts:
    `c = d/2 + |grade₂|/(d-|SU(2)|) - 1/((|grade₂|+1)T(16))`. -/
def ewCoeffCandidateQ : ℚ :=
  ((2 ^ 4 : ℚ) / 2)
    + ((Nat.choose 4 2 : ℕ) : ℚ) / (((2 ^ 4) - magneticTriplet.card : ℕ) : ℚ)
    - 1 / (((Nat.choose 4 2 + 1 : ℕ) : ℚ) * (triangularNumber (2 ^ 4) : ℚ))

theorem ew_coeff_candidate_closed_form :
    ewCoeffCandidateQ = (8 : ℚ) + 6 / 13 - 1 / (7 * 136 : ℚ) := by
  unfold ewCoeffCandidateQ
  have hs : magneticTriplet.card = 3 := su2_dim
  have hpow : (2 ^ 4 : ℕ) = 16 := by native_decide
  have hT : triangularNumber (2 ^ 4) = 136 := by
    rw [hpow, T16_eq_136]
  rw [hs, hT]
  native_decide

/-- `p` candidate is within `2e-5` of frozen triangulation anchor. -/
theorem p_candidate_close_to_frozen :
    |pCandidateQ - pFrozenQ| < 1 / 50000 := by
  rw [p_candidate_eq_5749_over_420]
  native_decide

/-- `kappa` candidate is within `2e-5` of frozen triangulation anchor. -/
theorem kappa_candidate_close_to_frozen :
    |kappaCandidateQ - kappaFrozenQ| < 1 / 50000 := by
  rw [kappa_candidate_closed_form]
  native_decide

/-- EW coefficient candidate is within `1e-6` of frozen triangulation anchor. -/
theorem ew_coeff_candidate_close_to_frozen :
    |ewCoeffCandidateQ - ewCoeffFrozenQ| < 1 / 1000000 := by
  rw [ew_coeff_candidate_closed_form]
  native_decide

/-- Coupled EW prediction at M_Z:
    `sin²θ_W(M_Z) = 3/13 + α² * ewCoeffCandidateQ`. -/
def sin2ThetaWMzCoupledQ : ℚ :=
  (3 / 13 : ℚ) + (1 / (137 * 137 : ℚ)) * ewCoeffCandidateQ

theorem sin2_theta_w_mz_coupled_closed_form :
    sin2ThetaWMzCoupledQ =
      (3 / 13 : ℚ) + (1 / (137 * 137 : ℚ)) * ((8 : ℚ) + 6 / 13 - 1 / (7 * 136 : ℚ)) := by
  unfold sin2ThetaWMzCoupledQ
  rw [ew_coeff_candidate_closed_form]

/-- Coupled EW bridge lands on the 0.23122 target at sub-1e-9 absolute error. -/
theorem sin2_theta_w_mz_coupled_near_target :
    |sin2ThetaWMzCoupledQ - (23122 / 100000 : ℚ)| < 1 / 1000000000 := by
  rw [sin2_theta_w_mz_coupled_closed_form]
  native_decide

end Gutoe.TriangulatedConstants
