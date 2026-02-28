import Mathlib
import Gutoe.FineStructure
import Gutoe.DimensionalStructure
import Gutoe.DarkMatterSector
import Gutoe.NuclearFirstPrinciples

namespace Gutoe.AbiogenesisThreshold

open Gutoe.FineStructure
open Gutoe.DimensionalStructure
open Gutoe.DarkMatterSector
open Gutoe.NuclearFirstPrinciples

/-- Structural monomer count used by the abiogenesis lane:
    total Cl(1,3) basis states plus grade-1 channels, `16 + 4 = 20`. -/
def abiogenesisMonomerCountQ : ℚ :=
  cliffordDimQ + (grade1_4d.card : ℚ)

theorem abiogenesis_monomer_count_eq_20 :
    abiogenesisMonomerCountQ = 20 := by
  unfold abiogenesisMonomerCountQ
  rw [clifford_dim_q_eq_16, grade1_state_count_eq]
  norm_num

/-- Structural polymer-length scale from finite basis plus vacuum channel:
    `16 + 1 = 17`. -/
def abiogenesisPolymerScaleQ : ℚ :=
  cliffordDimQ + 1

theorem abiogenesis_polymer_scale_eq_17 :
    abiogenesisPolymerScaleQ = 17 := by
  unfold abiogenesisPolymerScaleQ
  rw [clifford_dim_q_eq_16]
  norm_num

/-- Conservative catalytic channel factor from the polymer scale minus
    grade-2 count (`C(4,2)=6`): `17 - 6 = 11`. -/
def abiogenesisCatalyticScaleMinQ : ℚ :=
  abiogenesisPolymerScaleQ - ((Nat.choose 4 2 : ℕ) : ℚ)

theorem abiogenesis_catalytic_scale_min_eq_11 :
    abiogenesisCatalyticScaleMinQ = 11 := by
  unfold abiogenesisCatalyticScaleMinQ
  rw [abiogenesis_polymer_scale_eq_17, grade2_dim_eq_6]
  norm_num

/-- Leading-order electromagnetic suppression used by the lane:
    `α_LO = 1/137`. -/
def abiogenesisAlphaLeadingQ : ℚ :=
  1 / (alphaInverse 4 : ℚ)

theorem abiogenesis_alpha_leading_eq_one_over_137 :
    abiogenesisAlphaLeadingQ = 1 / 137 := by
  unfold abiogenesisAlphaLeadingQ
  rw [alpha_inverse_d4]
  norm_num

/-- Geometric-contact factor from shared dark-sector geometry:
    `f_contact = 60/71`. -/
def abiogenesisContactFractionQ : ℚ :=
  geometricDarkFractionOfMatter

theorem abiogenesis_contact_fraction_eq_60_over_71 :
    abiogenesisContactFractionQ = 60 / 71 := by
  unfold abiogenesisContactFractionQ
  exact geometric_dark_fraction_of_matter_eq

/-- Conservative catalytic probability lower bound:
    `p_min = α_LO * 11 * (60/71)`. -/
def abiogenesisCatalyticProbabilityMinQ : ℚ :=
  abiogenesisAlphaLeadingQ * abiogenesisCatalyticScaleMinQ * abiogenesisContactFractionQ

theorem abiogenesis_catalytic_probability_min_eq :
    abiogenesisCatalyticProbabilityMinQ = 660 / 9727 := by
  unfold abiogenesisCatalyticProbabilityMinQ
  rw [abiogenesis_alpha_leading_eq_one_over_137,
    abiogenesis_catalytic_scale_min_eq_11,
    abiogenesis_contact_fraction_eq_60_over_71]
  norm_num

theorem abiogenesis_catalytic_probability_min_gt_five_percent :
    abiogenesisCatalyticProbabilityMinQ > 1 / 20 := by
  rw [abiogenesis_catalytic_probability_min_eq]
  norm_num

/-- Kauffman closure threshold in the structural lane. -/
def kauffmanClosureThresholdQ : ℚ := 1

/-- Structural closure control `N * p`. -/
def abiogenesisClosureControlQ : ℚ :=
  abiogenesisMonomerCountQ * abiogenesisCatalyticProbabilityMinQ

theorem abiogenesis_closure_control_eq :
    abiogenesisClosureControlQ = 13200 / 9727 := by
  unfold abiogenesisClosureControlQ
  rw [abiogenesis_monomer_count_eq_20, abiogenesis_catalytic_probability_min_eq]
  norm_num

/-- The Kauffman closure gate is crossed in the structural lane:
    `N * p > 1`. -/
theorem abiogenesis_kauffman_closure_exceeds_threshold :
    abiogenesisClosureControlQ > kauffmanClosureThresholdQ := by
  unfold kauffmanClosureThresholdQ
  rw [abiogenesis_closure_control_eq]
  norm_num

/-- Structural closure margin above threshold. -/
def abiogenesisClosureMarginQ : ℚ :=
  abiogenesisClosureControlQ - kauffmanClosureThresholdQ

theorem abiogenesis_closure_margin_eq :
    abiogenesisClosureMarginQ = 3473 / 9727 := by
  unfold abiogenesisClosureMarginQ kauffmanClosureThresholdQ
  rw [abiogenesis_closure_control_eq]
  norm_num

theorem abiogenesis_closure_margin_gt_quarter :
    abiogenesisClosureMarginQ > 1 / 4 := by
  rw [abiogenesis_closure_margin_eq]
  norm_num

end Gutoe.AbiogenesisThreshold
