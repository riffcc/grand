import Mathlib
import Gutoe.FineStructure
import Gutoe.DimensionalStructure
import Gutoe.Z3Uniqueness
import Gutoe.DarkMatterSector
import Gutoe.Inflation

namespace Gutoe.StrongCouplingCInfBridge

open Gutoe.FineStructure
open Gutoe.DimensionalStructure
open Gutoe.Z3Uniqueness
open Gutoe.DarkMatterSector
open Gutoe.Inflation

/-- Shared second-order correction from common finite Cl(1,3) counts:
    `1 + 1/(|grade-2| * |visible|) = 1 + 1/(6*11)`. -/
def sharedSecondOrderCorrectionQ : ℚ :=
  1 + 1 / ((grade2_4d.card : ℚ) * (visibleSectorStates.card : ℚ))

theorem shared_second_order_correction_eq_67_over_66 :
    sharedSecondOrderCorrectionQ = 67 / 66 := by
  unfold sharedSecondOrderCorrectionQ
  have h2 : grade2_4d.card = 6 := by decide
  rcases visible_dark_state_count_split with ⟨hVis, _, _, _⟩
  rw [h2, hVis]
  norm_num

/-- Leading structural strong coupling candidate at the Z-pole:
    `α_s = 2^4 / α_EM^{-1} = 16/137`. -/
def alphaSStructuralLeadingQ : ℚ :=
  ((2 ^ 4 : ℕ) : ℚ) / (alphaInverse 4 : ℚ)

theorem alpha_s_structural_leading_eq_16_over_137 :
    alphaSStructuralLeadingQ = 16 / 137 := by
  unfold alphaSStructuralLeadingQ
  rw [alpha_inverse_d4]
  norm_num [clifford_dim_eq_16]

/-- Cross-sector corrected strong coupling candidate:
    `(16/137) * (67/66)` using the same correction as inflation. -/
def alphaSStructuralCorrectedQ : ℚ :=
  alphaSStructuralLeadingQ * sharedSecondOrderCorrectionQ

theorem alpha_s_structural_corrected_eq :
    alphaSStructuralCorrectedQ = (16 / 137) * (67 / 66) := by
  unfold alphaSStructuralCorrectedQ
  rw [alpha_s_structural_leading_eq_16_over_137, shared_second_order_correction_eq_67_over_66]

/-- The correction factor used in the strong-coupling lane equals the inflation
    `C_inf` correction factor (`67/66`) in the shared finite-count lane. -/
theorem shared_correction_matches_inflation_cinf :
    (sharedSecondOrderCorrectionQ : ℝ) = inflationCorrectionCInf := by
  have hq : sharedSecondOrderCorrectionQ = 67 / 66 :=
    shared_second_order_correction_eq_67_over_66
  have hi : inflationCorrectionCInf = (67 : ℝ) / 66 :=
    inflation_cinf_eq_67_over_66
  rw [hq, hi]
  norm_num

end Gutoe.StrongCouplingCInfBridge
