import Mathlib
import Gutoe.FineStructure
import Gutoe.MassSpectrum
import Gutoe.GaugeConstants
import Gutoe.DarkMatterSector
import Gutoe.Z3Uniqueness

namespace Gutoe.HadronTransduction

open Gutoe.FineStructure
open Gutoe.MassSpectrum
open Gutoe.GaugeConstants
open Gutoe.DarkMatterSector
open Gutoe.Z3Uniqueness

/-- Grade-2 state count from the shared Cl(1,3) combinatorics: `C(4,2)=6`. -/
def grade2StateCount : ℕ := Nat.choose 4 2

/-- Shared denominator for the pion transduction factor:
    `dim Cl(1,3) * |grade2| = 16 * 6 = 96`. -/
def pionTransductionDenominator : ℕ := (2 ^ 4) * grade2StateCount

/-- Structural pion transduction factor used in GRAND-353:
    `(mp/me) / (16*6)`. -/
def pionTransductionFactorQ : ℚ :=
  (mpMeAlgebraic : ℚ) / (pionTransductionDenominator : ℚ)

theorem grade2_state_count_eq : grade2StateCount = 6 := by
  native_decide

theorem pion_transduction_denominator_eq : pionTransductionDenominator = 96 := by
  unfold pionTransductionDenominator grade2StateCount
  native_decide

theorem pion_transduction_factor_closed_form :
    pionTransductionFactorQ = 153 / 8 := by
  unfold pionTransductionFactorQ
  rw [mp_me_eq_1836, pion_transduction_denominator_eq]
  norm_num

/-- Structural neutron-proton split factor from the pion lane:
    `α * (mZ/mW)^2 = (1/137)*(13/10)`. -/
def neutronSplitFromPionFactorQ : ℚ :=
  ((1 : ℚ) / (alphaInverse 4 : ℚ)) *
    ((1 : ℚ) / (1 - (magneticTriplet.card : ℚ) / (2 ^ 4 - magneticTriplet.card : ℚ)))

theorem neutron_split_from_pion_factor_closed_form :
    neutronSplitFromPionFactorQ = 13 / 1370 := by
  unfold neutronSplitFromPionFactorQ
  rw [alpha_inverse_d4, mZ_over_mW_sq_from_z3]
  norm_num

/-- Structural damping factor applied to the QCD scale:
    `(visible-1)/16 = (11-1)/16 = 5/8`. -/
def qcdVisibilityDampingFactorQ : ℚ :=
  ((visibleSectorStates.card : ℚ) - 1) / (2 ^ 4 : ℚ)

theorem qcd_visibility_damping_factor_closed_form :
    qcdVisibilityDampingFactorQ = 5 / 8 := by
  unfold qcdVisibilityDampingFactorQ
  rcases visible_dark_state_count_split with ⟨hVis, _, _, _⟩
  rw [hVis]
  norm_num

/-- Structural kaon/pion factor from corrected dark/visible ratio:
    `1 + (115/22)/2 = 159/44`. -/
def kaonToPionFactorQ : ℚ :=
  1 + correctedUnifiedBudgetDarkToVisibleRatio / 2

theorem kaon_to_pion_factor_closed_form :
    kaonToPionFactorQ = 159 / 44 := by
  unfold kaonToPionFactorQ
  rw [corrected_unified_budget_dark_to_visible_ratio_eq]
  norm_num

/-- GRAND-353 Lean parity spine for hadron transduction factors. -/
theorem hadron_transduction_structural_spine :
    pionTransductionFactorQ = 153 / 8 ∧
    neutronSplitFromPionFactorQ = 13 / 1370 ∧
    qcdVisibilityDampingFactorQ = 5 / 8 ∧
    kaonToPionFactorQ = 159 / 44 := by
  exact ⟨pion_transduction_factor_closed_form,
    neutron_split_from_pion_factor_closed_form,
    qcd_visibility_damping_factor_closed_form,
    kaon_to_pion_factor_closed_form⟩

end Gutoe.HadronTransduction
