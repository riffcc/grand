import Mathlib
import Gutoe.FineStructure
import Gutoe.MassSpectrum
import Gutoe.DarkMatterSector
import Gutoe.Inflation
import Gutoe.GravityMetric

namespace Gutoe.ElectronScaleTransduction

open Gutoe.FineStructure
open Gutoe.MassSpectrum
open Gutoe.DarkMatterSector
open Gutoe.Inflation
open Gutoe.GravityMetric

/-- Structural leading EM factor used in the electron transduction lane:
    `α = 1/137` from `alphaInverse 4`. -/
noncomputable def alphaStructural : ℝ :=
  1 / (alphaInverse 4 : ℝ)

/-- Structural corrected dark-to-visible ratio used after GRAND-356:
    `R = 115/22`. -/
noncomputable def correctedDarkVisibleRatio : ℝ :=
  (correctedUnifiedBudgetDarkToVisibleRatio : ℝ)

/-- Structural count split ratio:
    `d = 5/11`. -/
noncomputable def darkVisibleCountRatio : ℝ :=
  (darkToVisibleCountRatio : ℝ)

/-- Gauge-generator count from the SM decomposition:
    `N_g = 12`. -/
def gaugeGeneratorCount : ℕ := nLayers

/-- Flagged candidate electron-scale transduction factor:
    `F = α^13 * R^3 * C_inf * λ_QG^-3`. -/
noncomputable def electronScaleFactorFlagged : ℝ :=
  alphaStructural ^ 13
    * correctedDarkVisibleRatio ^ 3
    * inflationCorrectionCInf
    * lambda_qg ^ (-3 : ℤ)

/-- Equivalent candidate with explicit gauge-cube amplification:
    `F = α^13 * R^3 * C_inf * 12^3`. -/
noncomputable def electronScaleFactorGaugeCube : ℝ :=
  alphaStructural ^ 13
    * correctedDarkVisibleRatio ^ 3
    * inflationCorrectionCInf
    * ((gaugeGeneratorCount : ℝ) ^ 3)

theorem alpha_structural_eq :
    alphaStructural = (1 : ℝ) / 137 := by
  unfold alphaStructural
  rw [alpha_inverse_d4]
  norm_num

theorem corrected_dark_visible_ratio_eq :
    correctedDarkVisibleRatio = (115 : ℝ) / 22 := by
  unfold correctedDarkVisibleRatio
  norm_num [corrected_unified_budget_dark_to_visible_ratio_eq]

theorem dark_visible_count_ratio_eq :
    darkVisibleCountRatio = (5 : ℝ) / 11 := by
  unfold darkVisibleCountRatio
  norm_num [dark_to_visible_count_ratio_eq]

theorem gauge_generator_count_eq_12 :
    gaugeGeneratorCount = 12 := by
  rfl

theorem lambda_qg_inv_cube_eq_12_cube :
    lambda_qg ^ (-3 : ℤ) = (((12 : ℕ) : ℝ) ^ 3) := by
  unfold lambda_qg
  norm_num

theorem electron_scale_flagged_eq_gauge_cube :
    electronScaleFactorFlagged = electronScaleFactorGaugeCube := by
  unfold electronScaleFactorFlagged electronScaleFactorGaugeCube
  rw [lambda_qg_inv_cube_eq_12_cube]
  rw [gauge_generator_count_eq_12]

end Gutoe.ElectronScaleTransduction
