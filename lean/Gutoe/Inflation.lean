import Mathlib
import Gutoe.DarkMatterSector
import Gutoe.GravityMetric
import Gutoe.FineStructure
import Gutoe.Baryogenesis
import Gutoe.Z3Uniqueness

namespace Gutoe.Inflation

open Gutoe.DarkMatterSector
open Gutoe.GravityMetric
open Gutoe.FineStructure
open Gutoe.Baryogenesis
open Gutoe.Z3Uniqueness

/-- Structural inflation e-fold count from shared Clifford dark-sector counts:
    `N = 12 * 5 = 60`. -/
noncomputable def inflationEfoldCount : ℝ :=
  (geometricDarkAmplificationQ : ℝ) * (darkSectorCandidates.card : ℝ)

/-- Slow-roll epsilon in the structural plateau lane. -/
noncomputable def slowRollEpsilon (N : ℝ) : ℝ :=
  3 / (4 * N ^ 2)

/-- Slow-roll eta in the structural plateau lane. -/
noncomputable def slowRollEta (N : ℝ) : ℝ :=
  -1 / N

/-- Scalar spectral index from first-order slow-roll observables. -/
noncomputable def scalarSpectralIndex (N : ℝ) : ℝ :=
  1 - 6 * slowRollEpsilon N + 2 * slowRollEta N

/-- Tensor-to-scalar ratio from epsilon. -/
noncomputable def tensorToScalarRatio (N : ℝ) : ℝ :=
  16 * slowRollEpsilon N

/-- End-of-inflation proxy where `ε = 1` in this lane. -/
noncomputable def inflationEndN : ℝ :=
  Real.sqrt 3 / 2

/-- Total expansion factor `exp(N)`. -/
noncomputable def inflationExpansionFactor : ℝ :=
  Real.exp inflationEfoldCount

/-- Structural inflation Hubble ratio `H/M_pl` from shared primitives:
    `α_LO^2 * (60/11) * (1-λ_QG) * (3/6) * 1/sqrt(486)`. -/
noncomputable def inflationHubbleRatio : ℝ :=
  alphaLeadingOrder ^ 2
    * (geometricDarkToVisibleRatio : ℝ)
    * (1 - lambda_qg)
    * ((magneticTriplet.card : ℝ) / ((Nat.choose 4 2 : ℕ) : ℝ))
    * (1 / Real.sqrt (baryoMicroModeCount : ℝ))

/-- Structural inflation micro-correction candidate from shared finite counts:
    `C_inf = 1 + 1/(|grade-2| * |visible|) = 1 + 1/(6*11)`. -/
noncomputable def inflationCorrectionCInf : ℝ :=
  1 + 1 / (((Nat.choose 4 2 : ℕ) : ℝ) * (visibleSectorStates.card : ℝ))

/-- Scalar amplitude from slow-roll relation:
    `A_s = (H/M_pl)^2 / (8π² ε)`. -/
noncomputable def scalarAmplitude : ℝ :=
  inflationHubbleRatio ^ 2 / (8 * Real.pi ^ 2 * slowRollEpsilon inflationEfoldCount)

/-- Reheating equation-of-state proxy from shared dark fraction. -/
noncomputable def reheatingW : ℝ :=
  (darkFractionOfTotalStates : ℝ)

/-- Structural reheating e-fold count: `N_reh = N / 12 = 5`. -/
noncomputable def reheatingEfoldCount : ℝ :=
  inflationEfoldCount / (geometricDarkAmplificationQ : ℝ)

theorem dark_sector_card_eq_five :
    darkSectorCandidates.card = 5 := by
  rcases visible_dark_state_count_split with ⟨_, hDark, _, _⟩
  exact hDark

theorem inflation_efolds_eq_60 :
    inflationEfoldCount = 60 := by
  unfold inflationEfoldCount
  rw [geometric_dark_amplification_eq, dark_sector_card_eq_five]
  norm_num

theorem epsilon_structural_eq :
    slowRollEpsilon inflationEfoldCount = (1 : ℝ) / 4800 := by
  rw [inflation_efolds_eq_60]
  norm_num [slowRollEpsilon]

theorem eta_structural_eq :
    slowRollEta inflationEfoldCount = -(1 : ℝ) / 60 := by
  rw [inflation_efolds_eq_60]
  norm_num [slowRollEta]

theorem ns_structural_eq :
    scalarSpectralIndex inflationEfoldCount = (2317 : ℝ) / 2400 := by
  rw [inflation_efolds_eq_60]
  norm_num [scalarSpectralIndex, slowRollEpsilon, slowRollEta]

theorem r_structural_eq :
    tensorToScalarRatio inflationEfoldCount = (1 : ℝ) / 300 := by
  rw [inflation_efolds_eq_60]
  norm_num [tensorToScalarRatio, slowRollEpsilon]

theorem inflation_hubble_ratio_eq :
    inflationHubbleRatio
      = ((1 : ℝ) / 137) ^ 2
          * ((60 : ℝ) / 11)
          * ((11 : ℝ) / 12)
          * ((3 : ℝ) / 6)
          * (1 / Real.sqrt 486) := by
  unfold inflationHubbleRatio
  rw [alpha_leading_order_eq, geometric_dark_to_visible_ratio_eq]
  have hLam : lambda_qg = (1 : ℝ) / 12 := by
    unfold lambda_qg
    norm_num
  rw [hLam]
  have hs : magneticTriplet.card = 3 := su2_dim
  have h20 : (Nat.choose 4 2 : ℕ) = 6 := by native_decide
  rw [hs]
  rw [h20]
  norm_num [baryo_micro_mode_count_eq_486]

/-- Inflation Hubble ratio with the corrected shared dark budget
    `115/22 = (60 - 5/2)/11`, keeping all other structural factors unchanged. -/
theorem inflation_hubble_ratio_with_corrected_budget_eq :
    alphaLeadingOrder ^ 2
      * (correctedUnifiedBudgetDarkToVisibleRatio : ℝ)
      * (1 - lambda_qg)
      * ((magneticTriplet.card : ℝ) / ((Nat.choose 4 2 : ℕ) : ℝ))
      * (1 / Real.sqrt (baryoMicroModeCount : ℝ))
      = ((1 : ℝ) / 137) ^ 2
          * ((115 : ℝ) / 22)
          * ((11 : ℝ) / 12)
          * ((3 : ℝ) / 6)
          * (1 / Real.sqrt 486) := by
  rw [alpha_leading_order_eq]
  have hcorr : (correctedUnifiedBudgetDarkToVisibleRatio : ℝ) = (115 : ℝ) / 22 := by
    norm_num [corrected_unified_budget_dark_to_visible_ratio_eq]
  rw [hcorr]
  have hLam : lambda_qg = (1 : ℝ) / 12 := by
    unfold lambda_qg
    norm_num
  rw [hLam]
  have hs : magneticTriplet.card = 3 := su2_dim
  have h20 : (Nat.choose 4 2 : ℕ) = 6 := by native_decide
  rw [hs, h20]
  norm_num [baryo_micro_mode_count_eq_486]

/-- Exact structural inflation correction from shared counts:
    `C_inf = 1 + 1/(6*11) = 67/66`. -/
theorem inflation_cinf_eq :
    inflationCorrectionCInf
      = 1 + 1 / ((6 : ℝ) * 11) := by
  unfold inflationCorrectionCInf
  have hvis : visibleSectorStates.card = 11 := by
    rcases visible_dark_state_count_split with ⟨hVis, _, _, _⟩
    exact hVis
  have hchoose : (Nat.choose 4 2 : ℕ) = 6 := by
    native_decide
  rw [hchoose]
  norm_num [hvis]

theorem inflation_cinf_eq_67_over_66 :
    inflationCorrectionCInf = (67 : ℝ) / 66 := by
  rw [inflation_cinf_eq]
  norm_num

theorem inflation_hubble_ratio_pos :
    0 < inflationHubbleRatio := by
  rw [inflation_hubble_ratio_eq]
  positivity

theorem scalar_amplitude_pos :
    0 < scalarAmplitude := by
  unfold scalarAmplitude
  have heps : 0 < slowRollEpsilon inflationEfoldCount := by
    rw [epsilon_structural_eq]
    norm_num
  have hnum : 0 < inflationHubbleRatio ^ 2 := by
    exact sq_pos_of_pos inflation_hubble_ratio_pos
  have hden : 0 < 8 * Real.pi ^ 2 * slowRollEpsilon inflationEfoldCount := by
    have hpi : 0 < Real.pi := Real.pi_pos
    positivity
  exact div_pos hnum hden

theorem reheating_w_eq :
    reheatingW = (5 : ℝ) / 16 := by
  unfold reheatingW
  norm_num [dark_fraction_of_total_states_eq]

theorem reheating_efolds_eq :
    reheatingEfoldCount = 5 := by
  unfold reheatingEfoldCount
  rw [inflation_efolds_eq_60, geometric_dark_amplification_eq]
  norm_num

theorem ns_in_observational_window :
    (0.955 : ℝ) ≤ scalarSpectralIndex inflationEfoldCount ∧
    scalarSpectralIndex inflationEfoldCount ≤ 0.975 := by
  rw [ns_structural_eq]
  constructor <;> norm_num

theorem r_below_current_upper_bound :
    tensorToScalarRatio inflationEfoldCount < 0.06 := by
  rw [r_structural_eq]
  norm_num

theorem graceful_exit_condition :
    slowRollEpsilon inflationEndN = 1 := by
  unfold slowRollEpsilon inflationEndN
  have hsq : (Real.sqrt 3) ^ 2 = 3 := by
    have h3 : (0 : ℝ) ≤ 3 := by positivity
    simpa using Real.sq_sqrt h3
  calc
    3 / (4 * (Real.sqrt 3 / 2) ^ 2)
        = 3 / (4 * ((Real.sqrt 3) ^ 2 / 4)) := by ring
    _ = 3 / ((Real.sqrt 3) ^ 2) := by ring
    _ = 3 / 3 := by rw [hsq]
    _ = 1 := by norm_num

theorem expansion_factor_gt_zero :
    0 < inflationExpansionFactor := by
  unfold inflationExpansionFactor
  positivity

end Gutoe.Inflation
