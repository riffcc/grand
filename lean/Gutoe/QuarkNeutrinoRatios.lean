import Mathlib
import Gutoe.FlavorMixing
import Gutoe.DarkMatterSector
import Gutoe.StellarFusion
import Gutoe.FineStructure

namespace Gutoe.QuarkNeutrinoRatios

open Gutoe.Z3Uniqueness
open Gutoe.DarkMatterSector
open Gutoe.StellarFusion
open Gutoe.FineStructure

/-- Shared CKM suppression inverse-square:
    `λ⁻² = (2^4 + |SU(2)|) = 19`. -/
def ckmLambdaInvSqQ : ℚ :=
  (((2 ^ 4) + magneticTriplet.card : ℕ) : ℚ)

theorem ckm_lambda_inv_sq_eq_19 :
    ckmLambdaInvSqQ = 19 := by
  native_decide

/-- Shared finite-count correction `C_inf = 1 + 1/(|grade₂|*|visible|)`. -/
def cInfQ : ℚ :=
  1 + 1 / (((Nat.choose 4 2 : ℕ) : ℚ) * (visibleSectorStates.card : ℚ))

theorem c_inf_eq_67_over_66 :
    cInfQ = 67 / 66 := by
  unfold cInfQ
  have h2 : (Nat.choose 4 2 : ℕ) = 6 := by native_decide
  rcases visible_dark_state_count_split with ⟨hVis, _, _, _⟩
  rw [h2, hVis]
  norm_num

/-- Ticket GRAND-82 lane: first-generation quark ratio. -/
def muOverMdQ : ℚ :=
  ((magneticTriplet.card ^ 2 - 1 : ℕ) : ℚ) / (((2 ^ 4) + 1 : ℕ) : ℚ)

theorem mu_over_md_eq_8_over_17 :
    muOverMdQ = 8 / 17 := by
  native_decide

/-- Ticket GRAND-83 lane: second/third generation split ratio `m_c/m_s`. -/
def mcOverMsQ : ℚ :=
  (((2 ^ 4) - magneticTriplet.card : ℕ) : ℚ) /
      ((magneticTriplet.card * ((Nat.choose 4 2) + 1) : ℕ) : ℚ)
    * ckmLambdaInvSqQ * cInfQ

/-- Ticket GRAND-83 lane: third-generation heavy split ratio `m_t/m_b`. -/
def mtOverMbQ : ℚ :=
  (((2 ^ 4) - magneticTriplet.card : ℕ) : ℚ) / ((Nat.choose 4 2 : ℕ) : ℚ)
    * ckmLambdaInvSqQ * cInfQ

/-- Ticket GRAND-83 lane: cross-generation ratio `m_c/m_u`. -/
def mcOverMuQ : ℚ :=
  ((magneticTriplet.card ^ 2 - 1 : ℕ) : ℚ) /
      (((Nat.choose 4 1) + 1 : ℕ) : ℚ)
    * (ckmLambdaInvSqQ ^ 2) * cInfQ

/-- Ticket GRAND-83 lane: cross-generation ratio `m_t/m_c`. -/
def mtOverMcQ : ℚ :=
  ((magneticTriplet.card ^ 2 - 1 : ℕ) : ℚ) * (((2 ^ 4) + 1 : ℕ) : ℚ)

/-- Ticket GRAND-83 lane: down-sector split `m_s/m_d`. -/
def msOverMdQ : ℚ := ckmLambdaInvSqQ

/-- Ticket GRAND-83 lane: down-sector heavy split `m_b/m_s`. -/
def mbOverMsQ : ℚ :=
  ((magneticTriplet.card ^ 2 - 1 : ℕ) : ℚ) / (magneticTriplet.card : ℚ)
    * ckmLambdaInvSqQ * cInfQ

theorem mt_over_mc_eq_136 :
    mtOverMcQ = 136 := by
  unfold mtOverMcQ
  native_decide

/-- Legacy GRAND ticket windows are approximate (order-level), so we gate at 10%
    relative tolerance against the historical central values. -/
theorem quark_ratio_within_legacy_ten_percent :
    |muOverMdQ - (47 / 100 : ℚ)| ≤ (47 / 100 : ℚ) / 10 ∧
    |mcOverMsQ - (117 / 10 : ℚ)| ≤ (117 / 10 : ℚ) / 10 ∧
    |mtOverMbQ - (413 / 10 : ℚ)| ≤ (413 / 10 : ℚ) / 10 ∧
    |mcOverMuQ - (580 : ℚ)| ≤ (580 : ℚ) / 10 ∧
    |mtOverMcQ - (136 : ℚ)| ≤ (136 : ℚ) / 10 ∧
    |msOverMdQ - (20 : ℚ)| ≤ (20 : ℚ) / 10 ∧
    |mbOverMsQ - (51 : ℚ)| ≤ (51 : ℚ) / 10 := by
  unfold muOverMdQ mcOverMsQ mtOverMbQ mcOverMuQ mtOverMcQ msOverMdQ mbOverMsQ
  rw [ckm_lambda_inv_sq_eq_19, c_inf_eq_67_over_66]
  native_decide

/-- Structural neutrino/electron scale suppression used by the tiny-mass lane:
    `α^4 * (60/11)` from shared finite Cl(1,3) counts. -/
def neutrinoScaleOverElectronQ : ℚ :=
  ((1 : ℚ) / (alphaInverse 4 : ℚ)) ^ 4 * geometricDarkToVisibleRatio

theorem neutrino_scale_over_electron_eq :
    neutrinoScaleOverElectronQ = ((1 : ℚ) / 137) ^ 4 * (60 / 11 : ℚ) := by
  unfold neutrinoScaleOverElectronQ
  rw [alpha_inverse_d4, geometric_dark_to_visible_ratio_eq]
  ring

/-- Electron mass anchor in eV from shared fusion mass table. -/
def electronMassAnchorEvQ : ℚ :=
  electronRestMassMeV * 1000000

theorem electron_anchor_ev_eq :
    electronMassAnchorEvQ = 511000 := by
  unfold electronMassAnchorEvQ
  norm_num [electronRestMassMeV]

/-- Neutrino absolute scale in eV used in the Rust tiny-mass report. -/
def neutrinoScaleEvQ : ℚ :=
  electronMassAnchorEvQ * neutrinoScaleOverElectronQ

theorem neutrino_scale_ev_pos :
    0 < neutrinoScaleEvQ := by
  unfold neutrinoScaleEvQ
  rw [electron_anchor_ev_eq, neutrino_scale_over_electron_eq]
  positivity

/-- KATRIN-scale cap witness (`m_ν < 0.8 eV`) for the structural scale. -/
theorem neutrino_scale_ev_below_katrin_cap :
    neutrinoScaleEvQ < (4 / 5 : ℚ) := by
  unfold neutrinoScaleEvQ
  rw [electron_anchor_ev_eq, neutrino_scale_over_electron_eq]
  native_decide

/-- Three-mode sum cap witness (`Σm_ν < 0.12 eV`) at structural scale upper bound. -/
theorem neutrino_three_mode_sum_below_cosmology_cap :
    3 * neutrinoScaleEvQ < (3 / 25 : ℚ) := by
  unfold neutrinoScaleEvQ
  rw [electron_anchor_ev_eq, neutrino_scale_over_electron_eq]
  native_decide

/-- Any normalized positive mode scaled by the structural neutrino factor is
    strictly nonzero and below the KATRIN cap. -/
theorem normalized_mode_scaled_nonzero_tiny
    (r : ℚ) (hr_pos : 0 < r) (hr_le : r ≤ 1) :
    0 < r * neutrinoScaleEvQ ∧ r * neutrinoScaleEvQ < (4 / 5 : ℚ) := by
  constructor
  · have hs : 0 < neutrinoScaleEvQ := neutrino_scale_ev_pos
    positivity
  · have hs : neutrinoScaleEvQ < (4 / 5 : ℚ) := neutrino_scale_ev_below_katrin_cap
    have hs_pos : 0 ≤ neutrinoScaleEvQ := le_of_lt neutrino_scale_ev_pos
    have hmul : r * neutrinoScaleEvQ ≤ 1 * neutrinoScaleEvQ := by
      nlinarith
    have h1 : 1 * neutrinoScaleEvQ = neutrinoScaleEvQ := by ring
    rw [h1] at hmul
    exact lt_of_le_of_lt hmul hs

end Gutoe.QuarkNeutrinoRatios
