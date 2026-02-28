import Mathlib
import Gutoe.FineStructure
import Gutoe.GravityMetric
import Gutoe.Z3Uniqueness
import Gutoe.DarkMatterSector
import Gutoe.Baryogenesis
import Gutoe.Inflation

namespace Gutoe.BBN

open Gutoe.FineStructure
open Gutoe.GravityMetric
open Gutoe.Z3Uniqueness
open Gutoe.DarkMatterSector
open Gutoe.Baryogenesis
open Gutoe.Inflation

/-- Primordial He-4 mass-fraction anchor used by the BBN lane. -/
def ypTargetQ : ℚ := 245 / 1000

/-- Primordial D/H anchor used by the BBN lane. -/
def dOverHTargetQ : ℚ := 2547 / 100000000

/-- Primordial ³He/H anchor used by the BBN lane. -/
def he3OverHTargetQ : ℚ := 11 / 1000000

/-- Primordial ⁷Li/H observational anchor used by the BBN lane. -/
def li7OverHObservedQ : ℚ := 16 / 100000000000

/-- Structural `η₁₀` reference from shared inflation/Clifford counts:
    `(12*5)/(4+6) = 6`. -/
noncomputable def eta10Ref : ℝ :=
  ((geometricDarkAmplificationQ : ℝ) * (darkSectorCandidates.card : ℝ)) /
    (((Nat.choose 4 1 : ℕ) : ℝ) + ((Nat.choose 4 2 : ℕ) : ℝ))

/-- Structural deuterium exponent `(6+2)/(4+1) = 8/5`. -/
noncomputable def deuteriumEtaExponent : ℝ :=
  ((((Nat.choose 4 2 : ℕ) : ℝ) + 2) / (((Nat.choose 4 1 : ℕ) : ℝ) + 1))

/-- Structural helium-3 exponent `|SU(2)|/(4+1) = 3/5`. -/
noncomputable def helium3EtaExponent : ℝ :=
  ((magneticTriplet.card : ℝ) / (((Nat.choose 4 1 : ℕ) : ℝ) + 1))

/-- Structural lithium tension amplification `12/4 = 3`. -/
noncomputable def lithium7TensionAmplification : ℝ :=
  (geometricDarkAmplificationQ : ℝ) / ((Nat.choose 4 1 : ℕ) : ℝ)

/-- Structural finite-mode void correction for lithium:
    `(6*11)/(6*11+1) = 66/67`. -/
noncomputable def lithium7VoidCorrection : ℝ :=
  (((Nat.choose 4 2 : ℕ) : ℝ) * (visibleSectorStates.card : ℝ)) /
    ((((Nat.choose 4 2 : ℕ) : ℝ) * (visibleSectorStates.card : ℝ)) + 1)

/-- Structural BBN source from baryogenesis: `η₁₀ = 10¹⁰ η_B`. -/
noncomputable def eta10FromBaryogenesis : ℝ :=
  (10 : ℝ) ^ 10 * etaBaryonStructural

/-- He-4 mass fraction lane used by GRAND-349. -/
noncomputable def primordialHelium4MassFraction (eta10 : ℝ) : ℝ :=
  (ypTargetQ : ℝ) + (lambda_qg / 50) * (eta10 - eta10Ref)

/-- D/H lane used by GRAND-349. -/
noncomputable def primordialDeuteriumRatio (eta10 : ℝ) : ℝ :=
  (dOverHTargetQ : ℝ) * Real.rpow (eta10Ref / eta10) deuteriumEtaExponent

/-- ³He/H lane used by GRAND-349. -/
noncomputable def primordialHelium3Ratio (eta10 : ℝ) : ℝ :=
  (he3OverHTargetQ : ℝ) * Real.rpow (eta10Ref / eta10) helium3EtaExponent

/-- Observed-anchored ⁷Li/H diagnostics lane (legacy comparator). -/
noncomputable def primordialLithium7RatioObservedAnchored (eta10 : ℝ) : ℝ :=
  (li7OverHObservedQ : ℝ) * (eta10 / eta10Ref) ^ (2 : ℕ) * lithium7TensionAmplification

/-- Structural direct Li-7 branch fraction:
    one direct identity-like channel out of total finite states (`1/16`). -/
noncomputable def lithium7DirectChannelFraction : ℝ :=
  (1 : ℝ) / (totalFiniteStateCount : ℝ)

/-- Structural Be-7 mediated Li-7 branch fraction (`15/16`). -/
noncomputable def lithium7Be7ChannelFraction : ℝ :=
  1 - lithium7DirectChannelFraction

/-- Be-7 branch dark suppression from shared occupancy + void factors:
    `(5/16) * (66/67) = 165/536`. -/
noncomputable def lithium7Be7DarkSuppression : ℝ :=
  (darkFractionOfTotalStates : ℝ) * lithium7VoidCorrection

/-- Visible-lane occupancy fraction entering the Li-7 source (`11/16`). -/
noncomputable def lithium7VisibleFraction : ℝ :=
  (visibleSectorStates.card : ℝ) / (totalFiniteStateCount : ℝ)

/-- Structural Li-7 reaction-network gain `(12/4)*(11/16) = 33/16`. -/
noncomputable def lithium7ReactionNetworkGain : ℝ :=
  lithium7TensionAmplification * lithium7VisibleFraction

/-- Absolute Li-7 source from reaction channels:
`(D/H) * (³He/H) * (Yp / Yp_target) * gain`. -/
noncomputable def lithium7ReactionNetworkSource (eta10 : ℝ) : ℝ :=
  primordialDeuteriumRatio eta10
    * primordialHelium3Ratio eta10
    * (primordialHelium4MassFraction eta10 / (ypTargetQ : ℝ))
    * lithium7ReactionNetworkGain

/-- Channel-coupled Li-7 closure factor:
    direct branch unchanged + Be-7 branch dark-suppressed. -/
noncomputable def lithium7ChannelCoupledFactor : ℝ :=
  lithium7DirectChannelFraction + lithium7Be7ChannelFraction * lithium7Be7DarkSuppression

/-- Direct Li-7 component from the absolute reaction source. -/
noncomputable def lithium7DirectComponent (eta10 : ℝ) : ℝ :=
  lithium7ReactionNetworkSource eta10 * lithium7DirectChannelFraction

/-- Be-7 precursor component before dark coupling. -/
noncomputable def lithium7Be7ComponentUnsuppressed (eta10 : ℝ) : ℝ :=
  lithium7ReactionNetworkSource eta10 * lithium7Be7ChannelFraction

/-- Be-7 precursor component after dark coupling. -/
noncomputable def lithium7Be7ComponentDarkCoupled (eta10 : ℝ) : ℝ :=
  lithium7Be7ComponentUnsuppressed eta10 * lithium7Be7DarkSuppression

/-- Channel-coupled Li-7 abundance lane (absolute predictive lane). -/
noncomputable def primordialLithium7RatioChannelCoupled (eta10 : ℝ) : ℝ :=
  lithium7DirectComponent eta10 + lithium7Be7ComponentDarkCoupled eta10

/-- Predictive ⁷Li/H lane used by the active BBN gate. -/
noncomputable def primordialLithium7Ratio (eta10 : ℝ) : ℝ :=
  primordialLithium7RatioChannelCoupled eta10

/-- Corrected lane compatibility alias (routes to channel-coupled predictive lane). -/
noncomputable def primordialLithium7RatioCorrected (eta10 : ℝ) : ℝ :=
  primordialLithium7RatioChannelCoupled eta10

/-- ⁷Li tension ratio lane (`pred/observed`). -/
noncomputable def lithium7TensionRatio (eta10 : ℝ) : ℝ :=
  primordialLithium7Ratio eta10 / (li7OverHObservedQ : ℝ)

/-- Corrected ⁷Li tension ratio lane (`pred_corrected/observed`). -/
noncomputable def lithium7TensionRatioCorrected (eta10 : ℝ) : ℝ :=
  primordialLithium7RatioCorrected eta10 / (li7OverHObservedQ : ℝ)

/-- Channel-coupled Li-7 tension ratio lane (`pred_channel/observed`). -/
noncomputable def lithium7TensionRatioChannelCoupled (eta10 : ℝ) : ℝ :=
  primordialLithium7RatioChannelCoupled eta10 / (li7OverHObservedQ : ℝ)

theorem eta10_ref_eq_6 : eta10Ref = 6 := by
  unfold eta10Ref
  have h41 : (Nat.choose 4 1 : ℕ) = 4 := by native_decide
  have h42 : (Nat.choose 4 2 : ℕ) = 6 := by native_decide
  rw [geometric_dark_amplification_eq, dark_sector_card_eq_five]
  rw [h41, h42]
  norm_num

theorem eta10_ref_pos : 0 < eta10Ref := by
  rw [eta10_ref_eq_6]
  norm_num

theorem deuterium_eta_exponent_eq : deuteriumEtaExponent = (8 : ℝ) / 5 := by
  unfold deuteriumEtaExponent
  have h41 : (Nat.choose 4 1 : ℕ) = 4 := by native_decide
  have h42 : (Nat.choose 4 2 : ℕ) = 6 := by native_decide
  rw [h41, h42]
  norm_num

theorem helium3_eta_exponent_eq : helium3EtaExponent = (3 : ℝ) / 5 := by
  unfold helium3EtaExponent
  have h41 : (Nat.choose 4 1 : ℕ) = 4 := by native_decide
  rw [su2_dim]
  rw [h41]
  norm_num

theorem lithium7_tension_amplification_eq :
    lithium7TensionAmplification = 3 := by
  unfold lithium7TensionAmplification
  have h41 : (Nat.choose 4 1 : ℕ) = 4 := by native_decide
  rw [geometric_dark_amplification_eq]
  rw [h41]
  norm_num

theorem lithium7_void_correction_eq :
    lithium7VoidCorrection = (66 : ℝ) / 67 := by
  unfold lithium7VoidCorrection
  have h42 : (Nat.choose 4 2 : ℕ) = 6 := by native_decide
  rcases visible_dark_state_count_split with ⟨hVis, _, _, _⟩
  rw [h42, hVis]
  norm_num

theorem lithium7_corrected_closure_factor_eq :
    lithium7TensionAmplification * (darkFractionOfTotalStates : ℝ) * lithium7VoidCorrection
      = (495 : ℝ) / 536 := by
  rw [lithium7_tension_amplification_eq, dark_fraction_of_total_states_eq, lithium7_void_correction_eq]
  norm_num

theorem eta10_from_baryogenesis_pos : 0 < eta10FromBaryogenesis := by
  unfold eta10FromBaryogenesis
  have hpow : 0 < (10 : ℝ) ^ 10 := by positivity
  exact mul_pos hpow eta_baryon_structural_pos

theorem primordial_helium4_at_reference :
    primordialHelium4MassFraction eta10Ref = (ypTargetQ : ℝ) := by
  unfold primordialHelium4MassFraction
  ring

theorem primordial_deuterium_at_reference :
    primordialDeuteriumRatio eta10Ref = (dOverHTargetQ : ℝ) := by
  unfold primordialDeuteriumRatio
  have hdiv : eta10Ref / eta10Ref = (1 : ℝ) := by
    field_simp [eta10_ref_pos.ne']
  simp [hdiv]

theorem primordial_helium3_at_reference :
    primordialHelium3Ratio eta10Ref = (he3OverHTargetQ : ℝ) := by
  unfold primordialHelium3Ratio
  have hdiv : eta10Ref / eta10Ref = (1 : ℝ) := by
    field_simp [eta10_ref_pos.ne']
  simp [hdiv]

theorem lithium7_observed_anchored_ratio_at_reference :
    primordialLithium7RatioObservedAnchored eta10Ref = (li7OverHObservedQ : ℝ) * 3 := by
  unfold primordialLithium7RatioObservedAnchored
  have hdiv : eta10Ref / eta10Ref = (1 : ℝ) := by
    field_simp [eta10_ref_pos.ne']
  rw [hdiv]
  rw [one_pow, lithium7_tension_amplification_eq]
  ring

theorem lithium7_direct_channel_fraction_eq :
    lithium7DirectChannelFraction = (1 : ℝ) / 16 := by
  unfold lithium7DirectChannelFraction
  rw [total_finite_state_count_eq]
  norm_num

theorem lithium7_be7_channel_fraction_eq :
    lithium7Be7ChannelFraction = (15 : ℝ) / 16 := by
  unfold lithium7Be7ChannelFraction
  rw [lithium7_direct_channel_fraction_eq]
  norm_num

theorem lithium7_branch_fractions_sum_unity :
    lithium7DirectChannelFraction + lithium7Be7ChannelFraction = 1 := by
  rw [lithium7_direct_channel_fraction_eq, lithium7_be7_channel_fraction_eq]
  norm_num

theorem lithium7_be7_dark_suppression_eq :
    lithium7Be7DarkSuppression = (165 : ℝ) / 536 := by
  unfold lithium7Be7DarkSuppression
  rw [dark_fraction_of_total_states_eq, lithium7_void_correction_eq]
  norm_num

theorem lithium7_channel_coupled_factor_eq :
    lithium7ChannelCoupledFactor = (3011 : ℝ) / 8576 := by
  unfold lithium7ChannelCoupledFactor
  rw [lithium7_direct_channel_fraction_eq, lithium7_be7_channel_fraction_eq, lithium7_be7_dark_suppression_eq]
  norm_num

theorem lithium7_visible_fraction_eq :
    lithium7VisibleFraction = (11 : ℝ) / 16 := by
  unfold lithium7VisibleFraction
  rcases visible_dark_state_count_split with ⟨hVis, _, _, _⟩
  rw [hVis, total_finite_state_count_eq]
  norm_num

theorem lithium7_reaction_network_gain_eq :
    lithium7ReactionNetworkGain = (33 : ℝ) / 16 := by
  unfold lithium7ReactionNetworkGain
  rw [lithium7_tension_amplification_eq, lithium7_visible_fraction_eq]
  norm_num

theorem lithium7_reaction_network_source_at_reference :
    lithium7ReactionNetworkSource eta10Ref = (924561 : ℝ) / 1600000000000000 := by
  unfold lithium7ReactionNetworkSource
  rw [primordial_deuterium_at_reference, primordial_helium3_at_reference,
    primordial_helium4_at_reference, lithium7_reaction_network_gain_eq]
  norm_num [dOverHTargetQ, he3OverHTargetQ, ypTargetQ]

theorem lithium7_channel_coupled_factorization (eta10 : ℝ) :
    primordialLithium7RatioChannelCoupled eta10 =
      lithium7ReactionNetworkSource eta10 * lithium7ChannelCoupledFactor := by
  unfold primordialLithium7RatioChannelCoupled lithium7DirectComponent
    lithium7Be7ComponentDarkCoupled lithium7Be7ComponentUnsuppressed
    lithium7ChannelCoupledFactor
  ring

theorem lithium7_ratio_at_reference :
    primordialLithium7Ratio eta10Ref = (2783853171 : ℝ) / 13721600000000000000 := by
  unfold primordialLithium7Ratio
  rw [lithium7_channel_coupled_factorization, lithium7_reaction_network_source_at_reference,
    lithium7_channel_coupled_factor_eq]
  norm_num

theorem lithium7_tension_ratio_at_reference :
    lithium7TensionRatio eta10Ref = (2783853171 : ℝ) / 2195456000 := by
  unfold lithium7TensionRatio
  rw [lithium7_ratio_at_reference]
  norm_num [li7OverHObservedQ]

theorem lithium7_tension_ratio_reference_window :
    (4 / 5 : ℝ) ≤ lithium7TensionRatio eta10Ref ∧
      lithium7TensionRatio eta10Ref ≤ (7 / 5 : ℝ) := by
  rw [lithium7_tension_ratio_at_reference]
  constructor <;> norm_num

theorem lithium7_corrected_tension_ratio_at_reference :
    lithium7TensionRatioCorrected eta10Ref = (2783853171 : ℝ) / 2195456000 := by
  simpa [lithium7TensionRatioCorrected, primordialLithium7RatioCorrected, primordialLithium7Ratio]
    using lithium7_tension_ratio_at_reference

theorem lithium7_corrected_ratio_reference_window :
    (4 / 5 : ℝ) ≤ lithium7TensionRatioCorrected eta10Ref ∧
      lithium7TensionRatioCorrected eta10Ref ≤ (7 / 5 : ℝ) := by
  rw [lithium7_corrected_tension_ratio_at_reference]
  constructor <;> norm_num

theorem lithium7_channel_coupled_tension_ratio_at_reference :
    lithium7TensionRatioChannelCoupled eta10Ref = (2783853171 : ℝ) / 2195456000 := by
  unfold lithium7TensionRatioChannelCoupled primordialLithium7RatioChannelCoupled
    lithium7DirectComponent lithium7Be7ComponentDarkCoupled
    lithium7Be7ComponentUnsuppressed lithium7ReactionNetworkSource
  rw [primordial_deuterium_at_reference, primordial_helium3_at_reference,
    primordial_helium4_at_reference, lithium7_direct_channel_fraction_eq,
    lithium7_be7_channel_fraction_eq, lithium7_be7_dark_suppression_eq,
    lithium7_reaction_network_gain_eq]
  norm_num [dOverHTargetQ, he3OverHTargetQ, ypTargetQ, li7OverHObservedQ]

theorem lithium7_channel_coupled_ratio_reference_window :
    (4 / 5 : ℝ) ≤ lithium7TensionRatioChannelCoupled eta10Ref ∧
      lithium7TensionRatioChannelCoupled eta10Ref ≤ (7 / 5 : ℝ) := by
  rw [lithium7_channel_coupled_tension_ratio_at_reference]
  constructor <;> norm_num

end Gutoe.BBN
