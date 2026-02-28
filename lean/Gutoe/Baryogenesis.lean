import Mathlib
import Gutoe.FineStructure
import Gutoe.FlavorMixing
import Gutoe.DarkMatterSector
import Gutoe.GravityMetric
import Gutoe.Z3Uniqueness

namespace Gutoe.Baryogenesis

open Gutoe.FineStructure
open Gutoe.FlavorMixing
open Gutoe.DarkMatterSector
open Gutoe.GravityMetric
open Gutoe.Z3Uniqueness

/-- Leading-order electromagnetic suppression from the shared `α⁻¹ = 137` theorem. -/
noncomputable def alphaLeadingOrder : ℝ :=
  1 / (alphaInverse 4 : ℝ)

/-- GRAND-348 micro-mode count from shared `|SU(2)| = 3`: `2 * 3^5 = 486`. -/
def baryoMicroModeCount : ℕ :=
  2 * magneticTriplet.card ^ 5

/-- GRAND-348 finite-mode rescale: `N/(N-1)` with `N = 486`. -/
noncomputable def baryoMicroFiniteRescale : ℝ :=
  (baryoMicroModeCount : ℝ) / ((baryoMicroModeCount - 1 : ℕ) : ℝ)

/-- Departure-from-equilibrium survival factor from shared lattice primitives:
    `(1 - λ_QG) * (486/485)`. -/
noncomputable def baryoNonequilibriumSurvival : ℝ :=
  (1 - lambda_qg) * baryoMicroFiniteRescale

/-- Structural baryogenesis prefactor:
    `α² * (5/11) * f_neq`. -/
noncomputable def baryogenesisPrefactor : ℝ :=
  alphaLeadingOrder ^ 2 * (darkToVisibleCountRatio : ℝ) * baryoNonequilibriumSurvival

/-- Quantitative structural prediction for baryon-to-photon ratio. -/
noncomputable def etaBaryonStructural : ℝ :=
  jarlskog ckmSin12 ckmSin23 ckmSin13 ckmDelta * baryogenesisPrefactor

/-- PMNS θ23 void scalar imported from the flavor α² correction lane:
    `sin²θ23_direct - sin²θ23_corrected`. -/
noncomputable def pmnsTheta23VoidScalar : ℝ :=
  (pmnsSin23SqDirectQ - pmnsSin23SqCorrectedQ : ℚ)

/-- Default PMNS-linked leptogenesis multiplier used by runtime:
    `1 + (sin²θ23_direct - sin²θ23_corrected)`. -/
noncomputable def leptogenesisPmnsMultiplier : ℝ :=
  1 + pmnsTheta23VoidScalar

/-- PMNS-corrected baryogenesis prediction for the default structural gain. -/
noncomputable def etaBaryonStructuralWithPmns : ℝ :=
  etaBaryonStructural * leptogenesisPmnsMultiplier

theorem alpha_leading_order_eq :
    alphaLeadingOrder = (1 : ℝ) / 137 := by
  unfold alphaLeadingOrder
  norm_num [alpha_inverse_d4]

theorem baryo_micro_mode_count_eq_486 :
    baryoMicroModeCount = 486 := by
  unfold baryoMicroModeCount
  have hs : magneticTriplet.card = 3 := su2_dim
  norm_num [hs]

theorem baryo_micro_finite_rescale_eq :
    baryoMicroFiniteRescale = (486 : ℝ) / 485 := by
  unfold baryoMicroFiniteRescale
  norm_num [baryo_micro_mode_count_eq_486]

theorem baryo_nonequilibrium_survival_eq :
    baryoNonequilibriumSurvival = ((11 : ℝ) / 12) * ((486 : ℝ) / 485) := by
  norm_num [baryoNonequilibriumSurvival, baryoMicroFiniteRescale, baryoMicroModeCount, su2_dim, lambda_qg]

theorem baryo_nonequilibrium_survival_bounds :
    0 < baryoNonequilibriumSurvival ∧ baryoNonequilibriumSurvival < 1 := by
  rw [baryo_nonequilibrium_survival_eq]
  constructor
  · positivity
  · norm_num

theorem dark_ratio_real_eq :
    (darkToVisibleCountRatio : ℝ) = (5 : ℝ) / 11 := by
  norm_num [dark_to_visible_count_ratio_eq]

theorem baryogenesis_prefactor_eq :
    baryogenesisPrefactor
      = ((1 : ℝ) / 137) ^ 2 * ((5 : ℝ) / 11) *
        (((11 : ℝ) / 12) * ((486 : ℝ) / 485)) := by
  unfold baryogenesisPrefactor
  rw [alpha_leading_order_eq, dark_ratio_real_eq, baryo_nonequilibrium_survival_eq]

theorem baryogenesis_prefactor_pos :
    0 < baryogenesisPrefactor := by
  rw [baryogenesis_prefactor_eq]
  positivity

theorem eta_baryon_structural_pos :
    0 < etaBaryonStructural := by
  unfold etaBaryonStructural
  exact mul_pos ckm_jarlskog_positive baryogenesis_prefactor_pos

theorem eta_baryon_structural_from_shared_primitives :
    etaBaryonStructural
      = jarlskog ckmSin12 ckmSin23 ckmSin13 ckmDelta *
        (((1 : ℝ) / 137) ^ 2 * ((5 : ℝ) / 11) *
         (((11 : ℝ) / 12) * ((486 : ℝ) / 485))) := by
  unfold etaBaryonStructural
  rw [baryogenesis_prefactor_eq]

theorem pmns_theta23_void_scalar_eq :
    pmnsTheta23VoidScalar = (1 : ℝ) / 548 := by
  norm_num [pmnsTheta23VoidScalar, pmns_theta23_void_term]

theorem leptogenesis_pmns_multiplier_eq :
    leptogenesisPmnsMultiplier = (549 : ℝ) / 548 := by
  rw [leptogenesisPmnsMultiplier, pmns_theta23_void_scalar_eq]
  norm_num

theorem eta_baryon_structural_with_pmns_ratio :
    etaBaryonStructuralWithPmns = etaBaryonStructural * ((549 : ℝ) / 548) := by
  unfold etaBaryonStructuralWithPmns
  rw [leptogenesis_pmns_multiplier_eq]

theorem eta_baryon_structural_with_pmns_pos :
    0 < etaBaryonStructuralWithPmns := by
  unfold etaBaryonStructuralWithPmns
  have hmul : 0 < leptogenesisPmnsMultiplier := by
    rw [leptogenesis_pmns_multiplier_eq]
    positivity
  exact mul_pos eta_baryon_structural_pos hmul

end Gutoe.Baryogenesis
