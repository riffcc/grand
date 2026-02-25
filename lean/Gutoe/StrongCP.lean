/-
 * GUTOE — Strong CP Structural Gate (GRAND-125)
 *
 * This module formalizes the current Cl(1,3) structural claim:
 * - CP-odd support tracks the grade-2 rotation/boost imbalance.
 * - Cl(1,3) fixes a 3+3 Lorentz bivector split, so the imbalance is zero.
 * - Therefore the structural θ_QCD proxy is zero.
 *
 * Scope note:
 * This is a structural theorem chain. A full nonperturbative QCD vacuum-angle
 * derivation remains a separate milestone.
 -/

import Mathlib
import Gutoe.LorentzInvariance

namespace Gutoe.StrongCP

open Gutoe.DimensionalStructure
open Gutoe.Z3Uniqueness
open Gutoe.LorentzInvariance

/-- CP-odd sector source proxy: bivector rotation/boost imbalance in Cl(1,3). -/
def cpOddSectorImbalance : ℤ :=
  (magneticTriplet.card : ℤ) - (emTriplet.card : ℤ)

/-- Cl(1,3) Lorentz split forces zero CP-odd sector imbalance. -/
theorem cp_odd_sector_imbalance_zero : cpOddSectorImbalance = 0 := by
  unfold cpOddSectorImbalance
  rcases lorentz_algebra_decomposition with ⟨_, _, hmag, hem⟩
  rw [hmag, hem]
  norm_num

/-- Structural θ_QCD proxy carried by the CP-odd imbalance source. -/
def thetaQcdStructural : ℝ := (cpOddSectorImbalance : ℝ)

/-- The structural Strong-CP theorem in this model: θ_QCD = 0. -/
theorem theta_qcd_structural_zero : thetaQcdStructural = 0 := by
  unfold thetaQcdStructural
  rw [cp_odd_sector_imbalance_zero]
  norm_num

/-- Standard bridge coefficient used for neutron EDM estimate:
    |d_n| ≈ 2.4e-16 * |θ_QCD| e·cm. -/
def neutronEdmBridgeCoeff : ℝ := 2.4e-16

/-- EDM bridge map from θ_QCD to a neutron EDM estimate. -/
def neutronEdmFromTheta (thetaQcd : ℝ) : ℝ :=
  neutronEdmBridgeCoeff * |thetaQcd|

/-- If θ_QCD is structurally zero, EDM estimate is structurally zero. -/
theorem neutron_edm_from_structural_theta_zero :
    neutronEdmFromTheta thetaQcdStructural = 0 := by
  unfold neutronEdmFromTheta
  rw [theta_qcd_structural_zero]
  norm_num [neutronEdmBridgeCoeff]

/-- Structural EDM estimate stays below the catalog gate bound. -/
theorem neutron_edm_structural_within_catalog_bound :
    neutronEdmFromTheta thetaQcdStructural ≤ (1e-26 : ℝ) := by
  rw [neutron_edm_from_structural_theta_zero]
  norm_num

/-- Bridge theorem: if an effective θ is proportional to the structural source,
    Cl(1,3) balance forces that effective θ to vanish. -/
theorem theta_zero_of_proportional_cp_odd_source
    (thetaEff k : ℝ)
    (hprop : thetaEff = k * thetaQcdStructural) :
    thetaEff = 0 := by
  rw [hprop, theta_qcd_structural_zero]
  ring

end Gutoe.StrongCP
