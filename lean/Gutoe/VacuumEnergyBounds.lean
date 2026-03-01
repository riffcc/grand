import Mathlib
import Gutoe.EWSBHiggs
import Gutoe.DimensionalStructure
import Gutoe.Z3Uniqueness
import Gutoe.VoidRearFace

/-!
GUTOE — Vacuum Energy Bounds (Research Lane)

This lane formalizes the constraint skeleton for local negative-energy proposals:

1. Structural EW barrier proxy from Cl(1,3)-derived constants (`13/100`, `40/153`, `3/16`).
2. Ford-Roman-style quantum inequality as a proposition and its direct consequences.
3. Casimir `a⁻⁴` scaling as a hard geometric suppression mechanism.

This module is intentionally conservative: it proves bound logic, not an engine.
-/

namespace Gutoe.VacuumEnergyBounds

open Gutoe.EWSBHiggs
open Gutoe.DimensionalStructure
open Gutoe.Z3Uniqueness
open Gutoe.VoidRearFace

/-- Structural void fraction from the Cl(1,3) split:
    `f_void = |Z₃| / 2⁴ = 3/16`. -/
def voidFractionQ : ℚ := (magneticTriplet.card : ℚ) / (2 ^ 4 : ℚ)

theorem void_fraction_eq_3_16 : voidFractionQ = (3 : ℚ) / 16 := by
  unfold voidFractionQ
  rw [su2_dim]
  norm_num

/-- Structural electroweak barrier proxy (dimensionless in this lane's unit convention):
    `ΔV_proxy = λ_H * (v/mp)^4 / 4`. -/
def ewsbBarrierProxyQ : ℚ := higgsQuartic * (vevOverProton ^ 4) / 4

theorem ewsb_barrier_proxy_closed_form :
    ewsbBarrierProxyQ = (13 : ℚ) / 100 * (((40 : ℚ) / 153) ^ 4) / 4 := by
  unfold ewsbBarrierProxyQ
  rw [higgs_quartic_eq_13_100, vev_over_proton_eq_40_153]

theorem ewsb_barrier_proxy_pos : 0 < ewsbBarrierProxyQ := by
  rw [ewsb_barrier_proxy_closed_form]
  positivity

/-- Ford-Roman-style quantum inequality schema:
    higher negative energy magnitude forces shorter allowable duration. -/
def FordRomanBound (rhoNeg tau qeiK : ℝ) : Prop :=
  |rhoNeg| * tau ^ 4 ≤ qeiK

/-- Any attempt to sustain a minimum negative-energy magnitude for at least a minimum
duration is impossible once it exceeds the Ford-Roman budget. -/
theorem ford_roman_no_durable_window
    {rhoNeg tau qeiK rhoTarget tauMin : ℝ}
    (hFR : FordRomanBound rhoNeg tau qeiK)
    (hRhoTargetPos : 0 < rhoTarget)
    (hTauMinPos : 0 < tauMin)
    (hMag : rhoTarget ≤ |rhoNeg|)
    (hTau : tauMin ≤ tau)
    (hBudgetFail : qeiK < rhoTarget * tauMin ^ 4) :
    False := by
  have hTauPow : tauMin ^ 4 ≤ tau ^ 4 := by
    exact pow_le_pow_left₀ hTauMinPos.le hTau 4
  have hTargetStep1 : rhoTarget * tauMin ^ 4 ≤ rhoTarget * tau ^ 4 := by
    exact mul_le_mul_of_nonneg_left hTauPow hRhoTargetPos.le
  have hTauPowNonneg : 0 ≤ tau ^ 4 := by positivity
  have hTargetStep2 : rhoTarget * tau ^ 4 ≤ |rhoNeg| * tau ^ 4 := by
    exact mul_le_mul_of_nonneg_right hMag hTauPowNonneg
  have hTargetLeQei : rhoTarget * tauMin ^ 4 ≤ qeiK := by
    exact le_trans hTargetStep1 (le_trans hTargetStep2 hFR)
  exact not_le_of_gt hBudgetFail hTargetLeQei

/-- Idealized Casimir-magnitude model with geometric `a⁻⁴` suppression. -/
noncomputable def casimirMagnitude (kappa a : ℝ) : ℝ :=
  kappa / a ^ 4

/-- For nonnegative prefactor `κ`, enlarging the gap weakens Casimir magnitude. -/
theorem casimir_magnitude_antitone_in_gap
    {kappa aMin a : ℝ}
    (hKappa : 0 ≤ kappa)
    (hAMinPos : 0 < aMin)
    (hGap : aMin ≤ a) :
    casimirMagnitude kappa a ≤ casimirMagnitude kappa aMin := by
  unfold casimirMagnitude
  have hAPos : 0 < a := lt_of_lt_of_le hAMinPos hGap
  have hPow : aMin ^ 4 ≤ a ^ 4 := by
    exact pow_le_pow_left₀ hAMinPos.le hGap 4
  have hInv :
      (1 / a ^ 4) ≤ (1 / aMin ^ 4) := by
    exact one_div_le_one_div_of_le (pow_pos hAMinPos 4) hPow
  have hMul :
      kappa * (1 / a ^ 4) ≤ kappa * (1 / aMin ^ 4) := by
    exact mul_le_mul_of_nonneg_left hInv hKappa
  simpa [div_eq_mul_inv, mul_comm, mul_left_comm, mul_assoc] using hMul

/-- If a minimum achievable gap already undershoots a target density, any larger gap
also undershoots: geometric no-go from `a⁻⁴` scaling. -/
theorem casimir_no_go_from_min_gap
    {kappa aMin a rhoTarget : ℝ}
    (hKappa : 0 ≤ kappa)
    (hAMinPos : 0 < aMin)
    (hGap : aMin ≤ a)
    (hUnder : casimirMagnitude kappa aMin < rhoTarget) :
    casimirMagnitude kappa a < rhoTarget := by
  have hLe : casimirMagnitude kappa a ≤ casimirMagnitude kappa aMin :=
    casimir_magnitude_antitone_in_gap hKappa hAMinPos hGap
  exact lt_of_le_of_lt hLe hUnder

/-- Higgs-orientation gradient-energy proxy for a wall-like configuration:
    `E_grad = (f_c * v^2 * A * (Δθ)^2) / (2 * L)`.
    This captures the ungauged orientation-mode scaling skeleton. -/
noncomputable def higgsOrientationGradientEnergy
    (fc v area dTheta thickness : ℝ) : ℝ :=
  (fc * v ^ 2 * area * dTheta ^ 2) / (2 * thickness)

/-- Equivalent inverse-thickness form: `E_grad = K * (1/L)`. -/
theorem higgs_orientation_gradient_energy_inverse_thickness_form
    (fc v area dTheta thickness : ℝ) (hL : thickness ≠ 0) :
    higgsOrientationGradientEnergy fc v area dTheta thickness =
      ((fc * v ^ 2 * area * dTheta ^ 2) / 2) * (1 / thickness) := by
  unfold higgsOrientationGradientEnergy
  field_simp [hL]

/-- For fixed `fc, v, Δθ, L`, gradient energy is linear in wall area. -/
theorem higgs_orientation_gradient_energy_linear_in_area
    (fc v area dTheta thickness : ℝ) (hL : thickness ≠ 0) :
    higgsOrientationGradientEnergy fc v area dTheta thickness =
      area * ((fc * v ^ 2 * dTheta ^ 2) / (2 * thickness)) := by
  unfold higgsOrientationGradientEnergy
  field_simp [hL]

/-- For fixed prefactor and positive thicknesses, gradient energy is antitone in
wall thickness (`L↑ => E↓`). -/
theorem higgs_orientation_gradient_energy_antitone_in_thickness
    {fc v area dTheta L1 L2 : ℝ}
    (hfc : 0 ≤ fc)
    (harea : 0 ≤ area)
    (hL1 : 0 < L1)
    (hL2 : 0 < L2)
    (hL : L1 ≤ L2) :
    higgsOrientationGradientEnergy fc v area dTheta L2 ≤
      higgsOrientationGradientEnergy fc v area dTheta L1 := by
  have hk : 0 ≤ ((fc * v ^ 2 * area * dTheta ^ 2) / 2) := by
    have hv2 : 0 ≤ v ^ 2 := by positivity
    have hd2 : 0 ≤ dTheta ^ 2 := by positivity
    have hnum : 0 ≤ fc * v ^ 2 * area * dTheta ^ 2 := by
      exact mul_nonneg (mul_nonneg (mul_nonneg hfc hv2) harea) hd2
    exact div_nonneg hnum (by norm_num)
  have hInv : 1 / L2 ≤ 1 / L1 := one_div_le_one_div_of_le hL1 hL
  have hL1ne : L1 ≠ 0 := ne_of_gt hL1
  have hL2ne : L2 ≠ 0 := ne_of_gt hL2
  calc
    higgsOrientationGradientEnergy fc v area dTheta L2
        = ((fc * v ^ 2 * area * dTheta ^ 2) / 2) * (1 / L2) := by
            exact higgs_orientation_gradient_energy_inverse_thickness_form fc v area dTheta L2 hL2ne
    _ ≤ ((fc * v ^ 2 * area * dTheta ^ 2) / 2) * (1 / L1) :=
          mul_le_mul_of_nonneg_left hInv hk
    _ = higgsOrientationGradientEnergy fc v area dTheta L1 := by
          symm
          exact higgs_orientation_gradient_energy_inverse_thickness_form fc v area dTheta L1 hL1ne

/-- Structural specialization with the Cl(1,3) void fraction `f_c = 3/16`. -/
noncomputable def higgsOrientationGradientEnergyStructural
    (v area dTheta thickness : ℝ) : ℝ :=
  higgsOrientationGradientEnergy ((3 : ℝ) / 16) v area dTheta thickness

theorem higgs_orientation_gradient_energy_structural_eq
    (v area dTheta thickness : ℝ) :
    higgsOrientationGradientEnergyStructural v area dTheta thickness =
      higgsOrientationGradientEnergy ((3 : ℝ) / 16) v area dTheta thickness := by
  rfl

/-- Rear-face suppression factor imported from `VoidRearFace` as a real scalar. -/
noncomputable def rearFaceSuppressionR : ℝ := (rearCostFactor : ℝ)

theorem rear_face_suppression_eq_one_tenth :
    rearFaceSuppressionR = (1 : ℝ) / 10 := by
  norm_num [rearFaceSuppressionR, rear_cost_factor_eq_one_tenth]

/-- Relativistic wall-Lorentz factor under areal drive `E/A` and wall tension `σ`:
    `γ = 1 + (E/A)/σ`.
    (This lane tracks kinematics under the wall-surfing hypothesis.) -/
noncomputable def wallLorentzGamma (arealDrive sigma : ℝ) : ℝ :=
  1 + arealDrive / sigma

/-- Rear-face wall tension under the imported `1/10` suppression. -/
noncomputable def rearFaceWallTension (sigmaFront : ℝ) : ℝ :=
  rearFaceSuppressionR * sigmaFront

theorem rear_face_wall_tension_eq_one_tenth
    (sigmaFront : ℝ) :
    rearFaceWallTension sigmaFront = sigmaFront / 10 := by
  unfold rearFaceWallTension
  rw [rear_face_suppression_eq_one_tenth]
  ring

/-- Closed-form rear-face Lorentz factor:
    `γ_rear = 1 + 10 * (E/A)/σ_front`. -/
theorem wall_gamma_rear_face_closed_form
    (arealDrive sigmaFront : ℝ) (hσ : sigmaFront ≠ 0) :
    wallLorentzGamma arealDrive (rearFaceWallTension sigmaFront) =
      1 + 10 * arealDrive / sigmaFront := by
  unfold wallLorentzGamma
  rw [rear_face_wall_tension_eq_one_tenth]
  field_simp [hσ]

/-- Under nonnegative drive and positive front tension, rear-face suppression gives
an equal-or-higher Lorentz factor than front-face propagation. -/
theorem wall_gamma_rear_ge_front
    {arealDrive sigmaFront : ℝ}
    (hE : 0 ≤ arealDrive)
    (hσ : 0 < sigmaFront) :
    wallLorentzGamma arealDrive (rearFaceWallTension sigmaFront) ≥
      wallLorentzGamma arealDrive sigmaFront := by
  have hσ0 : sigmaFront ≠ 0 := ne_of_gt hσ
  rw [wall_gamma_rear_face_closed_form arealDrive sigmaFront hσ0]
  unfold wallLorentzGamma
  have hfrac : 0 ≤ arealDrive / sigmaFront := div_nonneg hE hσ.le
  have hscale : arealDrive / sigmaFront ≤ 10 * (arealDrive / sigmaFront) := by
    nlinarith [hfrac]
  have hten : 10 * arealDrive / sigmaFront = 10 * (arealDrive / sigmaFront) := by
    field_simp [hσ0]
  have hsum : 1 + 10 * arealDrive / sigmaFront ≥ 1 + arealDrive / sigmaFront := by
    calc
      1 + 10 * arealDrive / sigmaFront = 1 + 10 * (arealDrive / sigmaFront) := by rw [hten]
      _ ≥ 1 + arealDrive / sigmaFront := by linarith [hscale]
  exact hsum

end Gutoe.VacuumEnergyBounds
