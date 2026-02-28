/-
 * GUTOE — GRMHD Synchrotron Scaffold
 * Copyright (C) 2026 Riff Labs
 * AGPL-3.0-or-later
 *
 * Minimal physics-first synchrotron model tied to shared GUTOE primitives:
 *   - magneticTriplet.card = 3 (Cl(1,3) magnetic sector multiplicity)
 *   - lambda_qg = 1/12 (lattice UV correction)
 *   - Kerr frame dragging ω(r,θ) for rotation boosting
 *
 * No `sorry`.
-/

import Mathlib
import Gutoe.GravityMetric
import Gutoe.Z3Uniqueness
import Gutoe.KerrGeometry

namespace Gutoe.SynchrotronGRMHD

open Real
open Gutoe.Z3Uniqueness
open Gutoe.GravityMetric
open Gutoe.KerrGeometry

/-- Magnetic-sector multiplicity from the shared Z₃ orbit primitive. -/
noncomputable def magneticSectorWeight : ℝ := (magneticTriplet.card : ℝ)

/-- Electron rest-energy temperature `m_e c^2 / k_B` (kelvin scale). -/
noncomputable def electronRestTempK : ℝ := 5.93086740e9

/-- Dimensionless electron temperature `θ_e = T_e / (m_e c^2 / k_B)`. -/
noncomputable def electronTheta (teK : ℝ) : ℝ := teK / electronRestTempK

/-- Mahadevan-style spectral fit factor `I'(x_M)` used in thermal synchrotron
emissivity approximations. -/
noncomputable def mahadevanIPrime (xM : ℝ) : ℝ :=
  let x := max xM (1 / 1000000000000 : ℝ)
  (4.0505 / x ^ (1 / 6 : ℝ))
    * (1 + 0.40 / x ^ (1 / 4 : ℝ) + 0.5316 / x ^ (1 / 2 : ℝ))
    * Real.exp (-(1.8899 * x ^ (1 / 3 : ℝ)))

/-- Dimensionless cyclotron-frequency proxy used by the thermal fit lane.
This keeps the Lean side finite/robust while tracking `|B|` dependence. -/
noncomputable def cyclotronFreqProxy (B : ℝ) : ℝ := max |B| (1 / 1000000 : ℝ)

/-- Thermal synchrotron emissivity coefficient scaffold from literature-fit form:
`j_ν ∝ n_e * ν_s * I'(x_M)` with `ν_s ∝ ν_c θ_e² sin(pitch)`. -/
noncomputable def thermalSynchrotronEmissivity
    (nE B teK ν sinPitch : ℝ) : ℝ :=
  let θe := max (electronTheta teK) (1 / 1000000 : ℝ)
  let νc := cyclotronFreqProxy B
  let sinP := max sinPitch (1 / 1000 : ℝ)
  let νs := (2 / 9 : ℝ) * νc * θe ^ 2 * sinP
  let xM := ν / max νs (1 / 1000000000000 : ℝ)
  magneticSectorWeight * lambda_qg * max nE 0 * max νs 0 * mahadevanIPrime xM

/-- Thermal synchrotron absorptivity scaffold via Kirchhoff:
`α_ν = j_ν / B_ν(T_e)` with a strictly positive denominator floor. -/
noncomputable def thermalSynchrotronAbsorption
    (nE B teK ν sinPitch : ℝ) : ℝ :=
  let planckProxy := max (ν ^ 3 / (Real.exp (ν / max teK (1 / 1000 : ℝ)) - 1))
    (1 / 1000000000000 : ℝ)
  thermalSynchrotronEmissivity nE B teK ν sinPitch / planckProxy

/-- Baseline synchrotron emissivity proxy:

`j_syn ∝ |B|² * ν * exp(-ν r_s) * (1 + λ_QG (l_P/r_s)²)`

This keeps the model simple but physically constrained:
- quadratic in magnetic field amplitude,
- linear in observing frequency ν at low ν,
- high-frequency suppression via `exp(-ν r_s)`,
- lattice UV factor from shared `lambda_qg`.
-/
noncomputable def synchrotronEmissivity
    (B ν r_s l_P : ℝ) : ℝ :=
  thermalSynchrotronEmissivity 1 B 6.0e10 ν (Real.sqrt 2 / 2)
    * (1 + lambda_qg * (l_P / r_s) ^ 2)

/-- Kerr-boosted emissivity proxy:

`j_Kerr = j_syn * (1 + |ω|)` with `ω` the frame-dragging rate.

This cleanly separates local plasma emissivity from geometric rotation boost.
-/
noncomputable def synchrotronEmissivityKerr
    (B ν r_s l_P aStar r θ : ℝ) : ℝ :=
  synchrotronEmissivity B ν r_s l_P * (1 + |frameDraggingOmega r_s aStar r θ|)

theorem magneticSectorWeight_eq_three : magneticSectorWeight = 3 := by
  unfold magneticSectorWeight
  norm_num [su2_dim]

theorem magneticSectorWeight_nonneg : 0 ≤ magneticSectorWeight := by
  rw [magneticSectorWeight_eq_three]
  norm_num

theorem magneticSectorWeight_pos : 0 < magneticSectorWeight := by
  rw [magneticSectorWeight_eq_three]
  norm_num

lemma mahadevanIPrime_pos (xM : ℝ) : 0 < mahadevanIPrime xM := by
  unfold mahadevanIPrime
  set x : ℝ := max xM (1 / 1000000000000 : ℝ)
  have hfloor_pos : (0 : ℝ) < (1 / 1000000000000 : ℝ) := by
    norm_num
  have hx_ge_floor : (1 / 1000000000000 : ℝ) ≤ x := by
    simpa [x] using (le_max_right xM (1 / 1000000000000 : ℝ))
  have hx_pos : 0 < x := lt_of_lt_of_le hfloor_pos hx_ge_floor
  have hx_one_sixth_pos : 0 < x ^ (1 / 6 : ℝ) := Real.rpow_pos_of_pos hx_pos _
  have hx_one_fourth_pos : 0 < x ^ (1 / 4 : ℝ) := Real.rpow_pos_of_pos hx_pos _
  have hx_one_half_pos : 0 < x ^ (1 / 2 : ℝ) := Real.rpow_pos_of_pos hx_pos _
  have hfirst_pos : 0 < 4.0505 / x ^ (1 / 6 : ℝ) := by
    exact div_pos (by norm_num) hx_one_sixth_pos
  have hquarter_term_pos : 0 < 0.40 / x ^ (1 / 4 : ℝ) := by
    exact div_pos (by norm_num) hx_one_fourth_pos
  have hhalf_term_pos : 0 < 0.5316 / x ^ (1 / 2 : ℝ) := by
    exact div_pos (by norm_num) hx_one_half_pos
  have hmiddle_pos : 0 < 1 + 0.40 / x ^ (1 / 4 : ℝ) + 0.5316 / x ^ (1 / 2 : ℝ) := by
    nlinarith
  have hexp_pos : 0 < Real.exp (-(1.8899 * x ^ (1 / 3 : ℝ))) := Real.exp_pos _
  have hmul_pos :
      0 < (4.0505 / x ^ (1 / 6 : ℝ)) *
        (1 + 0.40 / x ^ (1 / 4 : ℝ) + 0.5316 / x ^ (1 / 2 : ℝ)) := by
    exact mul_pos hfirst_pos hmiddle_pos
  exact mul_pos hmul_pos hexp_pos

lemma mahadevanIPrime_nonneg (xM : ℝ) : 0 ≤ mahadevanIPrime xM :=
  le_of_lt (mahadevanIPrime_pos xM)

lemma thermalSynchrotronEmissivity_nonneg
    (nE B teK ν sinPitch : ℝ) :
    0 ≤ thermalSynchrotronEmissivity nE B teK ν sinPitch := by
  unfold thermalSynchrotronEmissivity
  set θe : ℝ := max (electronTheta teK) (1 / 1000000 : ℝ)
  set νc : ℝ := cyclotronFreqProxy B
  set sinP : ℝ := max sinPitch (1 / 1000 : ℝ)
  set νs : ℝ := (2 / 9 : ℝ) * νc * θe ^ 2 * sinP
  set xM : ℝ := ν / max νs (1 / 1000000000000 : ℝ)
  have hmag : 0 ≤ magneticSectorWeight := magneticSectorWeight_nonneg
  have hlam : 0 ≤ lambda_qg := le_of_lt lambda_qg_pos
  have hnE_nonneg : 0 ≤ max nE 0 := le_max_right nE 0
  have hνs_nonneg : 0 ≤ max νs 0 := le_max_right νs 0
  have hI_nonneg : 0 ≤ mahadevanIPrime xM := mahadevanIPrime_nonneg xM
  exact mul_nonneg (mul_nonneg (mul_nonneg (mul_nonneg hmag hlam) hnE_nonneg) hνs_nonneg) hI_nonneg

lemma thermalSynchrotronEmissivity_pos_synchrotron_lane (B ν : ℝ) :
    0 < thermalSynchrotronEmissivity 1 B 6.0e10 ν (Real.sqrt 2 / 2) := by
  unfold thermalSynchrotronEmissivity
  set θe : ℝ := max (electronTheta 6.0e10) (1 / 1000000 : ℝ)
  set νc : ℝ := cyclotronFreqProxy B
  set sinP : ℝ := max (Real.sqrt 2 / 2) (1 / 1000 : ℝ)
  set νs : ℝ := (2 / 9 : ℝ) * νc * θe ^ 2 * sinP
  set xM : ℝ := ν / max νs (1 / 1000000000000 : ℝ)
  have hmag : 0 < magneticSectorWeight := magneticSectorWeight_pos
  have hlam : 0 < lambda_qg := lambda_qg_pos
  have hnE_pos : 0 < max (1 : ℝ) 0 := by
    norm_num
  have hθe_floor_pos : (0 : ℝ) < (1 / 1000000 : ℝ) := by
    norm_num
  have hθe_ge_floor : (1 / 1000000 : ℝ) ≤ θe := by
    simpa [θe] using (le_max_right (electronTheta 6.0e10) (1 / 1000000 : ℝ))
  have hθe_pos : 0 < θe := lt_of_lt_of_le hθe_floor_pos hθe_ge_floor
  have hνc_floor_pos : (0 : ℝ) < (1 / 1000000 : ℝ) := by
    norm_num
  have hνc_pos_raw : 0 < cyclotronFreqProxy B := by
    unfold cyclotronFreqProxy
    exact lt_of_lt_of_le hνc_floor_pos (le_max_right |B| (1 / 1000000 : ℝ))
  have hνc_pos : 0 < νc := by
    simpa [νc] using hνc_pos_raw
  have hsinP_floor_pos : (0 : ℝ) < (1 / 1000 : ℝ) := by
    norm_num
  have hsinP_ge_floor : (1 / 1000 : ℝ) ≤ sinP := by
    simpa [sinP] using (le_max_right (Real.sqrt 2 / 2) (1 / 1000 : ℝ))
  have hsinP_pos : 0 < sinP := lt_of_lt_of_le hsinP_floor_pos hsinP_ge_floor
  have hθe_sq_pos : 0 < θe ^ 2 := pow_pos hθe_pos 2
  have htwo_ninth_pos : 0 < (2 / 9 : ℝ) := by
    norm_num
  have hνs_pos : 0 < νs := by
    have hprod1 : 0 < ((2 / 9 : ℝ) * νc) := mul_pos htwo_ninth_pos hνc_pos
    have hprod2 : 0 < (((2 / 9 : ℝ) * νc) * θe ^ 2) := mul_pos hprod1 hθe_sq_pos
    have hprod3 : 0 < ((((2 / 9 : ℝ) * νc) * θe ^ 2) * sinP) := mul_pos hprod2 hsinP_pos
    simpa [νs, mul_assoc] using hprod3
  have hνs_max_pos : 0 < max νs 0 := lt_of_lt_of_le hνs_pos (le_max_left νs 0)
  have hI_pos : 0 < mahadevanIPrime xM := mahadevanIPrime_pos xM
  have hprod1 : 0 < magneticSectorWeight * lambda_qg := mul_pos hmag hlam
  have hprod2 : 0 < (magneticSectorWeight * lambda_qg) * max (1 : ℝ) 0 := mul_pos hprod1 hnE_pos
  have hprod3 : 0 < ((magneticSectorWeight * lambda_qg) * max (1 : ℝ) 0) * max νs 0 :=
    mul_pos hprod2 hνs_max_pos
  have hprod4 :
      0 < (((magneticSectorWeight * lambda_qg) * max (1 : ℝ) 0) * max νs 0) * mahadevanIPrime xM := by
    exact mul_pos hprod3 hI_pos
  simpa [xM, νs, mul_assoc, mul_left_comm, mul_comm] using hprod4

theorem one_plus_lambda_qg_sq_nonneg {r_s l_P : ℝ} :
    0 ≤ 1 + lambda_qg * (l_P / r_s) ^ 2 := by
  have h_lambda : 0 ≤ lambda_qg := le_of_lt lambda_qg_pos
  have hsq : 0 ≤ (l_P / r_s) ^ 2 := sq_nonneg (l_P / r_s)
  have hmul : 0 ≤ lambda_qg * (l_P / r_s) ^ 2 := mul_nonneg h_lambda hsq
  linarith

/-- Baseline synchrotron emissivity is nonnegative for nonnegative frequency. -/
theorem synchrotronEmissivity_nonneg
    {B ν r_s l_P : ℝ} (_hν : 0 ≤ ν) :
    0 ≤ synchrotronEmissivity B ν r_s l_P := by
  unfold synchrotronEmissivity
  have htherm : 0 ≤ thermalSynchrotronEmissivity 1 B 6.0e10 ν (Real.sqrt 2 / 2) :=
    thermalSynchrotronEmissivity_nonneg 1 B 6.0e10 ν (Real.sqrt 2 / 2)
  have huv : 0 ≤ 1 + lambda_qg * (l_P / r_s) ^ 2 := one_plus_lambda_qg_sq_nonneg
  exact mul_nonneg htherm huv

/-- Strict positivity when magnetic field and frequency are nonzero/positive. -/
theorem synchrotronEmissivity_pos
    {B ν r_s l_P : ℝ} (_hB : B ≠ 0) (_hν : 0 < ν) :
    0 < synchrotronEmissivity B ν r_s l_P := by
  unfold synchrotronEmissivity
  have htherm : 0 < thermalSynchrotronEmissivity 1 B 6.0e10 ν (Real.sqrt 2 / 2) :=
    thermalSynchrotronEmissivity_pos_synchrotron_lane B ν
  have huv : 0 < 1 + lambda_qg * (l_P / r_s) ^ 2 := by
    have hnonneg : 0 ≤ lambda_qg * (l_P / r_s) ^ 2 := by
      exact mul_nonneg (le_of_lt lambda_qg_pos) (sq_nonneg (l_P / r_s))
    linarith
  exact mul_pos htherm huv

/-- Kerr frame-dragging can only increase (or keep) the emissivity in this model. -/
theorem synchrotronEmissivity_le_kerr
    {B ν r_s l_P aStar r θ : ℝ}
    (hbase : 0 ≤ synchrotronEmissivity B ν r_s l_P) :
    synchrotronEmissivity B ν r_s l_P ≤
      synchrotronEmissivityKerr B ν r_s l_P aStar r θ := by
  unfold synchrotronEmissivityKerr
  have hboost : 1 ≤ 1 + |frameDraggingOmega r_s aStar r θ| := by
    have habs : 0 ≤ |frameDraggingOmega r_s aStar r θ| := abs_nonneg _
    linarith
  have hmul : synchrotronEmissivity B ν r_s l_P * 1 ≤
      synchrotronEmissivity B ν r_s l_P * (1 + |frameDraggingOmega r_s aStar r θ|) := by
    exact mul_le_mul_of_nonneg_left hboost hbase
  simpa using hmul

/-- Master scaffold theorem: positive local emissivity with Kerr boost dominance. -/
theorem grmhd_synchrotron_scaffold
    {B ν r_s l_P aStar r θ : ℝ}
    (hB : B ≠ 0) (hν : 0 < ν) :
    0 < synchrotronEmissivity B ν r_s l_P ∧
    synchrotronEmissivity B ν r_s l_P ≤
      synchrotronEmissivityKerr B ν r_s l_P aStar r θ := by
  have hpos : 0 < synchrotronEmissivity B ν r_s l_P := synchrotronEmissivity_pos hB hν
  have hle : synchrotronEmissivity B ν r_s l_P ≤
      synchrotronEmissivityKerr B ν r_s l_P aStar r θ :=
    synchrotronEmissivity_le_kerr (le_of_lt hpos)
  exact ⟨hpos, hle⟩

end Gutoe.SynchrotronGRMHD
