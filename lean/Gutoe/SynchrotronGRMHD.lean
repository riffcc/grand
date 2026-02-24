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
  magneticSectorWeight
    * lambda_qg
    * B ^ 2
    * ν
    * Real.exp (-(ν * r_s))
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

theorem one_plus_lambda_qg_sq_nonneg {r_s l_P : ℝ} :
    0 ≤ 1 + lambda_qg * (l_P / r_s) ^ 2 := by
  have h_lambda : 0 ≤ lambda_qg := le_of_lt lambda_qg_pos
  have hsq : 0 ≤ (l_P / r_s) ^ 2 := sq_nonneg (l_P / r_s)
  have hmul : 0 ≤ lambda_qg * (l_P / r_s) ^ 2 := mul_nonneg h_lambda hsq
  linarith

/-- Baseline synchrotron emissivity is nonnegative for nonnegative frequency. -/
theorem synchrotronEmissivity_nonneg
    {B ν r_s l_P : ℝ} (hν : 0 ≤ ν) :
    0 ≤ synchrotronEmissivity B ν r_s l_P := by
  unfold synchrotronEmissivity
  have hmag : 0 ≤ magneticSectorWeight := by
    rw [magneticSectorWeight_eq_three]
    norm_num
  have h_lambda : 0 ≤ lambda_qg := le_of_lt lambda_qg_pos
  have hB : 0 ≤ B ^ 2 := sq_nonneg B
  have hexp : 0 ≤ Real.exp (-(ν * r_s)) := le_of_lt (Real.exp_pos _)
  have huv : 0 ≤ 1 + lambda_qg * (l_P / r_s) ^ 2 := one_plus_lambda_qg_sq_nonneg
  exact mul_nonneg
    (mul_nonneg
      (mul_nonneg
        (mul_nonneg
          (mul_nonneg hmag h_lambda)
          hB)
        hν)
      hexp)
    huv

/-- Strict positivity when magnetic field and frequency are nonzero/positive. -/
theorem synchrotronEmissivity_pos
    {B ν r_s l_P : ℝ} (hB : B ≠ 0) (hν : 0 < ν) :
    0 < synchrotronEmissivity B ν r_s l_P := by
  unfold synchrotronEmissivity
  have hmag : 0 < magneticSectorWeight := by
    rw [magneticSectorWeight_eq_three]
    norm_num
  have h_lambda : 0 < lambda_qg := lambda_qg_pos
  have hB2 : 0 < B ^ 2 := by
    have habs : 0 < |B| := abs_pos.mpr hB
    have hsq : 0 < |B| ^ 2 := sq_pos_of_pos habs
    simpa [sq_abs] using hsq
  have hexp : 0 < Real.exp (-(ν * r_s)) := Real.exp_pos _
  have huv : 0 < 1 + lambda_qg * (l_P / r_s) ^ 2 := by
    have hnonneg : 0 ≤ lambda_qg * (l_P / r_s) ^ 2 := by
      exact mul_nonneg (le_of_lt lambda_qg_pos) (sq_nonneg (l_P / r_s))
    linarith
  exact mul_pos
    (mul_pos
      (mul_pos
        (mul_pos
          (mul_pos hmag h_lambda)
          hB2)
        hν)
      hexp)
    huv

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
