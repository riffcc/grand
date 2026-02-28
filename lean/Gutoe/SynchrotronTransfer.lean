/-
 * GUTOE — Covariant Synchrotron Transfer Scaffold
 * Copyright (C) 2026 Riff Labs
 * AGPL-3.0-or-later
 *
 * First formal transfer step beyond emissivity-only rendering:
 *   - covariant redshift scaling for emissivity/absorption coefficients
 *   - single-step radiative transfer map
 *   - positivity and interval bounds (no vacuous theorems)
 *
 * No `sorry`.
-/

import Mathlib
import Gutoe.SynchrotronGRMHD

namespace Gutoe.SynchrotronTransfer

open Real

/-- Covariant emissivity scaling from the GR invariant `j_ν / ν^2`:
`j_obs = j_em * g^2`, with `g = ν_obs/ν_em`. -/
noncomputable def covariantEmissivity (jLocal g : ℝ) : ℝ :=
  jLocal * g ^ 2

/-- Covariant absorption scaling from the GR invariant `α_ν * ν`:
`α_obs = α_em / g`. -/
noncomputable def covariantAbsorption (alphaLocal g : ℝ) : ℝ :=
  alphaLocal / g

/-- One-zone transfer update:
`I_out = I_in * exp(-τ) + S * (1 - exp(-τ))` where `S = j/α` (source function). -/
noncomputable def transferStep (iIn source tau : ℝ) : ℝ :=
  iIn * Real.exp (-tau) + source * (1 - Real.exp (-tau))

/-- Equivalent affine form:
`I_out = I_in + (S - I_in) * (1 - exp(-τ))`. -/
theorem transferStep_affine (iIn source tau : ℝ) :
    transferStep iIn source tau =
      iIn + (source - iIn) * (1 - Real.exp (-tau)) := by
  unfold transferStep
  ring

lemma exp_neg_tau_le_one {tau : ℝ} (htau : 0 ≤ tau) : Real.exp (-tau) ≤ 1 := by
  have hneg : -tau ≤ 0 := by linarith
  simpa using (Real.exp_le_one_iff.mpr hneg)

lemma one_minus_exp_neg_tau_nonneg {tau : ℝ} (htau : 0 ≤ tau) :
    0 ≤ 1 - Real.exp (-tau) := by
  linarith [exp_neg_tau_le_one htau]

/-- Covariant emissivity stays nonnegative when local emissivity and redshift are nonnegative. -/
theorem covariantEmissivity_nonneg {jLocal g : ℝ}
    (hj : 0 ≤ jLocal) (hg : 0 ≤ g) :
    0 ≤ covariantEmissivity jLocal g := by
  unfold covariantEmissivity
  exact mul_nonneg hj (pow_nonneg hg 2)

/-- Covariant absorption stays nonnegative when local absorption and redshift are nonnegative. -/
theorem covariantAbsorption_nonneg {alphaLocal g : ℝ}
    (ha : 0 ≤ alphaLocal) (hg : 0 < g) :
    0 ≤ covariantAbsorption alphaLocal g := by
  unfold covariantAbsorption
  exact div_nonneg ha (le_of_lt hg)

/-- Transfer step preserves nonnegativity for nonnegative inputs and optical depth. -/
theorem transferStep_nonneg {iIn source tau : ℝ}
    (hi : 0 ≤ iIn) (hs : 0 ≤ source) (htau : 0 ≤ tau) :
    0 ≤ transferStep iIn source tau := by
  unfold transferStep
  have h_exp_nonneg : 0 ≤ Real.exp (-tau) := le_of_lt (Real.exp_pos _)
  have h_one_minus_nonneg : 0 ≤ 1 - Real.exp (-tau) := one_minus_exp_neg_tau_nonneg htau
  exact add_nonneg
    (mul_nonneg hi h_exp_nonneg)
    (mul_nonneg hs h_one_minus_nonneg)

/-- If `I_in ≤ S`, transfer moves intensity upward (or equal). -/
theorem transferStep_ge_input {iIn source tau : ℝ}
    (hIS : iIn ≤ source) (htau : 0 ≤ tau) :
    iIn ≤ transferStep iIn source tau := by
  rw [transferStep_affine]
  have hfactor : 0 ≤ 1 - Real.exp (-tau) := one_minus_exp_neg_tau_nonneg htau
  have hdiff : 0 ≤ source - iIn := sub_nonneg.mpr hIS
  have hterm : 0 ≤ (source - iIn) * (1 - Real.exp (-tau)) := mul_nonneg hdiff hfactor
  linarith

/-- If `I_in ≤ S`, transfer never overshoots the source function. -/
theorem transferStep_le_source {iIn source tau : ℝ}
    (hIS : iIn ≤ source) (htau : 0 ≤ tau) :
    transferStep iIn source tau ≤ source := by
  rw [transferStep_affine]
  have hfactor : 0 ≤ 1 - Real.exp (-tau) := one_minus_exp_neg_tau_nonneg htau
  have hfactor_le : 1 - Real.exp (-tau) ≤ 1 := by
    have h_exp_nonneg : 0 ≤ Real.exp (-tau) := le_of_lt (Real.exp_pos _)
    linarith
  have hdiff : 0 ≤ source - iIn := sub_nonneg.mpr hIS
  have hscaled : (source - iIn) * (1 - Real.exp (-tau)) ≤ source - iIn := by
    calc
      (source - iIn) * (1 - Real.exp (-tau))
          ≤ (source - iIn) * 1 := by
            exact mul_le_mul_of_nonneg_left hfactor_le hdiff
      _ = source - iIn := by ring
  linarith

/-- Master transfer interval theorem:
for `I_in ≤ S` and `τ ≥ 0`, `I_out` stays in `[I_in, S]`. -/
theorem transferStep_in_interval {iIn source tau : ℝ}
    (hIS : iIn ≤ source) (htau : 0 ≤ tau) :
    iIn ≤ transferStep iIn source tau ∧ transferStep iIn source tau ≤ source := by
  exact ⟨transferStep_ge_input hIS htau, transferStep_le_source hIS htau⟩

end Gutoe.SynchrotronTransfer
