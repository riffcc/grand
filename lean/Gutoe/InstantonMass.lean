/-
 * GUTOE — Z₃ Instanton Action and Mass Ratio Threshold
 * Copyright (C) 2026  Riff Labs
 * AGPL-3.0-or-later
 *
 * Three qualitative theorems about the Z₃ instanton action S_inst(t).
 *
 * Physical setup:
 *   cycle_prob_rg(t) = cp · (1 − b · ln(t+1))   one-loop RG fugacity
 *   where b = (β₀/2π) · α_UV, Landau pole t_* = exp(1/b) − 1.
 *
 *   S_inst(t) = −ln(fugacity(t))
 *
 * The mass ratio threshold t* is where S_inst first crosses ln(m_p/m_e) ≈ 7.515.
 * Numerically (default LatticeConfig): t_Landau ≈ 148, t* = 141.
 *
 * Theorem A (s_inst_strictMono):
 *   S_inst is strictly increasing on [0, t_Landau b).
 *
 * Theorem B (s_inst_unbounded):
 *   For any M : ℝ, ∃ t ∈ [0, t_Landau b) with S_inst(t) ≥ M.
 *
 * Theorem C (s_inst_threshold):
 *   For any S_target ≥ S_inst(0) = −ln(cp), ∃ x ∈ [0, t_Landau b) with
 *   S_inst(x) = S_target.  (IVT: the action hits every value above its minimum.)
 *
 * Corollary (mass_ratio_threshold_exists):
 *   For cp ≥ 1/1836.15, ∃ x* with S_inst(x*) = ln(m_p/m_e).
 *
 * All theorems no sorry.
 -/

import Mathlib
import Gutoe.Z3Uniqueness

namespace Gutoe.InstantonMass

open Gutoe.Z3Uniqueness

open Real Set

-- ══════════════════════════════════════════════════════════════════════════════
-- Definitions
-- ══════════════════════════════════════════════════════════════════════════════

/-- The Z₃ instanton fugacity at real time t ∈ [0, t_Landau).
    One-loop RG: fugacity(t, b, cp) = cp · (1 − b · ln(t+1)). -/
noncomputable def fugacity (t b cp : ℝ) : ℝ :=
  cp * (1 - b * Real.log (t + 1))

/-- The Z₃ instanton action: S_inst(t) = −ln(fugacity(t)).
    Strictly increasing on [0, t_Landau); diverges to +∞ at the Landau pole. -/
noncomputable def s_inst (t b cp : ℝ) : ℝ :=
  -Real.log (fugacity t b cp)

/-- The Landau pole: t_Landau(b) = exp(1/b) − 1. -/
noncomputable def t_landau (b : ℝ) : ℝ := Real.exp (1 / b) - 1

-- ══════════════════════════════════════════════════════════════════════════════
-- Auxiliary lemmas
-- ══════════════════════════════════════════════════════════════════════════════

lemma t_plus_one_pos {t : ℝ} (ht : 0 ≤ t) : 0 < t + 1 := by linarith

/-- The Landau pole is strictly positive when b > 0. -/
lemma t_landau_pos {b : ℝ} (hb : 0 < b) : 0 < t_landau b := by
  unfold t_landau
  have h : (1 : ℝ) < Real.exp (1 / b) := by
    apply Real.one_lt_exp_iff.mpr
    positivity
  linarith

/-- At t = 0, fugacity(0, b, cp) = cp. -/
lemma fugacity_at_zero (b cp : ℝ) : fugacity 0 b cp = cp := by
  simp [fugacity, Real.log_one]

/-- S_inst(0) = −ln(cp). -/
lemma s_inst_at_zero (b cp : ℝ) : s_inst 0 b cp = -Real.log cp := by
  simp [s_inst, fugacity_at_zero]

/-- For t ∈ [0, t_Landau b), the fugacity is strictly positive. -/
lemma fugacity_pos {t b cp : ℝ} (ht : 0 ≤ t) (hb : 0 < b) (hcp : 0 < cp)
    (ht_lt : t < t_landau b) : 0 < fugacity t b cp := by
  unfold fugacity
  apply mul_pos hcp
  rw [sub_pos]
  have ht1 : 0 < t + 1 := t_plus_one_pos ht
  have hlt : t + 1 < Real.exp (1 / b) := by
    unfold t_landau at ht_lt; linarith
  have hlog : Real.log (t + 1) < 1 / b := by
    rwa [← Real.log_exp (1 / b), Real.log_lt_log_iff ht1 (Real.exp_pos _)]
  -- b * log(t+1) < b * (1/b) = 1
  have hmul : 0 < b * (1 / b - Real.log (t + 1)) := mul_pos hb (by linarith)
  have heq : b * (1 / b - Real.log (t + 1)) = b * (1 / b) - b * Real.log (t + 1) := by ring
  have hone : b * (1 / b) = 1 := mul_one_div_cancel (ne_of_gt hb)
  linarith

/-- fugacity is continuous on [0, ∞) since log(t+1) is smooth there. -/
lemma fugacity_continuousOn_Ici (b cp : ℝ) :
    ContinuousOn (fun t => fugacity t b cp) (Set.Ici 0) := by
  unfold fugacity
  apply ContinuousOn.mul continuousOn_const
  apply ContinuousOn.sub continuousOn_const
  apply ContinuousOn.mul continuousOn_const
  -- log(t+1) continuous on Ici 0 since t+1 ≠ 0 there
  apply ContinuousOn.comp Real.continuousOn_log
    ((continuous_id.add continuous_const).continuousOn)
  intro t ht
  simp only [Set.mem_compl_iff, Set.mem_singleton_iff]
  exact ne_of_gt (t_plus_one_pos (Set.mem_Ici.mp ht))

/-- s_inst is continuous on [0, t_Landau b) where fugacity > 0. -/
lemma s_inst_continuousOn (b cp : ℝ) (hb : 0 < b) (hcp : 0 < cp) :
    ContinuousOn (fun t => s_inst t b cp) (Set.Ico 0 (t_landau b)) := by
  unfold s_inst
  apply ContinuousOn.neg
  apply ContinuousOn.log
  · exact (fugacity_continuousOn_Ici b cp).mono Set.Ico_subset_Ici_self
  · intro t ht
    exact ne_of_gt (fugacity_pos (Set.mem_Ico.mp ht).1 hb hcp (Set.mem_Ico.mp ht).2)

-- ══════════════════════════════════════════════════════════════════════════════
-- Theorem A: S_inst is strictly monotone on [0, t_Landau b)
-- ══════════════════════════════════════════════════════════════════════════════

/-- **Theorem A**: The Z₃ instanton action is strictly increasing on [0, t_Landau b).

    Physical meaning: as the colour coupling α_s runs toward the Landau pole,
    the Z₃ tunnelling fugacity decreases, so the action increases. The lepton
    mass scale is set at the unique t* where the action first reaches ln(m_p/m_e). -/
theorem s_inst_strictMono (b cp : ℝ) (hb : 0 < b) (hcp : 0 < cp)
    {t₁ t₂ : ℝ} (ht₁ : 0 ≤ t₁) (h12 : t₁ < t₂) (ht₂_lt : t₂ < t_landau b) :
    s_inst t₁ b cp < s_inst t₂ b cp := by
  unfold s_inst
  apply neg_lt_neg
  apply Real.log_lt_log
  · exact fugacity_pos (le_of_lt (lt_of_le_of_lt ht₁ h12)) hb hcp ht₂_lt
  · -- fugacity t₂ < fugacity t₁ : cp*(1-b*log(t₂+1)) < cp*(1-b*log(t₁+1))
    unfold fugacity
    have hlog12 : Real.log (t₁ + 1) < Real.log (t₂ + 1) :=
      Real.log_lt_log (t_plus_one_pos ht₁) (by linarith)
    have hmul : 0 < b * (Real.log (t₂ + 1) - Real.log (t₁ + 1)) :=
      mul_pos hb (by linarith)
    have hcp_mul : 0 < cp * (b * (Real.log (t₂ + 1) - Real.log (t₁ + 1))) :=
      mul_pos hcp hmul
    nlinarith [hcp_mul]

-- ══════════════════════════════════════════════════════════════════════════════
-- Theorem B: S_inst is unbounded above on [0, t_Landau b)
-- ══════════════════════════════════════════════════════════════════════════════

/-- **Theorem B**: For any bound M, there exists t ∈ [0, t_Landau b) with S_inst(t) ≥ M.

    Constructive proof: given M, choose t₀ = exp((1 − exp(−M)/cp) / b) − 1.
    Then fugacity(t₀) = exp(−M), so S_inst(t₀) = −ln(exp(−M)) = M. -/
theorem s_inst_unbounded (b cp : ℝ) (hb : 0 < b) (hcp : 0 < cp) (M : ℝ) :
    ∃ t : ℝ, 0 ≤ t ∧ t < t_landau b ∧ M ≤ s_inst t b cp := by
  by_cases hM : M ≤ -Real.log cp
  · -- Case 1: s_inst(0) = −log(cp) ≥ M — use t = 0
    exact ⟨0, le_refl 0, t_landau_pos hb, by simp [s_inst_at_zero]; exact hM⟩
  · -- Case 2: M > −log(cp) — construct t₀ explicitly
    push_neg at hM
    -- Let e = exp(−M)/cp ∈ (0, 1) since M > −log(cp)
    have he_pos : 0 < Real.exp (-M) / cp := div_pos (Real.exp_pos _) hcp
    have he_lt1 : Real.exp (-M) / cp < 1 := by
      rw [div_lt_one hcp, ← Real.exp_log hcp]
      exact Real.exp_lt_exp.mpr (by linarith)
    -- t₀ = exp((1 − e)/b) − 1 ≥ 0  since (1−e) > 0 → exp > 1
    have ht₀_exp_pos : 0 < (1 - Real.exp (-M) / cp) / b :=
      div_pos (by linarith) hb
    -- Witness
    refine ⟨Real.exp ((1 - Real.exp (-M) / cp) / b) - 1, ?_, ?_, ?_⟩
    · -- t₀ ≥ 0
      have : (1 : ℝ) ≤ Real.exp ((1 - Real.exp (-M) / cp) / b) :=
        Real.one_le_exp_iff.mpr (le_of_lt ht₀_exp_pos)
      linarith
    · -- t₀ < t_landau b
      unfold t_landau
      apply sub_lt_sub_right
      rw [Real.exp_lt_exp]
      -- (1−e)/b < 1/b since e > 0 and b > 0
      have hdiff : 1 / b - (1 - Real.exp (-M) / cp) / b = Real.exp (-M) / cp / b := by
        ring
      linarith [div_pos he_pos hb]
    · -- M ≤ s_inst(t₀): equality holds, S_inst(t₀) = M
      apply le_of_eq
      unfold s_inst fugacity
      -- log(t₀ + 1) = (1 − e)/b  since t₀ + 1 = exp((1-e)/b)
      have hlog : Real.log (Real.exp ((1 - Real.exp (-M) / cp) / b) - 1 + 1) =
                  (1 - Real.exp (-M) / cp) / b := by
        have hrw : Real.exp ((1 - Real.exp (-M) / cp) / b) - 1 + 1 =
                   Real.exp ((1 - Real.exp (-M) / cp) / b) := by ring
        rw [hrw, Real.log_exp]
      rw [hlog]
      -- cp * (1 − b · (1−e)/b) = cp * e = exp(−M)
      have hfug : cp * (1 - b * ((1 - Real.exp (-M) / cp) / b)) = Real.exp (-M) := by
        have hcancel : b * ((1 - Real.exp (-M) / cp) / b) = 1 - Real.exp (-M) / cp := by
          field_simp
        rw [hcancel]
        field_simp
        ring
      rw [hfug, Real.log_exp]
      ring

-- ══════════════════════════════════════════════════════════════════════════════
-- Theorem C: The mass ratio threshold exists (Intermediate Value Theorem)
-- ══════════════════════════════════════════════════════════════════════════════

/-- **Theorem C**: For any S_target ≥ S_inst(0) = −ln(cp), there exists
    x ∈ [0, t_Landau b) with S_inst(x) = S_target.

    Proof: S_inst is continuous on [0, t_Landau b), equals −ln(cp) at t=0,
    and is unbounded above (Theorem B). IVT gives the existence of x. -/
theorem s_inst_threshold (b cp S_target : ℝ) (hb : 0 < b) (hcp : 0 < cp)
    (hS : s_inst 0 b cp ≤ S_target) :
    ∃ x : ℝ, 0 ≤ x ∧ x < t_landau b ∧ s_inst x b cp = S_target := by
  -- Step 1: Get t₀ ∈ [0, t_Landau) with s_inst(t₀) ≥ S_target (Theorem B)
  obtain ⟨t₀, ht₀_nn, ht₀_lt, ht₀_ge⟩ := s_inst_unbounded b cp hb hcp S_target
  -- Step 2: s_inst is continuous on [0, t₀] ⊆ [0, t_Landau)
  have hcont : ContinuousOn (fun t => s_inst t b cp) (Set.Icc 0 t₀) :=
    (s_inst_continuousOn b cp hb hcp).mono
      (fun x hx => Set.mem_Ico.mpr ⟨(Set.mem_Icc.mp hx).1,
                                     lt_of_le_of_lt (Set.mem_Icc.mp hx).2 ht₀_lt⟩)
  -- Step 3: Apply IVT on the preconnected set [0, t₀]
  -- IsPreconnected.intermediate_value₂ : f(a) ≤ g(a) ∧ g(b) ≤ f(b) → ∃ c, f(c) = g(c)
  obtain ⟨c, hc_mem, hc_eq⟩ := isPreconnected_Icc.intermediate_value₂
    (Set.left_mem_Icc.mpr ht₀_nn)
    (Set.right_mem_Icc.mpr ht₀_nn)
    hcont
    continuousOn_const
    hS
    ht₀_ge
  -- Step 4: c ∈ [0, t₀] ⊆ [0, t_Landau)
  exact ⟨c, (Set.mem_Icc.mp hc_mem).1,
         lt_of_le_of_lt (Set.mem_Icc.mp hc_mem).2 ht₀_lt,
         hc_eq⟩

-- ══════════════════════════════════════════════════════════════════════════════
-- Corollary: The proton-to-electron mass ratio threshold exists
-- ══════════════════════════════════════════════════════════════════════════════

/-- **Corollary**: For parameters b > 0 and cp ≥ 1/1836.15 > 0, there exists
    x* ∈ [0, t_Landau b) with S_inst(x*) = ln(m_p/m_e) ≈ 7.515.

    The hypothesis cp ≥ 1/1836.15 ensures S_inst(0) = −ln(cp) ≤ ln(1836.15),
    i.e., the initial action is already below the mass-ratio target.

    For the default LatticeConfig (cp = 0.05 ≈ 91.8 × 1/1836.15), x* ≈ 141
    and t_Landau ≈ 148 (7 steps before confinement). -/
theorem mass_ratio_threshold_exists (b cp : ℝ) (hb : 0 < b) (hcp : 0 < cp)
    (hcp_lb : 1 / 1836.15 ≤ cp) :
    ∃ x : ℝ, 0 ≤ x ∧ x < t_landau b ∧ s_inst x b cp = Real.log 1836.15 := by
  apply s_inst_threshold b cp (Real.log 1836.15) hb hcp
  -- Need: s_inst(0) = −log(cp) ≤ log(1836.15)
  rw [s_inst_at_zero]
  -- −log(cp) ≤ log(1836.15)
  -- ↔ 0 ≤ log(1836.15) + log(cp)
  -- ↔ 0 ≤ log(1836.15 * cp)   since both > 0
  -- ↔ 1 ≤ 1836.15 * cp        since log x ≥ 0 ↔ x ≥ 1
  -- ↔ cp ≥ 1/1836.15 ✓
  have h1836 : (0 : ℝ) < 1836.15 := by norm_num
  have hprod : (1 : ℝ) ≤ 1836.15 * cp := by
    rw [show (1 : ℝ) = 1836.15 * (1 / 1836.15) from by norm_num]
    exact mul_le_mul_of_nonneg_left hcp_lb (le_of_lt h1836)
  have hlog_nonneg : (0 : ℝ) ≤ Real.log (1836.15 * cp) :=
    Real.log_nonneg (by linarith)
  rw [Real.log_mul (ne_of_gt h1836) (ne_of_gt hcp)] at hlog_nonneg
  linarith

-- ══════════════════════════════════════════════════════════════════════════════
-- Z₃ UV instanton action = |magneticTriplet|
-- ══════════════════════════════════════════════════════════════════════════════

/-- **Z₃ Instanton UV Action**: When cp = exp(−|magneticTriplet|), the initial
    instanton action S_inst(0) equals exactly |magneticTriplet| = 3.

    Physical meaning: the UV fugacity cp = exp(−3) encodes "one unit of action
    per Z₃ quark-colour corner" — the instanton traverses all three vacua
    {γ¹², γ¹³, γ²³} and pays one unit of action at each.

    This links the algebraic input cp to the Z₃ structure of Cl(1,3):
      cp = exp(−|{γ¹², γ¹³, γ²³}|) = exp(−3)  →  S_inst(0) = 3. -/
theorem z3_instanton_initial_action (b : ℝ) :
    s_inst 0 b (Real.exp (-(magneticTriplet.card : ℝ))) = (magneticTriplet.card : ℝ) := by
  rw [s_inst_at_zero, Real.log_exp]
  ring

end Gutoe.InstantonMass
