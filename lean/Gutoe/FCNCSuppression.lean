/-
 * GUTOE — FCNC Suppression / GIM Mechanism (GRAND-132)
 *
 * Structural lane:
 *   - CKM off-diagonal suppressions from Cl(1,3) counts (`1/24`, `1/272`)
 *   - GIM cancellation at constant kernel from unitarity sum-rule algebra
 *   - Mass-difference form showing only flavor splittings survive
 *
 * No `sorry`.
 -/

import Mathlib
import Gutoe.FlavorMixing

namespace Gutoe.FCNCSuppression

open Gutoe.FlavorMixing

/-- Structural FCNC loop suppressor from CKM off-diagonal mixings. -/
noncomputable def ckmLoopSuppressor : ℝ :=
  ckmSin23 * ckmSin13

/-- Squared loop suppressor proxy. -/
noncomputable def ckmLoopSuppressorSq : ℝ :=
  ckmLoopSuppressor ^ 2

theorem ckm_loop_suppressor_eq :
    ckmLoopSuppressor = (1 : ℝ) / 6528 := by
  unfold ckmLoopSuppressor
  rw [ckm_structural_values.1, ckm_structural_values.2]
  norm_num

theorem ckm_loop_suppressor_sq_eq :
    ckmLoopSuppressorSq = (1 : ℝ) / 42614784 := by
  unfold ckmLoopSuppressorSq
  rw [ckm_loop_suppressor_eq]
  norm_num

theorem ckm_loop_suppressor_sq_lt_micro :
    ckmLoopSuppressorSq < (1 : ℝ) / 10000000 := by
  rw [ckm_loop_suppressor_sq_eq]
  norm_num

/-- GIM cancellation at tree/degenerate-kernel level:
    if `λ_u + λ_c + λ_t = 0`, then `Σ λ_i κ = 0`. -/
theorem gim_constant_kernel_cancel
    {lu lc lt κ : ℂ}
    (hsum : lu + lc + lt = 0) :
    lu * κ + lc * κ + lt * κ = 0 := by
  calc
    lu * κ + lc * κ + lt * κ = (lu + lc + lt) * κ := by ring
    _ = 0 := by simp [hsum]

/-- GIM mass-difference form:
    under `λ_u + λ_c + λ_t = 0`, constant pieces cancel and
    the loop amplitude keeps only flavor differences. -/
theorem gim_mass_difference_form
    {lu lc lt fu fc ft : ℂ}
    (hsum : lu + lc + lt = 0) :
    lu * fu + lc * fc + lt * ft
      = lc * (fc - fu) + lt * (ft - fu) := by
  calc
    lu * fu + lc * fc + lt * ft
        = (lu + lc + lt) * fu + lc * (fc - fu) + lt * (ft - fu) := by ring
    _ = lc * (fc - fu) + lt * (ft - fu) := by simp [hsum]

/-- GRAND-132 closure gate:
    structural CKM suppressor is tiny, and the GIM cancellation identities hold. -/
theorem fcnc_suppression_gate :
    ckmLoopSuppressor = (1 : ℝ) / 6528 ∧
    ckmLoopSuppressorSq < (1 : ℝ) / 10000000 ∧
    (∀ {lu lc lt κ : ℂ}, lu + lc + lt = 0 →
      lu * κ + lc * κ + lt * κ = 0) ∧
    (∀ {lu lc lt fu fc ft : ℂ}, lu + lc + lt = 0 →
      lu * fu + lc * fc + lt * ft = lc * (fc - fu) + lt * (ft - fu)) := by
  refine ⟨ckm_loop_suppressor_eq, ckm_loop_suppressor_sq_lt_micro, ?_, ?_⟩
  · intro lu lc lt κ hsum
    exact gim_constant_kernel_cancel hsum
  · intro lu lc lt fu fc ft hsum
    exact gim_mass_difference_form hsum

end Gutoe.FCNCSuppression
