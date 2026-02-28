/-
 * GUTOE — Singularity Resolution (GRAND-128)
 *
 * Bridge lane:
 *   - Black-hole core regularization from `Gutoe.GravityMetric`
 *   - Big-Bang replacement via finite-density bounce kernel
 *
 * No `sorry`.
 -/

import Mathlib
import Gutoe.GravityMetric

namespace Gutoe.SingularityResolution

open Gutoe.GravityMetric

/-- Structural critical density from shared `C_inf` and an input Planck-density scale. -/
noncomputable def rhoCrit (rhoPlanck : ℝ) : ℝ := C_inf * rhoPlanck

/-- Bounce kernel used in the regularized Friedmann lane:
    `κ(ρ) = ρ * (1 - ρ/ρ_crit)`. -/
noncomputable def bounceKernel (rho rhoCrit : ℝ) : ℝ :=
  rho * (1 - rho / rhoCrit)

/-- Lattice-regularized Friedmann RHS proxy:
    `H² = max(0, (8π/3) * κ(ρ))`. -/
noncomputable def hubbleSq (rho rhoCrit : ℝ) : ℝ :=
  max 0 (((8 : ℝ) * Real.pi / 3) * bounceKernel rho rhoCrit)

theorem rho_crit_pos {rhoPlanck : ℝ} (h : 0 < rhoPlanck) :
    0 < rhoCrit rhoPlanck := by
  unfold rhoCrit
  exact mul_pos C_inf_pos h

theorem bounce_kernel_zero_at_zero {rhoCrit : ℝ} :
    bounceKernel 0 rhoCrit = 0 := by
  simp [bounceKernel]

theorem bounce_kernel_zero_at_critical {rhoCrit : ℝ} (h : rhoCrit ≠ 0) :
    bounceKernel rhoCrit rhoCrit = 0 := by
  unfold bounceKernel
  have hdiv : rhoCrit / rhoCrit = (1 : ℝ) := by
    field_simp [h]
  rw [hdiv]
  ring

theorem bounce_kernel_midpoint_eq {rhoCrit : ℝ} (h : rhoCrit ≠ 0) :
    bounceKernel (rhoCrit / 2) rhoCrit = rhoCrit / 4 := by
  unfold bounceKernel
  field_simp [h]
  ring

theorem bounce_kernel_midpoint_pos {rhoCrit : ℝ} (h : 0 < rhoCrit) :
    0 < bounceKernel (rhoCrit / 2) rhoCrit := by
  rw [bounce_kernel_midpoint_eq (ne_of_gt h)]
  nlinarith

theorem bounce_kernel_nonpos_of_ge {rho rhoCrit : ℝ}
    (hCrit : 0 < rhoCrit) (hge : rhoCrit ≤ rho) :
    bounceKernel rho rhoCrit ≤ 0 := by
  unfold bounceKernel
  have hrho_nonneg : 0 ≤ rho := le_trans (le_of_lt hCrit) hge
  have hsub_nonpos : rhoCrit - rho ≤ 0 := sub_nonpos.mpr hge
  have hdiv_nonpos : (rhoCrit - rho) / rhoCrit ≤ 0 := by
    exact div_nonpos_of_nonpos_of_nonneg hsub_nonpos (le_of_lt hCrit)
  have hrewrite : 1 - rho / rhoCrit = (rhoCrit - rho) / rhoCrit := by
    field_simp [ne_of_gt hCrit]
  rw [hrewrite]
  exact mul_nonpos_of_nonneg_of_nonpos hrho_nonneg hdiv_nonpos

theorem hubble_sq_nonneg {rho rhoCrit : ℝ} :
    0 ≤ hubbleSq rho rhoCrit := by
  unfold hubbleSq
  exact le_max_left 0 (((8 : ℝ) * Real.pi / 3) * bounceKernel rho rhoCrit)

theorem hubble_sq_at_critical {rhoCrit : ℝ} (h : rhoCrit ≠ 0) :
    hubbleSq rhoCrit rhoCrit = 0 := by
  unfold hubbleSq
  rw [bounce_kernel_zero_at_critical h]
  norm_num

theorem hubble_sq_midpoint_pos {rhoCrit : ℝ} (hCrit : 0 < rhoCrit) :
    0 < hubbleSq (rhoCrit / 2) rhoCrit := by
  unfold hubbleSq
  have hconst_pos : 0 < ((8 : ℝ) * Real.pi / 3) := by positivity
  have hkernel_pos : 0 < bounceKernel (rhoCrit / 2) rhoCrit :=
    bounce_kernel_midpoint_pos hCrit
  have hterm_pos : 0 < ((8 : ℝ) * Real.pi / 3) * bounceKernel (rhoCrit / 2) rhoCrit :=
    mul_pos hconst_pos hkernel_pos
  rw [max_eq_right (le_of_lt hterm_pos)]
  exact hterm_pos

theorem hubble_sq_above_critical_eq_zero {rho rhoCrit : ℝ}
    (hCrit : 0 < rhoCrit) (hge : rhoCrit ≤ rho) :
    hubbleSq rho rhoCrit = 0 := by
  unfold hubbleSq
  have hconst_nonneg : 0 ≤ ((8 : ℝ) * Real.pi / 3) := by positivity
  have hkernel_nonpos : bounceKernel rho rhoCrit ≤ 0 :=
    bounce_kernel_nonpos_of_ge hCrit hge
  have hterm_nonpos : ((8 : ℝ) * Real.pi / 3) * bounceKernel rho rhoCrit ≤ 0 :=
    mul_nonpos_of_nonneg_of_nonpos hconst_nonneg hkernel_nonpos
  exact max_eq_left hterm_nonpos

/-- GRAND-128 closure gate:
    - BH origin remains finite via `r_eff` core regularization.
    - Big-Bang singularity is replaced by a finite-density bounce with
      `H²(ρ_crit)=0`, positive pre-bounce branch, and no continuation above `ρ_crit`. -/
theorem singularity_resolution_gate
    {l_P r_s rhoPlanck : ℝ}
    (hlP : 0 < l_P)
    (hrhoPlanck : 0 < rhoPlanck) :
    0 < r_eff 0 l_P ∧
    g_tt 0 r_s l_P = -(1 - r_s / r_core l_P) ∧
    0 < rhoCrit rhoPlanck ∧
    hubbleSq (rhoCrit rhoPlanck) (rhoCrit rhoPlanck) = 0 ∧
    0 < hubbleSq (rhoCrit rhoPlanck / 2) (rhoCrit rhoPlanck) ∧
    hubbleSq (((11 : ℝ) / 10) * rhoCrit rhoPlanck) (rhoCrit rhoPlanck) = 0 := by
  have hCritPos : 0 < rhoCrit rhoPlanck := rho_crit_pos hrhoPlanck
  have hCritNe : rhoCrit rhoPlanck ≠ 0 := ne_of_gt hCritPos
  refine ⟨r_eff_pos hlP, g_tt_at_origin hlP, hCritPos, ?_, ?_, ?_⟩
  · exact hubble_sq_at_critical hCritNe
  · exact hubble_sq_midpoint_pos hCritPos
  · have hge : rhoCrit rhoPlanck ≤ ((11 : ℝ) / 10) * rhoCrit rhoPlanck := by
      nlinarith [hCritPos]
    exact hubble_sq_above_critical_eq_zero hCritPos hge

end Gutoe.SingularityResolution
