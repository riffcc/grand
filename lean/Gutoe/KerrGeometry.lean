/- 
 * GUTOE — Kerr Geometry Scaffold (GR baseline for GRAND-159)
 * Copyright (C) 2026 Riff Labs
 * AGPL-3.0-or-later
 *
 * This module is intentionally a lightweight, fully-proven scaffold:
 * it formalizes the rotating (Kerr) horizon/ergosphere algebra used by
 * the renderer roadmap, without claiming a full GUTOE-Kerr derivation yet.
 *
 * No `sorry`.
-/

import Mathlib
import Gutoe.GravityMetric

namespace Gutoe.KerrGeometry

open Real

/-- Kerr mass parameter in geometric units, using repo convention `r_s = 2M`. -/
noncomputable def mass (r_s : ℝ) : ℝ := r_s / 2

/-- Kerr spin length parameter `a = a_* M`. -/
noncomputable def spinLength (r_s aStar : ℝ) : ℝ := aStar * mass r_s

/-- Kerr horizon discriminant `M² - a² = M² (1 - a_*²)`. -/
noncomputable def horizonDiscriminant (r_s aStar : ℝ) : ℝ :=
  (mass r_s) ^ 2 - (spinLength r_s aStar) ^ 2

/-- Outer horizon radius `r_+ = M + √(M² - a²)`. -/
noncomputable def rPlus (r_s aStar : ℝ) : ℝ :=
  mass r_s + Real.sqrt (horizonDiscriminant r_s aStar)

/-- Inner horizon radius `r_- = M - √(M² - a²)`. -/
noncomputable def rMinus (r_s aStar : ℝ) : ℝ :=
  mass r_s - Real.sqrt (horizonDiscriminant r_s aStar)

/-- Static-limit (ergosphere) radius `r_erg(θ) = M + √(M² - a² cos²θ)`. -/
noncomputable def ergosphereRadius (r_s aStar θ : ℝ) : ℝ :=
  mass r_s + Real.sqrt ((mass r_s) ^ 2 - (spinLength r_s aStar) ^ 2 * (Real.cos θ) ^ 2)

/-- Kerr horizon angular velocity `Ω_H = a / (r_+² + a²)`. -/
noncomputable def horizonAngularVelocity (r_s aStar : ℝ) : ℝ :=
  spinLength r_s aStar / ((rPlus r_s aStar) ^ 2 + (spinLength r_s aStar) ^ 2)

/-- Kerr radial polynomial `Δ(r) = r² - r_s r + a²` (with `r_s = 2M`). -/
noncomputable def kerrDelta (r_s aStar r : ℝ) : ℝ :=
  r ^ 2 - (2 * mass r_s) * r + (spinLength r_s aStar) ^ 2

/-- Kerr frame-dragging denominator
    `A = (r²+a²)² - a² Δ sin²θ`. -/
noncomputable def frameDraggingDenom (r_s aStar r θ : ℝ) : ℝ :=
  (r ^ 2 + (spinLength r_s aStar) ^ 2) ^ 2
    - (spinLength r_s aStar) ^ 2 * (kerrDelta r_s aStar r) * (Real.sin θ) ^ 2

/-- ZAMO frame-dragging angular velocity
    `ω = (2 M a r) / A` (with an exact zero-guard at `A = 0`). -/
noncomputable def frameDraggingOmega (r_s aStar r θ : ℝ) : ℝ :=
  if frameDraggingDenom r_s aStar r θ = 0 then 0
  else (2 * mass r_s * spinLength r_s aStar * r) / frameDraggingDenom r_s aStar r θ

lemma mass_nonneg {r_s : ℝ} (hrs : 0 ≤ r_s) : 0 ≤ mass r_s := by
  unfold mass
  nlinarith

lemma abs_aStar_sq_le_one {aStar : ℝ} (ha : |aStar| ≤ 1) : aStar ^ 2 ≤ 1 := by
  have hsq : aStar ^ 2 ≤ (1 : ℝ) ^ 2 := by
    exact sq_le_sq.mpr (by simpa using ha)
  simpa using hsq

/-- For physical spin `|a_*| ≤ 1`, the Kerr horizon discriminant is nonnegative. -/
theorem horizonDiscriminant_nonneg {r_s aStar : ℝ}
    (ha : |aStar| ≤ 1) :
    0 ≤ horizonDiscriminant r_s aStar := by
  have hm2_nonneg : 0 ≤ (mass r_s) ^ 2 := sq_nonneg (mass r_s)
  have hfac_nonneg : 0 ≤ 1 - aStar ^ 2 := by
    have hs : aStar ^ 2 ≤ 1 := abs_aStar_sq_le_one ha
    linarith
  have hmul : 0 ≤ (mass r_s) ^ 2 * (1 - aStar ^ 2) := mul_nonneg hm2_nonneg hfac_nonneg
  have hrewrite : (mass r_s) ^ 2 * (1 - aStar ^ 2) = horizonDiscriminant r_s aStar := by
    unfold horizonDiscriminant spinLength
    ring
  simpa [hrewrite] using hmul

/-- Horizon ordering is always `r_- ≤ r_+` when horizons are real. -/
theorem rMinus_le_rPlus {r_s aStar : ℝ}
    (ha : |aStar| ≤ 1) :
    rMinus r_s aStar ≤ rPlus r_s aStar := by
  unfold rMinus rPlus
  have hdisc_nonneg : 0 ≤ horizonDiscriminant r_s aStar :=
    horizonDiscriminant_nonneg ha
  nlinarith [Real.sqrt_nonneg (horizonDiscriminant r_s aStar)]

/-- Schwarzschild limit (`a_* = 0`): `r_+ = r_s` and `r_- = 0`. -/
theorem schwarzschild_limit_horizons {r_s : ℝ} (hrs : 0 ≤ r_s) :
    rPlus r_s 0 = r_s ∧ rMinus r_s 0 = 0 := by
  have hm_nonneg : 0 ≤ mass r_s := mass_nonneg hrs
  have hsqrt : Real.sqrt ((mass r_s) ^ 2) = mass r_s := by
    simpa [abs_of_nonneg hm_nonneg] using (Real.sqrt_sq_eq_abs (mass r_s))
  constructor
  · unfold rPlus horizonDiscriminant spinLength
    simp
    rw [hsqrt]
    unfold mass
    ring
  · unfold rMinus horizonDiscriminant spinLength
    simp [hsqrt]

/-- Extremal spin (`|a_*| = 1`) gives coincident horizons `r_+ = r_- = M`. -/
theorem extremal_horizons_coincide {r_s aStar : ℝ} (ha : |aStar| = 1) :
    rPlus r_s aStar = mass r_s ∧ rMinus r_s aStar = mass r_s := by
  have habs_sq : |aStar| ^ 2 = (1 : ℝ) ^ 2 := by
    exact congrArg (fun x : ℝ => x ^ 2) ha
  have hsq : aStar ^ 2 = 1 := by
    have hsq' : aStar ^ 2 = (1 : ℝ) ^ 2 := by
      simpa [sq_abs] using habs_sq
    nlinarith [hsq']
  have hdisc : horizonDiscriminant r_s aStar = 0 := by
    unfold horizonDiscriminant spinLength
    nlinarith [hsq]
  constructor
  · unfold rPlus
    simp [hdisc]
  · unfold rMinus
    simp [hdisc]

/-- At the equator `θ = π/2`, the static limit radius is exactly `r_s`. -/
theorem ergosphere_equator_eq_r_s {r_s aStar : ℝ} (hrs : 0 ≤ r_s) :
    ergosphereRadius r_s aStar (Real.pi / 2) = r_s := by
  have hm_nonneg : 0 ≤ mass r_s := mass_nonneg hrs
  unfold ergosphereRadius spinLength
  have hcos : Real.cos (Real.pi / 2) = 0 := by simp
  rw [hcos]
  simp
  have hsqrt : Real.sqrt ((mass r_s) ^ 2) = mass r_s := by
    simpa [abs_of_nonneg hm_nonneg] using (Real.sqrt_sq_eq_abs (mass r_s))
  calc
    mass r_s + Real.sqrt ((mass r_s) ^ 2) = mass r_s + mass r_s := by rw [hsqrt]
    _ = r_s := by
      unfold mass
      ring

/-- At the pole `θ = 0`, static limit and outer horizon coincide. -/
theorem ergosphere_pole_eq_rPlus {r_s aStar : ℝ} :
    ergosphereRadius r_s aStar 0 = rPlus r_s aStar := by
  unfold ergosphereRadius rPlus horizonDiscriminant
  simp [Real.cos_zero]

/-- Non-spinning case has zero horizon angular velocity (`Ω_H = 0`). -/
theorem horizonAngularVelocity_zero_spin {r_s : ℝ} :
    horizonAngularVelocity r_s 0 = 0 := by
  unfold horizonAngularVelocity spinLength
  simp

/-- `r_+` is a root of the Kerr radial polynomial `Δ(r)`. -/
theorem kerrDelta_rPlus_eq_zero {r_s aStar : ℝ}
    (ha : |aStar| ≤ 1) :
    kerrDelta r_s aStar (rPlus r_s aStar) = 0 := by
  have hdisc_nonneg : 0 ≤ horizonDiscriminant r_s aStar :=
    horizonDiscriminant_nonneg ha
  have hdisc_nonneg' : 0 ≤ (mass r_s) ^ 2 - (spinLength r_s aStar) ^ 2 := by
    simpa [horizonDiscriminant] using hdisc_nonneg
  have hsq :
      (Real.sqrt ((mass r_s) ^ 2 - (spinLength r_s aStar) ^ 2)) ^ 2 =
      (mass r_s) ^ 2 - (spinLength r_s aStar) ^ 2 := by
    exact Real.sq_sqrt hdisc_nonneg'
  have hmul :
      Real.sqrt ((mass r_s) ^ 2 - (spinLength r_s aStar) ^ 2) *
      Real.sqrt ((mass r_s) ^ 2 - (spinLength r_s aStar) ^ 2) =
      (mass r_s) ^ 2 - (spinLength r_s aStar) ^ 2 := by
    exact Real.mul_self_sqrt hdisc_nonneg'
  unfold kerrDelta rPlus horizonDiscriminant
  ring_nf
  nlinarith [hsq, hmul]

/-- `r_-` is also a root of the Kerr radial polynomial `Δ(r)`. -/
theorem kerrDelta_rMinus_eq_zero {r_s aStar : ℝ}
    (ha : |aStar| ≤ 1) :
    kerrDelta r_s aStar (rMinus r_s aStar) = 0 := by
  have hdisc_nonneg : 0 ≤ horizonDiscriminant r_s aStar :=
    horizonDiscriminant_nonneg ha
  have hdisc_nonneg' : 0 ≤ (mass r_s) ^ 2 - (spinLength r_s aStar) ^ 2 := by
    simpa [horizonDiscriminant] using hdisc_nonneg
  have hsq :
      (Real.sqrt ((mass r_s) ^ 2 - (spinLength r_s aStar) ^ 2)) ^ 2 =
      (mass r_s) ^ 2 - (spinLength r_s aStar) ^ 2 := by
    exact Real.sq_sqrt hdisc_nonneg'
  have hmul :
      Real.sqrt ((mass r_s) ^ 2 - (spinLength r_s aStar) ^ 2) *
      Real.sqrt ((mass r_s) ^ 2 - (spinLength r_s aStar) ^ 2) =
      (mass r_s) ^ 2 - (spinLength r_s aStar) ^ 2 := by
    exact Real.mul_self_sqrt hdisc_nonneg'
  unfold kerrDelta rMinus horizonDiscriminant
  ring_nf
  nlinarith [hsq, hmul]

/-- Non-spinning case has vanishing frame-dragging (`ω = 0`) everywhere. -/
theorem frameDraggingOmega_zero_spin {r_s r θ : ℝ} :
    frameDraggingOmega r_s 0 r θ = 0 := by
  unfold frameDraggingOmega frameDraggingDenom spinLength kerrDelta
  simp

/-- Master scaffold: physical spin guarantees real, ordered horizons. -/
theorem kerr_horizon_scaffold {r_s aStar : ℝ}
    (ha : |aStar| ≤ 1) :
    0 ≤ horizonDiscriminant r_s aStar ∧ rMinus r_s aStar ≤ rPlus r_s aStar := by
  exact ⟨horizonDiscriminant_nonneg ha, rMinus_le_rPlus ha⟩

end Gutoe.KerrGeometry
