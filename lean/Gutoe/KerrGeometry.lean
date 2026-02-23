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

lemma mass_nonneg {r_s : ℝ} (hrs : 0 ≤ r_s) : 0 ≤ mass r_s := by
  unfold mass
  nlinarith

lemma abs_aStar_sq_le_one {aStar : ℝ} (ha : |aStar| ≤ 1) : aStar ^ 2 ≤ 1 := by
  nlinarith [sq_le_sq.mpr ha]

/-- For physical spin `|a_*| ≤ 1`, the Kerr horizon discriminant is nonnegative. -/
theorem horizonDiscriminant_nonneg {r_s aStar : ℝ}
    (hrs : 0 ≤ r_s) (ha : |aStar| ≤ 1) :
    0 ≤ horizonDiscriminant r_s aStar := by
  unfold horizonDiscriminant spinLength
  have hm2_nonneg : 0 ≤ (mass r_s) ^ 2 := sq_nonneg (mass r_s)
  have hfac_nonneg : 0 ≤ 1 - aStar ^ 2 := by
    have hs : aStar ^ 2 ≤ 1 := abs_aStar_sq_le_one ha
    linarith
  ring_nf
  exact mul_nonneg hm2_nonneg hfac_nonneg

/-- Horizon ordering is always `r_- ≤ r_+` when horizons are real. -/
theorem rMinus_le_rPlus {r_s aStar : ℝ}
    (hrs : 0 ≤ r_s) (ha : |aStar| ≤ 1) :
    rMinus r_s aStar ≤ rPlus r_s aStar := by
  unfold rMinus rPlus
  have hdisc_nonneg : 0 ≤ horizonDiscriminant r_s aStar :=
    horizonDiscriminant_nonneg hrs ha
  nlinarith [Real.sqrt_nonneg (horizonDiscriminant r_s aStar)]

/-- Schwarzschild limit (`a_* = 0`): `r_+ = r_s` and `r_- = 0`. -/
theorem schwarzschild_limit_horizons {r_s : ℝ} (hrs : 0 ≤ r_s) :
    rPlus r_s 0 = r_s ∧ rMinus r_s 0 = 0 := by
  have hm_nonneg : 0 ≤ mass r_s := mass_nonneg hrs
  unfold rPlus rMinus horizonDiscriminant spinLength
  have hsqrt : Real.sqrt ((mass r_s) ^ 2) = mass r_s := by
    simpa [abs_of_nonneg hm_nonneg] using (Real.sqrt_sq_eq_abs (mass r_s))
  constructor
  · calc
      mass r_s + Real.sqrt ((mass r_s) ^ 2 - 0) = mass r_s + Real.sqrt ((mass r_s) ^ 2) := by ring
      _ = mass r_s + mass r_s := by rw [hsqrt]
      _ = r_s := by
        unfold mass
        ring
  · calc
      mass r_s - Real.sqrt ((mass r_s) ^ 2 - 0) = mass r_s - Real.sqrt ((mass r_s) ^ 2) := by ring
      _ = mass r_s - mass r_s := by rw [hsqrt]
      _ = 0 := by ring

/-- Extremal spin (`|a_*| = 1`) gives coincident horizons `r_+ = r_- = M`. -/
theorem extremal_horizons_coincide {r_s aStar : ℝ} (ha : |aStar| = 1) :
    rPlus r_s aStar = mass r_s ∧ rMinus r_s aStar = mass r_s := by
  have hsq : aStar ^ 2 = 1 := by
    nlinarith [congrArg (fun x : ℝ => x ^ 2) ha]
  unfold rPlus rMinus horizonDiscriminant spinLength
  have hdisc : (mass r_s) ^ 2 - (aStar * mass r_s) ^ 2 = 0 := by
    nlinarith [hsq]
  constructor <;> simp [hdisc]

/-- At the equator `θ = π/2`, the static limit radius is exactly `r_s`. -/
theorem ergosphere_equator_eq_r_s {r_s aStar : ℝ} (hrs : 0 ≤ r_s) :
    ergosphereRadius r_s aStar (Real.pi / 2) = r_s := by
  have hm_nonneg : 0 ≤ mass r_s := mass_nonneg hrs
  unfold ergosphereRadius spinLength
  have hcos : Real.cos (Real.pi / 2) = 0 := by simpa using Real.cos_pi_div_two
  rw [hcos]
  simp
  have hsqrt : Real.sqrt ((mass r_s) ^ 2) = mass r_s := by
    simpa [abs_of_nonneg hm_nonneg] using (Real.sqrt_sq_eq_abs (mass r_s))
  calc
    mass r_s + Real.sqrt ((mass r_s) ^ 2 - (aStar * mass r_s) ^ 2 * 0 ^ 2) =
        mass r_s + Real.sqrt ((mass r_s) ^ 2) := by ring
    _ = mass r_s + mass r_s := by rw [hsqrt]
    _ = r_s := by
      unfold mass
      ring

/-- Master scaffold: physical spin guarantees real, ordered horizons. -/
theorem kerr_horizon_scaffold {r_s aStar : ℝ}
    (hrs : 0 ≤ r_s) (ha : |aStar| ≤ 1) :
    0 ≤ horizonDiscriminant r_s aStar ∧ rMinus r_s aStar ≤ rPlus r_s aStar := by
  exact ⟨horizonDiscriminant_nonneg hrs ha, rMinus_le_rPlus hrs ha⟩

end Gutoe.KerrGeometry
