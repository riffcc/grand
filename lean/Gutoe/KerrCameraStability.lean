/-
 * GUTOE — Kerr Camera Stability Constraints
 * Copyright (C) 2026 Riff Labs
 * AGPL-3.0-or-later
 *
 * Formal constraints for the Kerr image-plane mapping so renderer behavior is
 * stable and physically interpretable.
 *
 * No `sorry`.
-/

import Mathlib
import Gutoe.KerrGeometry

namespace Gutoe.KerrCameraStability

open Real

/-- Kerr camera constants (Bardeen/Carter form):
`ξ = -α sin θ_obs`, `η = β² + (α² - a²) cos² θ_obs`. -/
noncomputable def kerrXi (alpha thetaObs : ℝ) : ℝ :=
  -alpha * Real.sin thetaObs

noncomputable def kerrEta (alpha beta a thetaObs : ℝ) : ℝ :=
  beta ^ 2 + (alpha ^ 2 - a ^ 2) * (Real.cos thetaObs) ^ 2

/-- If screen `β` was already projected by `sin i`, applying `θ_obs = i` again
introduces an extra suppression in the `β` contribution (`sin² i` factor). -/
noncomputable def preprojectedBeta (beta inc : ℝ) : ℝ :=
  beta * Real.sin inc

/-- Equatorial-observer limit (`θ_obs = π/2`): `ξ = -α`. -/
theorem xi_equatorial (alpha : ℝ) :
    kerrXi alpha (Real.pi / 2) = -alpha := by
  unfold kerrXi
  simp

/-- Equatorial-observer limit (`θ_obs = π/2`): `η = β²`. -/
theorem eta_equatorial (alpha beta a : ℝ) :
    kerrEta alpha beta a (Real.pi / 2) = beta ^ 2 := by
  unfold kerrEta
  simp

/-- `η` is always nonnegative in the equatorial-observer limit. -/
theorem eta_equatorial_nonneg (alpha beta a : ℝ) :
    0 ≤ kerrEta alpha beta a (Real.pi / 2) := by
  rw [eta_equatorial]
  exact sq_nonneg beta

/-- `η` is even in `β` (camera up/down symmetry in constants). -/
theorem eta_even_in_beta (alpha beta a thetaObs : ℝ) :
    kerrEta alpha (-beta) a thetaObs = kerrEta alpha beta a thetaObs := by
  unfold kerrEta
  ring_nf

/-- Double-tilt expansion:
if `β` is already preprojected by `sin i`, and `θ_obs = i` is also used in the
Kerr constants, the `β` term is suppressed by `sin² i`.

This is the algebraic signature of over-application of inclination.
-/
theorem eta_double_tilt_expansion (alpha beta a inc : ℝ) :
    kerrEta alpha (preprojectedBeta beta inc) a inc
      = beta ^ 2 * (Real.sin inc) ^ 2
        + (alpha ^ 2 - a ^ 2) * (Real.cos inc) ^ 2 := by
  unfold kerrEta preprojectedBeta
  ring

/-- In the same double-tilt setup, `ξ = -α sin i` rather than `-α`.
For small `i`, horizontal impact constants are also compressed.
-/
theorem xi_double_tilt (alpha inc : ℝ) :
    kerrXi alpha inc = -alpha * Real.sin inc := by
  unfold kerrXi
  ring

/-- Face-on limit (`i = 0`): `ξ = 0`. Horizontal impact constants collapse. -/
theorem xi_face_on (alpha : ℝ) :
    kerrXi alpha 0 = 0 := by
  unfold kerrXi
  simp

/-- Face-on double-tilt limit (`i = 0`): the `β` contribution is removed
because `sin 0 = 0`, leaving only the `(α,a)` term in `η`. -/
theorem eta_double_tilt_face_on (alpha beta a : ℝ) :
    kerrEta alpha (preprojectedBeta beta 0) a 0 = alpha ^ 2 - a ^ 2 := by
  unfold kerrEta preprojectedBeta
  simp

/-- Preprojecting `β` by `sin i` never increases its squared magnitude. -/
theorem preprojected_beta_sq_le (beta inc : ℝ) :
    (preprojectedBeta beta inc) ^ 2 ≤ beta ^ 2 := by
  have hsin_id : Real.sin inc ^ 2 + Real.cos inc ^ 2 = 1 := Real.sin_sq_add_cos_sq inc
  have hcos_nonneg : 0 ≤ Real.cos inc ^ 2 := sq_nonneg (Real.cos inc)
  have hsin_le_one : Real.sin inc ^ 2 ≤ 1 := by
    nlinarith [hsin_id, hcos_nonneg]
  have hbeta_nonneg : 0 ≤ beta ^ 2 := sq_nonneg beta
  calc
    (preprojectedBeta beta inc) ^ 2 = beta ^ 2 * (Real.sin inc) ^ 2 := by
      unfold preprojectedBeta
      ring
    _ ≤ beta ^ 2 * 1 := by
      exact mul_le_mul_of_nonneg_left hsin_le_one hbeta_nonneg
    _ = beta ^ 2 := by ring

/-- Stable compatibility constraint used by renderer plumbing:
if screen coordinates are already inclination-projected, use equatorial observer
constants to avoid re-applying inclination in Kerr mapping. -/
theorem projected_screen_equatorial_constants (alpha beta a inc : ℝ) :
    kerrXi alpha (Real.pi / 2) = -alpha ∧
    kerrEta alpha (preprojectedBeta beta inc) a (Real.pi / 2)
      = beta ^ 2 * (Real.sin inc) ^ 2 := by
  constructor
  · exact xi_equatorial alpha
  · unfold kerrEta preprojectedBeta
    have hcos : Real.cos (Real.pi / 2) = 0 := by simp
    simp [hcos]
    ring_nf

end Gutoe.KerrCameraStability
