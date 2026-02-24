/- 
 * GUTOE — 3D Camera Ray Projection Invariants
 * Copyright (C) 2026 Riff Labs
 * AGPL-3.0-or-later
 *
 * This module formalizes core 3D pinhole-camera invariants used by the live
 * renderer. It is the Lean scaffold for moving from 2D impact-parameter
 * shortcuts to true 3D ray construction.
 *
 * No `sorry`.
 -/

import Mathlib
import Gutoe.KerrCameraStability

namespace Gutoe.Geodesic3DProjection

open Real

/-- 3D vector represented by coordinates `(x,y,z)`. -/
abbrev R3 := Fin 3 → ℝ

/-- Coordinate accessors. -/
noncomputable def xComp (v : R3) : ℝ := v 0
noncomputable def yComp (v : R3) : ℝ := v 1
noncomputable def zComp (v : R3) : ℝ := v 2

/-- Unnormalized pinhole-camera ray through image-plane coordinates `(α, β)`. -/
noncomputable def rayVec (α β : ℝ) : R3
| 0 => α
| 1 => β
| _ => 1

/-- Squared image-plane impact radius. -/
noncomputable def impactRadiusSq (α β : ℝ) : ℝ := α ^ 2 + β ^ 2

/-- Image-plane impact radius (nonnegative by construction). -/
noncomputable def impactRadius (α β : ℝ) : ℝ := Real.sqrt (impactRadiusSq α β)

/-- Euclidean squared norm specialized to `R3`. -/
noncomputable def normSq3 (v : R3) : ℝ := (xComp v) ^ 2 + (yComp v) ^ 2 + (zComp v) ^ 2

/-- Ray squared norm before normalization. -/
noncomputable def rayNormSq (α β : ℝ) : ℝ := normSq3 (rayVec α β)

/-- Unit-ray direction from `(α, β)` and focal length `1`. -/
noncomputable def rayDir (α β : ℝ) : R3 := fun i => rayVec α β i / Real.sqrt (rayNormSq α β)

theorem rayNormSq_eval (α β : ℝ) :
    rayNormSq α β = impactRadiusSq α β + 1 := by
  unfold rayNormSq normSq3 impactRadiusSq xComp yComp zComp rayVec
  ring

theorem rayNormSq_pos (α β : ℝ) : 0 < rayNormSq α β := by
  rw [rayNormSq_eval]
  have hnonneg : 0 ≤ α ^ 2 + β ^ 2 := by nlinarith [sq_nonneg α, sq_nonneg β]
  exact add_pos_of_nonneg_of_pos hnonneg zero_lt_one

theorem impactRadiusSq_nonneg (α β : ℝ) : 0 ≤ impactRadiusSq α β := by
  unfold impactRadiusSq
  nlinarith [sq_nonneg α, sq_nonneg β]

theorem impactRadius_nonneg (α β : ℝ) : 0 ≤ impactRadius α β := by
  unfold impactRadius
  exact Real.sqrt_nonneg _

theorem impactRadius_sq (α β : ℝ) :
    impactRadius α β ^ 2 = impactRadiusSq α β := by
  unfold impactRadius
  exact Real.sq_sqrt (impactRadiusSq_nonneg α β)

/-- Reflection across the horizontal axis leaves impact radius unchanged. -/
theorem impactRadius_even_beta (α β : ℝ) :
    impactRadius α (-β) = impactRadius α β := by
  unfold impactRadius impactRadiusSq
  ring_nf

/-- `rayDir` is unit-normalized in `normSq3`. -/
theorem rayDir_unit_normSq (α β : ℝ) :
    normSq3 (rayDir α β) = 1 := by
  have hpos : 0 ≤ rayNormSq α β := le_of_lt (rayNormSq_pos α β)
  have hsqrt_nonzero : Real.sqrt (rayNormSq α β) ≠ 0 := by
    exact ne_of_gt (Real.sqrt_pos.mpr (rayNormSq_pos α β))
  have hsqrt_sq : (Real.sqrt (rayNormSq α β)) ^ 2 = rayNormSq α β := by
    exact Real.sq_sqrt hpos
  unfold normSq3 xComp yComp zComp rayDir
  field_simp [hsqrt_nonzero, hsqrt_sq]
  rw [hsqrt_sq]
  rw [rayNormSq_eval]
  unfold rayVec
  unfold impactRadiusSq
  ring

/-- Positive forward component: camera rays point into `+z` half-space. -/
theorem rayDir_z_positive (α β : ℝ) : 0 < zComp (rayDir α β) := by
  unfold zComp rayDir rayVec
  exact one_div_pos.mpr (Real.sqrt_pos.mpr (rayNormSq_pos α β))

/-- `kerrXi` depends only on the horizontal screen coordinate `α`, not `β`. -/
theorem kerrXi_beta_invariant (α β₁ β₂ θObs : ℝ) :
    Gutoe.KerrCameraStability.kerrXi (xComp (rayVec α β₁)) θObs
      = Gutoe.KerrCameraStability.kerrXi (xComp (rayVec α β₂)) θObs := by
  unfold xComp rayVec
  simp

/-- Equatorial-observer Kerr limit for 3D pinhole rays:
`η = β²` when `θ_obs = π/2`. -/
theorem kerrEta_equatorial_from_ray (α β : ℝ) :
    Gutoe.KerrCameraStability.kerrEta (xComp (rayVec α β)) (yComp (rayVec α β)) 0 (Real.pi / 2)
      = β ^ 2 := by
  unfold xComp yComp rayVec
  simpa using (Gutoe.KerrCameraStability.eta_equatorial α β 0)

end Gutoe.Geodesic3DProjection
