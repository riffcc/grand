/- 
 * GUTOE — Kerr Tracer Equation Parity
 * Copyright (C) 2026 Riff Labs
 * AGPL-3.0-or-later
 *
 * Lean-side algebraic parity for the Rust/CUDA Kerr tracer core equations.
 * No `sorry`.
-/

import Mathlib
import Gutoe.KerrGeometry
import Gutoe.KerrCameraStability

namespace Gutoe.KerrTracerEquations

open Real
open Gutoe.KerrGeometry
open Gutoe.KerrCameraStability

/-- Kerr image-plane constants, parity aliases for tracer notation. -/
noncomputable def xi (alpha thetaObs : ℝ) : ℝ := kerrXi alpha thetaObs
noncomputable def eta (alpha beta a thetaObs : ℝ) : ℝ := kerrEta alpha beta a thetaObs

/-- Kerr radial potential used by the tracer (`E = 1` Carter form). -/
noncomputable def radialPotential (r xi eta r_s aStar : ℝ) : ℝ :=
  let a := spinLength r_s aStar
  let delta := kerrDelta r_s aStar r
  let t := (r ^ 2 + a ^ 2) - a * xi
  t ^ 2 - delta * ((xi - a) ^ 2 + eta)

/-- Kerr polar potential used by the tracer (`E = 1`, cot² form). -/
noncomputable def polarPotential (theta xi eta r_s aStar : ℝ) : ℝ :=
  let a := spinLength r_s aStar
  eta + a ^ 2 * (Real.cos theta) ^ 2 - xi ^ 2 * (Real.cos theta / Real.sin theta) ^ 2

theorem xi_eq_impl (alpha thetaObs : ℝ) :
    xi alpha thetaObs = -alpha * Real.sin thetaObs := by
  rfl

theorem eta_eq_impl (alpha beta a thetaObs : ℝ) :
    eta alpha beta a thetaObs = beta ^ 2 + (alpha ^ 2 - a ^ 2) * (Real.cos thetaObs) ^ 2 := by
  rfl

theorem xi_beta_invariant (alpha _beta₁ _beta₂ thetaObs : ℝ) :
    xi alpha thetaObs = xi alpha thetaObs := by
  rfl

theorem eta_equatorial (alpha beta a : ℝ) :
    eta alpha beta a (Real.pi / 2) = beta ^ 2 := by
  simpa [eta] using Gutoe.KerrCameraStability.eta_equatorial alpha beta a

/-- Schwarzschild limit of the radial potential (`a* = 0`) used in code parity checks. -/
theorem radialPotential_schwarzschild_limit (r xi eta r_s : ℝ) :
    radialPotential r xi eta r_s 0
      = r ^ 4 - (r ^ 2 - r_s * r) * (xi ^ 2 + eta) := by
  unfold radialPotential kerrDelta
  simp [spinLength, mass]
  ring

/-- Equatorial observer + zero spin gives `eta = beta²`, matching renderer mapping. -/
theorem eta_equatorial_zero_spin (alpha beta thetaObs : ℝ)
    (hθ : thetaObs = Real.pi / 2) :
    eta alpha beta 0 thetaObs = beta ^ 2 := by
  subst hθ
  simpa [eta] using Gutoe.KerrCameraStability.eta_equatorial alpha beta 0

end Gutoe.KerrTracerEquations
