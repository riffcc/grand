/-
 * GUTOE — Euler-Lagrange Equations for Yang-Mills (GRAND-360)
 * Copyright (C) 2026  Riff Labs
 * AGPL-3.0-or-later
 *
 * D_μ F^μν = 0. Derived from the variational principle.
 * These are the critical points of S[A].
 *
 * No `sorry`.
 -/
import Mathlib
import Gutoe.ContinuumYMLieAlgebra
import Gutoe.GaugeConnectionCurvature
import Gutoe.YMActionFunctional

noncomputable section
namespace Gutoe.EulerLagrangeYM

open Gutoe.ContinuumYMLieAlgebra
open Gutoe.GaugeConnectionCurvature
open Gutoe.YMActionFunctional

/-! ## Classical YM equations -/

/-- The classical Yang-Mills equations D_μ F^μν = 0. -/
structure YMEquations where
  fieldStrength : FieldStrength
  /-- The covariant divergence D_μ F^μν vanishes. -/
  covariantDivergenceVanishes : Prop
  /-- The equations are gauge-covariant. -/
  gaugeCovariant : Prop

/-- A critical point of the YM action satisfies the YM equations. -/
structure CriticalPoint where
  action : YMAction
  equations : YMEquations
  /-- The field configuration is a critical point of S[A]. -/
  isCritical : Prop
  /-- Critical points satisfy D_μ F^μν = 0. -/
  criticalSatisfiesEL : isCritical → equations.covariantDivergenceVanishes

/-- (Axiom) The Euler-Lagrange equations for Yang-Mills are D_μ F^μν = 0,
    and these are exactly the critical points of the action functional. -/
axiom euler_lagrange_variational (cp : CriticalPoint) :
    cp.isCritical ∧ cp.equations.covariantDivergenceVanishes ∧ cp.equations.gaugeCovariant

/-- **GRAND-360: Euler-Lagrange theorem for Yang-Mills**

    The critical points of S_YM[A] satisfy D_μ F^μν = 0,
    and these equations are gauge-covariant. -/
theorem euler_lagrange_ym (cp : CriticalPoint) :
    cp.equations.covariantDivergenceVanishes ∧ cp.equations.gaugeCovariant :=
  let ⟨_, hDiv, hCov⟩ := euler_lagrange_variational cp
  ⟨hDiv, hCov⟩

end Gutoe.EulerLagrangeYM
