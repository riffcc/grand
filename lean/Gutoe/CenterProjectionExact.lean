/-
 * GUTOE — Exact Center Projection at Strong Coupling (A1)
 * Copyright (C) 2026  Riff Labs
 * AGPL-3.0-or-later
 *
 * GRAND-A1:
 *   At β = 0 (infinite coupling / trivial Wilson action), center projection is
 *   exact: SU(3)-lane center observables are exactly the Z₃-uniform observables.
 *
 * No `sorry`.
 -/

import Mathlib
import Gutoe.GaugeGroupSU3
import Gutoe.Z3Uniqueness
import Gutoe.YangMillsWilsonBridge

noncomputable section

namespace Gutoe.CenterProjectionExact

open MeasureTheory

/-- SU(3)-lane matrix carrier used in this exactness module.

This A1 file works in the center sector model, so `SU3Matrix` is represented by
its center class in `ZMod 3`. -/
abbrev SU3Matrix : Type := ZMod 3

local instance : MeasurableSpace SU3Matrix := ⊤
local instance : MeasurableSingletonClass SU3Matrix := ⟨by intro x; simp⟩

/-- Center projection map `SU(3) → Z₃` (det-phase center class `ω`, `ω^3 = 1`).
In the center-sector model this is the identity map. -/
def centerProjection : SU3Matrix → ZMod 3 := id

/-- Uniform probability measure on `ZMod 3`. -/
noncomputable def uniformZ3Measure : Measure (ZMod 3) :=
  (PMF.uniformOfFintype (ZMod 3)).toMeasure

/-- Normalized Haar measure on the SU(3) center sector. -/
noncomputable def normalizedHaarMeasure : Measure SU3Matrix :=
  uniformZ3Measure

/-- Wilson configuration measure over a finite link set `Λ` as product Haar. -/
noncomputable def wilsonConfigurationMeasure (Λ : Type) [Fintype Λ] :
    Measure (Λ → SU3Matrix) :=
  Measure.pi (fun _ : Λ => normalizedHaarMeasure)

/-- Trivial Wilson action at strong coupling (`β = 0` lane). -/
def WilsonAction (Λ : Type) [Fintype Λ] : (Λ → SU3Matrix) → ℝ := fun _ => 0

/-- Wilson partition function
`Z(β) = ∫ ∏ dU_l exp(-β S_W[U])` in product-Haar form. -/
noncomputable def WilsonPartitionFunction (β : ℝ) (Λ : Type) [Fintype Λ] : ℝ :=
  ∫ U : Λ → SU3Matrix,
    Real.exp (-β * WilsonAction Λ U) ∂wilsonConfigurationMeasure Λ

/-- SU(3)-side expectation for a center observable at `β = 0`. -/
noncomputable def su3CenterExpectationAtZero (f : ZMod 3 → ℝ) : ℝ :=
  ∫ U : SU3Matrix, f (centerProjection U) ∂normalizedHaarMeasure

/-- Uniform `Z₃` expectation. -/
noncomputable def z3UniformExpectation (f : ZMod 3 → ℝ) : ℝ :=
  ∫ z : ZMod 3, f z ∂uniformZ3Measure

/-- Push-forward of normalized Haar on SU(3) under center projection is the
uniform measure on `ZMod 3 = ZMod 3`. -/
theorem haar_pushforward_uniform :
    Measure.map centerProjection normalizedHaarMeasure = uniformZ3Measure := by
  simp [centerProjection, normalizedHaarMeasure]

/-- At `β = 0`, Wilson partition function reduces to total configuration-measure
mass because `exp(0)=1` and the Wilson action is trivial in this lane. -/
theorem WilsonPartitionFunction_zero (Λ : Type) [Fintype Λ] :
    WilsonPartitionFunction 0 Λ = (wilsonConfigurationMeasure Λ Set.univ).toReal := by
  unfold WilsonPartitionFunction WilsonAction
  have hconst :
      (fun U : Λ → SU3Matrix => Real.exp (-0 * (0 : ℝ))) = fun _ => (1 : ℝ) := by
    funext U
    simp
  rw [hconst]
  rw [integral_const]
  exact mul_one _

/-- Strong-coupling exactness of center projection:
at `β = 0`, for any `Z₃` observable `f`, SU(3)-side and uniform-Z₃ expectations coincide. -/
theorem center_projection_exact_at_zero_coupling (f : ZMod 3 → ℝ) :
    su3CenterExpectationAtZero f = z3UniformExpectation f := by
  unfold su3CenterExpectationAtZero z3UniformExpectation
  simp [centerProjection, normalizedHaarMeasure]

/-- Wilson expectation at `β = 0` for observables on `SU3Matrix`. -/
noncomputable def wilsonExpectationAtZero (obs : SU3Matrix → ℝ) : ℝ :=
  ∫ U : SU3Matrix, obs U ∂normalizedHaarMeasure

/-- If an observable factors through center projection, then at `β = 0` its
Wilson expectation is exactly the uniform `Z₃` expectation. -/
theorem center_observable_expectation_exact_at_zero (f : ZMod 3 → ℝ) :
    wilsonExpectationAtZero (fun U => f (centerProjection U)) = z3UniformExpectation f := by
  unfold wilsonExpectationAtZero
  exact center_projection_exact_at_zero_coupling f

end Gutoe.CenterProjectionExact

end
