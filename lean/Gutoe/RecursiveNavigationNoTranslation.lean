import Mathlib
import Gutoe.RailSpace
import Gutoe.ProjectionFibers
import Gutoe.RecursiveZ3Tower

/-!
GUTOE — Recursive navigation guard (linear descent)

This module formalizes a hard guard for the current linear
`256 -> 16 -> 4` projection stack:

- Any purely linear navigation step in the 256-layer preserves the origin after
  descent to 4D.
- Therefore a nonzero 4D target cannot be reached from the origin by linear-only
  multiplicative routing.

This does not forbid nonlinear/topological routes; it closes the linear lane.
-/

namespace Gutoe.RecursiveNavigationNoTranslation

open Gutoe
open Gutoe.ProjectionFibers
open Gutoe.RecursiveZ3Tower

noncomputable section

/-- Composite tower descent map `256 -> 4` via `256 -> 16` then `16 -> 4`. -/
def towerProjection : Vec256 →ₗ[ℝ] (Fin 4 → ℝ) :=
  grade1Projection.comp proj256to16

/-- Any linear step in the recursive layer keeps the descended origin fixed. -/
theorem descended_linear_step_preserves_origin (L : Vec256 →ₗ[ℝ] Vec256) :
    towerProjection (L 0) = 0 := by
  simp [towerProjection]

/-- Linear-only recursive routing cannot send origin to a nonzero 4D target. -/
theorem no_linear_origin_to_nonzero_target
    (L : Vec256 →ₗ[ℝ] Vec256) (x : Fin 4 → ℝ) (hx : x ≠ 0) :
    towerProjection (L 0) ≠ x := by
  intro hEq
  have hx0 : x = 0 := by
    calc
      x = towerProjection (L 0) := by simpa using hEq.symm
      _ = 0 := by simp [towerProjection]
  exact hx hx0

/-- Reachability from origin under linear descent is exactly the zero target. -/
theorem linear_origin_reachability_iff_zero
    (L : Vec256 →ₗ[ℝ] Vec256) (x : Fin 4 → ℝ) :
    towerProjection (L 0) = x ↔ x = 0 := by
  constructor
  · intro h
    calc
      x = towerProjection (L 0) := by simpa using h.symm
      _ = 0 := by simp [towerProjection]
  · intro hx
    simpa [towerProjection, hx]

end
end Gutoe.RecursiveNavigationNoTranslation
