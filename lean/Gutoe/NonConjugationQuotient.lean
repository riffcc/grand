import Mathlib
import Gutoe.RecursiveNavigationNoTranslation

/-!
GUTOE — Non-conjugation quotient lane scaffold

This lane separates two classes:

1. Homogeneous non-conjugation products (`x ↦ A x`) remain origin-preserving.
2. Inhomogeneous extensions (`x ↦ A x + b`) can move origin.

It does not assert a specific Clifford product realization; it keeps the
mechanism classes explicit for falsification.
-/

namespace Gutoe.NonConjugationQuotient

open Gutoe
open Gutoe.RecursiveNavigationNoTranslation
open Gutoe.RecursiveZ3Tower

noncomputable section

/-- Homogeneous non-conjugation action proxy on the recursive layer. -/
def nonConjugationHomogeneous (A : Vec256 →ₗ[ℝ] Vec256) (x : Vec256) : Vec256 :=
  A x

/-- Inhomogeneous non-conjugation action proxy. -/
def nonConjugationAffine
    (A : Vec256 →ₗ[ℝ] Vec256) (b : Vec256) (x : Vec256) : Vec256 :=
  A x + b

/-- Homogeneous non-conjugation remains origin-preserving under descent. -/
theorem homogeneous_nonconj_preserves_origin
    (A : Vec256 →ₗ[ℝ] Vec256) :
    towerProjection (nonConjugationHomogeneous A 0) = 0 := by
  simp [nonConjugationHomogeneous, towerProjection]

/-- Inhomogeneous non-conjugation reaches any descended target that matches
the descended offset. -/
theorem affine_nonconj_reaches_target_if_offset_projects
    (A : Vec256 →ₗ[ℝ] Vec256) (b : Vec256) (x : Fin 4 → ℝ)
    (hx : towerProjection b = x) :
    towerProjection (nonConjugationAffine A b 0) = x := by
  calc
    towerProjection (nonConjugationAffine A b 0) = towerProjection b := by
      simp [nonConjugationAffine, towerProjection]
    _ = x := hx

end
end Gutoe.NonConjugationQuotient
