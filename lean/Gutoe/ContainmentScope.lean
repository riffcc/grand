import Mathlib
import Gutoe.ProjectionFibers

/-!
GUTOE — Containment Scope (law-level vs state-level)

This module sharpens the containment claim into provable statements:

1. Every fiber is an affine translate of a single shared kernel.
2. All fibers are nonempty and share the same kernel rank (`12`).
3. Full global state reconstruction from the 4D projection is impossible
   (the projection is non-injective, so no left inverse exists).
-/

namespace Gutoe.ContainmentScope

open Gutoe
open Gutoe.ProjectionFibers

noncomputable section

/-- Canonical basepoint in the fiber over `x`. -/
def fiberBase (x : Fin 4 → ℝ) : Vec16 := grade1Section x

/-- The canonical fiber basepoint projects back to `x`. -/
theorem grade1Projection_fiberBase (x : Fin 4 → ℝ) :
    grade1Projection (fiberBase x) = x := by
  -- Reuse the explicit section computation from the surjectivity witness.
  funext i
  have hsum :
      (∑ x' : Fin 4, (if (↑i : ℕ) = ↑x' then x x' else 0)) = x i := by
    have hsum' : (∑ x' : Fin 4, (if i = x' then x x' else 0)) = x i := by
      simp
    simpa [Fin.ext_iff] using hsum'
  simpa [fiberBase, grade1Projection, grade1Section, railBasisVec, railIndex4] using hsum

/-- Fiber membership is exactly "subtract the canonical basepoint and land in
the shared kernel". -/
theorem fiber_membership_iff_sub_base_in_kernel
    (x : Fin 4 → ℝ) (v : Vec16) :
    v ∈ fiberAt x ↔ v - fiberBase x ∈ LinearMap.ker grade1Projection := by
  constructor
  · intro hv
    change grade1Projection (v - fiberBase x) = 0
    have hv' : grade1Projection v = x := hv
    calc
      grade1Projection (v - fiberBase x)
          = grade1Projection v - grade1Projection (fiberBase x) := by simp
      _ = x - x := by simpa [hv', grade1Projection_fiberBase]
      _ = 0 := by simp
  · intro hk
    change grade1Projection v = x
    have hk' : grade1Projection (v - fiberBase x) = 0 := hk
    have hbase : grade1Projection (fiberBase x) = x := grade1Projection_fiberBase x
    have hcalc : grade1Projection v - x = 0 := by
      calc
        grade1Projection v - x
            = grade1Projection v - grade1Projection (fiberBase x) := by simpa [hbase]
        _ = grade1Projection (v - fiberBase x) := by simp
        _ = 0 := hk'
    exact sub_eq_zero.mp hcalc

/-- Set-level form: each fiber is an affine translate of the same kernel. -/
theorem fiber_eq_base_plus_kernel (x : Fin 4 → ℝ) :
    fiberAt x = {v : Vec16 | v - fiberBase x ∈ LinearMap.ker grade1Projection} := by
  ext v
  exact fiber_membership_iff_sub_base_in_kernel x v

/-- Every fiber is nonempty (there is always a canonical basepoint). -/
theorem every_fiber_nonempty (x : Fin 4 → ℝ) : (fiberAt x).Nonempty := by
  refine ⟨fiberBase x, ?_⟩
  exact grade1Projection_fiberBase x

/-- Translation by basepoint delta identifies fibers over `x` and `y`. -/
theorem fibers_translate_by_base_delta
    (x y : Fin 4 → ℝ) (v : Vec16) :
    v ∈ fiberAt x ↔ v + (fiberBase y - fiberBase x) ∈ fiberAt y := by
  rw [fiber_membership_iff_sub_base_in_kernel x v,
      fiber_membership_iff_sub_base_in_kernel y (v + (fiberBase y - fiberBase x))]
  have hEq : v + (fiberBase y - fiberBase x) - fiberBase y = v - fiberBase x := by
    simp [sub_eq_add_neg, add_assoc, add_comm, add_left_comm]
  simpa [hEq]

/-- All fibers share the same kernel rank (`12`). -/
theorem fibers_share_kernel_rank_12 (x : Fin 4 → ℝ) :
    Module.finrank ℝ (LinearMap.ker grade1Projection) = 12 := by
  exact grade1Projection_kernel_finrank

/-- No full global state reconstructor can exist from `grade1Projection` alone.
This is the formal state-level limit behind non-injectivity. -/
theorem no_global_state_reconstructor :
    ¬ ∃ rec : (Fin 4 → ℝ) → Vec16, Function.LeftInverse rec grade1Projection := by
  intro h
  rcases h with ⟨rec, hleft⟩
  have hinj : Function.Injective grade1Projection := Function.LeftInverse.injective hleft
  exact grade1Projection_not_injective hinj

end
end Gutoe.ContainmentScope
