import Mathlib
import Gutoe.CreationLanes
import Gutoe.DynamicTopologyCreation
import Gutoe.RecursiveNavigationNoTranslation
import Gutoe.RecursiveZ3Tower
import Gutoe.ProjectionFibers

/-!
GUTOE — Topological Access Controllability

This module answers a scoped question:

Can adding a non-linear/topological access primitive (affine offset lane) evade
the existing linear no-go for nonzero descended targets?

Result:
- Linear-only lane: no (already proven).
- Affine-offset lane: yes for 4D descended targets (surjective).
- This does not assert full 16D/256D state reconstruction.
-/

namespace Gutoe.TopologicalAccessControllability

open Gutoe
open Gutoe.CreationLanes
open Gutoe.DynamicTopologyCreation
open Gutoe.RecursiveNavigationNoTranslation
open Gutoe.RecursiveZ3Tower
open Gutoe.ProjectionFibers

noncomputable section

/-- `towerProjection : Vec256 -> (Fin 4 -> ℝ)` is surjective. -/
theorem tower_projection_surjective :
    Function.Surjective towerProjection := by
  intro x
  rcases grade1Projection_surjective x with ⟨w, hw⟩
  rcases proj256to16_surjective w with ⟨v, hv⟩
  refine ⟨v, ?_⟩
  unfold towerProjection
  calc
    grade1Projection (proj256to16 v)
        = grade1Projection w := by simpa [hv]
    _ = x := hw

/-- Affine-offset lane can reach any descended 4D target from origin
in one step (with `L = 0` and suitable offset). -/
theorem affine_origin_reaches_any_descended_target
    (x : Fin 4 → ℝ) :
    ∃ t : Vec256, towerProjection (affineStep (0 : Vec256 →ₗ[ℝ] Vec256) t 0) = x := by
  rcases tower_projection_surjective x with ⟨t, ht⟩
  refine ⟨t, ?_⟩
  exact affine_origin_reaches_target_if_offset_projects
      (L := (0 : Vec256 →ₗ[ℝ] Vec256)) (t := t) (x := x) ht

/-- Linear-only lane cannot reach nonzero targets; affine lane can.
This is the explicit bypass theorem for descended 4D coordinates. -/
theorem affine_bypasses_linear_no_go
    (x : Fin 4 → ℝ) (hx : x ≠ 0) :
    (∀ L : Vec256 →ₗ[ℝ] Vec256, towerProjection (L 0) ≠ x) ∧
    (∃ t : Vec256, towerProjection (affineStep (0 : Vec256 →ₗ[ℝ] Vec256) t 0) = x) := by
  refine ⟨?_, affine_origin_reaches_any_descended_target x⟩
  intro L
  exact no_linear_origin_to_nonzero_target L x hx

/-- If the dynamic local-creation gate is open, we can have both:
1) a local nontrivial identified time shift, and
2) affine controllability to an arbitrary descended 4D target.

This is a compositional statement; it does not claim full-state or
macroscopic engineering realizability. -/
theorem dynamic_gate_with_affine_targetability
    (budget radius period x_local : ℝ)
    (hGate : dynamicCreationGate budget radius period)
    (hxLocal : |x_local| ≤ radius)
    (target : Fin 4 → ℝ) :
    (∃ a b : Gutoe.CTCLegality.Event,
      sameOnLocalPatch period radius a b ∧ b.t ≠ a.t) ∧
    (∃ t : Vec256, towerProjection (affineStep (0 : Vec256 →ₗ[ℝ] Vec256) t 0) = target) := by
  refine ⟨?_, affine_origin_reaches_any_descended_target target⟩
  exact dynamic_gate_implies_local_shift budget radius period x_local hGate hxLocal

end
end Gutoe.TopologicalAccessControllability
