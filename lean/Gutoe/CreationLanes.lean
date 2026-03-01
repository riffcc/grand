import Mathlib
import Gutoe.RecursiveNavigationNoTranslation
import Gutoe.CTCLegality

/-!
GUTOE — Creation Lanes (toy closure)

This module captures two explicit "creation" lanes that remain after homogeneous
linear no-go results:

1. Affine lane: adding an inhomogeneous offset can move the descended origin.
2. Local-topology lane: a compactly-supported identification patch can admit
   nontrivial local shifts while forbidding them outside support.

These are structural statements for lane classification, not engineering claims.
-/

namespace Gutoe.CreationLanes

open Gutoe
open Gutoe.RecursiveNavigationNoTranslation
open Gutoe.RecursiveZ3Tower
open Gutoe.CTCLegality

noncomputable section

-- ── Affine creation lane on the recursive tower ──────────────────────────────

/-- Affine update on the 256-layer state space. -/
def affineStep
    (L : RecursiveZ3Tower.Vec256 →ₗ[ℝ] RecursiveZ3Tower.Vec256)
    (t : RecursiveZ3Tower.Vec256)
    (v : RecursiveZ3Tower.Vec256) : RecursiveZ3Tower.Vec256 :=
  L v + t

/-- Descended image of the affine-origin state is exactly the descended offset. -/
theorem affine_origin_descends_as_offset
    (L : RecursiveZ3Tower.Vec256 →ₗ[ℝ] RecursiveZ3Tower.Vec256)
    (t : RecursiveZ3Tower.Vec256) :
    towerProjection (affineStep L t 0) = towerProjection t := by
  unfold affineStep
  simp [towerProjection]

/-- Any target that is realized as the descent of an affine offset is reachable
from origin in one affine step. -/
theorem affine_origin_reaches_target_if_offset_projects
    (L : RecursiveZ3Tower.Vec256 →ₗ[ℝ] RecursiveZ3Tower.Vec256)
    (t : RecursiveZ3Tower.Vec256) (x : Fin 4 → ℝ)
    (hx : towerProjection t = x) :
    towerProjection (affineStep L t 0) = x := by
  rw [affine_origin_descends_as_offset, hx]

-- ── Compact local identification patch lane ──────────────────────────────────

/-- Time-cylinder identification localized to a compact spatial patch `|x| ≤ R`. -/
def sameOnLocalPatch (T R : ℝ) (a b : Event) : Prop :=
  |a.x| ≤ R ∧ |b.x| ≤ R ∧ sameOnTimeCylinder T a b

/-- Inside the patch, a nontrivial identified shift exists (for `T > 0`). -/
theorem local_patch_nontrivial_shift_exists
    (T R x : ℝ) (hT : T > 0) (hx : |x| ≤ R) :
    ∃ a b : Event, sameOnLocalPatch T R a b ∧ b.t ≠ a.t := by
  refine ⟨⟨0, x⟩, ⟨T, x⟩, ?_, ?_⟩
  · refine ⟨hx, hx, ?_⟩
    unfold sameOnTimeCylinder
    refine ⟨rfl, ?_⟩
    refine ⟨1, ?_⟩
    ring
  · have hne : T ≠ 0 := by linarith
    simpa using hne

/-- Outside the compact support, the localized relation cannot hold on fixed-`x`
events carrying a nontrivial time shift. -/
theorem local_patch_no_nontrivial_shift_outside
    (T R x : ℝ) (hx : R < |x|) :
    ¬ ∃ a b : Event,
      sameOnLocalPatch T R a b ∧ a.x = x ∧ b.x = x ∧ b.t ≠ a.t := by
  intro h
  rcases h with ⟨a, b, hPatch, hax, hbx, _hne⟩
  rcases hPatch with ⟨haR, hbR, _hid⟩
  have hxa : |x| ≤ R := by simpa [hax] using haR
  exact not_le_of_gt hx hxa

end
end Gutoe.CreationLanes
