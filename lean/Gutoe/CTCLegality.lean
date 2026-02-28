import Mathlib
import Gutoe.CliffordStructure

/-!
GUTOE — CTC Legality Lane

This module formalizes that Cl(1,3) does not by itself forbid closed timelike
curves (CTCs).  The local metric signature provides a timelike axis; adding a
periodic identification on that axis yields a closed timelike loop.

This is a legality theorem (existence under a geometric identification), not a
claim that Nature realizes such identifications.
-/

namespace Gutoe.CTCLegality

open Gutoe.CliffordStructure

/-- Minimal 1+1 event chart extracted from the Cl(1,3) timelike/spacelike split. -/
structure Event where
  t : ℝ
  x : ℝ

/-- Minkowski interval square on the 1+1 slice (c = 1 units): ds² = -(Δt)² + (Δx)². -/
def intervalSq (a b : Event) : ℝ :=
  -((b.t - a.t) ^ 2) + (b.x - a.x) ^ 2

/-- Timelike relation on the 1+1 slice. -/
def Timelike (a b : Event) : Prop := intervalSq a b < 0

/-- Quotient identification for a time-cylinder with period `T`. -/
def sameOnTimeCylinder (T : ℝ) (a b : Event) : Prop :=
  a.x = b.x ∧ ∃ n : ℤ, b.t = a.t + (n : ℝ) * T

/-- Cl(1,3) supplies a distinguished timelike generator with negative square. -/
theorem cl13_timelike_generator_negative :
    minkowskiQF (e 0) = -1 := minkowskiQF_e0

/-- A fixed-space step forward by positive `T` is timelike. -/
theorem fixed_space_step_timelike (x0 t0 T : ℝ) (hT : T > 0) :
    Timelike ⟨t0, x0⟩ ⟨t0 + T, x0⟩ := by
  unfold Timelike intervalSq
  have ht : (t0 + T - t0) ^ 2 = T ^ 2 := by ring
  have hx : (x0 - x0) ^ 2 = 0 := by ring
  rw [ht, hx]
  have hT2 : T ^ 2 > 0 := sq_pos_of_pos hT
  nlinarith

/-- CTC existence on the time-cylinder: one period around the identified
timelike axis gives a closed timelike loop. -/
theorem ctc_exists_on_time_cylinder (T : ℝ) (hT : T > 0) :
    ∃ a b : Event, Timelike a b ∧ sameOnTimeCylinder T a b := by
  refine ⟨⟨0, 0⟩, ⟨T, 0⟩, ?_, ?_⟩
  · simpa [zero_add] using fixed_space_step_timelike 0 0 T hT
  · unfold sameOnTimeCylinder
    refine ⟨rfl, ?_⟩
    refine ⟨1, ?_⟩
    ring

/-- Cl(1,3)-anchored legality statement:
with the timelike sign fixed by the algebra, periodic timelike identification
admits closed timelike curves. -/
theorem ctc_legal_if_periodic_timelike_identification (T : ℝ) (hT : T > 0) :
    minkowskiQF (e 0) = -1 ∧
      ∃ a b : Event, Timelike a b ∧ sameOnTimeCylinder T a b := by
  refine ⟨cl13_timelike_generator_negative, ?_⟩
  exact ctc_exists_on_time_cylinder T hT

end Gutoe.CTCLegality
