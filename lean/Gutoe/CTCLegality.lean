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

/-- A fixed-space step across `n` periods is timelike for `n > 0`. -/
theorem fixed_space_n_steps_timelike (x0 t0 T : ℝ) (n : ℕ) (hn : 0 < n) (hT : T > 0) :
    Timelike ⟨t0, x0⟩ ⟨t0 + (n : ℝ) * T, x0⟩ := by
  unfold Timelike intervalSq
  have ht : (t0 + (n : ℝ) * T - t0) ^ 2 = ((n : ℝ) * T) ^ 2 := by ring
  have hx : (x0 - x0) ^ 2 = 0 := by ring
  rw [ht, hx]
  have hnR : (n : ℝ) > 0 := by exact_mod_cast hn
  have hmul : (n : ℝ) * T > 0 := mul_pos hnR hT
  have hsq : (((n : ℝ) * T) ^ 2) > 0 := sq_pos_of_pos hmul
  nlinarith

/-- Exact interval on a fixed-space step of size `T`:
`ds² = -T²`. -/
theorem fixed_space_step_interval_exact (x0 t0 T : ℝ) :
    intervalSq ⟨t0, x0⟩ ⟨t0 + T, x0⟩ = -(T ^ 2) := by
  unfold intervalSq
  ring

/-- Exact interval on a fixed-space `n`-period step:
`ds² = -((n*T)²)`. -/
theorem fixed_space_n_steps_interval_exact (x0 t0 T : ℝ) (n : ℕ) :
    intervalSq ⟨t0, x0⟩ ⟨t0 + (n : ℝ) * T, x0⟩ = -(((n : ℝ) * T) ^ 2) := by
  unfold intervalSq
  ring

/-- On the time-cylinder, `n` periods returns to the same identified event class. -/
theorem same_on_time_cylinder_after_n_periods (x0 t0 T : ℝ) (n : ℕ) :
    sameOnTimeCylinder T ⟨t0, x0⟩ ⟨t0 + (n : ℝ) * T, x0⟩ := by
  unfold sameOnTimeCylinder
  refine ⟨rfl, ?_⟩
  refine ⟨(n : ℤ), ?_⟩
  change t0 + (n : ℝ) * T = t0 + ((n : ℤ) : ℝ) * T
  simp

/-- Escher-stair form: each lap is locally timelike and finite, while global
covering-time after `n` laps is `n*T`. -/
theorem escher_stair_n_laps (T : ℝ) (hT : T > 0) :
    ∀ n : ℕ, 0 < n →
      ∃ a b : Event,
        Timelike a b ∧ sameOnTimeCylinder T a b ∧ b.t - a.t = (n : ℝ) * T := by
  intro n hn
  refine ⟨⟨0, 0⟩, ⟨(n : ℝ) * T, 0⟩, ?_, ?_, ?_⟩
  · simpa [zero_add] using fixed_space_n_steps_timelike 0 0 T n hn hT
  · simpa [zero_add] using same_on_time_cylinder_after_n_periods 0 0 T n
  · ring

/-- Unbounded covering-time growth on the identified timelike axis:
for any target `M`, enough laps exceed it while every lap remains local-timelike. -/
theorem escher_stair_cover_time_unbounded (T : ℝ) (hT : T > 0) :
    ∀ M : ℝ, ∃ n : ℕ, M < (n : ℝ) * T := by
  intro M
  rcases exists_nat_gt (M / T) with ⟨n, hn⟩
  exact ⟨n, (div_lt_iff₀ hT).1 hn⟩

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

/-- "Escher-stair" loop statement:
there exists a locally timelike step with strictly positive local coordinate
time that closes only via periodic identification. -/
theorem ctc_step_forward_but_identified (T : ℝ) (hT : T > 0) :
    ∃ a b : Event, Timelike a b ∧ sameOnTimeCylinder T a b ∧ b.t > a.t := by
  refine ⟨⟨0, 0⟩, ⟨T, 0⟩, ?_, ?_, ?_⟩
  · simpa [zero_add] using fixed_space_step_timelike 0 0 T hT
  · unfold sameOnTimeCylinder
    refine ⟨rfl, ?_⟩
    refine ⟨1, ?_⟩
    ring
  · simpa using hT

/-- There exists a nontrivial (`Δt ≠ 0`) identified time-shift at fixed `x`. -/
def nontrivialTimeShiftAtX (T x : ℝ) : Prop :=
  ∃ a b : Event, a.x = x ∧ b.x = x ∧ sameOnTimeCylinder T a b ∧ b.t ≠ a.t

/-- In the current time-cylinder model, nontrivial identifications exist at
every spatial `x` (global boundary-condition style). -/
theorem nontrivial_shift_at_every_x (T : ℝ) (hT : T > 0) :
    ∀ x : ℝ, nontrivialTimeShiftAtX T x := by
  intro x
  refine ⟨⟨0, x⟩, ⟨T, x⟩, rfl, rfl, ?_, ?_⟩
  · unfold sameOnTimeCylinder
    refine ⟨rfl, ?_⟩
    refine ⟨1, ?_⟩
    ring
  · have hne : T ≠ 0 := by linarith
    simpa using hne

/-- The current time-cylinder identification cannot be compactly supported in
space: there is no finite `R` beyond which nontrivial identifications vanish. -/
theorem sameOnTimeCylinder_not_compactly_supported (T : ℝ) (hT : T > 0) :
    ¬ ∃ R : ℝ, 0 ≤ R ∧ ∀ x : ℝ, |x| > R → ¬ nontrivialTimeShiftAtX T x := by
  intro h
  rcases h with ⟨R, hR, houtside⟩
  let x0 : ℝ := R + 1
  have hx0_nonneg : 0 ≤ x0 := by
    dsimp [x0]
    linarith
  have hx0_abs_gt : |x0| > R := by
    rw [abs_of_nonneg hx0_nonneg]
    dsimp [x0]
    linarith
  have hnone : ¬ nontrivialTimeShiftAtX T x0 := houtside x0 hx0_abs_gt
  have hsome : nontrivialTimeShiftAtX T x0 := (nontrivial_shift_at_every_x T hT) x0
  exact hnone hsome

/-- Classification theorem for the current CTC legality model:
it formalizes a global boundary condition (Path A shape), not a localized
"create-the-staircase-here" compact-support mechanism. -/
theorem current_model_global_boundary_condition (T : ℝ) (hT : T > 0) :
    (∀ x : ℝ, nontrivialTimeShiftAtX T x) ∧
    (¬ ∃ R : ℝ, 0 ≤ R ∧ ∀ x : ℝ, |x| > R → ¬ nontrivialTimeShiftAtX T x) := by
  exact ⟨nontrivial_shift_at_every_x T hT, sameOnTimeCylinder_not_compactly_supported T hT⟩

end Gutoe.CTCLegality
