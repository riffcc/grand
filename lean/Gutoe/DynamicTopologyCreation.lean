import Mathlib
import Gutoe.CreationLanes
import Gutoe.VacuumEnergyBounds

/-!
GUTOE — Dynamic Topology Creation Lane (Path B scaffold)

This lane formalizes the *creation gate* skeleton: a compact local identification
is considered dynamically creatable only if an explicit budget exceeds a
structural threshold tied to Cl(1,3) constants.

This does not prove physical realizability; it keeps the decision surface explicit.
-/

namespace Gutoe.DynamicTopologyCreation

open Gutoe
open Gutoe.CreationLanes
open Gutoe.VacuumEnergyBounds

/-- Structural threshold proxy from Cl(1,3) lane constants.
`f_void = 3/16` controls the scale; radius and period provide geometric size. -/
def structuralCreationThreshold (radius period : ℝ) : ℝ :=
  (voidFractionQ : ℝ) * |radius| * |period|

/-- Dynamic creation gate (Path B): finite local patch creation requires
positive period and enough budget against the structural threshold. -/
def dynamicCreationGate (budget radius period : ℝ) : Prop :=
  0 < period ∧ budget ≥ structuralCreationThreshold radius period

/-- The structural threshold is nonnegative. -/
theorem structural_creation_threshold_nonneg (radius period : ℝ) :
    0 ≤ structuralCreationThreshold radius period := by
  unfold structuralCreationThreshold
  have hvoid : (0 : ℝ) ≤ (voidFractionQ : ℝ) := by
    norm_num [void_fraction_eq_3_16]
  exact mul_nonneg (mul_nonneg hvoid (abs_nonneg radius)) (abs_nonneg period)

/-- If the dynamic gate passes, a nontrivial local identified shift exists
inside the declared support radius. -/
theorem dynamic_gate_implies_local_shift
    (budget radius period x : ℝ)
    (hGate : dynamicCreationGate budget radius period)
    (hx : |x| ≤ radius) :
    ∃ a b : Gutoe.CTCLegality.Event,
      sameOnLocalPatch period radius a b ∧ b.t ≠ a.t := by
  exact local_patch_nontrivial_shift_exists period radius x hGate.1 hx

/-- Failing the budget side blocks the dynamic creation gate. -/
theorem insufficient_budget_blocks_dynamic_gate
    (budget radius period : ℝ)
    (hBudget : budget < structuralCreationThreshold radius period) :
    ¬ dynamicCreationGate budget radius period := by
  intro hGate
  exact not_le_of_gt hBudget hGate.2

end Gutoe.DynamicTopologyCreation
