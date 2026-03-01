import Mathlib
import Gutoe.DynamicTopologyCreation
import Gutoe.CTCPathAEffectiveArrival

/-!
GUTOE — Path-A loop gain closure from creation budget

This module removes the free `q` knob in the Path-A lane by defining loop gain
as the ratio of explicit creation budget to structural creation threshold.

`q_eff = budget / ((3/16)|R||T|)`

So:
- `q_eff = 1` is break-even,
- `q_eff > 1` is over-threshold.
-/

namespace Gutoe.CTCPathAQClosure

open Gutoe.DynamicTopologyCreation

noncomputable section

/-- Loop-gain closure: gain is budget divided by structural threshold. -/
def loopGainFromBudget (budget radius period : ℝ) : ℝ :=
  budget / structuralCreationThreshold radius period

/-- Real-valued effective arrival using budget-closed loop gain. -/
def effectiveArrivalFromBudget
    (dtIn dtOut T : ℝ) (n : ℕ)
    (budget radius period : ℝ) : ℝ :=
  dtIn + dtOut + (n : ℝ) * (1 - loopGainFromBudget budget radius period) * T

/-- Positive radius and positive period imply a positive structural threshold. -/
theorem structural_threshold_pos_of_pos
    (radius period : ℝ)
    (hR : 0 < radius) (hP : 0 < period) :
    0 < structuralCreationThreshold radius period := by
  unfold structuralCreationThreshold
  have hvoid : (0 : ℝ) < (Gutoe.VacuumEnergyBounds.voidFractionQ : ℝ) := by
    norm_num [Gutoe.VacuumEnergyBounds.void_fraction_eq_3_16]
  have hRabs : 0 < |radius| := by simpa [abs_of_pos hR] using hR
  have hPabs : 0 < |period| := by simpa [abs_of_pos hP] using hP
  exact mul_pos (mul_pos hvoid hRabs) hPabs

/-- At exact threshold budget, loop gain is exactly one. -/
theorem loop_gain_eq_one_of_budget_eq_threshold
    (budget radius period : ℝ)
    (hEq : budget = structuralCreationThreshold radius period)
    (hThr : structuralCreationThreshold radius period ≠ 0) :
    loopGainFromBudget budget radius period = 1 := by
  unfold loopGainFromBudget
  rw [hEq]
  field_simp [hThr]

/-- Over-threshold budget implies loop gain strictly above one. -/
theorem loop_gain_gt_one_of_budget_gt_threshold
    (budget radius period : ℝ)
    (hGt : structuralCreationThreshold radius period < budget)
    (hThrPos : 0 < structuralCreationThreshold radius period) :
    1 < loopGainFromBudget budget radius period := by
  let thr := structuralCreationThreshold radius period
  have hCore : 1 < budget / thr := by
    have hEq : (1 : ℝ) < budget / thr ↔ thr < budget := by
      exact one_lt_div (by simpa [thr] using hThrPos)
    exact hEq.2 (by simpa [thr] using hGt)
  simpa [loopGainFromBudget, thr] using hCore

/-- Dynamic creation gate plus positive radius implies non-lossy loop gain (`q_eff ≥ 1`). -/
theorem dynamic_gate_implies_loop_gain_ge_one
    (budget radius period : ℝ)
    (hGate : dynamicCreationGate budget radius period)
    (hR : 0 < radius) :
    1 ≤ loopGainFromBudget budget radius period := by
  let thr := structuralCreationThreshold radius period
  have hThrPos : 0 < structuralCreationThreshold radius period :=
    structural_threshold_pos_of_pos radius period hR hGate.1
  have hCore : 1 ≤ budget / thr := by
    have hEq : (1 : ℝ) ≤ budget / thr ↔ thr ≤ budget := by
      exact one_le_div (by simpa [thr] using hThrPos)
    exact hEq.2 (by simpa [thr] using hGate.2)
  simpa [loopGainFromBudget, thr] using hCore

/-- Break-even budget implies no loop-based coordinate gain in the Path-A formula. -/
theorem threshold_budget_no_coordinate_gain
    (dtIn dtOut T : ℝ) (n : ℕ)
    (radius period : ℝ)
    (hThr : structuralCreationThreshold radius period ≠ 0) :
    effectiveArrivalFromBudget dtIn dtOut T n
      (structuralCreationThreshold radius period) radius period = dtIn + dtOut := by
  unfold effectiveArrivalFromBudget
  have hq : loopGainFromBudget (structuralCreationThreshold radius period) radius period = 1 := by
    apply loop_gain_eq_one_of_budget_eq_threshold
    · rfl
    · exact hThr
  rw [hq]
  ring

/-- If budget is strictly above the structural threshold, then the closed-gain
`q_eff` is strictly above one; with positive loop period this guarantees
pre-departure effective coordinate arrival for sufficiently many loops. -/
theorem over_threshold_budget_predeparture_possible
    (dtIn dtOut T : ℝ)
    (budget radius period : ℝ)
    (hT : 0 < T)
    (hR : 0 < radius)
    (hP : 0 < period)
    (hBudget : structuralCreationThreshold radius period < budget) :
    ∃ n : ℕ, effectiveArrivalFromBudget dtIn dtOut T n budget radius period < 0 := by
  have hThrPos : 0 < structuralCreationThreshold radius period :=
    structural_threshold_pos_of_pos radius period hR hP
  have hq_gt_one : 1 < loopGainFromBudget budget radius period :=
    loop_gain_gt_one_of_budget_gt_threshold budget radius period hBudget hThrPos
  let k : ℝ := (loopGainFromBudget budget radius period - 1) * T
  have hk : 0 < k := by
    have hkg : 0 < loopGainFromBudget budget radius period - 1 := by linarith
    exact mul_pos hkg hT
  rcases exists_nat_gt ((dtIn + dtOut) / k) with ⟨n, hn⟩
  refine ⟨n, ?_⟩
  have hnk : dtIn + dtOut < (n : ℝ) * k := by
    exact (div_lt_iff₀ hk).1 hn
  calc
    effectiveArrivalFromBudget dtIn dtOut T n budget radius period
        = dtIn + dtOut + (n : ℝ) *
            (1 - loopGainFromBudget budget radius period) * T := by
            rfl
    _ = dtIn + dtOut - (n : ℝ) * k := by
          dsimp [k]
          ring
    _ < 0 := by linarith [hnk]

end

end Gutoe.CTCPathAQClosure
