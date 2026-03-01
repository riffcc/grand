import Mathlib
import Gutoe.CTCDoorEnergyBookkeeping

/-!
GUTOE — Closed-loop temporal-loan lane

This module resolves a common confusion:

- "Closed loop forbids positive export" is true only under *no drawdown* and
  nonnegative loss (already proven in `CTCDoorEnergyBookkeeping`).
- A closed loop can still exhibit a **transient positive local window** if one
  phase draws down door-state energy and a later phase repays it.

So conservation can hold globally while local windows show positive export.
-/

namespace Gutoe.CTCClosedLoopTemporalLoan

open Gutoe.CTCDoorEnergyBookkeeping

/-- In a closed packet cycle, export equals drawdown minus loss. -/
theorem closed_cycle_export_eq_drawdown_minus_loss
    (Ein Eprev Eout Enext Export Loss : ℝ)
    (hCon : LoopConservation Ein Eprev Eout Enext Export Loss)
    (hClosed : ClosedPacketCycle Ein Eout) :
    Export = (Eprev - Enext) - Loss := by
  have hExp := export_equals_inflow_plus_drawdown_minus_loss
    Ein Eprev Eout Enext Export Loss hCon
  unfold ClosedPacketCycle at hClosed
  rw [hClosed] at hExp
  linarith

/-- In a closed packet cycle, strictly positive export is equivalent to
drawdown exceeding losses. -/
theorem closed_cycle_positive_export_iff_drawdown_gt_loss
    (Ein Eprev Eout Enext Export Loss : ℝ)
    (hCon : LoopConservation Ein Eprev Eout Enext Export Loss)
    (hClosed : ClosedPacketCycle Ein Eout) :
    (0 < Export) ↔ (Loss < Eprev - Enext) := by
  have hEq : Export = (Eprev - Enext) - Loss :=
    closed_cycle_export_eq_drawdown_minus_loss
      Ein Eprev Eout Enext Export Loss hCon hClosed
  constructor <;> intro h
  · rw [hEq] at h
    linarith
  · rw [hEq]
    linarith

/-- A closed cycle can have a transient positive-export phase (loan draw) and a
later negative-export phase (repayment), with zero net export and restored door
state at cycle end. -/
theorem two_phase_temporal_loan_exists :
    ∃ E0 E1 ExportA ExportB : ℝ,
      -- Phase A (loan draw)
      LoopConservation 0 E0 0 E1 ExportA 0 ∧
      ClosedPacketCycle 0 0 ∧
      ExportA > 0 ∧
      E1 < E0 ∧
      -- Phase B (repay)
      LoopConservation 0 E1 0 E0 ExportB 0 ∧
      ClosedPacketCycle 0 0 ∧
      ExportB < 0 ∧
      -- End-to-end closure
      ExportA + ExportB = 0 ∧
      E0 = E0 := by
  refine ⟨10, 9, 1, -1, ?_⟩
  refine ⟨?_, ?_⟩
  · unfold LoopConservation
    norm_num
  · refine ⟨?_, ?_⟩
    · rfl
    · refine ⟨by norm_num, ?_⟩
      refine ⟨by norm_num, ?_⟩
      refine ⟨?_, ?_⟩
      · unfold LoopConservation
        norm_num
      · refine ⟨?_, ?_⟩
        · rfl
        · refine ⟨by norm_num, ?_⟩
          refine ⟨by norm_num, ?_⟩
          · rfl

/-- Persistent per-step positive export is impossible in a strict closed cycle
with no drawdown and nonnegative loss. (Restates the hard guard explicitly.) -/
theorem persistent_positive_export_blocked_under_nodraw
    (Ein Eprev Eout Enext Export Loss : ℝ)
    (hCon : LoopConservation Ein Eprev Eout Enext Export Loss)
    (hClosed : ClosedPacketCycle Ein Eout)
    (hNoDraw : NoDoorDrawdown Eprev Enext)
    (hLoss : 0 ≤ Loss) :
    ¬ (Export > 0) := by
  intro hPos
  have hLe : Export ≤ 0 := closed_cycle_no_positive_export
    Ein Eprev Eout Enext Export Loss hCon hClosed hNoDraw hLoss
  exact (not_lt_of_ge hLe) hPos

end Gutoe.CTCClosedLoopTemporalLoan
