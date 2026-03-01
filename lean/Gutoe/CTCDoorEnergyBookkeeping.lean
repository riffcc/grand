import Mathlib
import Gutoe.DynamicTopologyCreation

/-!
GUTOE — CTC door energy bookkeeping

This module formalizes the per-loop conservation identity for a dynamic
"door-state" model and proves the key contradiction condition:

Positive exported work cannot occur in a closed packet-energy cycle unless
either:
1) there is net packet-energy inflow (`Ein > Eout`), or
2) the door energy is drawn down (`Enext < Eprev`).
-/

namespace Gutoe.CTCDoorEnergyBookkeeping

open Gutoe.DynamicTopologyCreation

/-- Per-loop conservation identity (all terms in Joules):
`Ein + Eprev = Eout + Enext + Export + Loss`.
-/
def LoopConservation
    (Ein Eprev Eout Enext Export Loss : ℝ) : Prop :=
  Ein + Eprev = Eout + Enext + Export + Loss

/-- Closed packet-energy cycle: no net packet-energy transfer across the loop. -/
def ClosedPacketCycle (Ein Eout : ℝ) : Prop := Ein = Eout

/-- Non-drawdown door condition: next-cycle door energy is at least prior cycle. -/
def NoDoorDrawdown (Eprev Enext : ℝ) : Prop := Eprev ≤ Enext

/-- In a closed packet cycle with nonnegative losses and no door drawdown,
exported work is nonpositive. -/
theorem closed_cycle_no_positive_export
    (Ein Eprev Eout Enext Export Loss : ℝ)
    (hCon : LoopConservation Ein Eprev Eout Enext Export Loss)
    (hClosed : ClosedPacketCycle Ein Eout)
    (hNoDraw : NoDoorDrawdown Eprev Enext)
    (hLoss : 0 ≤ Loss) :
    Export ≤ 0 := by
  unfold LoopConservation ClosedPacketCycle NoDoorDrawdown at *
  linarith

/-- Rearranged identity for per-loop export. -/
theorem export_equals_inflow_plus_drawdown_minus_loss
    (Ein Eprev Eout Enext Export Loss : ℝ)
    (hCon : LoopConservation Ein Eprev Eout Enext Export Loss) :
    Export = (Ein - Eout) + (Eprev - Enext) - Loss := by
  unfold LoopConservation at hCon
  linarith

/-- If export is strictly positive with nonnegative losses, then either there is
net packet-energy inflow or the door is being drawn down. -/
theorem positive_export_requires_inflow_or_drawdown
    (Ein Eprev Eout Enext Export Loss : ℝ)
    (hCon : LoopConservation Ein Eprev Eout Enext Export Loss)
    (hLoss : 0 ≤ Loss)
    (hExp : 0 < Export) :
    (Ein > Eout) ∨ (Enext < Eprev) := by
  by_cases hIn : Ein > Eout
  · exact Or.inl hIn
  · right
    have hInLe : Ein ≤ Eout := le_of_not_gt hIn
    unfold LoopConservation at hCon
    linarith

/-- If a door starts above structural threshold and is not drawn down, it stays
above threshold after one loop. -/
theorem threshold_preserved_under_nondrawdown
    (radius period Eprev Enext : ℝ)
    (hPrev : structuralCreationThreshold radius period ≤ Eprev)
    (hNoDraw : Eprev ≤ Enext) :
    structuralCreationThreshold radius period ≤ Enext := by
  exact le_trans hPrev hNoDraw

end Gutoe.CTCDoorEnergyBookkeeping
