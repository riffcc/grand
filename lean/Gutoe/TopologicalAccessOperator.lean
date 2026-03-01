import Mathlib
import Gutoe.TopologicalAccessControllability
import Gutoe.CTCDoorEnergyBookkeeping
import Gutoe.TopologicalDefectBundle

/-!
GUTOE — Topological Access Operator (entry/exit + conservation)

This module composes three previously separated lanes:

1. Entry witness: local nontrivial identification shift (Path-B gate).
2. Exit witness: controllable descended target via affine/topological offset.
3. Bookkeeping: closed-cycle conservation forbids positive net export.

It does not claim full-state reconstruction from the 4D projection.
-/

namespace Gutoe.TopologicalAccessOperator

open Gutoe
open Gutoe.CreationLanes
open Gutoe.DynamicTopologyCreation
open Gutoe.RecursiveNavigationNoTranslation
open Gutoe.RecursiveZ3Tower
open Gutoe.TopologicalAccessControllability
open Gutoe.CTCDoorEnergyBookkeeping
open Gutoe.TopologicalDefectBundle

noncomputable section

/-- Combined controllable operator witness:
local entry shift + chosen descended target exit. -/
theorem controllable_entry_exit_if_gate
    (budget radius period x_local : ℝ)
    (hGate : dynamicCreationGate budget radius period)
    (hxLocal : |x_local| ≤ radius)
    (target : Fin 4 → ℝ) :
    ∃ a b : Gutoe.CTCLegality.Event, ∃ t : Vec256,
      sameOnLocalPatch period radius a b ∧
      b.t ≠ a.t ∧
      towerProjection (affineStep (0 : Vec256 →ₗ[ℝ] Vec256) t 0) = target := by
  rcases dynamic_gate_with_affine_targetability budget radius period x_local hGate hxLocal target with
    ⟨hEntry, hExit⟩
  rcases hEntry with ⟨a, b, hab, hneq⟩
  rcases hExit with ⟨t, ht⟩
  exact ⟨a, b, t, hab, hneq, ht⟩

/-- Linear lane remains blocked for nonzero targets, even in operator context. -/
theorem linear_lane_still_blocks_nonzero
    (x : Fin 4 → ℝ) (hx : x ≠ 0) :
    ∀ L : Vec256 →ₗ[ℝ] Vec256, towerProjection (L 0) ≠ x := by
  intro L
  exact no_linear_origin_to_nonzero_target L x hx

/-- Explicit bypass statement:
for nonzero descended targets, linear-only fails but topological/affine route
has a controllable witness when the dynamic gate is open. -/
theorem topological_operator_bypasses_linear_coordinate_nogo
    (budget radius period x_local : ℝ)
    (hGate : dynamicCreationGate budget radius period)
    (hxLocal : |x_local| ≤ radius)
    (x : Fin 4 → ℝ) (hx : x ≠ 0) :
    (∀ L : Vec256 →ₗ[ℝ] Vec256, towerProjection (L 0) ≠ x) ∧
    (∃ a b : Gutoe.CTCLegality.Event, ∃ t : Vec256,
      sameOnLocalPatch period radius a b ∧
      b.t ≠ a.t ∧
      towerProjection (affineStep (0 : Vec256 →ₗ[ℝ] Vec256) t 0) = x) := by
  refine ⟨linear_lane_still_blocks_nonzero x hx, ?_⟩
  exact controllable_entry_exit_if_gate budget radius period x_local hGate hxLocal x

/-- Conservation-compliant operator use in a closed packet cycle:
no positive net export is allowed under non-drawdown and nonnegative loss. -/
theorem operator_closed_cycle_export_nonpositive
    (Ein Eprev Eout Enext Export Loss : ℝ)
    (hCon : LoopConservation Ein Eprev Eout Enext Export Loss)
    (hClosed : ClosedPacketCycle Ein Eout)
    (hNoDraw : NoDoorDrawdown Eprev Enext)
    (hLoss : 0 ≤ Loss) :
    Export ≤ 0 := by
  exact closed_cycle_no_positive_export Ein Eprev Eout Enext Export Loss
    hCon hClosed hNoDraw hLoss

/-- 1D quotient-bridge kinematic shortcut witness (endpoint identification). -/
theorem quotient_bridge_endpoint_shortcut
    (l r : ℝ) (hneq : l ≠ r) :
    defectDistance l r l r < baseDistance l r := by
  exact bridge_endpoints_strict_shortcut l r hneq

end
end Gutoe.TopologicalAccessOperator
