/-
 * GUTOE - Void Rear Face Correspondence
 * Copyright (C) 2026  Riff Labs
 *
 * AGPL-3.0-or-later
 *
 * Formal lane for the dual-layer "front face / rear face" experiments:
 * each front state has a unique rear counterpart and the two faces
 * have matched Clifford-traced cardinality.
 -/

import Mathlib
import Gutoe.HexFermions
import Gutoe.Z3Uniqueness

namespace Gutoe.VoidRearFace

open HexState
open Gutoe.HexFermions
open Gutoe.Z3Uniqueness

/-!
`FrontState` and `RearState` are the two 6-state faces in the dual-layer hex
model.  Rear-face correspondence is the `negate` involution from `HexState`.
-/

/-- States on the front (positive) face. -/
def FrontState := { s : HexState // s.isPos = true }

/-- States on the rear (negative/void-complement) face. -/
def RearState := { s : HexState // s.isNeg = true }

/-- Rear-face correspondence map: front state ↦ its negated rear counterpart. -/
def toRear (s : FrontState) : RearState :=
  ⟨s.1.negate, negate_pos_to_neg s.1 s.2⟩

/-- Inverse correspondence: rear state ↦ its negated front counterpart. -/
def toFront (s : RearState) : FrontState :=
  ⟨s.1.negate, negate_neg_to_pos s.1 s.2⟩

/-- `toFront` is a left inverse of `toRear`. -/
theorem toFront_leftInverse_toRear : Function.LeftInverse toFront toRear := by
  intro s
  apply Subtype.ext
  simp [toFront, toRear, HexState.negate_involutive]

/-- `toFront` is a right inverse of `toRear`. -/
theorem toFront_rightInverse_toRear : Function.RightInverse toFront toRear := by
  intro s
  apply Subtype.ext
  simp [toFront, toRear, HexState.negate_involutive]

/-- Rear-face correspondence is bijective: every front state has one unique rear
counterpart and every rear state comes from exactly one front state. -/
theorem rear_face_correspondence_bijective : Function.Bijective toRear := by
  refine ⟨toFront_leftInverse_toRear.injective, ?_⟩
  exact toFront_rightInverse_toRear.surjective

/-- Shared Cl(1,3) structural count: the grade-2 index set has cardinality 6. -/
theorem grade2_4d_card_eq_six : grade2_4d.card = 6 := by
  native_decide

/-- Each face has cardinality equal to the shared grade-2 structural count. -/
theorem face_cardinality_matches_grade2 :
    posFace.length = grade2_4d.card ∧ negFace.length = grade2_4d.card := by
  constructor <;> simpa [grade2_4d_card_eq_six, posFace, negFace]

/-- Total dual-layer state count is exactly `2 * grade2_4d.card = 12`. -/
theorem dual_layer_count_eq_two_times_grade2 :
    HexState.all.length = 2 * grade2_4d.card := by
  simpa [grade2_4d_card_eq_six, HexState.all]

/-!
Experiment-aligned transfer-cost lane:
- normal/front hop weight = 1
- rear/void-channel hop weight = 1/10
For equal 3-hop routes, rear-channel route is strictly shorter.
-/

/-- Front-face (normal) edge weight from the transfer experiment lane. -/
def normalHopCost : ℚ := 1

/-- Rear-face (void-channel) edge weight from the transfer experiment lane. -/
def rearHopCost : ℚ := (1 : ℚ) / 10

/-- A representative 3-hop front-face route cost. -/
def frontRouteCost : ℚ := 3 * normalHopCost

/-- A representative 3-hop rear-face route cost. -/
def rearRouteCost : ℚ := 3 * rearHopCost

/-- Under shared hop-count, rear-face channel routing is strictly cheaper. -/
theorem rear_route_strictly_shorter : rearRouteCost < frontRouteCost := by
  norm_num [rearRouteCost, frontRouteCost, rearHopCost, normalHopCost]

/-- Rear/front route-cost factor under the experiment lane. -/
def rearCostFactor : ℚ := rearHopCost / normalHopCost

/-- Rear cost factor is exactly one tenth. -/
theorem rear_cost_factor_eq_one_tenth : rearCostFactor = (1 : ℚ) / 10 := by
  norm_num [rearCostFactor, rearHopCost, normalHopCost]

/-- Linear energy scaling law: rear-channel budget is one tenth of front budget. -/
theorem rear_energy_linear_scaling (Efront : ℚ) :
    Efront * rearCostFactor = Efront / 10 := by
  rw [rear_cost_factor_eq_one_tenth]
  ring

/-- Under positive front-channel budget, rear-channel budget is strictly smaller. -/
theorem rear_energy_strictly_smaller_of_positive (Efront : ℚ) (hE : Efront > 0) :
    Efront * rearCostFactor < Efront := by
  rw [rear_energy_linear_scaling]
  nlinarith

/-- Concrete wall translation in the linear lane:
front ratio `10^32` maps to rear ratio `10^31`. -/
theorem wall_ratio_reduces_one_order_linear :
    ((10 : ℚ) ^ 32) * rearCostFactor = (10 : ℚ) ^ 31 := by
  rw [rear_cost_factor_eq_one_tenth]
  have hpow : (10 : ℚ) ^ 32 = (10 : ℚ) ^ 31 * (10 : ℚ) := by
    simpa [pow_succ]
  rw [hpow]
  ring

end Gutoe.VoidRearFace
