import Mathlib
import Gutoe.FTLGeometry
import Gutoe.VacuumEnergyBounds

/-!
GUTOE — Finite Mechanism Exhaustion Lane for FTL Claims

This module does not claim a universal theorem over all conceivable physics.
It proves an exhaustion result over an explicit, finite mechanism universe used
by the current FTL board lane.
-/

namespace Gutoe.FTLExhaustion

/-- Enumerated mechanism universe for current FTL evaluation. -/
inductive FTLMechanism where
  | warpMetric
  | traversableWormhole
  | casimirNegativeEnergyDrive
  | higgsGradientWallSurf
  | rearFaceWallSurf
  | tachyonicSignal
  | entanglementMessaging
deriving DecidableEq, Repr

/-- Requires local superluminal signalling in the spacetime projection. -/
def requiresLocalSuperluminalSignal : FTLMechanism → Prop
  | .tachyonicSignal => True
  | _ => False

/-- Requires macroscopic negative-energy support. -/
def requiresMacroscopicNegativeEnergy : FTLMechanism → Prop
  | .warpMetric => True
  | .traversableWormhole => True
  | .casimirNegativeEnergyDrive => True
  | _ => False

/-- Requires controllable macroscopic wall-surf actuation in the Higgs/void lane. -/
def requiresWallSurfActuator : FTLMechanism → Prop
  | .higgsGradientWallSurf => True
  | .rearFaceWallSurf => True
  | _ => False

/-- Mechanisms that are non-traversable/non-signalling by construction. -/
def nonTraversableOrNonSignallingAsStated : FTLMechanism → Prop
  | .entanglementMessaging => True
  | _ => False

/-- Feasibility gate under the current finite mechanism lane. -/
def feasibleUnderCurrentLane (m : FTLMechanism) : Prop :=
  ¬ requiresLocalSuperluminalSignal m ∧
  ¬ requiresMacroscopicNegativeEnergy m ∧
  ¬ requiresWallSurfActuator m ∧
  ¬ nonTraversableOrNonSignallingAsStated m

/-- Every mechanism in the current finite universe triggers at least one blocker class. -/
theorem every_mechanism_has_a_blocker_class (m : FTLMechanism) :
    requiresLocalSuperluminalSignal m ∨
      requiresMacroscopicNegativeEnergy m ∨
      requiresWallSurfActuator m ∨
      nonTraversableOrNonSignallingAsStated m := by
  cases m <;> simp [requiresLocalSuperluminalSignal, requiresMacroscopicNegativeEnergy,
    requiresWallSurfActuator, nonTraversableOrNonSignallingAsStated]

/-- Exhaustion theorem over the declared finite mechanism universe:
every listed mechanism hits at least one blocker class, so no entry passes
`feasibleUnderCurrentLane`. -/
theorem exhaustive_no_go_declared_universe :
    ∀ m, ¬ feasibleUnderCurrentLane m := by
  intro m hFeasible
  rcases every_mechanism_has_a_blocker_class m with hL | hN | hW | hNT
  · exact (hFeasible.1) hL
  · exact (hFeasible.2.1) hN
  · exact (hFeasible.2.2.1) hW
  · exact (hFeasible.2.2.2) hNT

end Gutoe.FTLExhaustion
