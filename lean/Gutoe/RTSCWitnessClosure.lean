import Mathlib
import Gutoe.RTSCAdmissibility
import Gutoe.LatticeGeometry
import Gutoe.GaugeGroupSU3

namespace Gutoe.RTSCWitnessClosure

open Gutoe.RTSCAdmissibility
open Gutoe.LatticeGeometry
open Gutoe.GaugeGroupSU3

/-- Runtime-forced finite witness set from the RTSC Rust gate lane.
    This module closes structural arithmetic invariants on that set. -/
def rtscWitnessZ : Finset ℕ := {24, 30, 42, 48, 72, 78}

theorem rtsc_witness_card :
    rtscWitnessZ.card = 6 := by
  decide

theorem rtsc_witness_card_eq_coordination :
    rtscWitnessZ.card = coordinationNumber := by
  rw [rtsc_witness_card, coordination_number_is_6]

theorem rtsc_witness_triplet_residue_zero :
    ∀ z ∈ rtscWitnessZ, z % quarkOrbit.card = 0 := by
  intro z hz
  rw [quarkOrbit_card]
  fin_cases hz <;> decide

theorem rtsc_witness_pair_step_coordination :
    (30 - 24 = coordinationNumber) ∧
    (48 - 42 = coordinationNumber) ∧
    (78 - 72 = coordinationNumber) := by
  rw [coordination_number_is_6]
  norm_num

theorem rtsc_witness_lane_closed :
    rtscAdmissible ∧
    rtscWitnessZ.card = coordinationNumber ∧
    (∀ z ∈ rtscWitnessZ, z % quarkOrbit.card = 0) := by
  refine ⟨rtsc_gate_admissible, rtsc_witness_card_eq_coordination, rtsc_witness_triplet_residue_zero⟩

end Gutoe.RTSCWitnessClosure

