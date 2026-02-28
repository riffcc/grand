import Mathlib
import Gutoe.LatticeGeometry
import Gutoe.GaugeGroupSU3
import Gutoe.DarkMatterSector
import Gutoe.DimensionalStructure
import Gutoe.Z3Uniqueness

namespace Gutoe.RTSCAdmissibility

open Gutoe.LatticeGeometry
open Gutoe.GaugeGroupSU3
open Gutoe.DarkMatterSector
open Gutoe.DimensionalStructure
open Gutoe.Z3Uniqueness

/-- Forced lattice family candidates for the RTSC gate. -/
inductive LatticeFamily where
  | simpleCubic
  | rejected
deriving DecidableEq, Repr

/-- Forced lattice family from the Cl(1,3) spatial-bivector coordination. -/
def forcedLatticeFamily : LatticeFamily :=
  if coordinationNumber = 6 then LatticeFamily.simpleCubic else LatticeFamily.rejected

theorem forced_lattice_family_simple_cubic :
    forcedLatticeFamily = LatticeFamily.simpleCubic := by
  unfold forcedLatticeFamily
  simp [coordination_number_is_6]

/-- Forced filling multiplicity from the Z₃ quark orbit size. -/
def forcedFillingMultiplicity : ℕ := quarkOrbit.card

theorem forced_filling_triplet :
    forcedFillingMultiplicity = 3 := by
  unfold forcedFillingMultiplicity
  simpa using quarkOrbit_card

/-- Repulsive kernel term from half the inverse grade-2 multiplicity. -/
def pairingRepulsionQ : ℚ := 1 / (2 * (grade2_4d.card : ℚ))

theorem pairing_repulsion_eq_1_over_12 :
    pairingRepulsionQ = 1 / 12 := by
  unfold pairingRepulsionQ
  have h6 : grade2_4d.card = 6 := by decide
  rw [h6]
  norm_num

/-- Net pairing kernel for the RTSC gate:
    dark occupancy fraction minus short-range repulsion term. -/
def pairingKernelQ : ℚ := darkFractionOfTotalStates - pairingRepulsionQ

theorem pairing_kernel_eq_11_over_48 :
    pairingKernelQ = 11 / 48 := by
  unfold pairingKernelQ
  rw [dark_fraction_of_total_states_eq, pairing_repulsion_eq_1_over_12]
  norm_num

/-- Pairing kernel sign is attractive. -/
theorem pairing_kernel_attractive :
    0 < pairingKernelQ := by
  rw [pairing_kernel_eq_11_over_48]
  norm_num

/-- Structural Tc proxy in kelvin-equivalent units for the forced gate. -/
def tcStructuralQ : ℚ := 300 * (1 + pairingKernelQ)

theorem tc_structural_eq_1475_over_4 :
    tcStructuralQ = 1475 / 4 := by
  unfold tcStructuralQ
  rw [pairing_kernel_eq_11_over_48]
  norm_num

theorem tc_structural_ge_300 :
    300 ≤ tcStructuralQ := by
  rw [tc_structural_eq_1475_over_4]
  norm_num

/-- RTSC forced admissibility gate:
    SC lattice + triplet filling + attractive kernel + Tc>=300. -/
def rtscAdmissible : Prop :=
  forcedLatticeFamily = LatticeFamily.simpleCubic ∧
  forcedFillingMultiplicity = 3 ∧
  0 < pairingKernelQ ∧
  300 ≤ tcStructuralQ

theorem rtsc_gate_admissible :
    rtscAdmissible := by
  refine ⟨forced_lattice_family_simple_cubic, forced_filling_triplet, pairing_kernel_attractive, tc_structural_ge_300⟩

/-- Forced-gate witness theorem: one constrained family survives. -/
theorem exists_forced_rtsc_family :
    ∃ (L : LatticeFamily) (n : ℕ) (tc : ℚ),
      L = LatticeFamily.simpleCubic ∧
      n = 3 ∧
      tc = tcStructuralQ ∧
      300 ≤ tc := by
  refine ⟨LatticeFamily.simpleCubic, 3, tcStructuralQ, ?_⟩
  exact ⟨rfl, rfl, rfl, tc_structural_ge_300⟩

end Gutoe.RTSCAdmissibility
