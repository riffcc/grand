import Mathlib
import Gutoe.EvenSubalgebraSuppression
import Gutoe.VacuumEnergyBounds
import Gutoe.DimensionalStructure
import Gutoe.Z3Uniqueness

/-!
GUTOE — Overdetermined Topology Ratios

This lane promotes the "ratio closure" claim into hard constraints:

`G = branching * void * eta * infra = 1`

with each factor traced to Cl(1,3) counting invariants.
Any one factor is algebraically determined by the other three once `G=1`.
-/

namespace Gutoe.OverdeterminedTopologyRatios

open Gutoe.EvenSubalgebraSuppression
open Gutoe.VacuumEnergyBounds
open Gutoe.DimensionalStructure
open Gutoe.Z3Uniqueness

/-- Structural grade-1 count (`{γ⁰,γ¹,γ²,γ³}`). -/
def grade1CardQ : ℚ := (grade1_4d.card : ℚ)

/-- Structural grade-2 count (bivectors). -/
def grade2CardQ : ℚ := (grade2_4d.card : ℚ)

/-- Structural total basis count `2^4 = 16`. -/
def totalBasisQ : ℚ := (2 ^ 4 : ℚ)

theorem grade1_card_eq_4 : grade1CardQ = 4 := by
  unfold grade1CardQ
  native_decide

theorem grade2_card_eq_6 : grade2CardQ = 6 := by
  unfold grade2CardQ
  native_decide

theorem total_basis_eq_16 : totalBasisQ = 16 := by
  unfold totalBasisQ
  norm_num

/-- Structural branching ratio from `|Z₃| = 3`. -/
def branchingStructuralQ : ℚ := (magneticTriplet.card : ℚ)

/-- Structural void merge ratio from Cl(1,3) split. -/
def voidStructuralQ : ℚ := voidFractionQ

/-- Structural transport ratio: visible vectors over bivectors (`4/6`). -/
def etaStructuralQ : ℚ := grade1CardQ / grade2CardQ

/-- Structural infrastructure ratio: full basis over bivectors (`16/6`). -/
def infraStructuralQ : ℚ := totalBasisQ / grade2CardQ

theorem branching_structural_eq_3 : branchingStructuralQ = 3 := by
  unfold branchingStructuralQ
  norm_num [su2_dim]

theorem void_structural_eq_3_16 : voidStructuralQ = (3 : ℚ) / 16 := by
  exact void_fraction_eq_3_16

theorem eta_structural_eq_2_3 : etaStructuralQ = (2 : ℚ) / 3 := by
  unfold etaStructuralQ
  rw [grade1_card_eq_4, grade2_card_eq_6]
  norm_num

theorem infra_structural_eq_8_3 : infraStructuralQ = (8 : ℚ) / 3 := by
  unfold infraStructuralQ
  rw [total_basis_eq_16, grade2_card_eq_6]
  norm_num

/-- Ratio-closure gain functional. -/
def topologyGainQ : ℚ :=
  branchingStructuralQ * voidStructuralQ * etaStructuralQ * infraStructuralQ

/-- Full counting-closure identity: `3 * (3/16) * (2/3) * (8/3) = 1`. -/
theorem topology_gain_eq_one : topologyGainQ = 1 := by
  unfold topologyGainQ
  rw [branching_structural_eq_3, void_structural_eq_3_16,
      eta_structural_eq_2_3, infra_structural_eq_8_3]
  norm_num

/-- Equivalent `144/144 = 1` closed form. -/
theorem topology_gain_eq_144_over_144 :
    topologyGainQ = (144 : ℚ) / 144 := by
  rw [topology_gain_eq_one]
  norm_num

/-- This closure agrees with the existing Z3/void split gain lane. -/
theorem topology_gain_matches_existing_geff :
    topologyGainQ = geffZ3VoidSplitQ := by
  rw [topology_gain_eq_one, geff_z3_void_split_eq_one]

/-- This closure is compatible with the independent even/odd `1/2` closure lane. -/
theorem topology_closure_consistent_with_even_lane :
    topologyGainQ = 1 ∧ geffCanonicalQ = 1 := by
  exact ⟨topology_gain_eq_one, geff_canonical_eq_one⟩

/-- Overdetermination: if `G=1` and the first three ratios are fixed structurally,
the infrastructure ratio is uniquely forced. -/
theorem infra_forced_by_unit_gain
    (ξ : ℚ)
    (h : branchingStructuralQ * voidStructuralQ * etaStructuralQ * ξ = 1) :
    ξ = infraStructuralQ := by
  rw [branching_structural_eq_3, void_structural_eq_3_16, eta_structural_eq_2_3] at h
  ring_nf at h
  have hξ : ξ = (8 : ℚ) / 3 := by
    linarith
  rw [infra_structural_eq_8_3]
  exact hξ

/-- Overdetermination: if `G=1` and branching/void/infra are fixed, transport
ratio `eta` is uniquely forced. -/
theorem eta_forced_by_unit_gain
    (η : ℚ)
    (h : branchingStructuralQ * voidStructuralQ * η * infraStructuralQ = 1) :
    η = etaStructuralQ := by
  rw [branching_structural_eq_3, void_structural_eq_3_16, infra_structural_eq_8_3] at h
  ring_nf at h
  have hη : η = (2 : ℚ) / 3 := by
    linarith
  rw [eta_structural_eq_2_3]
  exact hη

/-- Overdetermination: if `G=1` and branching/eta/infra are fixed, void ratio
is uniquely forced. -/
theorem void_forced_by_unit_gain
    (ν : ℚ)
    (h : branchingStructuralQ * ν * etaStructuralQ * infraStructuralQ = 1) :
    ν = voidStructuralQ := by
  rw [branching_structural_eq_3, eta_structural_eq_2_3, infra_structural_eq_8_3] at h
  ring_nf at h
  have hν : ν = (3 : ℚ) / 16 := by linarith
  rw [void_structural_eq_3_16]
  exact hν

/-- Overdetermination: if `G=1` and void/eta/infra are fixed, branching is
uniquely forced. -/
theorem branching_forced_by_unit_gain
    (β : ℚ)
    (h : β * voidStructuralQ * etaStructuralQ * infraStructuralQ = 1) :
    β = branchingStructuralQ := by
  rw [void_structural_eq_3_16, eta_structural_eq_2_3, infra_structural_eq_8_3] at h
  ring_nf at h
  have hβ : β = 3 := by linarith
  rw [branching_structural_eq_3]
  exact hβ

end Gutoe.OverdeterminedTopologyRatios
