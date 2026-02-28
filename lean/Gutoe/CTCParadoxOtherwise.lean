import Mathlib
import Gutoe.CliffordStructure

/-!
GUTOE — CTC paradox "otherwise" lane.

This file formalizes two complementary statements:
1) In strict single-history closure, a forced traveler-kill loop is inconsistent.
2) In a branch-split closure model, a paradox-style assignment is consistent.

This is a logic theorem about closure axioms, not a claim about realized
cosmological topology.
-/

namespace Gutoe.CTCParadoxOtherwise

open Gutoe.CliffordStructure

/-- Single-history closure constraints for a grandfather-style loop.
`a0`: ancestor alive before intervention
`t`: traveler exists
`k`: kill action
`a1`: ancestor alive after intervention

Rules:
- traveler exists iff ancestor line exists
- traveler performs kill
- kill removes ancestor state
- loop closure enforces `a1 = a0`
-/
def singleHistoryConstraints (a0 t k a1 : Bool) : Prop :=
  t = a0 ∧
  k = t ∧
  a1 = (a0 && (!k)) ∧
  a1 = a0

/-- Forced traveler paradox branch is inconsistent in strict single-history closure. -/
theorem forced_traveler_unsat_single_history :
    ¬ ∃ a0 k a1 : Bool, singleHistoryConstraints a0 true k a1 := by
  intro h
  rcases h with ⟨a0, k, a1, ht, hk, ha1, hclosure⟩
  have ha0 : a0 = true := by simpa using ht.symm
  have hk' : k = true := by simpa [ha0] using hk
  have ha1false : a1 = false := by
    rw [ha1, ha0, hk']
    decide
  have ha1true : a1 = true := by simpa [ha0] using hclosure
  rw [ha1true] at ha1false
  cases ha1false

/-- Branch-split closure constraints:
origin branch `O`: traveler closure
target branch `T`: intervention dynamics
-/
def branchSplitConstraints
    (oA0 oT oA1 tA0 tK tA1 : Bool) : Prop :=
  oT = oA0 ∧
  oA1 = oA0 ∧
  tK = oT ∧
  tA1 = (tA0 && (!tK))

/-- Paradox-style assignment exists in branch-split closure:
origin retains traveler lineage, target ancestor is removed.
-/
theorem paradox_style_exists_branch_split :
    ∃ oA0 oT oA1 tA0 tK tA1 : Bool,
      branchSplitConstraints oA0 oT oA1 tA0 tK tA1 ∧
      oA0 = true ∧ oT = true ∧ tA0 = true ∧ tK = true ∧ tA1 = false := by
  refine ⟨true, true, true, true, true, false, ?_⟩
  constructor
  · unfold branchSplitConstraints
    decide
  · repeat constructor <;> decide

/-- Cl(1,3)-anchored "otherwise" theorem:
with the timelike sign fixed by the algebra, branch-split paradox-style
assignments are logically consistent while strict single-history forced loops
remain inconsistent. -/
theorem cl13_anchored_otherwise :
    minkowskiQF (e 0) = -1 ∧
      (¬ ∃ a0 k a1 : Bool, singleHistoryConstraints a0 true k a1) ∧
      (∃ oA0 oT oA1 tA0 tK tA1 : Bool,
          branchSplitConstraints oA0 oT oA1 tA0 tK tA1 ∧
          oA0 = true ∧ oT = true ∧ tA0 = true ∧ tK = true ∧ tA1 = false) := by
  refine ⟨minkowskiQF_e0, forced_traveler_unsat_single_history, ?_⟩
  exact paradox_style_exists_branch_split

-- ── Deutsch fixed-point lane and matter bookkeeping ─────────────────────────

/-- Deutsch NOT-map for the traveler branch weight. -/
def notMap (p : ℝ) : ℝ := 1 - p

/-- Fixed points of `notMap` are exactly `p = 1/2`. -/
theorem notMap_fixedpoint_iff_half (p : ℝ) :
    p = notMap p ↔ p = (1 : ℝ) / 2 := by
  unfold notMap
  constructor <;> intro h <;> linarith

/-- At the fixed point, both branch weights are exactly 1/2. -/
theorem deutsch_branch_weights_half (p : ℝ) (h : p = notMap p) :
    p = (1 : ℝ) / 2 ∧ notMap p = (1 : ℝ) / 2 := by
  have hp : p = (1 : ℝ) / 2 := (notMap_fixedpoint_iff_half p).mp h
  constructor
  · exact hp
  · rw [notMap, hp]
    norm_num

/-- Traveler packet mass present on the local slice. -/
def travelerMassPresent (mT : ℝ) (travelerPresent : Bool) : ℝ :=
  if travelerPresent then mT else 0

/-- Traveler packet mass in the loop channel (complement of local presence). -/
def travelerMassChannel (mT : ℝ) (travelerPresent : Bool) : ℝ :=
  if travelerPresent then 0 else mT

/-- Traveler packet mass is conserved across local+channel bookkeeping. -/
theorem traveler_packet_mass_conserved (mT : ℝ) (travelerPresent : Bool) :
    travelerMassPresent mT travelerPresent +
      travelerMassChannel mT travelerPresent = mT := by
  cases travelerPresent <;> simp [travelerMassPresent, travelerMassChannel]

/-- Total mass bookkeeping with ancestor packet `mA` stays constant:
ancestor mass plus traveler packet (local or channel) is invariant. -/
theorem total_mass_bookkeeping_conserved (mA mT : ℝ) (travelerPresent : Bool) :
    mA + travelerMassPresent mT travelerPresent +
      travelerMassChannel mT travelerPresent = mA + mT := by
  nlinarith [traveler_packet_mass_conserved mT travelerPresent]

end Gutoe.CTCParadoxOtherwise
