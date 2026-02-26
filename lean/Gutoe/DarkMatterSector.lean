/-
 * GUTOE — Dark Sector Candidates from Z₃ Orbit Structure
 * Copyright (C) 2026  Riff Labs
 *
 * AGPL-3.0-or-later
 *
 * GRAND-346:
 *   Isolate the Clifford/Z₃ sectors that are disjoint from the SM interaction
 *   carrier orbits used by the current gauge/matter lane.
 *
 * No `sorry`.
 -/

import Mathlib
import Gutoe.GaugeGroupSM
import Gutoe.Z3Uniqueness

namespace Gutoe.DarkMatterSector

open Gutoe.DimensionalStructure
open Gutoe.GaugeGroupSM
open Gutoe.Z3Uniqueness

/-- Low-grade Z₃ singlet pair used in the visible lane (`scalar + lepton`). -/
def lowSingletPair : Finset ℕ := {1, 2}

/-- High-grade Z₃ singlet pair (`γ¹²³`, `γ⁰¹²³`). -/
def highSingletPair : Finset ℕ := {15, 16}

/-- SM interaction carrier orbits in the current finite Cl(1,3) lane. -/
def smInteractionCarrier : Finset ℕ :=
  leptonState ∪ quarkTriplet ∪ emTriplet ∪ magneticTriplet

/-- Candidate dark sector: dual-EM triplet plus high-grade singlet pair. -/
def darkSectorCandidates : Finset ℕ :=
  dualEmTriplet ∪ highSingletPair

/-- Visible finite-state sector used by the current lattice lane. -/
def visibleSectorStates : Finset ℕ :=
  lowSingletPair ∪ quarkTriplet ∪ emTriplet ∪ magneticTriplet

/-- Candidate dark sector is exactly the dual-EM orbit plus high singlet pair. -/
theorem dark_sector_candidates_exact :
    darkSectorCandidates = dualEmTriplet ∪ ({15, 16} : Finset ℕ) := by
  decide

/-- High singlet pair is contained in the Z₃ singlet set. -/
theorem high_singlet_pair_in_z3_singlets :
    highSingletPair ⊆ z3_singlets := by
  simpa [highSingletPair] using (right_handed_singlet_pair).1

/-- Candidate dark sector is Z₃-invariant (closed under the Z₃ action). -/
theorem dark_sector_z3_closed :
    ∀ s ∈ darkSectorCandidates, z3_4d s ∈ darkSectorCandidates := by
  decide

/-- Candidate dark sector is disjoint from the SM interaction carrier. -/
theorem dark_sector_disjoint_from_sm_carrier :
    darkSectorCandidates ∩ smInteractionCarrier = ∅ := by
  decide

/-- Finite-state split: visible lane has 11 states and dark candidates 5 states. -/
theorem visible_dark_state_count_split :
    visibleSectorStates.card = 11 ∧
    darkSectorCandidates.card = 5 ∧
    visibleSectorStates ∩ darkSectorCandidates = ∅ ∧
    visibleSectorStates.card + darkSectorCandidates.card = 16 := by
  decide

/-- Count-level ratio from the structural 11/5 split. -/
def darkToVisibleCountRatio : ℚ :=
  (darkSectorCandidates.card : ℚ) / (visibleSectorStates.card : ℚ)

/-- Exact finite-state dark/visible ratio from the current orbit split. -/
theorem dark_to_visible_count_ratio_eq :
    darkToVisibleCountRatio = 5 / 11 := by
  unfold darkToVisibleCountRatio
  rcases visible_dark_state_count_split with ⟨hVis, hDark, _, _⟩
  rw [hVis, hDark]
  norm_num

end Gutoe.DarkMatterSector
