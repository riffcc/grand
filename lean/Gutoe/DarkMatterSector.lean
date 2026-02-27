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
import Gutoe.GravityMetric

namespace Gutoe.DarkMatterSector

open Gutoe.DimensionalStructure
open Gutoe.GaugeGroupSM
open Gutoe.Z3Uniqueness
open Gutoe.GravityMetric

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

/-- Fraction of dark-sector candidate states in the (visible + dark) split. -/
def darkFractionOfTotalStates : ℚ :=
  (darkSectorCandidates.card : ℚ) /
    ((darkSectorCandidates.card + visibleSectorStates.card : ℕ) : ℚ)

/-- Exact finite-state dark fraction from the structural 11/5 split. -/
theorem dark_fraction_of_total_states_eq :
    darkFractionOfTotalStates = 5 / 16 := by
  unfold darkFractionOfTotalStates
  rcases visible_dark_state_count_split with ⟨hVis, hDark, _, hTot⟩
  rw [hVis, hDark]
  norm_num

/-- Total finite state count in the visible+dark split. -/
def totalFiniteStateCount : ℕ :=
  darkSectorCandidates.card + visibleSectorStates.card

/-- The visible+dark split spans all 16 Clifford states. -/
theorem total_finite_state_count_eq :
    totalFiniteStateCount = 16 := by
  unfold totalFiniteStateCount
  rcases visible_dark_state_count_split with ⟨_, _, _, hTot⟩
  simpa [Nat.add_comm] using hTot

/-- Grade-1 has four states in Cl(1,3). -/
theorem grade1_state_count_eq : grade1_4d.card = 4 := by
  decide

/-- Geometric dark amplification from non-grade-1 Clifford channels:
    (total states) - (grade-1 states) = 16 - 4 = 12. -/
def geometricDarkAmplificationQ : ℚ :=
  ((totalFiniteStateCount - grade1_4d.card : ℕ) : ℚ)

/-- Exact geometric amplification from shared Clifford counts. -/
theorem geometric_dark_amplification_eq :
    geometricDarkAmplificationQ = 12 := by
  unfold geometricDarkAmplificationQ
  rw [total_finite_state_count_eq, grade1_state_count_eq]
  norm_num

/-- Geometric branch dark/visible ratio after structural amplification. -/
def geometricDarkToVisibleRatio : ℚ :=
  geometricDarkAmplificationQ * darkToVisibleCountRatio

/-- Exact geometric branch dark/visible ratio. -/
theorem geometric_dark_to_visible_ratio_eq :
    geometricDarkToVisibleRatio = 60 / 11 := by
  unfold geometricDarkToVisibleRatio
  rw [geometric_dark_amplification_eq, dark_to_visible_count_ratio_eq]
  norm_num

/-- Geometric branch dark fraction in total matter:
    f = (ρ_dark/ρ_visible) / (1 + ρ_dark/ρ_visible). -/
def geometricDarkFractionOfMatter : ℚ :=
  geometricDarkToVisibleRatio / (1 + geometricDarkToVisibleRatio)

/-- Exact geometric dark fraction from the structural amplified ratio. -/
theorem geometric_dark_fraction_of_matter_eq :
    geometricDarkFractionOfMatter = 60 / 71 := by
  unfold geometricDarkFractionOfMatter
  rw [geometric_dark_to_visible_ratio_eq]
  norm_num

/-- Unified branch budget ratio (how much): use the geometric amplified ratio. -/
def unifiedBudgetDarkToVisibleRatio : ℚ :=
  geometricDarkToVisibleRatio

/-- Unified branch budget ratio is exactly `60/11`. -/
theorem unified_budget_dark_to_visible_ratio_eq :
    unifiedBudgetDarkToVisibleRatio = 60 / 11 := by
  exact geometric_dark_to_visible_ratio_eq

/-- Unified local ratio (where): particle-like local clustering modulated by `κ`. -/
noncomputable def unifiedLocalDarkToVisibleRatio (kappa : ℝ) : ℝ :=
  (5 / 11 : ℝ) * kappa

/-- Unified local ratio is nonnegative when `κ` is nonnegative. -/
theorem unified_local_dark_to_visible_ratio_nonneg
    {kappa : ℝ}
    (hk : 0 ≤ kappa) :
    0 ≤ unifiedLocalDarkToVisibleRatio kappa := by
  unfold unifiedLocalDarkToVisibleRatio
  have hbase : 0 ≤ (5 / 11 : ℝ) := by norm_num
  exact mul_nonneg hbase hk

/-- Particle-branch effective dark density from the count ratio. -/
def effectiveDarkDensityQ (rhoVisible : ℚ) : ℚ :=
  darkToVisibleCountRatio * rhoVisible

/-- Effective particle-branch dark density is nonnegative for nonnegative visible
density. -/
theorem effective_dark_density_nonneg
    {rhoVisible : ℚ}
    (hρ : 0 ≤ rhoVisible) :
    0 ≤ effectiveDarkDensityQ rhoVisible := by
  unfold effectiveDarkDensityQ
  have hratio : 0 ≤ darkToVisibleCountRatio := by
    rw [dark_to_visible_count_ratio_eq]
    norm_num
  exact mul_nonneg hratio hρ

/-- Total (visible + effective-dark) density is nonnegative when visible density is
nonnegative. -/
theorem total_density_nonneg
    {rhoVisible : ℚ}
    (hρ : 0 ≤ rhoVisible) :
    0 ≤ rhoVisible + effectiveDarkDensityQ rhoVisible := by
  have hdark : 0 ≤ effectiveDarkDensityQ rhoVisible :=
    effective_dark_density_nonneg hρ
  linarith

/-- Vacuum-source curvature factor from Einstein source-term split:
    `1 + ρ_Λ / ρ_visible`. -/
noncomputable def vacuumCurvatureBoost (rhoVisible rhoVacuum : ℝ) : ℝ :=
  1 + rhoVacuum / rhoVisible

/-- UV curvature factor from the shared lattice correction:
    `1 + λ_QG (l_P / r)^2`. -/
noncomputable def uvCurvatureBoost (lP r : ℝ) : ℝ :=
  1 + lambda_qg * (lP / r) ^ 2

/-- Derived Einstein/cosmology curvature factor used by GRAND-346:
    `κ(r) = uv * vacuum-source`. -/
noncomputable def einsteinCosmologyKappa (rhoVisible rhoVacuum lP r : ℝ) : ℝ :=
  uvCurvatureBoost lP r * vacuumCurvatureBoost rhoVisible rhoVacuum

/-- The vacuum-source factor is ≥ 1 when visible density is positive and vacuum
    density is nonnegative. -/
theorem vacuum_curvature_boost_ge_one
    {rhoVisible rhoVacuum : ℝ}
    (hVis : 0 < rhoVisible)
    (hVac : 0 ≤ rhoVacuum) :
    1 ≤ vacuumCurvatureBoost rhoVisible rhoVacuum := by
  unfold vacuumCurvatureBoost
  have hdiv : 0 ≤ rhoVacuum / rhoVisible := by
    exact div_nonneg hVac (le_of_lt hVis)
  linarith

/-- The UV factor is ≥ 1 away from `r = 0`. -/
theorem uv_curvature_boost_ge_one
    {lP r : ℝ}
    (hr : r ≠ 0) :
    1 ≤ uvCurvatureBoost lP r := by
  have _hr := hr
  unfold uvCurvatureBoost
  have hLam : 0 ≤ lambda_qg := le_of_lt lambda_qg_pos
  have hsq : 0 ≤ (lP / r) ^ 2 := sq_nonneg (lP / r)
  have hprod : 0 ≤ lambda_qg * (lP / r) ^ 2 := mul_nonneg hLam hsq
  linarith

/-- The combined Einstein/cosmology curvature factor is ≥ 1 under the physical
    positivity assumptions. -/
theorem einstein_cosmology_kappa_ge_one
    {rhoVisible rhoVacuum lP r : ℝ}
    (hVis : 0 < rhoVisible)
    (hVac : 0 ≤ rhoVacuum)
    (hr : r ≠ 0) :
    1 ≤ einsteinCosmologyKappa rhoVisible rhoVacuum lP r := by
  unfold einsteinCosmologyKappa
  have huv : 1 ≤ uvCurvatureBoost lP r := uv_curvature_boost_ge_one hr
  have hvac : 1 ≤ vacuumCurvatureBoost rhoVisible rhoVacuum :=
    vacuum_curvature_boost_ge_one hVis hVac
  have huv_nonneg : 0 ≤ uvCurvatureBoost lP r := by linarith
  have hone_nonneg : 0 ≤ (1 : ℝ) := by norm_num
  calc
    1 = 1 * 1 := by ring
    _ ≤ uvCurvatureBoost lP r * vacuumCurvatureBoost rhoVisible rhoVacuum := by
      exact mul_le_mul huv hvac hone_nonneg huv_nonneg
    _ = einsteinCosmologyKappa rhoVisible rhoVacuum lP r := by rfl

/-- Under physical positivity assumptions, unified local dark ratio is at least
    the particle baseline `5/11`. -/
theorem unified_local_ratio_ge_particle_baseline
    {rhoVisible rhoVacuum lP r : ℝ}
    (hVis : 0 < rhoVisible)
    (hVac : 0 ≤ rhoVacuum)
    (hr : r ≠ 0) :
    (5 / 11 : ℝ) ≤ unifiedLocalDarkToVisibleRatio
        (einsteinCosmologyKappa rhoVisible rhoVacuum lP r) := by
  unfold unifiedLocalDarkToVisibleRatio
  have hk : 1 ≤ einsteinCosmologyKappa rhoVisible rhoVacuum lP r :=
    einstein_cosmology_kappa_ge_one hVis hVac hr
  have hfac_nonneg : 0 ≤ (5 / 11 : ℝ) := by norm_num
  have hmul :=
    mul_le_mul_of_nonneg_left hk hfac_nonneg
  simpa [mul_assoc, mul_comm, mul_left_comm] using hmul

end Gutoe.DarkMatterSector
