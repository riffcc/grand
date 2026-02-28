import Mathlib
import Gutoe.GaugeConstants
import Gutoe.StrongCouplingCInfBridge
import Gutoe.TriangulatedConstants
import Gutoe.TriangulatedClosureUniqueness
import Gutoe.KoideMasses
import Gutoe.MuonG2
import Gutoe.BBN

namespace Gutoe.FiveLockClosure

open Gutoe.GaugeConstants
open Gutoe.StrongCouplingCInfBridge
open Gutoe.TriangulatedConstants
open Gutoe.TriangulatedClosureUniqueness
open Gutoe.KoideMasses
open Gutoe.MuonG2
open Gutoe.BBN

/-- Single-RGE lock inputs from shared structural closures:
    EW bridge input and strong-coupling corrected input are both fixed. -/
theorem rge_input_lock :
    weinbergMZStructuralQ = (3 / 13 : ℚ) + 8 / (137 ^ 2 : ℚ) ∧
    alphaSStructuralCorrectedQ = (16 / 137 : ℚ) * (67 / 66 : ℚ) := by
  exact ⟨weinberg_mz_structural_eq, alpha_s_structural_corrected_eq⟩

/-- Neutrino triangulation lock at the parameter level:
    `p` and `κ` are near frozen anchors and grammar closure is unique. -/
theorem neutrino_triangulation_parameter_lock :
    |pCandidateQ - pFrozenQ| < (1 / 50000 : ℚ) ∧
    |kappaCandidateQ - kappaFrozenQ| < (1 / 50000 : ℚ) ∧
    pGoodSigns.card = 1 ∧
    kappaGoodTuples.card = 1 ∧
    ewGoodTuples.card = 1 := by
  refine ⟨p_candidate_close_to_frozen, kappa_candidate_close_to_frozen, ?_⟩
  exact constrained_grammar_uniqueness_closure

/-- Charged-lepton hierarchy lock in the Z₃/grade lane:
    Koide target and grade-ratio target coincide exactly at `2/3`. -/
theorem charged_lepton_hierarchy_lock :
    koideClifford = (2 / 3 : ℚ) ∧
    (leptonGradeDim : ℚ) / gaugeGradeDim = (2 / 3 : ℚ) := by
  exact ⟨koide_clifford_is_2_3, grade1_over_grade2_is_2_3⟩

/-- Dual g-2 lock (formal side): the structural muon unresolved-gap candidate
    is tightly bounded to the reference lane. -/
theorem g2_dual_lock :
    |muonG2GapCandidateQ - muonG2GapReferenceQ| < (1 / 100000000000 : ℚ) := by
  exact muon_g2_gap_candidate_within_1e11

/-- BBN three-isotope lock at the reference baryogenesis scale. -/
theorem bbn_three_isotope_lock :
    primordialHelium4MassFraction eta10Ref = (ypTargetQ : ℝ) ∧
    primordialDeuteriumRatio eta10Ref = (dOverHTargetQ : ℝ) ∧
    primordialHelium3Ratio eta10Ref = (he3OverHTargetQ : ℝ) := by
  exact ⟨primordial_helium4_at_reference, primordial_deuterium_at_reference,
    primordial_helium3_at_reference⟩

/-- Joint formal closure bundle for the five-lock CI package. -/
theorem five_lock_formal_closure :
    (weinbergMZStructuralQ = (3 / 13 : ℚ) + 8 / (137 ^ 2 : ℚ) ∧
      alphaSStructuralCorrectedQ = (16 / 137 : ℚ) * (67 / 66 : ℚ)) ∧
    (|pCandidateQ - pFrozenQ| < (1 / 50000 : ℚ) ∧
      |kappaCandidateQ - kappaFrozenQ| < (1 / 50000 : ℚ) ∧
      pGoodSigns.card = 1 ∧ kappaGoodTuples.card = 1 ∧ ewGoodTuples.card = 1) ∧
    (koideClifford = (2 / 3 : ℚ) ∧
      (leptonGradeDim : ℚ) / gaugeGradeDim = (2 / 3 : ℚ)) ∧
    (|muonG2GapCandidateQ - muonG2GapReferenceQ| < (1 / 100000000000 : ℚ)) ∧
    (primordialHelium4MassFraction eta10Ref = (ypTargetQ : ℝ) ∧
      primordialDeuteriumRatio eta10Ref = (dOverHTargetQ : ℝ) ∧
      primordialHelium3Ratio eta10Ref = (he3OverHTargetQ : ℝ)) := by
  exact ⟨rge_input_lock, neutrino_triangulation_parameter_lock,
    charged_lepton_hierarchy_lock, g2_dual_lock, bbn_three_isotope_lock⟩

end Gutoe.FiveLockClosure
