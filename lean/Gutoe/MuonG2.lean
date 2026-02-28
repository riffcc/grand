import Mathlib
import Gutoe.FineStructure
import Gutoe.GaugeConstants
import Gutoe.Z3Uniqueness

namespace Gutoe.MuonG2

open Gutoe.FineStructure
open Gutoe.GaugeConstants
open Gutoe.Z3Uniqueness

/-- Structural gauge-generator count from the SM algebra:
    `8 + 3 + 1 = 12`. -/
def gaugeGeneratorCountQ : ℚ :=
  (((3 ^ 2 - 1) + (2 ^ 2 - 1) + 1 : ℕ) : ℚ)

theorem gauge_generator_count_eq_12 :
    gaugeGeneratorCountQ = 12 := by
  unfold gaugeGeneratorCountQ
  have h : (3 ^ 2 - 1) + (2 ^ 2 - 1) + 1 = 12 := total_gauge_bosons
  norm_num [h]

/-- Clifford complement count from the SU(2) triplet split:
    `2^4 - |SU(2)| = 16 - 3 = 13`. -/
def cliffordComplementQ : ℚ :=
  (((2 ^ 4) - magneticTriplet.card : ℕ) : ℚ)

theorem clifford_complement_eq_13 :
    cliffordComplementQ = 13 := by
  unfold cliffordComplementQ
  rw [su2_dim]
  norm_num

/-- Leading-order structural alpha in the finite-count lane: `α = 1/137`. -/
def alphaLeadingQ : ℚ := 1 / (alphaInverse 4 : ℚ)

theorem alpha_leading_eq_1_over_137 :
    alphaLeadingQ = 1 / 137 := by
  unfold alphaLeadingQ
  rw [alpha_inverse_d4]
  norm_num

/-- Structural candidate for the unresolved muon `g-2` gap:
    `Δa_μ,cand = α^3 / (12 * 13)`. -/
def muonG2GapCandidateQ : ℚ :=
  alphaLeadingQ ^ 3 / (gaugeGeneratorCountQ * cliffordComplementQ)

theorem muon_g2_gap_candidate_closed_form :
    muonG2GapCandidateQ = 1 / 401131068 := by
  unfold muonG2GapCandidateQ
  rw [alpha_leading_eq_1_over_137, gauge_generator_count_eq_12, clifford_complement_eq_13]
  ring_nf

/-- Reference unresolved gap lane used in the runtime report:
    `Δa_μ,ref = 2.49×10^-9 = 249 / 10^11`. -/
def muonG2GapReferenceQ : ℚ := 249 / 100000000000

/-- The structural candidate is within `1e-11` absolute of the reference gap lane. -/
theorem muon_g2_gap_candidate_within_1e11 :
    |muonG2GapCandidateQ - muonG2GapReferenceQ| < 1 / 100000000000 := by
  rw [muon_g2_gap_candidate_closed_form]
  unfold muonG2GapReferenceQ
  native_decide

end Gutoe.MuonG2
