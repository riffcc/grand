import Mathlib
import Gutoe.RTSCWitnessClosure
import Gutoe.NuclearFirstPrinciples

namespace Gutoe.RTSCIsotopeStability

open Gutoe.RTSCWitnessClosure
open Gutoe.NuclearFirstPrinciples

/-- Witness isotope set (Z,N) used by the RTSC forced-candidate lane.
    These are concrete isotope data points attached to the six witness elements. -/
def rtscWitnessIsotopes : Finset (ℕ × ℕ) :=
  {(24, 28), (30, 34), (42, 56), (48, 66), (72, 108), (78, 117)}

theorem rtsc_witness_isotopes_card :
    rtscWitnessIsotopes.card = 6 := by
  decide

/-- Coarse structural binding proxy (MeV-scale up to global normalization):
    volume - Coulomb - asymmetry - pairing floor.
    Uses shared SEMF coefficients derived from Cl(1,3) primitives. -/
def boundScoreQ (z n : ℕ) : ℚ :=
  let a : ℚ := ((z + n : ℕ) : ℚ)
  semfAVQ * a
    - semfACQ * (((z * (z - 1) : ℕ) : ℚ) / a)
    - semfAAQ * (((a - 2 * z) ^ (2 : ℕ)) / a)
    - semfAPQ

/-- All witness isotopes project to the forced six-element witness-Z set. -/
theorem witness_isotopes_project_to_witness_z :
    ∀ zn ∈ rtscWitnessIsotopes, zn.1 ∈ rtscWitnessZ := by
  intro zn hzn
  fin_cases hzn <;> decide

/-- All witness isotopes are structurally bound in the coarse first-principles proxy. -/
theorem witness_isotopes_structurally_bound :
    ∀ zn ∈ rtscWitnessIsotopes, 0 < boundScoreQ zn.1 zn.2 := by
  intro zn hzn
  fin_cases hzn
  · -- Cr-52
    unfold boundScoreQ
    rw [semf_a_v_eq_95_over_6, semf_a_c_eq_2_over_3, semf_a_a_eq_23, semf_a_p_eq_12]
    norm_num
  · -- Zn-64
    unfold boundScoreQ
    rw [semf_a_v_eq_95_over_6, semf_a_c_eq_2_over_3, semf_a_a_eq_23, semf_a_p_eq_12]
    norm_num
  · -- Mo-98
    unfold boundScoreQ
    rw [semf_a_v_eq_95_over_6, semf_a_c_eq_2_over_3, semf_a_a_eq_23, semf_a_p_eq_12]
    norm_num
  · -- Cd-114
    unfold boundScoreQ
    rw [semf_a_v_eq_95_over_6, semf_a_c_eq_2_over_3, semf_a_a_eq_23, semf_a_p_eq_12]
    norm_num
  · -- Hf-180
    unfold boundScoreQ
    rw [semf_a_v_eq_95_over_6, semf_a_c_eq_2_over_3, semf_a_a_eq_23, semf_a_p_eq_12]
    norm_num
  · -- Pt-195
    unfold boundScoreQ
    rw [semf_a_v_eq_95_over_6, semf_a_c_eq_2_over_3, semf_a_a_eq_23, semf_a_p_eq_12]
    norm_num

/-- Bundle theorem for RTSC isotope/stability closure lane. -/
theorem rtsc_isotope_stability_bundle :
    rtscWitnessIsotopes.card = 6 ∧
    (∀ zn ∈ rtscWitnessIsotopes, zn.1 ∈ rtscWitnessZ) ∧
    (∀ zn ∈ rtscWitnessIsotopes, 0 < boundScoreQ zn.1 zn.2) := by
  exact ⟨rtsc_witness_isotopes_card, witness_isotopes_project_to_witness_z, witness_isotopes_structurally_bound⟩

end Gutoe.RTSCIsotopeStability

