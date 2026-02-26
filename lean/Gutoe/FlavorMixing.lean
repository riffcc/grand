import Mathlib
import Gutoe.FineStructure
import Gutoe.Z3Uniqueness

namespace Gutoe.FlavorMixing

open Gutoe.Z3Uniqueness

/-- CKM sine entry `s12` from Clifford dimension + SU(2) orbit mismatch. -/
noncomputable def ckmSin12 : ℝ :=
  1 / Real.sqrt (((2 ^ 4) + magneticTriplet.card : ℕ) : ℝ)

/-- CKM sine entry `s23` from grade-1 × grade-2 suppression. -/
noncomputable def ckmSin23 : ℝ :=
  1 / (((Nat.choose 4 1) * (Nat.choose 4 2) : ℕ) : ℝ)

/-- CKM sine entry `s13` from Clifford × augmented-dimension suppression. -/
noncomputable def ckmSin13 : ℝ :=
  1 / ((((2 ^ 4) * ((2 ^ 4) + 1)) : ℕ) : ℝ)

/-- CKM phase: Z₃ base phase π/3 plus a lattice shift arctan(1/7). -/
noncomputable def ckmDelta : ℝ :=
  Real.pi / 3 + Real.arctan (1 / (((Nat.choose 4 2) + 1 : ℕ) : ℝ))

/-- PMNS sine entry `s12` from `(grade1)/(Clifford−SU(2)) = 4/13`. -/
noncomputable def pmnsSin12 : ℝ :=
  Real.sqrt (((Nat.choose 4 1 : ℕ) : ℝ) / (((2 ^ 4) - magneticTriplet.card : ℕ) : ℝ))

/-- PMNS sine entry `s23` from `(grade1)/(grade2+1) = 4/7`. -/
noncomputable def pmnsSin23 : ℝ :=
  Real.sqrt (((Nat.choose 4 1 : ℕ) : ℝ) / (((Nat.choose 4 2) + 1 : ℕ) : ℝ))

/-- PMNS sine entry `s13` from `(grade2+1)⁻¹ = 1/7`. -/
noncomputable def pmnsSin13 : ℝ :=
  1 / ((((Nat.choose 4 2) + 1 : ℕ) : ℝ))

/-- PMNS phase: π plus a Z₃ correction arctan(1/3). -/
noncomputable def pmnsDelta : ℝ :=
  Real.pi + Real.arctan (1 / (magneticTriplet.card : ℝ))

/-- CKM texture coupling `M_d(1,2)` coefficient = `(4/5) * λ` with `λ = 1/√19`. -/
noncomputable def ckmMd12Coeff : ℝ :=
  ((Nat.choose 4 1 : ℕ) : ℝ) / (((Nat.choose 4 1) + 1 : ℕ) : ℝ) *
  (1 / Real.sqrt (((2 ^ 4) + magneticTriplet.card : ℕ) : ℝ))

/-- PMNS texture coupling `M_ν(1,2)` coefficient = `(3/4) * √(4/13)`. -/
noncomputable def pmnsMnu12Coeff : ℝ :=
  (magneticTriplet.card : ℝ) / ((Nat.choose 4 1 : ℕ) : ℝ) * pmnsSin12

/-- PMNS texture coupling `M_ν(2,3)` coefficient = `(2/3) * √(4/7)`. -/
noncomputable def pmnsMnu23Coeff : ℝ :=
  (2 : ℝ) / (magneticTriplet.card : ℝ) * pmnsSin23

/-- Standard Jarlskog invariant from mixing entries and phase. -/
noncomputable def jarlskog (s12 s23 s13 δ : ℝ) : ℝ :=
  let c12 := Real.sqrt (1 - s12 ^ 2)
  let c23 := Real.sqrt (1 - s23 ^ 2)
  let c13 := Real.sqrt (1 - s13 ^ 2)
  c12 * c23 * c13 ^ 2 * s12 * s23 * s13 * Real.sin δ

/-- Structural evaluations used for parity with Rust harness. -/
theorem ckm_structural_values :
    ckmSin23 = (1 : ℝ) / 24 ∧
    ckmSin13 = (1 : ℝ) / 272 := by
  have h10 : (Nat.choose 4 1 : ℕ) = 4 := by native_decide
  have h20 : (Nat.choose 4 2 : ℕ) = 6 := by native_decide
  constructor
  · simp [ckmSin23, h10, h20]
  · simp [ckmSin13]

/-- Texture coupling coefficients used by Rust `*_from_textures` path. -/
theorem texture_coeff_structural_values :
    ckmMd12Coeff = (4 : ℝ) / 5 * (1 / Real.sqrt 19) ∧
    pmnsMnu12Coeff = (3 : ℝ) / 4 * Real.sqrt ((4 : ℝ) / 13) ∧
    pmnsMnu23Coeff = (2 : ℝ) / 3 * Real.sqrt ((4 : ℝ) / 7) := by
  have hs : magneticTriplet.card = 3 := su2_dim
  have h10 : (Nat.choose 4 1 : ℕ) = 4 := by native_decide
  have h20 : (Nat.choose 4 2 : ℕ) = 6 := by native_decide
  constructor
  · simp [ckmMd12Coeff, h10, hs]
  · constructor
    · simp [pmnsMnu12Coeff, pmnsSin12, h10, hs]
    · simp [pmnsMnu23Coeff, pmnsSin23, h10, h20, hs]

/-- PMNS structural evaluations from the same Clifford counts. -/
theorem pmns_structural_values :
    pmnsSin12 ^ 2 = (4 : ℝ) / 13 ∧
    pmnsSin23 ^ 2 = (4 : ℝ) / 7 ∧
    pmnsSin13 = (1 : ℝ) / 7 := by
  have hs : magneticTriplet.card = 3 := su2_dim
  have h10 : (Nat.choose 4 1 : ℕ) = 4 := by native_decide
  have h20 : (Nat.choose 4 2 : ℕ) = 6 := by native_decide
  constructor
  · have harg :
        0 ≤ ((Nat.choose 4 1 : ℕ) : ℝ) / (((2 ^ 4) - magneticTriplet.card : ℕ) : ℝ) := by
        norm_num [h10, hs]
    calc
      pmnsSin12 ^ 2
          = (Real.sqrt (((Nat.choose 4 1 : ℕ) : ℝ) / (((2 ^ 4) - magneticTriplet.card : ℕ) : ℝ))) ^ 2 := by
              rfl
      _ = ((Nat.choose 4 1 : ℕ) : ℝ) / (((2 ^ 4) - magneticTriplet.card : ℕ) : ℝ) := by
            simpa using Real.sq_sqrt harg
      _ = (4 : ℝ) / 13 := by norm_num [h10, hs]
  · constructor
    · have harg :
          0 ≤ ((Nat.choose 4 1 : ℕ) : ℝ) / (((Nat.choose 4 2) + 1 : ℕ) : ℝ) := by
          norm_num [h10, h20]
      calc
        pmnsSin23 ^ 2
            = (Real.sqrt (((Nat.choose 4 1 : ℕ) : ℝ) / (((Nat.choose 4 2) + 1 : ℕ) : ℝ))) ^ 2 := by
                rfl
        _ = ((Nat.choose 4 1 : ℕ) : ℝ) / (((Nat.choose 4 2) + 1 : ℕ) : ℝ) := by
              simpa using Real.sq_sqrt harg
        _ = (4 : ℝ) / 7 := by norm_num [h10, h20]
    · simp [pmnsSin13, h20]

/-- CKM phase lies in (0, π), so CP term contributes with positive sine. -/
theorem ckm_delta_in_open_pi : 0 < ckmDelta ∧ ckmDelta < Real.pi := by
  have h20 : (Nat.choose 4 2 : ℕ) = 6 := by native_decide
  have harctan_pos : 0 < Real.arctan (1 / (((Nat.choose 4 2) + 1 : ℕ) : ℝ)) := by
    have : (0 : ℝ) < 1 / (((Nat.choose 4 2) + 1 : ℕ) : ℝ) := by
      norm_num [h20]
    exact Real.arctan_pos.mpr this
  have harctan_lt : Real.arctan (1 / (((Nat.choose 4 2) + 1 : ℕ) : ℝ)) < Real.pi / 2 := by
    exact Real.arctan_lt_pi_div_two _
  constructor
  · have h : 0 < Real.pi / 3 := by positivity
    simpa [ckmDelta] using add_pos h harctan_pos
  · have hsum :
        ckmDelta < Real.pi / 3 + Real.pi / 2 := by
      simpa [ckmDelta, add_comm, add_left_comm, add_assoc] using
        (add_lt_add_left harctan_lt (Real.pi / 3))
    have hbound : Real.pi / 3 + Real.pi / 2 < Real.pi := by
      nlinarith [Real.pi_pos]
    exact lt_trans hsum hbound

/-- CKM Jarlskog sign from structural phase window and positive suppressions. -/
theorem ckm_jarlskog_positive : 0 < jarlskog ckmSin12 ckmSin23 ckmSin13 ckmDelta := by
  have hdelta : 0 < ckmDelta ∧ ckmDelta < Real.pi := ckm_delta_in_open_pi
  have hsin : 0 < Real.sin ckmDelta := Real.sin_pos_of_pos_of_lt_pi hdelta.1 hdelta.2
  have hs12 : 0 < ckmSin12 := by
    have hroot : 0 < Real.sqrt ((((2 ^ 4) + magneticTriplet.card : ℕ) : ℝ)) := by
      have : (0 : ℝ) < (((2 ^ 4) + magneticTriplet.card : ℕ) : ℝ) := by
        have hs : magneticTriplet.card = 3 := su2_dim
        norm_num [hs]
      exact Real.sqrt_pos.mpr this
    exact one_div_pos.mpr hroot
  have hs23 : 0 < ckmSin23 := by
    have h : ckmSin23 = (1 : ℝ) / 24 := ckm_structural_values.1
    linarith [h]
  have hs13 : 0 < ckmSin13 := by
    have h : ckmSin13 = (1 : ℝ) / 272 := ckm_structural_values.2
    linarith [h]
  have hc12 : 0 < Real.sqrt (1 - ckmSin12 ^ 2) := by
    have hs : magneticTriplet.card = 3 := su2_dim
    have hsq : ckmSin12 ^ 2 = (1 : ℝ) / 19 := by
      have hsqrt_sq : (Real.sqrt (19 : ℝ)) ^ 2 = (19 : ℝ) := by
        have h19 : (0 : ℝ) ≤ 19 := by positivity
        simp [Real.sq_sqrt h19]
      have hfrac : (1 / Real.sqrt (19 : ℝ)) ^ 2 = (1 : ℝ) / 19 := by
        have hsqrt_ne : Real.sqrt (19 : ℝ) ≠ 0 := by positivity
        have hfrac_sq :
            (1 / Real.sqrt (19 : ℝ)) ^ 2 = (1 : ℝ) / ((Real.sqrt (19 : ℝ)) ^ 2) := by
          field_simp [hsqrt_ne]
        calc
          (1 / Real.sqrt (19 : ℝ)) ^ 2
              = (1 : ℝ) / ((Real.sqrt (19 : ℝ)) ^ 2) := hfrac_sq
          _ = (1 : ℝ) / 19 := by simp [hsqrt_sq]
      calc
        ckmSin12 ^ 2
            = (1 / Real.sqrt (19 : ℝ)) ^ 2 := by
                simp [ckmSin12, hs]
        _ = (1 : ℝ) / 19 := hfrac
    have harg : 0 < 1 - ckmSin12 ^ 2 := by nlinarith [hsq]
    exact Real.sqrt_pos.mpr harg
  have hc23 : 0 < Real.sqrt (1 - ckmSin23 ^ 2) := by
    have h23 : ckmSin23 = (1 : ℝ) / 24 := ckm_structural_values.1
    have hsq : ckmSin23 ^ 2 = (1 : ℝ) / 576 := by nlinarith [h23]
    have harg : 0 < 1 - ckmSin23 ^ 2 := by nlinarith [hsq]
    exact Real.sqrt_pos.mpr harg
  have hc13 : 0 < Real.sqrt (1 - ckmSin13 ^ 2) := by
    have h13 : ckmSin13 = (1 : ℝ) / 272 := ckm_structural_values.2
    have hsq : ckmSin13 ^ 2 = (1 : ℝ) / 73984 := by nlinarith [h13]
    have harg : 0 < 1 - ckmSin13 ^ 2 := by nlinarith [hsq]
    exact Real.sqrt_pos.mpr harg
  dsimp [jarlskog]
  positivity

/-- Diagnostic theorem: quark mixing suppressions are far below lepton sectors. -/
theorem quark_lepton_mixing_gap :
    ckmSin23 = (1 : ℝ) / 24 ∧
    ckmSin13 = (1 : ℝ) / 272 ∧
    pmnsSin12 ^ 2 = (4 : ℝ) / 13 ∧
    pmnsSin23 ^ 2 = (4 : ℝ) / 7 ∧
    pmnsSin13 = (1 : ℝ) / 7 := by
  refine ⟨ckm_structural_values.1, ckm_structural_values.2, ?_, ?_, ?_⟩
  · exact pmns_structural_values.1
  · exact pmns_structural_values.2.1
  · exact pmns_structural_values.2.2

end Gutoe.FlavorMixing
