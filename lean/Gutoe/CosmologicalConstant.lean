/-
 * GUTOE — Cosmological Constant Structural Suppression from Cl(1,3)
 * Copyright (C) 2026 Riff Labs
 * AGPL-3.0-or-later
 *
 * GRAND-92 slice:
 *   Λ_struct = λ_H^(α⁻¹_LO) / l_P^2
 * with λ_H = 13/100 and α⁻¹_LO = 137 from the shared proof chain.
 -/

import Mathlib
import Gutoe.EWSBHiggs
import Gutoe.FineStructure
import Gutoe.SignatureUniqueness
import Gutoe.Z3Uniqueness
import Gutoe.Chirality

namespace Gutoe.CosmologicalConstant

open Gutoe.EWSBHiggs
open Gutoe.FineStructure
open Gutoe.Z3Uniqueness
open Gutoe.Chirality

/-- Structural vacuum-energy suppression factor:
    s_Λ = λ_H^(α⁻¹_LO). -/
def lambdaSuppression : ℚ := higgsQuartic ^ (alphaInverse 4)

/-- Exact structural suppression value:
    s_Λ = (13/100)^137. -/
theorem lambda_suppression_eq_13_100_pow_137 :
    lambdaSuppression = ((13 : ℚ) / 100) ^ (137 : ℕ) := by
  unfold lambdaSuppression
  rw [higgs_quartic_eq_13_100, alpha_inverse_d4]

/-- Structural suppression is strictly positive. -/
theorem lambda_suppression_pos : 0 < lambdaSuppression := by
  rw [lambda_suppression_eq_13_100_pow_137]
  positivity

/-- Structural suppression is below unity. -/
theorem lambda_suppression_lt_one : lambdaSuppression < 1 := by
  rw [lambda_suppression_eq_13_100_pow_137]
  norm_num

/-- Canonical generator labels in `Fin 4`. -/
def f0 : Fin 4 := ⟨0, by decide⟩
def f1 : Fin 4 := ⟨1, by decide⟩
def f2 : Fin 4 := ⟨2, by decide⟩
def f3 : Fin 4 := ⟨3, by decide⟩

/-- The six independent grade-2 generator pairs in Cl(1,3). -/
def bivectorPairs13 : Finset (Fin 4 × Fin 4) :=
  {(f0, f1), (f0, f2), (f0, f3), (f1, f2), (f1, f3), (f2, f3)}

/-- Signature parity (`±1`) of each bivector pair. -/
def bivectorPairSignature (p : Fin 4 × Fin 4) : ℤ :=
  bivectorParity13 p.1 p.2

/-- Timelike-spacelike bivectors: signature `-1`. -/
def timelikeSpacelikePairs : Finset (Fin 4 × Fin 4) :=
  bivectorPairs13.filter (fun p => bivectorPairSignature p = -1)

/-- Spacelike-spacelike bivectors: signature `+1`. -/
def spacelikeSpacelikePairs : Finset (Fin 4 × Fin 4) :=
  bivectorPairs13.filter (fun p => bivectorPairSignature p = 1)

/-- Cl(1,3) bivector signature split: 3 mixed-sign and 3 same-sign pairs. -/
theorem bivector_signature_split_3_3 :
    timelikeSpacelikePairs.card = 3 ∧ spacelikeSpacelikePairs.card = 3 := by
  native_decide

/-- Lorentzian bivector split normalization candidate:
    k_sig = sqrt(|grade-2| / |temporal-bivector-orbit|). -/
noncomputable def lorentzSignatureNormalization : ℝ :=
  Real.sqrt ((grade2_4d.card : ℝ) / (emTriplet.card : ℝ))

/-- Exact structural value of the normalization candidate:
    sqrt(6/3) = sqrt(2). -/
theorem lorentz_signature_normalization_eq_sqrt2 :
    lorentzSignatureNormalization = Real.sqrt 2 := by
  unfold lorentzSignatureNormalization
  have hg2 : grade2_4d.card = 6 := by native_decide
  have hem : emTriplet.card = 3 := by native_decide
  rw [hg2, hem]
  norm_num
  rfl

/-- Equivalent signature normalization from explicit `(1,3)` bivector parity split:
    sqrt(total bivectors / timelike-spacelike bivectors) = sqrt(6/3) = sqrt(2). -/
noncomputable def lorentzSignatureNormalizationFromParity : ℝ :=
  Real.sqrt ((bivectorPairs13.card : ℝ) / (timelikeSpacelikePairs.card : ℝ))

/-- Exact value from the explicit parity split. -/
theorem lorentz_signature_normalization_from_parity_eq_sqrt2 :
    lorentzSignatureNormalizationFromParity = Real.sqrt 2 := by
  unfold lorentzSignatureNormalizationFromParity
  have hcard : bivectorPairs13.card = 6 := by native_decide
  have hsplit : timelikeSpacelikePairs.card = 3 := (bivector_signature_split_3_3).1
  rw [hcard, hsplit]
  norm_num
  rfl

/-- The two normalization views are equal:
    (orbit-count split) = (signature-parity split). -/
theorem lorentz_signature_normalization_eq_from_parity :
    lorentzSignatureNormalization = lorentzSignatureNormalizationFromParity := by
  rw [lorentz_signature_normalization_eq_sqrt2, lorentz_signature_normalization_from_parity_eq_sqrt2]

/-- Cosmological constant candidate from Planck curvature scaling:
    Λ_struct(l_P) = s_Λ / l_P². -/
noncomputable def lambdaCosmologicalFromPlanck (lP : ℝ) : ℝ :=
  ((lambdaSuppression : ℚ) : ℝ) / (lP ^ 2)

/-- GRAND-293 candidate normalization from Lorentzian signature:
    Λ_sig(l_P) = Λ_struct(l_P) / √2.

    This introduces no new continuous parameter: √2 is fixed and linked to
    the Minkowski-signature branch selected in `SignatureUniqueness`.
    (Conjectural bridge until the full bivector-normalization derivation closes.) -/
noncomputable def lambdaCosmologicalSignatureCandidate (lP : ℝ) : ℝ :=
  lambdaCosmologicalFromPlanck lP / Real.sqrt 2

/-- Same candidate written through the Lorentzian split normalization. -/
noncomputable def lambdaCosmologicalSignatureFromSplit (lP : ℝ) : ℝ :=
  lambdaCosmologicalFromPlanck lP / lorentzSignatureNormalization

/-- Exact real-form structural cosmological candidate:
    Λ_struct(l_P) = ((13/100)^137) / l_P². -/
theorem lambda_cosmological_from_planck_eq
    (lP : ℝ) :
    lambdaCosmologicalFromPlanck lP =
      ((((13 : ℚ) / 100) ^ (137 : ℕ) : ℚ) : ℝ) / (lP ^ 2) := by
  unfold lambdaCosmologicalFromPlanck
  rw [lambda_suppression_eq_13_100_pow_137]

/-- Signature-corrected candidate is exactly the structural term divided by √2. -/
theorem lambda_cosmological_signature_candidate_eq
    (lP : ℝ) :
    lambdaCosmologicalSignatureCandidate lP =
      ((((13 : ℚ) / 100) ^ (137 : ℕ) : ℚ) : ℝ) / (Real.sqrt 2 * (lP ^ 2)) := by
  unfold lambdaCosmologicalSignatureCandidate
  rw [lambda_cosmological_from_planck_eq]
  have hs2 : (Real.sqrt 2) ≠ 0 := by
    have hs2pos : 0 < Real.sqrt 2 := by positivity
    exact ne_of_gt hs2pos
  field_simp [hs2]

/-- The split-normalized and sqrt(2)-normalized candidates are definitionally equal. -/
theorem lambda_signature_from_split_eq_candidate (lP : ℝ) :
    lambdaCosmologicalSignatureFromSplit lP = lambdaCosmologicalSignatureCandidate lP := by
  unfold lambdaCosmologicalSignatureFromSplit lambdaCosmologicalSignatureCandidate
  rw [lorentz_signature_normalization_eq_sqrt2]

/-- For nonzero Planck length, the structural Λ candidate is positive. -/
theorem lambda_cosmological_from_planck_pos
    {lP : ℝ}
    (hlP : lP ≠ 0) :
    0 < lambdaCosmologicalFromPlanck lP := by
  unfold lambdaCosmologicalFromPlanck
  have hsuppQ : 0 < lambdaSuppression := lambda_suppression_pos
  have hsuppR : 0 < ((lambdaSuppression : ℚ) : ℝ) := by
    exact_mod_cast hsuppQ
  have hden : 0 < lP ^ 2 := by
    nlinarith [sq_pos_of_ne_zero hlP]
  exact div_pos hsuppR hden

/-- Positive signature-corrected candidate for nonzero Planck length. -/
theorem lambda_cosmological_signature_candidate_pos
    {lP : ℝ}
    (hlP : lP ≠ 0) :
    0 < lambdaCosmologicalSignatureCandidate lP := by
  unfold lambdaCosmologicalSignatureCandidate
  have hbase : 0 < lambdaCosmologicalFromPlanck lP := lambda_cosmological_from_planck_pos hlP
  have hs2 : 0 < Real.sqrt 2 := by
    exact Real.sqrt_pos.2 (by norm_num)
  exact div_pos hbase hs2

end Gutoe.CosmologicalConstant
