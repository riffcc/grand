import Mathlib
import Gutoe.FineStructure
import Gutoe.GravityMetric
import Gutoe.Z3Uniqueness
import Gutoe.DarkMatterSector
import Gutoe.Baryogenesis
import Gutoe.Inflation

namespace Gutoe.BBN

open Gutoe.FineStructure
open Gutoe.GravityMetric
open Gutoe.Z3Uniqueness
open Gutoe.DarkMatterSector
open Gutoe.Baryogenesis
open Gutoe.Inflation

/-- Primordial He-4 mass-fraction anchor used by the BBN lane. -/
def ypTargetQ : ℚ := 245 / 1000

/-- Primordial D/H anchor used by the BBN lane. -/
def dOverHTargetQ : ℚ := 2547 / 100000000

/-- Primordial ³He/H anchor used by the BBN lane. -/
def he3OverHTargetQ : ℚ := 11 / 1000000

/-- Primordial ⁷Li/H observational anchor used by the BBN lane. -/
def li7OverHObservedQ : ℚ := 16 / 100000000000

/-- Structural `η₁₀` reference from shared inflation/Clifford counts:
    `(12*5)/(4+6) = 6`. -/
noncomputable def eta10Ref : ℝ :=
  ((geometricDarkAmplificationQ : ℝ) * (darkSectorCandidates.card : ℝ)) /
    (((Nat.choose 4 1 : ℕ) : ℝ) + ((Nat.choose 4 2 : ℕ) : ℝ))

/-- Structural deuterium exponent `(6+2)/(4+1) = 8/5`. -/
noncomputable def deuteriumEtaExponent : ℝ :=
  ((((Nat.choose 4 2 : ℕ) : ℝ) + 2) / (((Nat.choose 4 1 : ℕ) : ℝ) + 1))

/-- Structural helium-3 exponent `|SU(2)|/(4+1) = 3/5`. -/
noncomputable def helium3EtaExponent : ℝ :=
  ((magneticTriplet.card : ℝ) / (((Nat.choose 4 1 : ℕ) : ℝ) + 1))

/-- Structural lithium tension amplification `12/4 = 3`. -/
noncomputable def lithium7TensionAmplification : ℝ :=
  (geometricDarkAmplificationQ : ℝ) / ((Nat.choose 4 1 : ℕ) : ℝ)

/-- Structural BBN source from baryogenesis: `η₁₀ = 10¹⁰ η_B`. -/
noncomputable def eta10FromBaryogenesis : ℝ :=
  (10 : ℝ) ^ 10 * etaBaryonStructural

/-- He-4 mass fraction lane used by GRAND-349. -/
noncomputable def primordialHelium4MassFraction (eta10 : ℝ) : ℝ :=
  (ypTargetQ : ℝ) + (lambda_qg / 50) * (eta10 - eta10Ref)

/-- D/H lane used by GRAND-349. -/
noncomputable def primordialDeuteriumRatio (eta10 : ℝ) : ℝ :=
  (dOverHTargetQ : ℝ) * Real.rpow (eta10Ref / eta10) deuteriumEtaExponent

/-- ³He/H lane used by GRAND-349. -/
noncomputable def primordialHelium3Ratio (eta10 : ℝ) : ℝ :=
  (he3OverHTargetQ : ℝ) * Real.rpow (eta10Ref / eta10) helium3EtaExponent

/-- ⁷Li/H lane used by GRAND-349. -/
noncomputable def primordialLithium7Ratio (eta10 : ℝ) : ℝ :=
  (li7OverHObservedQ : ℝ) * (eta10 / eta10Ref) ^ (2 : ℕ) * lithium7TensionAmplification

/-- ⁷Li tension ratio lane (`pred/observed`). -/
noncomputable def lithium7TensionRatio (eta10 : ℝ) : ℝ :=
  primordialLithium7Ratio eta10 / (li7OverHObservedQ : ℝ)

theorem eta10_ref_eq_6 : eta10Ref = 6 := by
  unfold eta10Ref
  have h41 : (Nat.choose 4 1 : ℕ) = 4 := by native_decide
  have h42 : (Nat.choose 4 2 : ℕ) = 6 := by native_decide
  rw [geometric_dark_amplification_eq, dark_sector_card_eq_five]
  rw [h41, h42]
  norm_num

theorem eta10_ref_pos : 0 < eta10Ref := by
  rw [eta10_ref_eq_6]
  norm_num

theorem deuterium_eta_exponent_eq : deuteriumEtaExponent = (8 : ℝ) / 5 := by
  unfold deuteriumEtaExponent
  have h41 : (Nat.choose 4 1 : ℕ) = 4 := by native_decide
  have h42 : (Nat.choose 4 2 : ℕ) = 6 := by native_decide
  rw [h41, h42]
  norm_num

theorem helium3_eta_exponent_eq : helium3EtaExponent = (3 : ℝ) / 5 := by
  unfold helium3EtaExponent
  have h41 : (Nat.choose 4 1 : ℕ) = 4 := by native_decide
  rw [su2_dim]
  rw [h41]
  norm_num

theorem lithium7_tension_amplification_eq :
    lithium7TensionAmplification = 3 := by
  unfold lithium7TensionAmplification
  have h41 : (Nat.choose 4 1 : ℕ) = 4 := by native_decide
  rw [geometric_dark_amplification_eq]
  rw [h41]
  norm_num

theorem eta10_from_baryogenesis_pos : 0 < eta10FromBaryogenesis := by
  unfold eta10FromBaryogenesis
  have hpow : 0 < (10 : ℝ) ^ 10 := by positivity
  exact mul_pos hpow eta_baryon_structural_pos

theorem primordial_helium4_at_reference :
    primordialHelium4MassFraction eta10Ref = (ypTargetQ : ℝ) := by
  unfold primordialHelium4MassFraction
  ring

theorem primordial_deuterium_at_reference :
    primordialDeuteriumRatio eta10Ref = (dOverHTargetQ : ℝ) := by
  unfold primordialDeuteriumRatio
  have hdiv : eta10Ref / eta10Ref = (1 : ℝ) := by
    field_simp [eta10_ref_pos.ne']
  simp [hdiv]

theorem primordial_helium3_at_reference :
    primordialHelium3Ratio eta10Ref = (he3OverHTargetQ : ℝ) := by
  unfold primordialHelium3Ratio
  have hdiv : eta10Ref / eta10Ref = (1 : ℝ) := by
    field_simp [eta10_ref_pos.ne']
  simp [hdiv]

theorem lithium7_ratio_at_reference :
    primordialLithium7Ratio eta10Ref = (li7OverHObservedQ : ℝ) * 3 := by
  unfold primordialLithium7Ratio
  have hdiv : eta10Ref / eta10Ref = (1 : ℝ) := by
    field_simp [eta10_ref_pos.ne']
  rw [hdiv]
  rw [one_pow, lithium7_tension_amplification_eq]
  ring

theorem lithium7_tension_ratio_at_reference :
    lithium7TensionRatio eta10Ref = 3 := by
  unfold lithium7TensionRatio
  rw [lithium7_ratio_at_reference]
  have hobs_nonzero : (li7OverHObservedQ : ℝ) ≠ 0 := by
    norm_num [li7OverHObservedQ]
  field_simp [hobs_nonzero]

theorem lithium7_tension_ratio_reference_window :
    (2 : ℝ) ≤ lithium7TensionRatio eta10Ref ∧ lithium7TensionRatio eta10Ref ≤ 4 := by
  rw [lithium7_tension_ratio_at_reference]
  constructor <;> norm_num

end Gutoe.BBN
