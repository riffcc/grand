import Mathlib
import Gutoe.FineStructure
import Gutoe.GaugeConstants
import Gutoe.GaugeGroupSM
import Gutoe.Chirality

namespace Gutoe.Homochirality

open Gutoe.FineStructure
open Gutoe.GaugeConstants
open Gutoe.GaugeGroupSM
open Gutoe.Chirality
open Gutoe.Z3Uniqueness

/-- Weak-sector share of the SM gauge algebra:
    `|SU(2)| / (8+3+1) = 3/12 = 1/4`. -/
def weakGaugeFractionQ : ℚ :=
  (magneticTriplet.card : ℚ) / (((3 ^ 2 - 1) + (2 ^ 2 - 1) + 1 : ℕ) : ℚ)

theorem weak_gauge_fraction_eq_one_quarter :
    weakGaugeFractionQ = 1 / 4 := by
  unfold weakGaugeFractionQ
  rw [su2_dim, total_gauge_bosons]
  norm_num

/-- Chiral projection from finite Cl(1,3) state count:
    `1 / 2^4 = 1/16`. -/
def chiralProjectionQ : ℚ := 1 / ((2 ^ 4 : ℕ) : ℚ)

theorem chiral_projection_eq_one_over_16 :
    chiralProjectionQ = 1 / 16 := by
  unfold chiralProjectionQ
  norm_num

/-- Weak nuclear charge lane:
    `Q_W = N - (1 - 4 sin²θ_W) Z`, with `sin²θ_W` from the shared Z₃ theorem. -/
def weakChargeQ (z n : ℚ) : ℚ :=
  n - (1 - 4 * ((magneticTriplet.card : ℚ) / (2 ^ 4 - magneticTriplet.card : ℚ))) * z

theorem weak_charge_q_simplifies (z n : ℚ) :
    weakChargeQ z n = n - z / 13 := by
  unfold weakChargeQ
  rw [weinberg_from_z3_orbits]
  ring

/-- Weak charge for ¹⁴N: `Q_W = 84/13`. -/
def nitrogenWeakChargeQ : ℚ := weakChargeQ 7 7

theorem nitrogen_weak_charge_eq :
    nitrogenWeakChargeQ = 84 / 13 := by
  unfold nitrogenWeakChargeQ
  rw [weak_charge_q_simplifies]
  norm_num

/-- Weak charge for ¹⁶O: `Q_W = 96/13`. -/
def oxygenWeakChargeQ : ℚ := weakChargeQ 8 8

theorem oxygen_weak_charge_eq :
    oxygenWeakChargeQ = 96 / 13 := by
  unfold oxygenWeakChargeQ
  rw [weak_charge_q_simplifies]
  norm_num

/-- Canonical amino-acid backbone weak/chiral source factor:
    one N + two O contributions weighted by `Z^3 Q_W`. -/
def aminoBackboneNuclearFactorQ : ℚ :=
  (7 : ℚ) ^ 3 * nitrogenWeakChargeQ + (2 : ℚ) * (8 : ℚ) ^ 3 * oxygenWeakChargeQ

theorem amino_backbone_nuclear_factor_eq :
    aminoBackboneNuclearFactorQ = 127116 / 13 := by
  unfold aminoBackboneNuclearFactorQ
  rw [nitrogen_weak_charge_eq, oxygen_weak_charge_eq]
  norm_num

/-- Structural parity factor before electromagnetic suppression. -/
def aminoBackboneParityFactorQ : ℚ :=
  weakGaugeFractionQ * chiralProjectionQ * aminoBackboneNuclearFactorQ

theorem amino_backbone_parity_factor_eq :
    aminoBackboneParityFactorQ = 31779 / 208 := by
  unfold aminoBackboneParityFactorQ
  rw [weak_gauge_fraction_eq_one_quarter, chiral_projection_eq_one_over_16,
    amino_backbone_nuclear_factor_eq]
  norm_num

/-- Leading-order structural fine-structure constant in rational form. -/
def alphaLeadingQ : ℚ := 1 / (alphaInverse 4 : ℚ)

theorem alpha_leading_eq_one_over_137 :
    alphaLeadingQ = 1 / 137 := by
  unfold alphaLeadingQ
  rw [alpha_inverse_d4]
  norm_num

/-- Alpha-suppressed backbone parity proxy used by the runtime lane:
    `proxy = parity_factor * α^4`. -/
def aminoBackboneParityProxyQ : ℚ :=
  aminoBackboneParityFactorQ * alphaLeadingQ ^ 4

theorem amino_backbone_parity_proxy_closed_form :
    aminoBackboneParityProxyQ = 31779 / 73273275088 := by
  unfold aminoBackboneParityProxyQ
  rw [amino_backbone_parity_factor_eq, alpha_leading_eq_one_over_137]
  norm_num

theorem amino_backbone_parity_proxy_positive :
    0 < aminoBackboneParityProxyQ := by
  rw [amino_backbone_parity_proxy_closed_form]
  norm_num

/-- Sign lane imported from the proved Cl(1,3) chirality theorem:
    SU(2)-quark coupling parity is negative (`-1`). -/
def aminoHandednessSign : ℤ :=
  bivectorParity13 ⟨1, by decide⟩ ⟨2, by decide⟩ * metricParity13 ⟨1, by decide⟩

theorem amino_handedness_sign_eq_neg_one :
    aminoHandednessSign = -1 := by
  unfold aminoHandednessSign
  exact su2_quark_coupling_parity

theorem amino_handedness_sign_negative :
    aminoHandednessSign < 0 := by
  rw [amino_handedness_sign_eq_neg_one]
  norm_num

end Gutoe.Homochirality

