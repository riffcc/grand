/-
 * GUTOE — Chiral Symmetry Breaking Structural Gate (GRAND-126)
 *
 * Structural lane:
 *   - nonzero quark condensate proxy from shared Cl(1,3) primitives
 *   - pseudo-Goldstone pion scaling
 *   - confinement-linked positive witness
 *
 * No `sorry`.
-/

import Mathlib
import Gutoe.FineStructure
import Gutoe.GravityMetric
import Gutoe.GaugeGroupSU3
import Gutoe.AsymptoticFreedomEntropy

namespace Gutoe.ChiralSymmetryBreaking

open Gutoe.FineStructure
open Gutoe.GravityMetric
open Gutoe.GaugeGroupSU3
open Gutoe.AsymptoticFreedomEntropy

/-- Leading-order structural alpha from the shared `α⁻¹ = 137` theorem. -/
noncomputable def alphaLeadingOrder : ℝ :=
  1 / (alphaInverse 4 : ℝ)

/-- Structural quark condensate proxy:
    `-(1 - λ_QG) * (|quarkOrbit| / dim Cl(1,3))`. -/
noncomputable def quarkCondensateProxy : ℝ :=
  -((1 - lambda_qg) * ((quarkOrbit.card : ℝ) / (((2 ^ 4 : ℕ) : ℝ))))

/-- Pion mass-squared proxy for the pseudo-Goldstone lane:
    `m_π² ∝ α * (-⟨qq⟩)`. -/
noncomputable def pionMassSqProxy : ℝ :=
  alphaLeadingOrder * (-quarkCondensateProxy)

/-- Pseudo-Goldstone ratio map:
    `m_π² / (-⟨qq⟩)` should recover the explicit-breaking scale. -/
noncomputable def pseudoGoldstoneRatio : ℝ :=
  pionMassSqProxy / (-quarkCondensateProxy)

/-- Chiral-limit map for explicit breaking parameter `ε`. -/
noncomputable def pionMassSqFromExplicitBreaking (ε : ℝ) : ℝ :=
  ε * (-quarkCondensateProxy)

/-- Positive witness linking confinement and chiral breaking:
    `β₀ * (-⟨qq⟩)`. -/
noncomputable def confinementChiralLinkStrength : ℝ :=
  beta0Clifford * (-quarkCondensateProxy)

theorem alpha_leading_order_eq :
    alphaLeadingOrder = (1 : ℝ) / 137 := by
  unfold alphaLeadingOrder
  norm_num [alpha_inverse_d4]

theorem quark_condensate_proxy_closed_form :
    quarkCondensateProxy = -(11 : ℝ) / 64 := by
  unfold quarkCondensateProxy
  norm_num [lambda_qg, quarkOrbit_card]

theorem quark_condensate_proxy_negative :
    quarkCondensateProxy < 0 := by
  rw [quark_condensate_proxy_closed_form]
  norm_num

theorem quark_condensate_proxy_nonzero :
    quarkCondensateProxy ≠ 0 := by
  exact ne_of_lt quark_condensate_proxy_negative

theorem pion_mass_sq_proxy_closed_form :
    pionMassSqProxy = (11 : ℝ) / 8768 := by
  unfold pionMassSqProxy
  rw [alpha_leading_order_eq, quark_condensate_proxy_closed_form]
  norm_num

theorem pion_mass_sq_proxy_positive :
    0 < pionMassSqProxy := by
  rw [pion_mass_sq_proxy_closed_form]
  positivity

theorem pseudo_goldstone_ratio_eq_alpha :
    pseudoGoldstoneRatio = alphaLeadingOrder := by
  unfold pseudoGoldstoneRatio pionMassSqProxy
  have hq : quarkCondensateProxy ≠ 0 := quark_condensate_proxy_nonzero
  field_simp [hq]

theorem pseudo_goldstone_ratio_closed_form :
    pseudoGoldstoneRatio = (1 : ℝ) / 137 := by
  rw [pseudo_goldstone_ratio_eq_alpha, alpha_leading_order_eq]

theorem pion_mass_sq_chiral_limit_zero :
    pionMassSqFromExplicitBreaking 0 = 0 := by
  unfold pionMassSqFromExplicitBreaking
  ring

theorem pion_mass_sq_from_positive_breaking_positive
    {ε : ℝ} (hε : 0 < ε) :
    0 < pionMassSqFromExplicitBreaking ε := by
  unfold pionMassSqFromExplicitBreaking
  have hcond : 0 < -quarkCondensateProxy := by
    nlinarith [quark_condensate_proxy_negative]
  exact mul_pos hε hcond

theorem confinement_chiral_link_closed_form :
    confinementChiralLinkStrength = (319 : ℝ) / 96 := by
  unfold confinementChiralLinkStrength
  rw [quark_condensate_proxy_closed_form]
  norm_num [beta0Clifford]

theorem confinement_chiral_link_positive :
    0 < confinementChiralLinkStrength := by
  rw [confinement_chiral_link_closed_form]
  positivity

/-- GRAND-126 structural gate:
    nonzero condensate, pseudo-Goldstone scaling, chiral-limit vanishing,
    and positive confinement-linked witness. -/
theorem chiral_symmetry_breaking_gate :
    quarkCondensateProxy < 0 ∧
    0 < pionMassSqProxy ∧
    pseudoGoldstoneRatio = (1 : ℝ) / 137 ∧
    pionMassSqFromExplicitBreaking 0 = 0 ∧
    0 < confinementChiralLinkStrength := by
  exact ⟨quark_condensate_proxy_negative,
         pion_mass_sq_proxy_positive,
         pseudo_goldstone_ratio_closed_form,
         pion_mass_sq_chiral_limit_zero,
         confinement_chiral_link_positive⟩

end Gutoe.ChiralSymmetryBreaking
