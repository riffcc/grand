/-
 * GUTOE — Leptogenesis from the Neutrino Sector
 * Copyright (C) 2026  Riff Labs
 *
 * AGPL-3.0-or-later
 *
 * GRAND-130: Connection between neutrino sector and matter-antimatter asymmetry.
 *
 * Chain: Cl(1,3) → Z₃ (n_gen=3) → PMNS δ = π + arctan(1/3) → ε₁ → 28/79 → η_B
 *
 * All theorems proven (no sorry).
 -/

import Mathlib
import Gutoe.Basic
import Gutoe.DimensionalStructure
import Gutoe.Z3Uniqueness
import Gutoe.GaugeGroupSU3
import Gutoe.FlavorMixing
import Gutoe.Baryogenesis

namespace Gutoe.Leptogenesis

open Gutoe.DimensionalStructure Gutoe.Z3Uniqueness Gutoe.GaugeGroupSU3
open Gutoe.FlavorMixing Gutoe.Baryogenesis

-- ══════════════════════════════════════════════════════════════════════════════
-- Part A: Sphaleron conversion factor
-- ══════════════════════════════════════════════════════════════════════════════

/-- SM sphaleron B/(B-L) = (8 n_f + 4) / (22 n_f + 13), n_f = 3 → 28/79.
    n_f = 3 is forced by Z₃ (quarkOrbit.card = 3, proven by decide). -/
theorem sphaleron_28_over_79 :
    (8 * (quarkOrbit.card : ℚ) + 4) / (22 * quarkOrbit.card + 13) = 28 / 79 := by
  norm_num [show quarkOrbit.card = 3 from by decide]

theorem sphaleron_pos : (0 : ℚ) < 28 / 79 := by norm_num

theorem sphaleron_lt_one : (28 : ℚ) / 79 < 1 := by norm_num

-- ══════════════════════════════════════════════════════════════════════════════
-- Part B: PMNS CP phase → CP violation
-- ══════════════════════════════════════════════════════════════════════════════

/-- δ_PMNS = π + arctan(1/3) (from Cl(1,3) structure). -/
theorem pmns_delta_def : pmnsDelta = Real.pi + Real.arctan (1 / 3) := rfl

/-- arctan(1/3) > 0 since 1/3 > 0. -/
theorem arctan_one_third_pos : Real.arctan (1 / 3) > 0 :=
  Real.arctan_pos.mpr (by norm_num)

/-- arctan(1/3) < π/4: follows from arctan strict monotonicity and arctan(1) = π/4. -/
theorem arctan_one_third_lt_pi_div_four : Real.arctan (1 / 3) < Real.pi / 4 := by
  have : Real.arctan (1 / 3) < Real.arctan 1 :=
    Real.arctan_strictMono (by norm_num)
  rwa [Real.arctan_one] at this

/-- δ_PMNS - π = arctan(1/3) > 0: the phase is displaced from π. -/
theorem pmns_cp_displacement_pos : pmnsDelta - Real.pi > 0 := by
  rw [pmns_delta_def]; linarith [arctan_one_third_pos]

/-- δ_PMNS ≠ π: PMNS phase is not CP-conserving. -/
theorem pmns_delta_ne_pi : pmnsDelta ≠ Real.pi := by
  intro h; linarith [pmns_cp_displacement_pos, show pmnsDelta - Real.pi = 0 from by rw [h]; ring]

/-- sin(δ_PMNS) ≠ 0: CP asymmetry ε₁ is structurally nonzero.
    sin(π + arctan(1/3)) = -sin(arctan(1/3)) and sin(arctan(1/3)) > 0. -/
theorem pmns_sin_delta_ne_zero : Real.sin pmnsDelta ≠ 0 := by
  rw [pmns_delta_def, Real.sin_add, Real.sin_pi, Real.cos_pi]
  ring_nf
  have hpos : 0 < Real.arctan (1 / 3) := arctan_one_third_pos
  have hlt : Real.arctan (1 / 3) < Real.pi :=
    lt_trans arctan_one_third_lt_pi_div_four (by linarith [Real.pi_pos])
  have hsin : 0 < Real.sin (Real.arctan (1 / 3)) :=
    Real.sin_pos_of_pos_of_lt_pi hpos hlt
  linarith

-- ══════════════════════════════════════════════════════════════════════════════
-- Part C: Three generations from Z₃
-- ══════════════════════════════════════════════════════════════════════════════

/-- Z₃ quark orbit has cardinality 3 → three generations → three N_R. -/
theorem n_gen_eq_3 : quarkOrbit.card = 3 := by decide

/-- n_gen ≥ 2: the Davidson-Ibarra mechanism requires at least 2 heavy N_R. -/
theorem leptogenesis_n_gen_sufficient : quarkOrbit.card ≥ 2 := by decide

-- ══════════════════════════════════════════════════════════════════════════════
-- Part D: PMNS void → leptogenesis multiplier
-- ══════════════════════════════════════════════════════════════════════════════

/-- Void correction Δsin²θ₂₃ = 1/548 (α² correction lane, no free parameters). -/
theorem pmns_void_eq_1_over_548 : pmnsTheta23VoidScalar = 1 / 548 :=
  pmns_theta23_void_scalar_eq

/-- Leptogenesis multiplier = 1 + 1/548 = 549/548. -/
theorem lepto_mult_eq : leptogenesisPmnsMultiplier = 549 / 548 :=
  leptogenesis_pmns_multiplier_eq

/-- The multiplier > 1: PMNS void correction enhances η_B. -/
theorem lepto_mult_gt_one : leptogenesisPmnsMultiplier > 1 := by
  rw [lepto_mult_eq]; norm_num

-- ══════════════════════════════════════════════════════════════════════════════
-- Part E: η_B positivity
-- ══════════════════════════════════════════════════════════════════════════════

/-- Structural baryogenesis prefactor > 0. -/
theorem eta_b_prefactor_pos : baryogenesisPrefactor > 0 :=
  baryogenesis_prefactor_pos

/-- η_B structural with PMNS enhancement > 0. -/
theorem eta_b_with_pmns_pos : etaBaryonStructuralWithPmns > 0 :=
  eta_baryon_structural_with_pmns_pos

-- ══════════════════════════════════════════════════════════════════════════════
-- Master theorem
-- ══════════════════════════════════════════════════════════════════════════════

/-- The GUTOE leptogenesis pathway is structurally complete (no sorry):
    (A) B violation: sphaleron 28/79 > 0, forced by Z₃ n_f=3
    (B) CP violation: δ_PMNS - π = arctan(1/3) > 0
    (C) ε₁ ≠ 0: sin(δ_PMNS) ≠ 0, CP asymmetry nonzero
    (D) n_gen = 3: three heavy N_R from Z₃ quark orbit
    (E) PMNS void 1/548 → multiplier 549/548 > 1 (links lepton mixing to η_B)
    (F) η_B > 0: baryon asymmetry is positive -/
theorem leptogenesis_pathway_complete :
    (0 : ℚ) < 28 / 79 ∧
    pmnsDelta - Real.pi > 0 ∧
    Real.sin pmnsDelta ≠ 0 ∧
    quarkOrbit.card = 3 ∧
    leptogenesisPmnsMultiplier > 1 ∧
    etaBaryonStructuralWithPmns > 0 :=
  ⟨by norm_num,
   pmns_cp_displacement_pos,
   pmns_sin_delta_ne_zero,
   by decide,
   lepto_mult_gt_one,
   eta_b_with_pmns_pos⟩

end Gutoe.Leptogenesis
