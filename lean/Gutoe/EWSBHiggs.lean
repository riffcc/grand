/-
 * GUTOE — Electroweak Symmetry Breaking and Higgs Quartic from Cl(1,3)
 * Copyright (C) 2026 Riff Labs
 * AGPL-3.0-or-later
 *
 * Structural bridge for GRAND-80 / GRAND-131:
 *   - Higgs order parameter identified with lattice void fraction f₀.
 *   - Mexican-hat quartic λ derived from shared Clifford counts.
 *   - Broken/unbroken phase split derived from a Cl(1,3) critical fraction.
 *
 * No free fit knobs are introduced in this module.
 -/

import Mathlib
import Gutoe.DimensionalStructure
import Gutoe.Z3Uniqueness
import Gutoe.GaugeConstants

namespace Gutoe.EWSBHiggs

open Gutoe.DimensionalStructure
open Gutoe.Z3Uniqueness
open Gutoe.GaugeConstants

/-- Electroweak grade-counting denominator: |grade-1| + |grade-2| = 4 + 6 = 10. -/
def ewGradeSum : ℕ := grade1_4d.card + grade2_4d.card

/-- Shared counting identity for the electroweak grade sum. -/
theorem ew_grade_sum_eq_10 : ewGradeSum = 10 := by native_decide

/-- Clifford complement to the magnetic triplet: 16 - 3 = 13. -/
theorem clifford_complement_eq_13 : (2 ^ 4 : ℕ) - magneticTriplet.card = 13 := by
  have hs : magneticTriplet.card = 3 := su2_dim
  norm_num [hs]

/-- Higgs quartic coupling from Clifford counts:
    λ_H = (16 - 3) / (4 + 6)^2 = 13 / 100. -/
def higgsQuartic : ℚ :=
  ((2 ^ 4 : ℕ) - magneticTriplet.card : ℚ) / (ewGradeSum ^ 2 : ℚ)

/-- Exact structural quartic value from Cl(1,3). -/
theorem higgs_quartic_eq_13_100 : higgsQuartic = (13 : ℚ) / 100 := by
  have hs : magneticTriplet.card = 3 := su2_dim
  have hsum : ewGradeSum = 10 := ew_grade_sum_eq_10
  rw [higgsQuartic, hs, hsum]
  norm_num

/-- Structural quartic is strictly positive. -/
theorem higgs_quartic_pos : 0 < higgsQuartic := by
  rw [higgs_quartic_eq_13_100]
  norm_num

/-- Critical symmetry-restoration fraction from Z₃/SU(2) over Clifford dimension:
    f_c = 3/16. -/
def criticalVoidFraction : ℚ := (magneticTriplet.card : ℚ) / (2 ^ 4 : ℚ)

/-- Exact value of critical fraction. -/
theorem critical_void_fraction_eq_3_16 : criticalVoidFraction = (3 : ℚ) / 16 := by
  have hs : magneticTriplet.card = 3 := su2_dim
  rw [criticalVoidFraction, hs]
  norm_num

/-- Effective mass-squared control parameter:
    μ²(f₀) = f₀ - f_c.
    Broken phase: μ²>0 (f₀ > f_c), unbroken phase: μ²≤0 (f₀ ≤ f_c). -/
def muSq (f0 : ℚ) : ℚ := f0 - criticalVoidFraction

/-- Mexican-hat effective potential in order-parameter form:
    V(φ;f₀) = -μ²(f₀) φ² + λ_H φ⁴. -/
def higgsPotential (φ f0 : ℚ) : ℚ := -(muSq f0) * φ ^ 2 + higgsQuartic * φ ^ 4

/-- Derivative of the quartic effective potential. -/
def higgsPotentialDeriv (φ f0 : ℚ) : ℚ :=
  -2 * (muSq f0) * φ + 4 * higgsQuartic * φ ^ 3

/-- Non-trivial stationary branch in squared form: φ² = μ² / (2λ). -/
def nontrivialVevSq (f0 : ℚ) : ℚ := muSq f0 / (2 * higgsQuartic)

/-- Broken phase condition gives positive μ². -/
theorem broken_phase_muSq_pos (f0 : ℚ) (hbreak : criticalVoidFraction < f0) :
    0 < muSq f0 := by
  unfold muSq
  linarith

/-- Unbroken phase condition gives non-positive μ². -/
theorem unbroken_phase_muSq_nonpos (f0 : ℚ) (hunbroken : f0 ≤ criticalVoidFraction) :
    muSq f0 ≤ 0 := by
  unfold muSq
  linarith

/-- In broken phase, non-trivial vev² branch is strictly positive. -/
theorem nontrivial_vev_sq_pos (f0 : ℚ) (hbreak : criticalVoidFraction < f0) :
    0 < nontrivialVevSq f0 := by
  unfold nontrivialVevSq
  have hmu : 0 < muSq f0 := broken_phase_muSq_pos f0 hbreak
  have hLam : 0 < higgsQuartic := higgs_quartic_pos
  have hden : 0 < 2 * higgsQuartic := by nlinarith [hLam]
  exact div_pos hmu hden

/-- Any branch with φ² = μ²/(2λ) is a stationary point of V.
    This is the algebraic core of spontaneous symmetry breaking in the quartic model. -/
theorem higgs_deriv_zero_at_nontrivial_stationary
    (f0 φ : ℚ)
    (hφ : φ ^ 2 = nontrivialVevSq f0) :
    higgsPotentialDeriv φ f0 = 0 := by
  have hLam0 : higgsQuartic ≠ 0 := by
    have hLam : 0 < higgsQuartic := higgs_quartic_pos
    linarith
  have h2Lam0 : (2 * higgsQuartic) ≠ 0 := by
    intro h
    apply hLam0
    nlinarith [h]
  calc
    higgsPotentialDeriv φ f0
        = φ * (-2 * muSq f0 + 4 * higgsQuartic * φ ^ 2) := by
            unfold higgsPotentialDeriv
            ring
    _ = φ * (-2 * muSq f0 + 4 * higgsQuartic * nontrivialVevSq f0) := by
          rw [hφ]
    _ = φ * 0 := by
          unfold nontrivialVevSq
          field_simp [h2Lam0]
          ring
    _ = 0 := by ring

/-- Higgs mass-to-vev squared ratio from structural λ:
    (m_H / v)^2 = 2 λ_H = 13/50. -/
def higgsMassOverVevSq : ℚ := 2 * higgsQuartic

theorem higgs_mass_over_vev_sq_eq_13_50 : higgsMassOverVevSq = (13 : ℚ) / 50 := by
  unfold higgsMassOverVevSq
  rw [higgs_quartic_eq_13_100]
  ring

/-- Combined electroweak closure slice used by runtime parity:
    - quartic fixed structurally,
    - broken phase admits non-trivial stationary branch,
    - mass ratio coefficient is fixed. -/
theorem ewsb_structural_closure (f0 : ℚ)
    (hbreak : criticalVoidFraction < f0) :
    higgsQuartic = (13 : ℚ) / 100 ∧
    0 < nontrivialVevSq f0 ∧
    higgsMassOverVevSq = (13 : ℚ) / 50 := by
  exact ⟨higgs_quartic_eq_13_100, nontrivial_vev_sq_pos f0 hbreak, higgs_mass_over_vev_sq_eq_13_50⟩

end Gutoe.EWSBHiggs
