/-
 * GUTOE — Falsifiability Catalog (GRAND-122/123/124)
 * Copyright (C) 2026  Riff Labs
 * AGPL-3.0-or-later
 *
 * This module formalizes three planning artifacts as Lean propositions:
 *   - GRAND-124: falsifiable core prediction gates
 *   - GRAND-122: null-result consistency gates
 *   - GRAND-123: Standard Model + GR limit recovery bundle
 *
 * No `sorry`.
-/

import Mathlib
import Gutoe.Z3Uniqueness
import Gutoe.GaugeConstants
import Gutoe.FineStructure
import Gutoe.GravityMetric
import Gutoe.KerrGeometry
import Gutoe.LorentzInvariance

namespace Gutoe.FalsifiabilityCatalog

open Gutoe.Z3Uniqueness
open Gutoe.GaugeConstants
open Gutoe.FineStructure
open Gutoe.GravityMetric
open Gutoe.KerrGeometry
open Gutoe.LorentzInvariance

/-! ### GRAND-124: Explicit falsifiable core predictions -/

/-- Minimal observation bundle for core algebraic targets. -/
structure CoreObservation where
  sin2ThetaW : ℚ
  gaugeBosons : ℕ
  alphaInv : ℕ

/-- Core GUTOE prediction gate: all three structural observables must match. -/
def corePredictionGate (obs : CoreObservation) : Prop :=
  obs.sin2ThetaW = (magneticTriplet.card : ℚ) / (2 ^ 4 - magneticTriplet.card : ℚ) ∧
  obs.gaugeBosons = 12 ∧
  obs.alphaInv = 137

/-- The algebraic core target point is internally consistent. -/
theorem core_prediction_target_is_consistent :
    corePredictionGate
      { sin2ThetaW := (3 / 13 : ℚ), gaugeBosons := 12, alphaInv := 137 } := by
  constructor
  · exact weinberg_from_z3_orbits.symm
  constructor
  · exact total_gauge_bosons
  · exact alpha_inverse_d4

/-- Any mismatch in `sin²θ_W` falsifies the core gate. -/
theorem core_falsified_of_sin2_mismatch
    (obs : CoreObservation)
    (h : obs.sin2ThetaW ≠ (magneticTriplet.card : ℚ) / (2 ^ 4 - magneticTriplet.card : ℚ)) :
    ¬ corePredictionGate obs := by
  intro hgate
  exact h hgate.1

/-- Any mismatch in total gauge boson count falsifies the core gate. -/
theorem core_falsified_of_gauge_count_mismatch
    (obs : CoreObservation)
    (h : obs.gaugeBosons ≠ 12) :
    ¬ corePredictionGate obs := by
  intro hgate
  exact h hgate.2.1

/-- Any mismatch in the leading-order alpha inverse target falsifies the core gate. -/
theorem core_falsified_of_alpha_inverse_mismatch
    (obs : CoreObservation)
    (h : obs.alphaInv ≠ 137) :
    ¬ corePredictionGate obs := by
  intro hgate
  exact h hgate.2.2

/-! ### GRAND-122: Null-result consistency gates -/

/-! ### GRAND-125: Strong-CP structural gate -/

/-- Structural θ_QCD proxy from Cl(1,3) bivector balance.
    In this model, CP-odd phase support tracks the rotation/boost asymmetry.
    Exact Lorentz-sector balance (3 rotations, 3 boosts) forces this to zero. -/
def thetaQcdStructural : ℚ :=
  (magneticTriplet.card : ℚ) - (emTriplet.card : ℚ)

/-- Cl(1,3) Lorentz decomposition forces `thetaQcdStructural = 0`. -/
theorem theta_qcd_structural_zero : thetaQcdStructural = 0 := by
  unfold thetaQcdStructural
  rcases lorentz_algebra_decomposition with ⟨_, _, hmag, hem⟩
  rw [hmag, hem]
  norm_num

/-- Runtime bridge constant for neutron EDM estimate:
    `|d_n| ≈ c * |θ_QCD|`, with `c = 2.4e-16 e*cm`. -/
def neutronEdmBridgeCoeff : ℝ := 2.4e-16

/-- Structural neutron EDM estimate from structural θ_QCD. -/
def neutronEdmStructural : ℝ :=
  neutronEdmBridgeCoeff * |(thetaQcdStructural : ℝ)|

/-- Structural θ=0 implies structural neutron EDM estimate is exactly zero. -/
theorem neutron_edm_structural_zero : neutronEdmStructural = 0 := by
  unfold neutronEdmStructural
  rw [theta_qcd_structural_zero]
  norm_num [neutronEdmBridgeCoeff]

/-- Structural neutron EDM gate passes the catalog bound `1e-26 e*cm`. -/
theorem neutron_edm_structural_within_bound :
    neutronEdmStructural ≤ (1e-26 : ℝ) := by
  rw [neutron_edm_structural_zero]
  norm_num

/-- Null-result bounds to be respected by any viable low-energy effective model. -/
structure NullBounds where
  minProtonLifetime : ℝ
  maxEDM : ℝ
  maxFifthForce : ℝ
  maxLorentzViolation : ℝ

/-- Measured observables corresponding to the null-result constraints. -/
structure NullObservation where
  protonLifetime : ℝ
  edmMagnitude : ℝ
  fifthForceStrength : ℝ
  lorentzViolationScale : ℝ

/-- Null-result consistency gate. -/
def respectsNullResults (b : NullBounds) (o : NullObservation) : Prop :=
  b.minProtonLifetime ≤ o.protonLifetime ∧
  o.edmMagnitude ≤ b.maxEDM ∧
  o.fifthForceStrength ≤ b.maxFifthForce ∧
  o.lorentzViolationScale ≤ b.maxLorentzViolation

/-- Proton lifetime below bound violates the null-result gate. -/
theorem null_falsified_of_short_proton_lifetime
    (b : NullBounds) (o : NullObservation)
    (h : o.protonLifetime < b.minProtonLifetime) :
    ¬ respectsNullResults b o := by
  intro hgate
  exact not_le_of_gt h hgate.1

/-- EDM above bound violates the null-result gate. -/
theorem null_falsified_of_edm_excess
    (b : NullBounds) (o : NullObservation)
    (h : b.maxEDM < o.edmMagnitude) :
    ¬ respectsNullResults b o := by
  intro hgate
  exact not_le_of_gt h hgate.2.1

/-- Fifth-force excess violates the null-result gate. -/
theorem null_falsified_of_fifth_force_excess
    (b : NullBounds) (o : NullObservation)
    (h : b.maxFifthForce < o.fifthForceStrength) :
    ¬ respectsNullResults b o := by
  intro hgate
  exact not_le_of_gt h hgate.2.2.1

/-- Lorentz-violation excess violates the null-result gate. -/
theorem null_falsified_of_lorentz_violation_excess
    (b : NullBounds) (o : NullObservation)
    (h : b.maxLorentzViolation < o.lorentzViolationScale) :
    ¬ respectsNullResults b o := by
  intro hgate
  exact not_le_of_gt h hgate.2.2.2

/-! ### GRAND-123: Standard Model + GR limit recovery bundle -/

/-- Formal bundle of SM and GR limit checks used by roadmap verification. -/
def smGrLimitBundle : Prop :=
  (1 - (magneticTriplet.card : ℚ) / (2 ^ 4 - magneticTriplet.card : ℚ) = 10 / 13) ∧
  ((3 ^ 2 - 1) + (2 ^ 2 - 1) + 1 = 12) ∧
  (thetaQcdStructural = 0) ∧
  (∀ r : ℝ, Real.sqrt (r ^ 2 + (r_core 0) ^ 2) = |r|) ∧
  (∀ r_s : ℝ, 0 ≤ r_s → rPlus r_s 0 = r_s ∧ rMinus r_s 0 = 0)

/-- Existing theorem chain implies the SM+GR limit bundle. -/
theorem sm_gr_limits_recovered : smGrLimitBundle := by
  constructor
  · exact cos_sq_theta_w_from_z3
  constructor
  · exact total_gauge_bosons
  constructor
  · exact theta_qcd_structural_zero
  constructor
  · intro r
    exact r_eff_classical_limit r
  · intro r_s hrs
    exact schwarzschild_limit_horizons hrs

end Gutoe.FalsifiabilityCatalog
