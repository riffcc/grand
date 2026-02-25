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
import Gutoe.StrongCP
import Gutoe.SMQCDUnification

namespace Gutoe.FalsifiabilityCatalog

open Gutoe.Z3Uniqueness
open Gutoe.GaugeConstants
open Gutoe.FineStructure
open Gutoe.GravityMetric
open Gutoe.KerrGeometry
open Gutoe.StrongCP
open Gutoe.SMQCDUnification
open Gutoe.StrongCPEmergence
open Gutoe.StrongCPVacuum

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

/-! ### SM × QCD acceptance gate (closure integration) -/

/-- Top-level acceptance gate joining SM+QCD unification and SM+GR limits. -/
def smQcdAcceptanceGate
    {X : Type}
    [TopologicalSpace X] [PreconnectedSpace X] [Nonempty X]
    (x0 : X)
    (qClass : HomotopyClass X Su3Matrix → ℤ) : Prop :=
  smQcdGeneralCaseBundle x0 qClass ∧ smGrLimitBundle

/-- Existing theorem chain satisfies the integrated SM×QCD acceptance gate. -/
theorem sm_qcd_acceptance_gate_holds
    {X : Type}
    [TopologicalSpace X] [PreconnectedSpace X] [Nonempty X]
    (x0 : X)
    (qClass : HomotopyClass X Su3Matrix → ℤ) :
    smQcdAcceptanceGate x0 qClass := by
  refine ⟨sm_qcd_general_case_bundle_holds x0 qClass, sm_gr_limits_recovered⟩

end Gutoe.FalsifiabilityCatalog
