/-
 * GUTOE - Real Gate Behavior (Formal Verification)
 * Copyright (C) 2026  Riff Labs
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 -/

import Mathlib
import Gutoe.GateProperties

namespace Gutoe

/-!
## Real Gate Behavior — Python Source vs Rust Port

The Rust codebase under formal scrutiny is a port of a Python original.
This module formalises the DIFFERENCES between the two and proves key
properties of the Python semantics that the port lost.

Key findings, all PROVEN:

1. **Real CNOT is Z₂ (self-inverse).**  The Python `TripartiteCNOT` swaps
   SINE↔COSINE when the control is SINE — a clean involution.  The Rust port
   replaced this with `cycle` (Z₃), which requires three applications to
   return to the original state.

2. **Real Hadamard is invertible when phase is tracked.**  The Python gate
   records whether the superposition came from SINE (+π/4) or COSINE (-π/4).
   A `PhaseState` type with an `OriginTag` captures this.  The state-only
   TriState Hadamard from `GateProperties` collapses both inputs to TANGENT
   and is irreversible; the phase-aware version is fully involutive.

3. **Quark mass is quadratic in λ_QG.**  Expanding the Python
   `particle_formation.py` formula shows `m ∝ λ_QG²`.  The Rust port uses
   λ_QG = 0.084372 (from the LQG curve fit below); the Python tuner converged
   on 0.120000.  This predicts more than a 2× mass difference — proven by `norm_num`.

4. **λ_QG = 0.084372 is the `b` parameter of a Loop QG velocity fit.**
   `experiment-28/metric_extraction.py` fit veracity-rail wave speeds around a
   simulated black hole to `v = v_max·(1 − exp(−(r − r_s)/a)) + b`, obtaining
   a = 29.366106, **b = 0.084372**, R² = 0.990642.  The `b` parameter is the
   residual velocity floor at the horizon.  This value was then adopted as λ_QG
   for the dispersion relation `ω² = v²k² − λ_QG l_P² k⁴` — an identification
   by analogy, since a velocity offset and a dispersion coupling are
   dimensionally distinct quantities.  The standard LQG Barbero–Immirzi
   parameter γ ≈ 0.2375 is notably closer to the tuned value (0.120) than to
   the fitted b (0.084372).
-/

-- ── Real CNOT: SINE↔COSINE swap (Z₂ involution) ─────────────────────────

/-- The Python TripartiteCNOT: swap SINE↔COSINE when control is SINE -/
def realCNOT (control target : TriState) : TriState :=
  if control = TriState.SINE then
    match target with
    | TriState.SINE      => TriState.COSINE
    | TriState.COSINE    => TriState.SINE
    | TriState.TANGENT => TriState.TANGENT
    | TriState.VOID      => TriState.VOID
  else
    target

/-- realCNOT with non-SINE control is the identity — REAL -/
theorem real_cnot_identity_when_not_sine
    (control target : TriState) (h : control ≠ TriState.SINE) :
    realCNOT control target = target := by
  simp [realCNOT, h]

/-- realCNOT is self-inverse (Z₂): applying it twice restores the target — REAL -/
theorem real_cnot_self_inverse (target : TriState) :
    realCNOT TriState.SINE (realCNOT TriState.SINE target) = target := by
  cases target <;> simp [realCNOT]

/-- The Rust port's cCycle on COSINE gives TANGENT; realCNOT gives SINE — REAL -/
theorem rust_cnot_differs_from_real_on_cosine :
    cCycle TriState.SINE TriState.COSINE ≠ realCNOT TriState.SINE TriState.COSINE := by
  decide

/-- The Rust port's cCycle on TANGENT gives SINE; realCNOT leaves it unchanged — REAL -/
theorem rust_cnot_differs_from_real_on_undefined :
    cCycle TriState.SINE TriState.TANGENT ≠ realCNOT TriState.SINE TriState.TANGENT := by
  decide

/-- Rust port is NOT self-inverse on COSINE (applying cCycle twice ≠ identity) — REAL -/
theorem rust_cnot_not_self_inverse_cosine :
    cCycle TriState.SINE (cCycle TriState.SINE TriState.COSINE) ≠ TriState.COSINE := by
  decide

-- ── Phase-aware Hadamard: involutive when origin tag is tracked ───────────

/-- Origin tag recording which basis state produced the superposition -/
inductive OriginTag : Type
  | fromSINE   : OriginTag
  | fromCOSINE : OriginTag
  | generic    : OriginTag
  deriving DecidableEq, Repr

/-- Extended state with phase metadata -/
inductive PhaseState : Type
  | VOID
  | SINE
  | COSINE
  | SUPER (tag : OriginTag)   -- superposition carrying its origin
  deriving DecidableEq, Repr

/-- Phase-aware Hadamard: invertible because the tag distinguishes origins -/
def phaseHadamard : PhaseState → PhaseState
  | .SINE              => .SUPER .fromSINE
  | .COSINE            => .SUPER .fromCOSINE
  | .SUPER .fromSINE   => .SINE
  | .SUPER .fromCOSINE => .COSINE
  | .SUPER .generic    => .SUPER .generic
  | .VOID              => .VOID

/-- phaseHadamard² = id on SINE — REAL -/
theorem phase_hadamard_involutive_sine :
    phaseHadamard (phaseHadamard PhaseState.SINE) = PhaseState.SINE := rfl

/-- phaseHadamard² = id on COSINE — REAL -/
theorem phase_hadamard_involutive_cosine :
    phaseHadamard (phaseHadamard PhaseState.COSINE) = PhaseState.COSINE := rfl

/-- phaseHadamard² = id on VOID — REAL -/
theorem phase_hadamard_involutive_void :
    phaseHadamard (phaseHadamard PhaseState.VOID) = PhaseState.VOID := rfl

/-- Phase-aware superpositions of SINE and COSINE are distinguishable — REAL -/
theorem phase_hadamard_distinguishes_sine_cosine :
    phaseHadamard PhaseState.SINE ≠ phaseHadamard PhaseState.COSINE := by decide

/-- In contrast, the state-only Hadamard (from GateProperties) still collapses
    both inputs to the same output — REAL -/
theorem state_only_hadamard_still_collapses :
    hadamard TriState.SINE = hadamard TriState.COSINE := rfl

-- ── Quark mass: quadratic in λ_QG ────────────────────────────────────────

/-- Quark mass formula from Python particle_formation.py:
    m = veracity × curvature × field_gradient / planck_length × λ_QG²
    Simplified to highlight the quadratic dependence on λ_QG. -/
noncomputable def quarkMass (v c g l lam : ℝ) : ℝ :=
  v * c * g / l * lam ^ 2

/-- Doubling λ_QG quadruples the predicted mass — REAL -/
theorem mass_quadratic_in_lambda (v c g l lam : ℝ) :
    quarkMass v c g l (2 * lam) = 4 * quarkMass v c g l lam := by
  simp [quarkMass]; ring

/-- The tuned Python λ_QG (0.120) vs Rust port (0.084372):
    mass ratio > 2 — the Rust port predicts less than half the mass — REAL -/
theorem rust_vs_tuned_lambda_mass_ratio :
    (0.12 : ℝ) ^ 2 / (0.084372 : ℝ) ^ 2 > 2 := by norm_num

/-- Any two λ values with ratio r give mass ratio r² — REAL -/
theorem mass_ratio_is_lambda_ratio_squared
    (v c g l lam1 lam2 : ℝ) (hl : l ≠ 0) (hlam1 : lam1 ≠ 0) (hv : v ≠ 0) (hc : c ≠ 0)
    (hg : g ≠ 0) :
    quarkMass v c g l lam2 / quarkMass v c g l lam1 = (lam2 / lam1) ^ 2 := by
  simp [quarkMass]
  field_simp

-- ── LQG curve-fit origin of λ_QG = 0.084372 ─────────────────────────────

/-!
### LQG velocity model (experiment-28)

`metric_extraction.py` fit the wave velocity profile around a simulated
black hole to the Loop Quantum Gravity model:

    v(r) = v_max · (1 − exp(−(r − r_s) / a)) + b     (r > r_s)

Fit result: a = 29.366106, **b = 0.084372**, R² = 0.990642.

The parameter `b` is the value of v(r) as r → r_s⁺ (the velocity floor at
the event horizon).  It was subsequently adopted as λ_QG throughout the
codebase, but the quantities have different physical roles:
- `b`: velocity offset in the LQG velocity model (simulation velocity units)
- `λ_QG`: coupling strength in the dispersion relation `ω² = v²k² − λ_QG l_P² k⁴`
-/

/-- With b = 0.084372 and v_max = 0.2449, the LQG velocity floor is ~34%
    of v_max — not a small perturbation — REAL -/
theorem lqg_b_is_large_fraction_of_vmax :
    (0.084372 : ℝ) / 0.2449 > 1 / 4 := by norm_num

/-- The tuned λ_QG (0.120) is closer to the Barbero–Immirzi γ proxy (0.2375)
    than the LQG-fitted b (0.084372) is — REAL
    (Using the rational proxy 19/80 ≈ 0.2375 for γ) -/
theorem tuned_lambda_closer_to_barbero_immirzi :
    |(0.12 : ℝ) - 0.2375| < |(0.084372 : ℝ) - 0.2375| := by norm_num

/-- The LQG velocity formula evaluated at the horizon (r = r_s) gives exactly b.
    v(r_s) = v_max * (1 - exp(0)) + b = v_max * 0 + b = b — REAL -/
theorem lqg_velocity_at_horizon_is_b (v_max a b : ℝ) :
    v_max * (1 - Real.exp (-(0 : ℝ) / a)) + b = b ↔ v_max * (1 - 1) = 0 := by
  simp [Real.exp_zero]

/-- The LQG velocity limit as r → r_s⁺ is b (horizon velocity floor) — REAL -/
theorem lqg_horizon_limit_is_b (v_max a b : ℝ) :
    Filter.Tendsto
      (fun r => v_max * (1 - Real.exp (-r / a)) + b)
      (nhds 0) (nhds b) := by
  have hc : Continuous (fun r : ℝ => v_max * (1 - Real.exp (-r / a)) + b) :=
    (continuous_const.mul
      (continuous_const.sub
        (Real.continuous_exp.comp (continuous_neg.div_const a)))).add continuous_const
  have h0 : (fun r : ℝ => v_max * (1 - Real.exp (-r / a)) + b) 0 = b := by
    simp [Real.exp_zero]
  have hca : Filter.Tendsto (fun r : ℝ => v_max * (1 - Real.exp (-r / a)) + b)
      (nhds 0) (nhds ((fun r : ℝ => v_max * (1 - Real.exp (-r / a)) + b) 0)) :=
    hc.continuousAt
  rw [h0] at hca
  exact hca

end Gutoe
