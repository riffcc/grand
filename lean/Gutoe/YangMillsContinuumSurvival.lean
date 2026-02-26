/-
 * GUTOE — Yang-Mills Continuum Survival Bridge
 * Copyright (C) 2026  Riff Labs
 * AGPL-3.0-or-later
 *
 * GRAND-299:
 *   Prove the structural Doeblin mass-gap lower bound remains non-vanishing
 *   under an explicit continuum-survival hypothesis package.
 *
 * No `sorry`.
-/

import Mathlib
import Gutoe.YangMillsMassGap

namespace Gutoe.YangMillsContinuumSurvival

open Real
open Filter
open scoped Topology
open Gutoe.YangMillsMassGap

/-- Explicit hypothesis package for a non-vanishing Doeblin gap lane over a
joint refinement schedule `n ↦ (a_t n, eps n)`.

`a_t` is the temporal lattice spacing used in the transfer step and `eps` is the
Doeblin minorization constant.

This package is intentionally explicit:
- positivity/range for each step,
- a uniform positive lower floor on `eps`,
- a uniform positive upper cap on `a_t`.
These are exactly the assumptions needed to force a strictly positive global
lower bound on the Doeblin gap estimator.
-/
def ContinuumSurvivalHypotheses (a_t eps : ℕ → ℝ) : Prop :=
  (∀ n, 0 < a_t n) ∧
  (∀ n, 0 < eps n ∧ eps n < 1) ∧
  (∃ epsFloor, 0 < epsFloor ∧ ∀ n, epsFloor ≤ eps n) ∧
  (∃ aCap, 0 < aCap ∧ ∀ n, a_t n ≤ aCap)

/-- Pointwise lower bound from explicit floor/cap assumptions. -/
theorem doeblin_gap_lower_bound_of_floor_cap
    (a_t eps : ℕ → ℝ)
    {epsFloor aCap : ℝ}
    (haPos : ∀ n, 0 < a_t n)
    (hepsRange : ∀ n, 0 < eps n ∧ eps n < 1)
    (hepsFloorLe : ∀ n, epsFloor ≤ eps n)
    (haCapPos : 0 < aCap)
    (haCap : ∀ n, a_t n ≤ aCap) :
    ∀ n, doeblinGapLowerBound aCap epsFloor ≤ doeblinGapLowerBound (a_t n) (eps n) := by
  intro n
  have hepsnPos : 0 < eps n := (hepsRange n).1
  have hepsnLt1 : eps n < 1 := (hepsRange n).2
  have honeMinusPos : 0 < 1 - eps n := by linarith
  have honeMinusLt : 1 - eps n ≤ 1 - epsFloor := by linarith [hepsFloorLe n]
  have hlogMon :
      Real.log (1 - eps n) ≤ Real.log (1 - epsFloor) := by
    exact Real.log_le_log honeMinusPos honeMinusLt
  have hnumMon :
      -Real.log (1 - epsFloor) ≤ -Real.log (1 - eps n) := by
    linarith
  have hnumNonneg : 0 ≤ -Real.log (1 - eps n) := by
    have hlogNeg : Real.log (1 - eps n) < 0 := by
      have hlt1 : 1 - eps n < 1 := by linarith
      have : Real.log (1 - eps n) < Real.log 1 :=
        Real.log_lt_log honeMinusPos hlt1
      simpa using this
    linarith
  have hstep1 :
      (-Real.log (1 - epsFloor)) / aCap ≤ (-Real.log (1 - eps n)) / aCap := by
    exact div_le_div_of_nonneg_right hnumMon (le_of_lt haCapPos)
  have hstep2 :
      (-Real.log (1 - eps n)) / aCap ≤ (-Real.log (1 - eps n)) / (a_t n) := by
    exact div_le_div_of_nonneg_left hnumNonneg (haPos n) (haCap n)
  simpa [doeblinGapLowerBound] using le_trans hstep1 hstep2

/-- Structural continuum-survival bridge:
if Doeblin floors/caps are uniform, the mass-gap lower bound is uniformly
strictly positive across the entire refinement schedule. -/
theorem continuum_survival_gap_nonvanishing
    (a_t eps : ℕ → ℝ)
    (h : ContinuumSurvivalHypotheses a_t eps) :
    ∃ c : ℝ, 0 < c ∧ ∀ n, c ≤ doeblinGapLowerBound (a_t n) (eps n) := by
  rcases h with ⟨haPos, hepsRange, ⟨epsFloor, hepsFloorPos, hepsFloorLe⟩, ⟨aCap, haCapPos, haCap⟩⟩
  refine ⟨doeblinGapLowerBound aCap epsFloor, ?_, ?_⟩
  · have hepsFloorLt1 : epsFloor < 1 := by
      exact lt_of_le_of_lt (hepsFloorLe 0) (hepsRange 0).2
    exact doeblin_bound_positive haCapPos hepsFloorPos hepsFloorLt1
  · exact doeblin_gap_lower_bound_of_floor_cap
      a_t eps haPos hepsRange hepsFloorLe haCapPos haCap

end Gutoe.YangMillsContinuumSurvival
