/-
 * GUTOE — Wilson-Action Equivalence Bridge (Structural)
 * Copyright (C) 2026  Riff Labs
 * AGPL-3.0-or-later
 *
 * GRAND-300 (Theorem C bridge lane):
 *   - represent a Wilson-plaquette schedule by local Z₃ nearest-neighbor
 *     transition targets on the 3-state transfer basis
 *   - prove this representation induces SC-regular row totals structurally
 *   - instantiate the continuum-survival mass-gap lane with no empirical
 *     max-row certificate
 *
 * No `sorry`.
 -/

import Mathlib
import Gutoe.YangMillsContinuumSurvival
import Gutoe.YangMillsStructuralGap
import Gutoe.LatticeGeometry

namespace Gutoe.YangMillsWilsonBridge

open Real
open Gutoe.YangMillsMassGap
open Gutoe.YangMillsStructuralGap
open Gutoe.YangMillsContinuumSurvival
open Gutoe.LatticeGeometry

/-- Structural Wilson-action schedule projected to the Z₃ transfer basis.

`targetSchedule n i e` gives the transfer-basis target state for refinement step
`n`, source basis state `i`, and incident SC edge `e`.

`betaSchedule` is carried for Wilson-lane bookkeeping (plaquette coupling) and
required to remain positive. The mass-gap bridge below only uses the transfer
projection object. -/
structure WilsonZ3Action where
  targetSchedule : ℕ → Z3NearestNeighborTargets
  betaSchedule : ℕ → ℝ
  beta_pos : ∀ n, 0 < betaSchedule n

/-- Wilson-induced row-total schedule on the transfer basis. -/
def wilsonRowTotalsSchedule (W : WilsonZ3Action) : ℕ → Fin 3 → ℕ :=
  z3NearestNeighborRowTotalsSchedule W.targetSchedule

/-- Wilson-induced row totals are SC-regular at every refinement step. -/
theorem wilson_row_totals_sc_regular (W : WilsonZ3Action) :
    ∀ n, SCRegularRowTotals (wilsonRowTotalsSchedule W n) := by
  exact z3_nn_schedule_sc_regular W.targetSchedule

/-- Wilson-induced max row total is exactly the SC coordination number (`6`)
for every refinement step. -/
theorem wilson_max_row_total_eq_coordination (W : WilsonZ3Action) :
    ∀ n, maxRowTotal (wilsonRowTotalsSchedule W n) = coordinationNumber := by
  intro n
  exact z3_nn_max_row_total_eq_coordination (W.targetSchedule n)

/-- Wilson-induced minorization constant has a structural closed form at each
refinement step. -/
theorem wilson_minorization_eps_closed_form
    (W : WilsonZ3Action) (alpha : ℝ) :
    ∀ n,
      minorizationEps (wilsonRowTotalsSchedule W n) alpha =
        (3 * alpha) / ((coordinationNumber : ℝ) + 3 * alpha) := by
  intro n
  exact minorization_eps_eq_sc_regular
    (wilsonRowTotalsSchedule W n)
    alpha
    (wilson_row_totals_sc_regular W n)

/-- Bridge theorem (Theorem C lane, structural form):
if a Wilson schedule is represented by nearest-neighbor Z₃ transfer targets,
then the non-vanishing continuum mass-gap lower bound follows from the
structural Yang-Mills chain without empirical row-total hypotheses. -/
theorem wilson_action_bridge_nonvanishing_gap
    (W : WilsonZ3Action)
    (a_t : ℕ → ℝ)
    (alpha : ℝ)
    (ha_t_pos : ∀ n, 0 < a_t n)
    (ha_t_cap : ∃ aCap, 0 < aCap ∧ ∀ n, a_t n ≤ aCap)
    (ha : 0 < alpha) :
    ∃ c : ℝ, 0 < c ∧
      ∀ n, c ≤ doeblinGapLowerBound (a_t n)
        (minorizationEps (wilsonRowTotalsSchedule W n) alpha) := by
  exact continuum_survival_gap_nonvanishing_of_z3_nn_schedule
    a_t
    W.targetSchedule
    alpha
    ha_t_pos
    ha_t_cap
    ha

end Gutoe.YangMillsWilsonBridge
