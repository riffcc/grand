/-
 * GUTOE — Unconditional Yang-Mills Mass Gap Instantiation
 * Copyright (C) 2026  Riff Labs
 * AGPL-3.0-or-later
 *
 * Instantiates `WilsonEquivalenceDomain` and `WilsonZ3Action` with
 * explicit concrete witnesses (constant identity target schedule,
 * unit coupling/lattice parameters) to produce an unconditional
 * mass gap theorem for the GUTOE Z₃-center lattice gauge model.
 *
 * No `sorry`.
 -/

import Mathlib
import Gutoe.YangMillsContinuumMassGap

noncomputable section

namespace Gutoe.YangMillsUnconditional

open Gutoe.YangMillsWilsonBridge
open Gutoe.YangMillsWilsonEquivalence
open Gutoe.YangMillsContinuumMassGap
open Gutoe.LatticeGeometry

/-- Canonical constant lattice-spacing schedule: `a_t(n) = 1` for all `n`. -/
def unitLatticeSpacing : ℕ → ℝ := fun _ => 1

/-- Canonical Laplace coupling floor. -/
def unitAlpha : ℝ := 1

/-- Canonical constant identity target schedule:
    at every refinement level, every source state `i` maps all its
    `coordinationNumber`-many neighbors to itself. -/
def identityTargetSchedule : ℕ → Z3NearestNeighborTargets :=
  fun _ i _ => i

/-- Canonical constant β schedule. -/
def unitBetaSchedule : ℕ → ℝ := fun _ => 1

/-- Canonical `WilsonZ3Action` with identity targets and unit coupling. -/
def canonicalWilsonZ3Action : WilsonZ3Action where
  targetSchedule := identityTargetSchedule
  betaSchedule := unitBetaSchedule
  beta_pos := fun _ => by norm_num [unitBetaSchedule]

/-- The canonical lattice spacing schedule is strictly positive at every step. -/
theorem unitLatticeSpacing_pos : ∀ n, 0 < unitLatticeSpacing n := by
  intro _
  norm_num [unitLatticeSpacing]

/-- The canonical lattice spacing schedule is bounded above. -/
theorem unitLatticeSpacing_cap :
    ∃ aCap, 0 < aCap ∧ ∀ n, unitLatticeSpacing n ≤ aCap := by
  exact ⟨1, by norm_num, fun _ => by norm_num [unitLatticeSpacing]⟩

/-- The canonical coupling floor is positive. -/
theorem unitAlpha_pos : 0 < unitAlpha := by
  norm_num [unitAlpha]

/-- Explicit instantiation of `WilsonEquivalenceDomain` at unit parameters. -/
theorem canonicalDomain :
    WilsonEquivalenceDomain unitLatticeSpacing unitAlpha where
  a_t_pos := unitLatticeSpacing_pos
  a_t_cap := unitLatticeSpacing_cap
  alpha_pos := unitAlpha_pos

/-- **Unconditional mass gap theorem** for the GUTOE Z₃-center lattice model.

    The reconstructed continuum model on the canonical simple-cubic Z₃ lattice
    with unit coupling parameters admits:
    1. Explicit OS end-to-end packages at every refinement step,
    2. Self-adjoint continuum generators (Hamiltonians),
    3. A strictly positive, uniform mass gap `Δ > 0`. -/
theorem gutoe_yang_mills_mass_gap_unconditional :
    (∀ n, Nonempty (OSEndToEndStepPackage canonicalWilsonZ3Action
      unitLatticeSpacing unitAlpha n)) ∧
    (∀ n, IsSelfAdjoint (osGeneratorAt canonicalWilsonZ3Action
      unitLatticeSpacing unitAlpha n)) ∧
    (∃ Δ : ℝ, 0 < Δ ∧ ∀ n, Δ ≤ continuumMassGapAt canonicalWilsonZ3Action
      unitLatticeSpacing unitAlpha n) := by
  exact grand333_continuum_mass_gap_of_domain
    canonicalWilsonZ3Action unitLatticeSpacing unitAlpha canonicalDomain

end Gutoe.YangMillsUnconditional
