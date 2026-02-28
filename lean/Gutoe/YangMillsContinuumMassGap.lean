/-
 * GUTOE — GRAND-333 Continuum Reconstructed Mass Gap
 * Copyright (C) 2026  Riff Labs
 * AGPL-3.0-or-later
 *
 * GRAND-333:
 *   Prove a strictly positive mass gap in the reconstructed continuum model
 *   itself (post-OS reconstruction), with a self-adjoint generator/Hamiltonian
 *   lane and uniform positive spectral floor.
 *
 * No `sorry`.
 -/

import Mathlib
import Gutoe.YangMillsOSEndToEnd

noncomputable section

namespace Gutoe.YangMillsContinuumMassGap

open Gutoe.YangMillsOSTextbook
open Gutoe.YangMillsOSCompletion
open Gutoe.YangMillsOSEndToEnd
open Gutoe.YangMillsWilsonBridge
open Gutoe.YangMillsWilsonEquivalence

/-- Continuum reconstructed mass-gap observable at refinement step `n`:
the positive Hamiltonian scale in the OS-completed continuum lane. -/
noncomputable def continuumMassGapAt
    (W : WilsonZ3Action)
    (a_t : ℕ → ℝ)
    (alpha : ℝ)
    (n : ℕ) : ℝ :=
  osHamiltonianAt W a_t alpha n

/-- GRAND-333 closure theorem:
in the reconstructed continuum lane (from GRAND-331/321), the model admits:
1. explicit OS end-to-end packages at every refinement step,
2. self-adjoint continuum generators,
3. a strict, uniform positive mass-gap floor `Δ > 0`. -/
theorem grand333_continuum_mass_gap_of_domain
    (W : WilsonZ3Action)
    (a_t : ℕ → ℝ)
    (alpha : ℝ)
    (hDom : WilsonEquivalenceDomain a_t alpha) :
    (∀ n, Nonempty (OSEndToEndStepPackage W a_t alpha n)) ∧
    (∀ n, IsSelfAdjoint (osGeneratorAt W a_t alpha n)) ∧
    (∃ Δ : ℝ, 0 < Δ ∧ ∀ n, Δ ≤ continuumMassGapAt W a_t alpha n) := by
  rcases grand331_end_to_end_os_reconstruction_of_domain W a_t alpha hDom with
    ⟨_, _, _, hPackages, _⟩
  refine ⟨hPackages, ?_, ?_⟩
  · intro n
    exact osGeneratorAt_selfAdjoint W a_t alpha n
  · rcases osGenerator_uniform_gap_floor_of_domain W a_t alpha hDom with ⟨Δ, hΔpos, hΔle⟩
    refine ⟨Δ, hΔpos, ?_⟩
    intro n
    simpa [continuumMassGapAt] using hΔle n

end Gutoe.YangMillsContinuumMassGap

