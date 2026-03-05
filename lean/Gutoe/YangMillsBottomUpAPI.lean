/-
 * GUTOE — Bottom-Up API Exports for Wilson Bridge
 * Copyright (C) 2026  Riff Labs
 * AGPL-3.0-or-later
 *
 * GRAND-482/495:
 *   Single clean export surface for the Phase 4 bridge lane.
 *
 * No `sorry`.
 -/

import Mathlib
import Gutoe.YangMillsWilsonBridge
import Gutoe.YangMillsWilsonEquivalence
import Gutoe.YangMillsOSEndToEnd
import Gutoe.YangMillsContinuumMassGap
import Gutoe.YangMillsUnconditional

noncomputable section

namespace Gutoe.YangMillsBottomUpAPI

/-- Re-export of structural nearest-neighbor transfer targets.

Type signature: `Fin 3 → Fin coordinationNumber → Fin 3`. -/
abbrev Z3NearestNeighborTargets : Type :=
  Gutoe.YangMillsStructuralGap.Z3NearestNeighborTargets

/-- Re-export of the Wilson action projected to the Z₃ transfer basis.

Type signature: `Type` (structure with target and β schedules). -/
abbrev WilsonZ3Action : Type :=
  Gutoe.YangMillsWilsonBridge.WilsonZ3Action

/-- Re-export of the explicit Theorem-C domain assumptions.

Type signature: `(a_t : ℕ → ℝ) → (alpha : ℝ) → Prop`. -/
abbrev WilsonEquivalenceDomain (a_t : ℕ → ℝ) (alpha : ℝ) : Prop :=
  Gutoe.YangMillsWilsonEquivalence.WilsonEquivalenceDomain a_t alpha

/-- Re-export of the GRAND-331 explicit per-step OS package.

Type signature:
`(W : WilsonZ3Action) → (a_t : ℕ → ℝ) → (alpha : ℝ) → (n : ℕ) → Type`. -/
abbrev OSEndToEndStepPackage
    (W : WilsonZ3Action)
    (a_t : ℕ → ℝ)
    (alpha : ℝ)
    (n : ℕ) : Type :=
  Gutoe.YangMillsOSEndToEnd.OSEndToEndStepPackage W a_t alpha n

/-- Re-export of GRAND-333 continuum mass-gap closure over a domain.

Type signature:
`(W : WilsonZ3Action) → (a_t : ℕ → ℝ) → (alpha : ℝ) →
  WilsonEquivalenceDomain a_t alpha →
  (∀ n, Nonempty (OSEndToEndStepPackage W a_t alpha n)) ∧
  (∀ n, IsSelfAdjoint (osGeneratorAt W a_t alpha n)) ∧
  (∃ Δ : ℝ, 0 < Δ ∧ ∀ n, Δ ≤ continuumMassGapAt W a_t alpha n)`. -/
abbrev grand333_continuum_mass_gap_of_domain
    (W : WilsonZ3Action)
    (a_t : ℕ → ℝ)
    (alpha : ℝ)
    (hDom : WilsonEquivalenceDomain a_t alpha) :
    (∀ n, Nonempty (OSEndToEndStepPackage W a_t alpha n)) ∧
    (∀ n, IsSelfAdjoint (Gutoe.YangMillsOSTextbook.osGeneratorAt W a_t alpha n)) ∧
    (∃ Δ : ℝ, 0 < Δ ∧
      ∀ n, Δ ≤ Gutoe.YangMillsContinuumMassGap.continuumMassGapAt W a_t alpha n) :=
  Gutoe.YangMillsContinuumMassGap.grand333_continuum_mass_gap_of_domain W a_t alpha hDom

/-- Re-export of the unconditional canonical Wilson-action witness.

Type signature: `WilsonZ3Action`. -/
abbrev canonicalWilsonZ3Action : WilsonZ3Action :=
  Gutoe.YangMillsUnconditional.canonicalWilsonZ3Action

/-- Re-export of the unconditional canonical domain witness.

Type signature:
`WilsonEquivalenceDomain unitLatticeSpacing unitAlpha`. -/
abbrev canonicalDomain :
    WilsonEquivalenceDomain
      Gutoe.YangMillsUnconditional.unitLatticeSpacing
      Gutoe.YangMillsUnconditional.unitAlpha :=
  Gutoe.YangMillsUnconditional.canonicalDomain

/-- Single Phase-4 bridge bundle for bottom-up consumers. -/
structure BottomUpAPIPackage where
  W : WilsonZ3Action
  a_t : ℕ → ℝ
  alpha : ℝ
  domain : WilsonEquivalenceDomain a_t alpha
  endToEnd :
    ∀ n, Nonempty (OSEndToEndStepPackage W a_t alpha n)
  selfAdjoint :
    ∀ n, IsSelfAdjoint (Gutoe.YangMillsOSTextbook.osGeneratorAt W a_t alpha n)
  massGap :
    ∃ Δ : ℝ, 0 < Δ ∧
      ∀ n, Δ ≤ Gutoe.YangMillsContinuumMassGap.continuumMassGapAt W a_t alpha n

/-- Constructor that packages all bridge obligations from `grand333`. -/
def mkBottomUpAPIPackage
    (W : WilsonZ3Action)
    (a_t : ℕ → ℝ)
    (alpha : ℝ)
    (hDom : WilsonEquivalenceDomain a_t alpha) :
    BottomUpAPIPackage := by
  rcases grand333_continuum_mass_gap_of_domain W a_t alpha hDom with
    ⟨hEndToEnd, hSelfAdjoint, hMassGap⟩
  exact
    { W := W
      a_t := a_t
      alpha := alpha
      domain := hDom
      endToEnd := hEndToEnd
      selfAdjoint := hSelfAdjoint
      massGap := hMassGap }

/-- Canonical unconditional Phase-4 bridge package. -/
def canonicalBottomUpAPIPackage : BottomUpAPIPackage :=
  mkBottomUpAPIPackage
    canonicalWilsonZ3Action
    Gutoe.YangMillsUnconditional.unitLatticeSpacing
    Gutoe.YangMillsUnconditional.unitAlpha
    canonicalDomain

end Gutoe.YangMillsBottomUpAPI
