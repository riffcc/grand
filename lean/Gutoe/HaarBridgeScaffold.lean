/-
 * GUTOE — Haar Bridge Scaffold (Path-2, structural layer)
 * Copyright (C) 2026  Riff Labs
 *
 * AGPL-3.0-or-later
 *
 * GRAND-308:
 *   Build the Lean group/quotient scaffold needed before full Haar-measure
 *   decomposition over SU(3)/Z3.
 *
 * This module intentionally stays at the group-theoretic layer:
 * - center subgroup and quotient object
 * - quotient lift/descent theorems
 * - fiber-invariance under center action for observables that factor through
 *   the quotient (or for homomorphisms killing center elements)
 *
 * No `sorry`.
 -/

import Mathlib
import Gutoe.GaugeGroupSU3
import Gutoe.YangMillsWilsonBridge

noncomputable section

namespace Gutoe.HaarBridgeScaffold

open QuotientGroup
open Gutoe.GaugeGroupSU3
open Gutoe.YangMillsWilsonBridge

/-- Cl(1,3)->Z3->SU(3) construction anchor reused by the Haar bridge lane. -/
theorem clifford_z3_su3_anchor :
    quarkOrbit.card = 3 ∧
    Nonempty ({s // s ∈ quarkOrbit} ≃ Fin 3) ∧
    quarkOrbit.card ^ 2 - 1 = 8 := by
  exact ⟨quarkOrbit_card, quarkOrbit_equiv_fin3, quarks_predict_gluon_count⟩

section GroupQuotient

variable {G M : Type*} [Group G] [Group M]

/-- The center subgroup used as the Path-2 quotient kernel prototype. -/
abbrev centerSubgroup : Subgroup G := Subgroup.center G

/-- The center is normal, so the quotient `G ⧸ centerSubgroup` is available. -/
theorem center_normal : (centerSubgroup (G := G)).Normal := inferInstance

/-- The quotient by center is inhabited (structural existence witness). -/
theorem center_quotient_nonempty : Nonempty (G ⧸ centerSubgroup (G := G)) :=
  ⟨(1 : G ⧸ centerSubgroup (G := G))⟩

/-- A homomorphism whose kernel contains center descends to the center quotient. -/
def descendHomThroughCenter
    (φ : G →* M)
    (hCenter : centerSubgroup (G := G) ≤ φ.ker) :
    G ⧸ centerSubgroup (G := G) →* M :=
  QuotientGroup.lift (centerSubgroup (G := G)) φ hCenter

/-- The descended homomorphism composes back to the original map. -/
theorem descendHomThroughCenter_comp_mk
    (φ : G →* M)
    (hCenter : centerSubgroup (G := G) ≤ φ.ker) :
    (descendHomThroughCenter (G := G) (M := M) φ hCenter).comp
      (QuotientGroup.mk' (centerSubgroup (G := G))) = φ := by
  exact QuotientGroup.lift_comp_mk' (centerSubgroup (G := G)) φ hCenter

/-- Uniqueness of the descended map from equality on the original carrier. -/
theorem descendHomThroughCenter_unique
    (φ : G →* M)
    (hCenter : centerSubgroup (G := G) ≤ φ.ker)
    (ψ : G ⧸ centerSubgroup (G := G) →* M)
    (hψ : ψ.comp (QuotientGroup.mk' (centerSubgroup (G := G))) = φ) :
    ψ = descendHomThroughCenter (G := G) (M := M) φ hCenter := by
  apply QuotientGroup.monoidHom_ext (N := centerSubgroup (G := G))
  calc
    ψ.comp (QuotientGroup.mk' (centerSubgroup (G := G))) = φ := hψ
    _ = (descendHomThroughCenter (G := G) (M := M) φ hCenter).comp
          (QuotientGroup.mk' (centerSubgroup (G := G))) := by
          simpa [descendHomThroughCenter] using
            (QuotientGroup.lift_comp_mk' (centerSubgroup (G := G)) φ hCenter).symm

/-- Fiber invariance predicate along center right-cosets. -/
def CenterFiberInvariant (f : G → M) : Prop :=
  ∀ g z, z ∈ centerSubgroup (G := G) → f (g * z) = f g

/-- If a homomorphism kills center elements, it is center-fiber invariant. -/
theorem hom_center_fiber_invariant
    (φ : G →* M)
    (hCenter : centerSubgroup (G := G) ≤ φ.ker) :
    CenterFiberInvariant (G := G) (M := M) φ := by
  intro g z hz
  have hzker : z ∈ φ.ker := hCenter hz
  have hzOne : φ z = 1 := hzker
  calc
    φ (g * z) = φ g * φ z := by simp
    _ = φ g * 1 := by rw [hzOne]
    _ = φ g := by simp

end GroupQuotient

end Gutoe.HaarBridgeScaffold

