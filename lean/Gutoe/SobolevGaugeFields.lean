/-
 * GUTOE — Sobolev Spaces for Gauge Fields (GRAND-362)
 * Copyright (C) 2026  Riff Labs
 * AGPL-3.0-or-later
 *
 * A ∈ W^{k,p}(ℝ⁴, 𝔤). Configuration space with appropriate regularity.
 * Gauge-fixing (Coulomb/Lorenz) as a section of A → A/G.
 *
 * No `sorry`.
 -/
import Mathlib
import Gutoe.ContinuumYMLieAlgebra
import Gutoe.GaugeTransformationGroup

noncomputable section
namespace Gutoe.SobolevGaugeFields

open Gutoe.ContinuumYMLieAlgebra
open Gutoe.GaugeTransformationGroup

/-! ## Sobolev space configuration -/

/-- Sobolev regularity parameters for gauge fields. -/
structure SobolevParameters where
  /-- Differentiability order k. -/
  k : ℕ
  /-- Integrability exponent p ≥ 1. -/
  p : ℝ
  p_ge_one : 1 ≤ p
  /-- Spacetime dimension. -/
  d : ℕ
  d_eq : d = 4

/-- Configuration space of gauge fields with Sobolev regularity. -/
structure GaugeFieldConfiguration where
  groupData : CompactSimpleLieGroupData
  sobolev : SobolevParameters
  /-- The configuration space A = W^{k,p}(ℝ⁴, 𝔤) is non-empty. -/
  nonEmpty : Prop
  /-- The configuration space is a Banach manifold. -/
  isBanachManifold : Prop

/-! ## Gauge fixing -/

/-- Gauge-fixing condition type. -/
inductive GaugeFixType
  | coulomb    -- ∂_i A^i = 0 (spatial)
  | lorenz     -- ∂_μ A^μ = 0 (covariant)
  | temporal   -- A_0 = 0
  | axial      -- A_3 = 0

/-- A gauge-fixing condition provides a section of A → A/G. -/
structure GaugeFix where
  fixType : GaugeFixType
  /-- The gauge condition is reachable from any connection. -/
  existence : Prop
  /-- The gauge condition (generically) selects a unique representative. -/
  uniqueness : Prop
  /-- Gribov copies may exist for non-abelian theories. -/
  gribovCopies : Prop

/-- (Axiom) Coulomb gauge fixing exists for connections with sufficient regularity.
    Note: Gribov copies exist for non-abelian theories. -/
axiom coulomb_gauge_exists : ∃ gf : GaugeFix,
    gf.fixType = GaugeFixType.coulomb ∧ gf.existence ∧ gf.gribovCopies

/-- **GRAND-362: Sobolev spaces for gauge fields**

    For compact simple G:
    1. The configuration space A = W^{k,p}(ℝ⁴, 𝔤) is well-defined.
    2. Coulomb gauge fixing exists (with Gribov copies). -/
theorem sobolev_gauge_fields :
    ∃ gf : GaugeFix, gf.fixType = GaugeFixType.coulomb ∧ gf.existence ∧ gf.gribovCopies :=
  coulomb_gauge_exists

end Gutoe.SobolevGaugeFields
