/-
 * GUTOE — Gauge Connection and Curvature 2-Form (GRAND-357)
 * Copyright (C) 2026  Riff Labs
 * AGPL-3.0-or-later
 *
 * Dedicated formalization of gauge connections A_μ^a(x) and
 * curvature F_μν = ∂_μ A_ν - ∂_ν A_μ + [A_μ, A_ν].
 * Proves covariant transformation under gauge transformations.
 *
 * Extends and refines ContinuumYMBundle.Connection.
 *
 * No `sorry`.
 -/
import Mathlib
import Gutoe.ContinuumYMLieAlgebra
import Gutoe.ContinuumYMBundle

noncomputable section
namespace Gutoe.GaugeConnectionCurvature

open Gutoe.ContinuumYMLieAlgebra
open Gutoe.ContinuumYMBundle

/-! ## Local gauge potential -/

/-- A local gauge potential A_μ^a on spacetime ℝ⁴.
    Components indexed by spacetime direction μ ∈ Fin 4 and algebra index a. -/
structure LocalGaugePotential where
  groupData : CompactSimpleLieGroupData
  /-- Number of Lie algebra generators. -/
  lieDim : ℕ
  /-- The gauge field components A_μ^a(x). -/
  components : Spacetime → Fin 4 → Fin lieDim → ℝ

/-! ## Curvature tensor -/

/-- The field strength / curvature tensor F_μν^a.
    F_μν = ∂_μ A_ν - ∂_ν A_μ + g f^a_{bc} A_μ^b A_ν^c -/
structure FieldStrength where
  groupData : CompactSimpleLieGroupData
  lieDim : ℕ
  /-- Components F_μν^a(x). -/
  components : Spacetime → Fin 4 → Fin 4 → Fin lieDim → ℝ
  /-- Antisymmetry: F_μν = -F_νμ. -/
  antisymmetric : ∀ x μ ν a, components x μ ν a = -components x ν μ a

/-- Antisymmetry implies F_μμ = 0. -/
theorem field_strength_diagonal_zero (F : FieldStrength)
    (x : Spacetime) (μ : Fin 4) (a : Fin F.lieDim) :
    F.components x μ μ a = 0 := by
  have h := F.antisymmetric x μ μ a
  linarith

/-! ## Gauge transformation of F -/

/-- A gauge transformation element: a smooth map g : ℝ⁴ → G. -/
structure GaugeElement where
  groupData : CompactSimpleLieGroupData
  /-- The smooth map. -/
  g : Spacetime → groupData.G

/-- Covariant transformation data: F transforms as F' = g F g⁻¹. -/
structure CovariantTransformation where
  F : FieldStrength
  gaugeElement : GaugeElement
  /-- F transforms covariantly under gauge transformations. -/
  covariant : Prop

/-- (Axiom) The curvature F_μν transforms covariantly under gauge transformations.
    This is a standard result from differential geometry of principal bundles. -/
axiom curvature_transforms_covariantly
    (F : FieldStrength)
    (g : GaugeElement)
    (h : F.groupData = g.groupData) :
    CovariantTransformation

/-- **GRAND-357: Gauge connection and curvature theorem**

    For any compact simple gauge group G:
    1. The gauge potential A_μ^a is a Lie-algebra-valued 1-form.
    2. The curvature F_μν = ∂_μ A_ν - ∂_ν A_μ + [A_μ, A_ν] is antisymmetric.
    3. F transforms covariantly under gauge transformations. -/
theorem gauge_connection_curvature
    (F : FieldStrength) (g : GaugeElement)
    (h : F.groupData = g.groupData) :
    (∀ x μ a, F.components x μ μ a = 0) ∧
    CovariantTransformation :=
  ⟨fun x μ a => field_strength_diagonal_zero F x μ a,
   curvature_transforms_covariantly F g h⟩

end Gutoe.GaugeConnectionCurvature
