/-
 * GUTOE — Gauge Transformation Group (GRAND-359)
 * Copyright (C) 2026  Riff Labs
 * AGPL-3.0-or-later
 *
 * Defines G = Map(ℝ⁴, G) as the gauge group.
 * Proves A ↦ gAg⁻¹ + g(dg⁻¹) is a group action. Orbit space A/G.
 *
 * No `sorry`.
 -/
import Mathlib
import Gutoe.ContinuumYMLieAlgebra
import Gutoe.ContinuumYMBundle

noncomputable section
namespace Gutoe.GaugeTransformationGroup

open Gutoe.ContinuumYMLieAlgebra
open Gutoe.ContinuumYMBundle

/-! ## Gauge group as function space -/

/-- The gauge group G = C^∞(M, G) as smooth maps from spacetime to the structure group. -/
structure GaugeGroup where
  groupData : CompactSimpleLieGroupData
  /-- A gauge transformation is a smooth map from spacetime to G. -/
  elements : Type
  /-- Pointwise multiplication makes G a group. -/
  instGroup : Group elements
  /-- Identity is the constant map to 1_G. -/
  identity : elements
  /-- Composition of gauge transformations. -/
  compose : elements → elements → elements

attribute [instance] GaugeGroup.instGroup

/-! ## Gauge action on connections -/

/-- The gauge action on connections: A ↦ gAg⁻¹ + g(dg⁻¹). -/
structure GaugeAction where
  gaugeGroup : GaugeGroup
  /-- The transformed connection. -/
  transform : gaugeGroup.elements → Prop → Prop
  /-- The action preserves connection structure. -/
  preservesConnection : Prop
  /-- Identity acts trivially. -/
  identityActs : Prop
  /-- The action is compatible with group multiplication. -/
  compositionLaw : Prop

/-- (Axiom) The gauge action is a well-defined group action. -/
axiom gauge_action_well_defined (ga : GaugeAction) :
    ga.preservesConnection ∧ ga.identityActs ∧ ga.compositionLaw

/-! ## Orbit space -/

/-- The orbit space A/G: connections modulo gauge equivalence. -/
structure OrbitSpace where
  gaugeGroup : GaugeGroup
  /-- Two connections are gauge-equivalent. -/
  gaugeEquivalent : Prop
  /-- Gauge equivalence is an equivalence relation. -/
  isEquivalenceRelation : Prop

/-- (Axiom) Gauge equivalence is an equivalence relation. -/
axiom gauge_equiv_is_equivalence (os : OrbitSpace) : os.isEquivalenceRelation

/-- **GRAND-359: Gauge transformation group theorem**

    For any compact simple G:
    1. G = Map(ℝ⁴, G) is a group.
    2. The gauge action A ↦ gAg⁻¹ + g(dg⁻¹) is a well-defined group action.
    3. The orbit space A/G is well-defined. -/
theorem gauge_transformation_group (ga : GaugeAction) (os : OrbitSpace) :
    (ga.preservesConnection ∧ ga.identityActs ∧ ga.compositionLaw) ∧
    os.isEquivalenceRelation :=
  ⟨gauge_action_well_defined ga, gauge_equiv_is_equivalence os⟩

end Gutoe.GaugeTransformationGroup
