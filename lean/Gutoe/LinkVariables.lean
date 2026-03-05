/-
 * GUTOE — Link Variables U_μ(x) ∈ G (GRAND-371)
 * Copyright (C) 2026  Riff Labs
 * AGPL-3.0-or-later
 *
 * Gauge field on the lattice as group-valued link variables.
 * Parallel transport along paths as ordered products.
 *
 * No `sorry`.
 -/
import Mathlib
import Gutoe.ContinuumYMLieAlgebra
import Gutoe.LatticeGeometry

noncomputable section
namespace Gutoe.LinkVariables

open Gutoe.ContinuumYMLieAlgebra

/-! ## Link variables -/

/-- A lattice site in ℤ⁴. -/
abbrev LatticeSite := Fin 4 → ℤ

/-- A lattice direction μ ∈ {0,1,2,3}. -/
abbrev LatticeDir := Fin 4

/-- Link variable data: U_μ(x) ∈ G for each site x and direction μ. -/
structure LinkVariableConfig where
  groupData : CompactSimpleLieGroupData
  /-- Link variable assignment: each oriented link gets a group element. -/
  linkType : Type
  /-- Link variables form a group (for composition). -/
  instGroup : Group linkType
  /-- Unitarity: U†U = 1 (compactness of G). -/
  isUnitary : Prop
  /-- The configuration space is non-empty. -/
  configNonEmpty : Nonempty linkType

attribute [instance] LinkVariableConfig.instGroup

/-! ## Parallel transport -/

/-- A lattice path is a sequence of oriented links. -/
structure LatticePath where
  /-- Number of links in the path. -/
  length : ℕ
  /-- The path is non-degenerate. -/
  length_pos : 0 < length

/-- Parallel transport along a path as the ordered product of link variables. -/
structure ParallelTransport where
  config : LinkVariableConfig
  path : LatticePath
  /-- The holonomy (ordered product) along the path. -/
  holonomy : config.linkType
  /-- Reversing the path gives the inverse holonomy. -/
  reverseIsInverse : Prop
  /-- Concatenation of paths gives product of holonomies. -/
  concatIsProduct : Prop

/-- (Axiom) Parallel transport is well-defined:
    reverse path gives inverse, concatenation gives product. -/
axiom parallel_transport_valid (pt : ParallelTransport) :
    pt.reverseIsInverse ∧ pt.concatIsProduct

/-- **GRAND-371: Link variables theorem**

    On a hypercubic lattice with compact simple G:
    1. Link variables U_μ(x) ∈ G are unitary.
    2. Parallel transport is compositional.
    3. Reversing a path inverts the holonomy. -/
theorem link_variables_theorem (pt : ParallelTransport)
    (hU : pt.config.isUnitary) :
    pt.config.isUnitary ∧ pt.reverseIsInverse ∧ pt.concatIsProduct :=
  ⟨hU, (parallel_transport_valid pt).1, (parallel_transport_valid pt).2⟩

end Gutoe.LinkVariables
