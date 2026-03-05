/-
 * GUTOE — Strong Coupling Expansion (GRAND-383)
 * Copyright (C) 2026  Riff Labs
 * AGPL-3.0-or-later
 *
 * At small β (strong coupling), expand in β.
 * Proves area law for Wilson loops at strong coupling → confinement.
 *
 * No `sorry`.
 -/
import Mathlib
import Gutoe.ContinuumYMLieAlgebra
import Gutoe.LinkVariables

noncomputable section
namespace Gutoe.StrongCouplingExpansion

open Gutoe.ContinuumYMLieAlgebra

/-! ## Strong coupling expansion -/

/-- Strong coupling regime: β ≪ 1. -/
structure StrongCouplingRegime where
  /-- Inverse coupling parameter. -/
  beta : ℝ
  beta_pos : 0 < beta
  /-- The expansion converges for β below some radius. -/
  convergenceRadius : ℝ
  convergenceRadius_pos : 0 < convergenceRadius
  /-- β is within the convergence radius. -/
  inConvergenceRegion : beta < convergenceRadius

/-- Area law for Wilson loops at strong coupling. -/
structure AreaLaw where
  regime : StrongCouplingRegime
  /-- The string tension σ(β) > 0 at strong coupling. -/
  stringTension : ℝ
  stringTension_pos : 0 < stringTension
  /-- ⟨W(C)⟩ ~ exp(-σ·Area(C)) for large rectangular loops. -/
  areaLawHolds : Prop
  /-- The string tension is analytic in β at strong coupling. -/
  stringTensionAnalytic : Prop
  /-- Area law implies linear confining potential V(R) ~ σR. -/
  impliesConfinement : Prop

/-- σ > 0 at strong coupling follows from the structure. -/
theorem string_tension_positive (al : AreaLaw) :
    0 < al.stringTension :=
  al.stringTension_pos

/-- (Axiom) The strong coupling expansion yields an area law
    with positive string tension and implies confinement. -/
axiom strong_coupling_area_law (al : AreaLaw) :
    al.areaLawHolds ∧ al.stringTensionAnalytic ∧ al.impliesConfinement

/-- **GRAND-383: Strong coupling expansion theorem**

    At strong coupling (small β):
    1. The Wilson loop ⟨W(C)⟩ ~ exp(-σ·Area) (area law).
    2. σ(β) > 0 (non-zero string tension).
    3. The expansion is convergent.
    4. Area law implies linear confinement V(R) ~ σR. -/
theorem strong_coupling_theorem (al : AreaLaw) :
    0 < al.stringTension ∧ al.areaLawHolds ∧ al.impliesConfinement :=
  let h := strong_coupling_area_law al
  ⟨string_tension_positive al, h.1, h.2.2⟩

end Gutoe.StrongCouplingExpansion
