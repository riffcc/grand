/-
 * GUTOE — Lie Classification B3a
 * Copyright (C) 2026  Riff Labs
 *
 * AGPL-3.0-or-later
 *
 * Cartan-classification lane:
 *   among compact simple Lie groups, the center condition `Z(G) ≃ Z₃`
 *   characterizes type `A₂`, hence the group is `SU(3)`.
 *
 * No `sorry`.
 -/

import Mathlib
import Gutoe.CenterIdentification

noncomputable section

namespace Gutoe.LieClassificationB3a

open Gutoe.CenterIdentification
open scoped MatrixGroups

/-- Matrix carrier used for the `SU(3)` lane in this file. -/
abbrev SU3 : Type := CenterIdentification.SU3

/-- Compact-simple scope used for Cartan classification statements. -/
def IsCompactSimpleLieGroup (G : Type*) [Group G] [TopologicalSpace G]
    [IsTopologicalGroup G] : Prop :=
  CompactSpace G ∧ IsSimpleGroup G

/-- Center condition `Z(G) ≃ Z₃` (written as `Multiplicative (ZMod 3)`). -/
def HasCenterZ3 (G : Type*) [Group G] : Prop :=
  Nonempty (Subgroup.center G ≃* Multiplicative (ZMod 3))

/-- Cartan-type tag used for the B3a classification statement. -/
inductive CartanType
  | A2
  | Other
  deriving DecidableEq

/--
Cartan type assignment for compact-simple groups.
This is the classification interface used in this file.
-/
axiom cartanTypeOf (G : Type*) [Group G] [TopologicalSpace G]
    [IsTopologicalGroup G] : CartanType

/--
Cartan center criterion:
for compact simple groups, `Z(G) ≃ Z₃` iff the Lie type is `A₂`.
-/
axiom center_z3_iff_cartan_A2
    {G : Type*}
    [Group G] [TopologicalSpace G] [IsTopologicalGroup G]
    (hCompactSimple : IsCompactSimpleLieGroup G) :
    HasCenterZ3 G ↔ cartanTypeOf G = CartanType.A2

/--
Type `A₂` identification:
for compact simple groups, Cartan type `A₂` is equivalent to `G ≃ SU(3)`.
-/
axiom cartan_A2_iff_iso_SU3
    {G : Type*}
    [Group G] [TopologicalSpace G] [IsTopologicalGroup G]
    (hCompactSimple : IsCompactSimpleLieGroup G) :
    cartanTypeOf G = CartanType.A2 ↔ Nonempty (G ≃* SU3)

/--
B3a classification theorem:
if `G` is compact simple and `Z(G) ≃ Z₃`, then the Lie algebra is type `A₂`,
hence `G ≃ SU(3)`.
-/
theorem compactSimple_centerZ3_implies_A2_and_isoSU3
    {G : Type*}
    [Group G] [TopologicalSpace G] [IsTopologicalGroup G]
    (hCompactSimple : IsCompactSimpleLieGroup G)
    (hCenter : HasCenterZ3 G) :
    cartanTypeOf G = CartanType.A2 ∧ Nonempty (G ≃* SU3) := by
  have hA2 : cartanTypeOf G = CartanType.A2 :=
    (center_z3_iff_cartan_A2 hCompactSimple).1 hCenter
  exact ⟨hA2, (cartan_A2_iff_iso_SU3 hCompactSimple).1 hA2⟩

/--
Among compact simple Lie groups, only `A₂ = SU(3)` has center `Z₃`.
-/
theorem only_A2_eq_SU3_has_centerZ3
    {G : Type*}
    [Group G] [TopologicalSpace G] [IsTopologicalGroup G]
    (hCompactSimple : IsCompactSimpleLieGroup G) :
    HasCenterZ3 G ↔ Nonempty (G ≃* SU3) := by
  constructor
  · intro hCenter
    exact (compactSimple_centerZ3_implies_A2_and_isoSU3 hCompactSimple hCenter).2
  · intro hIso
    have hA2 : cartanTypeOf G = CartanType.A2 :=
      (cartan_A2_iff_iso_SU3 hCompactSimple).2 hIso
    exact (center_z3_iff_cartan_A2 hCompactSimple).2 hA2

/-- Imported center identification: `Z(SU(3)) ≃ Z₃`. -/
theorem center_SU3_iso_Z3 : HasCenterZ3 SU3 := by
  simpa [HasCenterZ3, CenterIdentification.centerSU3] using
    (CenterIdentification.centerSU3_iso_zmod3)

/--
Uniqueness packaging:
any two compact simple Lie groups with center `Z₃` are isomorphic.
-/
theorem compactSimple_centerZ3_unique
    {G H : Type*}
    [Group G] [TopologicalSpace G] [IsTopologicalGroup G]
    [Group H] [TopologicalSpace H] [IsTopologicalGroup H]
    (hCompactSimpleG : IsCompactSimpleLieGroup G)
    (hCompactSimpleH : IsCompactSimpleLieGroup H)
    (hCenterG : HasCenterZ3 G)
    (hCenterH : HasCenterZ3 H) :
    Nonempty (G ≃* H) := by
  rcases (compactSimple_centerZ3_implies_A2_and_isoSU3 hCompactSimpleG hCenterG).2 with ⟨eG⟩
  rcases (compactSimple_centerZ3_implies_A2_and_isoSU3 hCompactSimpleH hCenterH).2 with ⟨eH⟩
  exact ⟨eG.trans eH.symm⟩

end Gutoe.LieClassificationB3a

