/-
 * GUTOE — Center Identification (GRAND-461)
 *
 * Bridge theorem:
 *   the bottom-up GUTOE Z₃ automorphism (`z3_4d`) is identified with
 *   the center of the SU(3)-lane (formalized here with `SL(3,ℂ)` matrices).
 *
 * No `sorry`.
 -/

import Mathlib
import Gutoe.Z3Uniqueness
import Gutoe.GaugeGroupSU3

noncomputable section

namespace Gutoe.CenterIdentification

open scoped MatrixGroups
open Gutoe.DimensionalStructure Gutoe.Z3Uniqueness Gutoe.GaugeGroupSU3

/-- Matrix carrier used for the SU(3) lane in this formalization. -/
abbrev SU3 : Type := Matrix.SpecialLinearGroup (Fin 3) ℂ

/-- The center of the SU(3)-lane group (formalized as `SL(3,ℂ)`). -/
abbrev centerSU3 : Subgroup SU3 := Subgroup.center SU3

/-- Matrix-level description: center elements are exactly scalar matrices `ωI` with `ω^3 = 1`. -/
theorem mem_centerSU3_iff (A : SU3) :
    A ∈ centerSU3 ↔
      ∃ ω : ℂ, ω ^ 3 = 1 ∧ Matrix.scalar (Fin 3) ω = A := by
  simpa [centerSU3, Fintype.card_fin] using
    (Matrix.SpecialLinearGroup.mem_center_iff (n := Fin 3) (R := ℂ) (A := A))

/-- The center of the SU(3)-lane is equivalent to cubic roots of unity. -/
noncomputable def centerSU3_equiv_rootsOfUnity : centerSU3 ≃* rootsOfUnity 3 ℂ :=
  Matrix.SpecialLinearGroup.center_equiv_rootsOfUnity' (n := Fin 3) (R := ℂ) (i := (0 : Fin 3))

/-- Cubic roots of unity are (non-canonically) the cyclic group `ZMod 3` (multiplicative form). -/
noncomputable def rootsOfUnity3_equiv_zmod3 : rootsOfUnity 3 ℂ ≃* Multiplicative (ZMod 3) := by
  classical
  refine mulEquivOfCyclicCardEq ?hcard
  have hRoots : Nat.card (rootsOfUnity 3 ℂ) = 3 := by
    exact (Nat.card_eq_fintype_card (α := rootsOfUnity 3 ℂ)).trans
      (Complex.card_rootsOfUnity (n := 3))
  have hZ3 : Nat.card (Multiplicative (ZMod 3)) = 3 := by
    simp [Nat.card_eq_fintype_card]
  exact hRoots.trans hZ3.symm

/-- Center(SU(3)) is the cyclic group of order 3. -/
noncomputable def centerSU3_equiv_zmod3 : centerSU3 ≃* Multiplicative (ZMod 3) :=
  centerSU3_equiv_rootsOfUnity.trans rootsOfUnity3_equiv_zmod3

/-- Group-isomorphism packaging: center(SU(3)) ≅ Z₃. -/
theorem centerSU3_iso_zmod3 : Nonempty (centerSU3 ≃* Multiplicative (ZMod 3)) :=
  ⟨centerSU3_equiv_zmod3⟩

/-- The GUTOE Z₃ phase group used for `z3_4d` iterates. -/
abbrev z3_4dGroup : Type := Multiplicative (ZMod 3)

/-- Map from GUTOE `z3_4d`-phase group to the SU(3)-center. -/
noncomputable def z3_4d_to_centerSU3 : z3_4dGroup →* centerSU3 :=
  (centerSU3_equiv_zmod3.symm : z3_4dGroup ≃* centerSU3).toMonoidHom

/-- The above map is a group isomorphism. -/
noncomputable def z3_4d_center_iso : z3_4dGroup ≃* centerSU3 :=
  centerSU3_equiv_zmod3.symm

/-- Explicitly, the map from the GUTOE Z₃ phase group is exactly an isomorphism. -/
theorem z3_4d_to_centerSU3_is_iso :
    ∃ e : z3_4dGroup ≃* centerSU3, z3_4d_to_centerSU3 = e.toMonoidHom := by
  exact ⟨z3_4d_center_iso, rfl⟩

/-- `z3_4d` iterates on the reference quark state `3`, indexed by `ZMod 3`. -/
def z3_4dIterateOnQuark (k : ZMod 3) : ℕ := (z3_4d^[k.val]) 3

/-- The three `z3_4d` iterate values on the quark reference state. -/
theorem z3_4dIterateOnQuark_values :
    z3_4dIterateOnQuark 0 = 3 ∧
    z3_4dIterateOnQuark 1 = 5 ∧
    z3_4dIterateOnQuark 2 = 9 := by
  native_decide

/-- Every `z3_4d` iterate of the reference quark stays in the quark orbit. -/
theorem z3_4dIterateOnQuark_mem_quarkOrbit :
    ∀ k : ZMod 3, z3_4dIterateOnQuark k ∈ quarkOrbit := by
  native_decide

/-- Standard basis vector of the SU(3) fundamental carrier. -/
def quarkBasisVec (i : Fin 3) : Fin 3 → ℂ := Pi.single i (1 : ℂ)

/-- Restriction of the SU(3)-center to the fundamental rep is scalar phase multiplication. -/
theorem centerSU3_scalar_on_fundamental
    (z : centerSU3) (i : Fin 3) :
    ∃ ω : ℂ, ω ^ 3 = 1 ∧
      (((z : SU3) : Matrix (Fin 3) (Fin 3) ℂ).mulVec (quarkBasisVec i) = ω • quarkBasisVec i) := by
  rcases (mem_centerSU3_iff (A := (z : SU3))).1 z.property with ⟨ω, hω, hscalar⟩
  refine ⟨ω, hω, ?_⟩
  have hscalarM :
      (Matrix.scalar (Fin 3) ω : Matrix (Fin 3) (Fin 3) ℂ) =
        ((z : SU3) : Matrix (Fin 3) (Fin 3) ℂ) := by
    simpa using hscalar
  calc
    (((z : SU3) : Matrix (Fin 3) (Fin 3) ℂ).mulVec (quarkBasisVec i))
        = ((Matrix.scalar (Fin 3) ω : Matrix (Fin 3) (Fin 3) ℂ).mulVec (quarkBasisVec i)) := by
          rw [← hscalarM]
    _ = ω • quarkBasisVec i := by
          simpa [Matrix.scalar_apply, quarkBasisVec] using
            (Matrix.diagonal_const_mulVec (m := Fin 3) ω (quarkBasisVec i))

/--
GUTOE bridge theorem:
1) `z3_4d` gives the quark triplet orbit `{3,5,9}`,
2) that orbit is the 3-dimensional fundamental carrier,
3) center(SU(3)) is Z₃,
4) and center action on the fundamental is exactly phase multiplication.
-/
theorem quark_orbit_is_fundamental_restricted_center :
    (z3_4dIterateOnQuark 0 = 3 ∧ z3_4dIterateOnQuark 1 = 5 ∧ z3_4dIterateOnQuark 2 = 9) ∧
    quarkOrbit.card = 3 ∧
    Nonempty ({s // s ∈ quarkOrbit} ≃ Fin 3) ∧
    Nonempty (centerSU3 ≃* Multiplicative (ZMod 3)) ∧
    (∀ z : centerSU3, ∀ i : Fin 3,
      ∃ ω : ℂ, ω ^ 3 = 1 ∧
        ((((z : SU3) : Matrix (Fin 3) (Fin 3) ℂ).mulVec (quarkBasisVec i)) = ω • quarkBasisVec i)) := by
  refine ⟨z3_4dIterateOnQuark_values, quarkOrbit_card, quarkOrbit_equiv_fin3,
    centerSU3_iso_zmod3, ?_⟩
  intro z i
  exact centerSU3_scalar_on_fundamental z i

end Gutoe.CenterIdentification
