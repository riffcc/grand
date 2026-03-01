import Mathlib
import Gutoe.DimensionalStructure
import Gutoe.RailSpace
import Gutoe.ProjectionFibers

/-!
GUTOE — Recursive Z₃ Tower Tests

This lane formalizes two test classes:

1) Upward recursion (one level): a lifted Z₃ action on a 256-state layer
   (`16 x 16`) and its order-3 closure.
2) Projection tower kernel profile: `256 -> 16 -> 4` with exact kernel ranks.

The goal is not a full dynamics claim; it is structural closure and rank accounting.
-/

namespace Gutoe.RecursiveZ3Tower

open Gutoe
open Gutoe.DimensionalStructure
open Gutoe.ProjectionFibers

noncomputable section

/-- Single recursive axis carrier (the 16-state Clifford layer index). -/
abbrev Axis16 := Fin 16

/-- One-level recursive index space (`16 x 16 = 256`). -/
abbrev Layer256 := Axis16 × Axis16

/-- One-level recursive vector space over the 256-state layer. -/
abbrev Vec256 := EuclideanSpace ℝ Layer256

/-- Distinguished zero index in the 16-state carrier. -/
def axis0 : Axis16 := ⟨0, by decide⟩

/-- Standard basis vectors of the 256-layer Euclidean space. -/
def layerBasisVec (p : Layer256) : Vec256 :=
  EuclideanSpace.basisFun Layer256 ℝ p

/-- `z3_4d` restricted/reindexed to the 16 Clifford basis states `{1..16}`. -/
def z3_16 (i : Axis16) : Axis16 :=
  match i.1 with
  | 0  => ⟨0, by decide⟩
  | 1  => ⟨1, by decide⟩
  | 2  => ⟨4, by decide⟩
  | 3  => ⟨5, by decide⟩
  | 4  => ⟨8, by decide⟩
  | 5  => ⟨9, by decide⟩
  | 6  => ⟨12, by decide⟩
  | 7  => ⟨13, by decide⟩
  | 8  => ⟨2, by decide⟩
  | 9  => ⟨3, by decide⟩
  | 10 => ⟨6, by decide⟩
  | 11 => ⟨7, by decide⟩
  | 12 => ⟨10, by decide⟩
  | 13 => ⟨11, by decide⟩
  | 14 => ⟨14, by decide⟩
  | 15 => ⟨15, by decide⟩
  | _  => ⟨0, by decide⟩

/-- The restricted 16-state Z₃ map has order 3. -/
theorem z3_16_order3 (i : Axis16) : z3_16 (z3_16 (z3_16 i)) = i := by
  fin_cases i <;> decide

/-- Lifted Z₃ action on one recursive layer (`16 x 16`). -/
def z3_256_index (p : Layer256) : Layer256 := (z3_16 p.1, z3_16 p.2)

/-- Lifted Z₃ still has order 3 on the recursive layer. -/
theorem z3_256_index_order3 (p : Layer256) :
    z3_256_index (z3_256_index (z3_256_index p)) = p := by
  rcases p with ⟨a, b⟩
  simp [z3_256_index, z3_16_order3]

/-- Upward/downward bridge: this restricted map exactly matches `z3_4d` on
states `1..16` after reindexing (`i = s-1`). -/
theorem z3_16_matches_z3_4d (i : Axis16) :
    (z3_16 i).1 + 1 = z3_4d (i.1 + 1) := by
  fin_cases i <;> decide

/-- Downward recursion witness already present in the base lane:
`grade1 = (fixed 1) + (orbit 3)`. -/
theorem downward_z3_split :
    (∃ s ∈ grade1_4d, z3_4d s = s) ∧
    (z3_4d 3 = 5 ∧ z3_4d 5 = 9 ∧ z3_4d 9 = 3) := by
  exact ⟨⟨2, by decide, z3_4d_gamma0_fixed⟩, z3_4d_quark_orbit⟩

-- ── Projection tower: 256 -> 16 -> 4 ─────────────────────────────────────────

/-- Projection `256 -> 16`: read one fixed second-index slice (`j = 0`). -/
def proj256to16 : Vec256 →ₗ[ℝ] Vec16 where
  toFun := fun v => ∑ i : Axis16, v (i, axis0) • railBasisVec i
  map_add' := by
    intro v w
    ext i
    simp [railBasisVec, add_smul, Finset.sum_add_distrib]
  map_smul' := by
    intro c v
    ext i
    simp [railBasisVec, Finset.smul_sum, mul_smul]

/-- Section witnessing surjectivity of `proj256to16`. -/
def section16to256 (w : Vec16) : Vec256 :=
  ∑ i : Axis16, w i • layerBasisVec (i, axis0)

/-- `proj256to16` is surjective. -/
theorem proj256to16_surjective : Function.Surjective proj256to16 := by
  intro w
  refine ⟨section16to256 w, ?_⟩
  have hsum : (∑ i : Axis16, w i • EuclideanSpace.single i (1 : ℝ)) = w := by
    simpa using
      (OrthonormalBasis.sum_repr (EuclideanSpace.basisFun Axis16 ℝ) w)
  simpa [proj256to16, section16to256, railBasisVec, layerBasisVec, axis0] using hsum

/-- Dimension of the recursive layer vector space. -/
theorem vec256_finrank : Module.finrank ℝ Vec256 = 256 := by
  simp [Vec256, finrank_euclideanSpace, Fintype.card_prod]

/-- Kernel rank for the first projection in the tower: `256 -> 16` gives `240`. -/
theorem proj256to16_kernel_finrank :
    Module.finrank ℝ (LinearMap.ker proj256to16) = 240 := by
  have hsum :
      Module.finrank ℝ (LinearMap.range proj256to16) +
        Module.finrank ℝ (LinearMap.ker proj256to16) =
      Module.finrank ℝ Vec256 := LinearMap.finrank_range_add_finrank_ker proj256to16
  have hrange :
      Module.finrank ℝ (LinearMap.range proj256to16) = 16 := by
    rw [LinearMap.range_eq_top.2 proj256to16_surjective]
    simpa using (vec16_dim : Module.finrank ℝ Vec16 = 16)
  have hdom : Module.finrank ℝ Vec256 = 256 := vec256_finrank
  omega

/-- Tower kernel profile lock:
`ker(256->16)=240`, `ker(16->4)=12`, total hidden dimensions `252`. -/
theorem projection_tower_kernel_profile :
    Module.finrank ℝ (LinearMap.ker proj256to16) = 240 ∧
    Module.finrank ℝ (LinearMap.ker grade1Projection) = 12 ∧
    240 + 12 = (252 : ℕ) := by
  refine ⟨proj256to16_kernel_finrank, grade1Projection_kernel_finrank, by decide⟩

/-- Index-level commutation: first projection of the lifted Z₃ layer equals
apply-Z₃-after-first-projection. -/
theorem proj_index_commutes_with_lifted_z3 (p : Layer256) :
    (z3_256_index p).1 = z3_16 p.1 := by
  rfl

end
end Gutoe.RecursiveZ3Tower
