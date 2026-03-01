import Mathlib
import Gutoe.DimensionalStructure
import Gutoe.RailSpace

/-!
GUTOE — Projection Fiber Structure (Cl(1,3) -> Minkowski axis carrier)

We model a concrete projection from the finite Cl(1,3) basis-state encoding
(`s ∈ {1..16}`, with mask `mi = s-1`) to a 4-axis Minkowski carrier label.

Interpretation:
- The map keeps only the first-present generator bit (timelike-first ordering).
- Higher-grade basis states therefore collapse onto grade-1 carriers.

This gives a testable fiber statement for the "many-to-one projection" idea:
the projection is non-injective on full basis states but injective on the
grade-1 subset `{2,3,5,9}`.
-/

namespace Gutoe.ProjectionFibers

open Gutoe.DimensionalStructure

noncomputable section

/-- Basis mask used throughout the finite-state Cl(1,3) encoding:
`mi = s - 1` for Clifford basis states `s ∈ {1..16}`. -/
def basisMask (s : ℕ) : ℕ := s - 1

/-- Projection from a Cl(1,3) basis-state label to one Minkowski axis carrier.
Priority is bit0, then bit1, then bit2, else bit3.

This is a finite surrogate for "project full multivector state onto one
spacetime carrier axis." -/
def cl13ToMinkowskiAxis (s : ℕ) : Fin 4 :=
  if Nat.testBit (basisMask s) 0 then ⟨0, by decide⟩
  else if Nat.testBit (basisMask s) 1 then ⟨1, by decide⟩
  else if Nat.testBit (basisMask s) 2 then ⟨2, by decide⟩
  else ⟨3, by decide⟩

/-- Grade-1 generators map to distinct carriers (time + three spatial axes). -/
theorem axis_of_grade1_generators :
    cl13ToMinkowskiAxis 2 = ⟨0, by decide⟩ ∧
    cl13ToMinkowskiAxis 3 = ⟨1, by decide⟩ ∧
    cl13ToMinkowskiAxis 5 = ⟨2, by decide⟩ ∧
    cl13ToMinkowskiAxis 9 = ⟨3, by decide⟩ := by
  native_decide

/-- Nontrivial fiber witness on valid Clifford basis states:
`3` (grade-1) and `7` (higher-grade) project to the same Minkowski carrier. -/
theorem projection_has_nontrivial_fiber :
    ∃ a b : ℕ,
      1 ≤ a ∧ a ≤ 16 ∧ 1 ≤ b ∧ b ≤ 16 ∧ a ≠ b ∧
      cl13ToMinkowskiAxis a = cl13ToMinkowskiAxis b := by
  refine ⟨3, 7, by decide, by decide, by decide, by decide, by decide, ?_⟩
  native_decide

/-- On the grade-1 subset, the projection is injective. -/
theorem projection_injective_on_grade1 :
    ∀ a ∈ grade1_4d, ∀ b ∈ grade1_4d,
      cl13ToMinkowskiAxis a = cl13ToMinkowskiAxis b → a = b := by
  intro a ha b hb h
  fin_cases ha <;> fin_cases hb <;> simp [axis_of_grade1_generators] at h ⊢

/-- Explicit mixed-grade fiber witness for axis-1:
one grade-1 state and one non-grade-1 state share the same projected axis. -/
theorem axis_one_fiber_grade1_plus_higher_grade :
    3 ∈ grade1_4d ∧ 7 ∉ grade1_4d ∧
    cl13ToMinkowskiAxis 3 = cl13ToMinkowskiAxis 7 := by
  refine ⟨by decide, by decide, ?_⟩
  native_decide

/-- Master fiber-structure statement for the finite Cl(1,3) basis projection:
non-injective globally, injective on grade-1, with explicit mixed-grade fiber. -/
theorem cl13_to_minkowski_projection_fiber_structure :
    (∃ a b : ℕ,
      1 ≤ a ∧ a ≤ 16 ∧ 1 ≤ b ∧ b ≤ 16 ∧ a ≠ b ∧
      cl13ToMinkowskiAxis a = cl13ToMinkowskiAxis b) ∧
    (∀ a ∈ grade1_4d, ∀ b ∈ grade1_4d,
      cl13ToMinkowskiAxis a = cl13ToMinkowskiAxis b → a = b) ∧
    (3 ∈ grade1_4d ∧ 7 ∉ grade1_4d ∧
      cl13ToMinkowskiAxis 3 = cl13ToMinkowskiAxis 7) := by
  exact ⟨projection_has_nontrivial_fiber,
    projection_injective_on_grade1,
    axis_one_fiber_grade1_plus_higher_grade⟩

-- ── Linear-algebra projection lane (Vec16 -> Minkowski 4-vector) ────────────

/-- Embed `Fin 4` coordinate indices into the first four rail indices of `Vec16`. -/
def railIndex4 (i : Fin 4) : Fin 16 :=
  ⟨i.1, Nat.lt_trans i.2 (by decide)⟩

/-- Grade-1 coordinate projection from the 16D rail state to a Minkowski 4-vector.
It keeps the first four rail coordinates and discards the remaining twelve. -/
def grade1Projection : Vec16 →ₗ[ℝ] (Fin 4 → ℝ) where
  toFun := fun v i => v (railIndex4 i)
  map_add' := by
    intro v w
    funext i
    simp [railIndex4]
  map_smul' := by
    intro c v
    funext i
    simp [railIndex4]

/-- A canonical section (right-inverse candidate) for `grade1Projection`:
fill first four coordinates from `w`, set the remaining twelve to zero. -/
def grade1Section (w : Fin 4 → ℝ) : Vec16 :=
  ∑ i : Fin 4, w i • railBasisVec (railIndex4 i)

/-- `grade1Projection` is surjective via the explicit section. -/
theorem grade1Projection_surjective : Function.Surjective grade1Projection := by
  classical
  intro w
  refine ⟨grade1Section w, ?_⟩
  funext i
  have hsum :
      (∑ x : Fin 4, (if (↑i : ℕ) = ↑x then w x else 0)) = w i := by
    have hsum' : (∑ x : Fin 4, (if i = x then w x else 0)) = w i := by
      simp
    simpa [Fin.ext_iff] using hsum'
  simpa [grade1Projection, grade1Section, railBasisVec, railIndex4] using hsum

/-- Range of `grade1Projection` is all of Minkowski 4-space. -/
theorem grade1Projection_range_top : LinearMap.range grade1Projection = ⊤ := by
  exact LinearMap.range_eq_top.2 grade1Projection_surjective

/-- The kernel of the 16→4 grade-1 projection has dimension 12. -/
theorem grade1Projection_kernel_finrank :
    Module.finrank ℝ (LinearMap.ker grade1Projection) = 12 := by
  have hsum :
      Module.finrank ℝ (LinearMap.range grade1Projection) +
        Module.finrank ℝ (LinearMap.ker grade1Projection) =
      Module.finrank ℝ Vec16 := LinearMap.finrank_range_add_finrank_ker grade1Projection
  have hrange :
      Module.finrank ℝ (LinearMap.range grade1Projection) = 4 := by
    rw [grade1Projection_range_top]
    simpa using (finrank_euclideanSpace_fin : Module.finrank ℝ (Fin 4 → ℝ) = 4)
  have hdom : Module.finrank ℝ Vec16 = 16 := vec16_dim
  omega

/-- Structural corollary: this projection cannot be injective (kernel has positive dimension). -/
theorem grade1Projection_not_injective : ¬ Function.Injective grade1Projection := by
  intro hinj
  have hker : LinearMap.ker grade1Projection = ⊥ := by
    exact LinearMap.ker_eq_bot.mpr hinj
  have hker0 : Module.finrank ℝ (LinearMap.ker grade1Projection) = 0 := by
    rw [hker]
    simp
  have hker12 : Module.finrank ℝ (LinearMap.ker grade1Projection) = 12 :=
    grade1Projection_kernel_finrank
  omega

/-- Fiber over base point `x` under the grade-1 projection. -/
def fiberAt (x : Fin 4 → ℝ) : Set Vec16 := { v | grade1Projection v = x }

/-- Coordinate-wise norm bound in `Vec16`: every coordinate magnitude is bounded by
the ambient Euclidean norm. -/
theorem abs_coord_le_norm (u : Vec16) (j : Fin 16) : |u j| ≤ ‖u‖ := by
  have hsum : ‖u j‖ ^ 2 ≤ ∑ i : Fin 16, ‖u i‖ ^ 2 := by
    refine Finset.single_le_sum (f := fun i : Fin 16 => ‖u i‖ ^ 2) ?hnonneg ?hmem
    · intro i hi
      positivity
    · simp
  have hsq : (‖u j‖) ^ 2 ≤ (Real.sqrt (∑ i : Fin 16, ‖u i‖ ^ 2)) ^ 2 := by
    nlinarith [hsum, Real.sq_sqrt (by positivity : 0 ≤ ∑ i : Fin 16, ‖u i‖ ^ 2)]
  have hroot : ‖u j‖ ≤ Real.sqrt (∑ i : Fin 16, ‖u i‖ ^ 2) := by
    have habs : |‖u j‖| ≤ |Real.sqrt (∑ i : Fin 16, ‖u i‖ ^ 2)| := (sq_le_sq).1 hsq
    simpa [abs_of_nonneg (norm_nonneg _), abs_of_nonneg (Real.sqrt_nonneg _)] using habs
  simpa [EuclideanSpace.norm_eq, Real.norm_eq_abs] using hroot

/-- Fiber-membership pins the projected axis-0 coordinate exactly. -/
theorem fiber_coord0 (x : Fin 4 → ℝ) {v : Vec16} (hv : v ∈ fiberAt x) :
    v (railIndex4 0) = x 0 := by
  change (grade1Projection v) 0 = x 0
  simpa [fiberAt] using congrArg (fun f : Fin 4 → ℝ => f 0) hv

/-- Any connector between fibers over `x` and `y` has norm at least the axis-0
base separation. -/
theorem fiber_connector_lower_bound_axis0
    (x y : Fin 4 → ℝ) {v w : Vec16}
    (hv : v ∈ fiberAt x) (hw : w ∈ fiberAt y) :
    |y 0 - x 0| ≤ ‖w - v‖ := by
  have hv0 : v (railIndex4 0) = x 0 := fiber_coord0 x hv
  have hw0 : w (railIndex4 0) = y 0 := fiber_coord0 y hw
  have hcoord : (w - v) (railIndex4 0) = y 0 - x 0 := by
    simp [hw0, hv0]
  have hnorm : |(w - v) (railIndex4 0)| ≤ ‖w - v‖ :=
    abs_coord_le_norm (w - v) (railIndex4 0)
  simpa [hcoord] using hnorm

/-- No uniform base-independent bound exists for connector lengths between all
fibers in this model. -/
theorem no_uniform_fiber_connector_bound :
    ∀ B : ℝ, ∃ x y : Fin 4 → ℝ,
      ∀ v ∈ fiberAt x, ∀ w ∈ fiberAt y, B < ‖w - v‖ := by
  intro B
  let x : Fin 4 → ℝ := fun _ => 0
  let y : Fin 4 → ℝ := fun i => if i = 0 then |B| + 1 else 0
  refine ⟨x, y, ?_⟩
  intro v hv w hw
  have hlow : |y 0 - x 0| ≤ ‖w - v‖ := fiber_connector_lower_bound_axis0 x y hv hw
  have hnonneg : 0 ≤ |B| + 1 := by positivity
  have hsep : |y 0 - x 0| = |B| + 1 := by
    simp [x, y, hnonneg]
  have hB : B < |B| + 1 := by
    have habs : B ≤ |B| := le_abs_self B
    linarith
  rw [hsep] at hlow
  linarith

end
end Gutoe.ProjectionFibers
