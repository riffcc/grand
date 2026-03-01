# 194 — Cl(1,3) Projection Fibers: Non-injective Map and 12D Kernel

## Scope
Formalize a testable core of the `#50` hypothesis:
- finite-state projection fibers from Cl(1,3) basis-state encoding,
- linear projection `Vec16 -> R^4` with kernel dimension theorem.

## Lean module
- `lean/Gutoe/ProjectionFibers.lean`
- Added to build roots in `lean/lakefile.lean` as `Gutoe.ProjectionFibers`.

## New formal content

### A) Finite-state fiber structure (basis-state encoding lane)
Definitions:
- `basisMask : ℕ -> ℕ` with `mi = s - 1`
- `cl13ToMinkowskiAxis : ℕ -> Fin 4`

Theorems:
- `projection_has_nontrivial_fiber`
- `projection_injective_on_grade1`
- `axis_one_fiber_grade1_plus_higher_grade`
- `cl13_to_minkowski_projection_fiber_structure`

Meaning:
- Projection is non-injective globally (explicit fiber collision witness).
- Projection remains injective on grade-1 subset `{2,3,5,9}`.

### B) Linear map kernel theorem (Vec16 -> Minkowski 4-vector)
Definitions:
- `railIndex4 : Fin 4 -> Fin 16`
- `grade1Projection : Vec16 ->ₗ[ℝ] (Fin 4 -> ℝ)`
- `grade1Section : (Fin 4 -> ℝ) -> Vec16`

Theorems:
- `grade1Projection_surjective`
- `grade1Projection_range_top`
- `grade1Projection_kernel_finrank` (`finrank ker = 12`)
- `grade1Projection_not_injective`

Meaning:
- The explicit 16->4 projection is surjective.
- Kernel has dimension exactly 12.
- Map is provably non-injective.

## Verification
- `lake build Gutoe.ProjectionFibers` ✅
- `lake build Gutoe` ✅ (`8140 jobs`)
- No errors; warning-only lint hints.

## Honest boundary
This proves fiber multiplicity/non-injectivity and kernel dimension claims.
It does **not** yet prove bounded-length connected paths between fibers over distant spacetime points (connection/dynamics/geodesic question remains open).
