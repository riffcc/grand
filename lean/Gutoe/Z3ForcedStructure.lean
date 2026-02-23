/-
 * GUTOE — Z₃ Forced Structure
 * Copyright (C) 2026  Riff Labs
 * AGPL-3.0-or-later
 *
 * THEOREM (z3_forces_hermitian_circulant): any 3×3 Hermitian matrix
 * commuting with the cyclic Z₃ permutation is a Hermitian circulant —
 * parametrised by exactly one diagonal entry (a : ℝ) and one off-diagonal
 * entry (ε : ℂ).  Zero free structure.
 *
 * This is the structural keystone: the Z₃ vacuum (PerturbativeSymmetry)
 * forces a mass matrix whose eigenvalues are exactly the Hermitian
 * circulant eigenvalues that give Koide = 2/3 (LeptonMass).
 *
 * All theorems no sorry.
 -/

import Mathlib
import Gutoe.LeptonMass

namespace Gutoe.Z3ForcedStructure

-- ══════════════════════════════════════════════════════════════════════════════
-- Definitions
-- ══════════════════════════════════════════════════════════════════════════════

/-- The Z₃ cyclic permutation matrix: generation 0→1→2→0. -/
def z3Perm : Matrix (Fin 3) (Fin 3) ℂ :=
  !![0, 1, 0; 0, 0, 1; 1, 0, 0]

/-- The Hermitian circulant mass matrix with real diagonal a and complex
    off-diagonal ε.  Rows are cyclic shifts of ⟨a, ε, ε★⟩. -/
noncomputable def instantonMassMatrix (a : ℝ) (ε : ℂ) : Matrix (Fin 3) (Fin 3) ℂ :=
  !![↑a, ε, star ε; star ε, ↑a, ε; ε, star ε, ↑a]

-- ══════════════════════════════════════════════════════════════════════════════
-- Hermitian property
-- ══════════════════════════════════════════════════════════════════════════════

/-- instantonMassMatrix is Hermitian. -/
theorem instMat_hermitian (a : ℝ) (ε : ℂ) :
    (instantonMassMatrix a ε).IsHermitian := by
  apply Matrix.IsHermitian.ext
  intro i j
  fin_cases i <;> fin_cases j <;>
    simp [instantonMassMatrix, Complex.conj_ofReal]

-- ══════════════════════════════════════════════════════════════════════════════
-- Z₃ Forcing Theorem
-- ══════════════════════════════════════════════════════════════════════════════

/-- **Z₃ Forcing Theorem**: any 3×3 Hermitian matrix commuting with the Z₃
    cyclic permutation equals instantonMassMatrix (M 0 0).re (M 0 1).

    Physical meaning: Z₃ symmetry + Hermiticity force the mass matrix to be
    a Hermitian circulant — two real parameters, no further freedom. -/
theorem z3_forces_hermitian_circulant
    (M : Matrix (Fin 3) (Fin 3) ℂ)
    (hH : M.IsHermitian)
    (hZ : z3Perm * M = M * z3Perm) :
    M = instantonMassMatrix (M 0 0).re (M 0 1) := by
  -- Diagonal entries are real: M i i = ↑(M i i).re
  have dr : ∀ i : Fin 3, M i i = ((M i i).re : ℂ) := fun i => by
    apply Complex.ext
    · simp
    · have hself : star (M i i) = M i i := hH.apply i i
      have him   : (star (M i i)).im = -(M i i).im := Complex.conj_im (M i i)
      simp only [Complex.ofReal_im]
      linarith [him.symm.trans (congrArg Complex.im hself)]
  -- Helper: extract one (z3Perm * M)[i,j] = (M * z3Perm)[i,j] at concrete indices,
  -- then fully reduce via simp (concrete indices → !! entries reduce to 0 or 1).
  let ze := fun i j => congrFun (congrFun hZ i) j
  -- Six commutativity equalities we need (each proved by reducing the matrix products)
  have z00 : M 1 0 = M 0 2 := by
    have h := ze 0 0
    simp [Matrix.mul_apply, z3Perm, Fin.sum_univ_three] at h; exact h
  have z01 : M 1 1 = M 0 0 := by
    have h := ze 0 1
    simp [Matrix.mul_apply, z3Perm, Fin.sum_univ_three] at h; exact h
  have z02 : M 1 2 = M 0 1 := by
    have h := ze 0 2
    simp [Matrix.mul_apply, z3Perm, Fin.sum_univ_three] at h; exact h
  have z10 : M 2 0 = M 1 2 := by
    have h := ze 1 0
    simp [Matrix.mul_apply, z3Perm, Fin.sum_univ_three] at h; exact h
  have z11 : M 2 1 = M 1 0 := by
    have h := ze 1 1
    simp [Matrix.mul_apply, z3Perm, Fin.sum_univ_three] at h; exact h
  have z12 : M 2 2 = M 1 1 := by
    have h := ze 1 2
    simp [Matrix.mul_apply, z3Perm, Fin.sum_univ_three] at h; exact h
  -- Prove M = instantonMassMatrix entry by entry.
  -- Each `show` relies on definitional reduction of !! matrix access at concrete indices.
  apply Matrix.ext; intro i j
  fin_cases i <;> fin_cases j
  · show M 0 0 = ((M 0 0).re : ℂ);         exact dr 0
  · show M 0 1 = M 0 1;                     rfl
  · show M 0 2 = star (M 0 1);              exact ((hH.apply 1 0).trans z00).symm
  · show M 1 0 = star (M 0 1);              exact (hH.apply 1 0).symm
  · show M 1 1 = ((M 0 0).re : ℂ);         exact z01.trans (dr 0)
  · show M 1 2 = M 0 1;                     exact z02
  · show M 2 0 = M 0 1;                     exact z10.trans z02
  · show M 2 1 = star (M 0 1);             exact z11.trans (hH.apply 1 0).symm
  · show M 2 2 = ((M 0 0).re : ℂ);         exact (z12.trans z01).trans (dr 0)

end Gutoe.Z3ForcedStructure
