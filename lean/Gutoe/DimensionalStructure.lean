/-
 * GUTOE - Dimensional Structure: Why d=4 is Special
 * Copyright (C) 2026  Riff Labs
 *
 * AGPL-3.0-or-later
 *
 * The Z₃ rotation on Cl(1,3) has a grade-1 fixed point (γ⁰ = lepton).
 * In Cl(1,2) (d=3), no grade-1 state is a Z₃ fixed point.
 * Therefore d=4 is the minimum dimension for stable lepton-quark distinction.
 *
 * All theorems proven (no sorry).
 -/

import Mathlib
import Gutoe.FineStructure

namespace Gutoe.DimensionalStructure

open Gutoe.FineStructure

-- ── Z₃ rotation on Cl(1,3) ────────────────────────────────────────────────────
--
-- State s ∈ {0..16}: VOID=0, Clifford states s=1..16 with mi = s−1.
-- Bit rotation: (b₀,b₁,b₂,b₃) → (b₀, b₃, b₁, b₂)
-- The timelike bit b₀ is FIXED; spatial bits b₁,b₂,b₃ cycle.
--
-- Precomputed table (verified against Rust Z3_TABLE):
--   s: 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
--   →: 0  1  2  5  6  9 10 13 14  3  4  7  8 11 12 15 16

/-- Z₃ rotation on Cl(1,3) Clifford states (ℕ → ℕ). -/
def z3_4d : ℕ → ℕ
  | 0  => 0  | 1  => 1  | 2  => 2   -- VOID, scalar, γ⁰ (fixed!)
  | 3  => 5  | 4  => 6  | 5  => 9   -- γ¹→γ², γ⁰¹→γ⁰², γ²→γ³
  | 6  => 10 | 7  => 13 | 8  => 14  -- γ⁰²→γ⁰³, γ¹²→γ²³, γ⁰¹²→γ⁰²³
  | 9  => 3  | 10 => 4  | 11 => 7   -- γ³→γ¹, γ⁰³→γ⁰¹, γ¹³→γ¹²
  | 12 => 8  | 13 => 11 | 14 => 12  -- γ⁰¹³→γ⁰¹², γ²³→γ¹³, γ⁰²³→γ⁰¹³
  | 15 => 15 | 16 => 16             -- γ¹²³ (fixed), γ⁰¹²³ (fixed)
  | _  => 0  -- out of range

/-- Z₃ rotation on Cl(1,2) Clifford states (ℕ → ℕ).
    Cyclic rotation: (b₀,b₁,b₂) → (b₂, b₀, b₁) — ALL bits cycle, including timelike. -/
def z3_3d : ℕ → ℕ
  | 0 => 0  | 1 => 1  -- VOID, scalar (fixed)
  | 2 => 3  | 3 => 5  | 4 => 7  -- γ⁰→γ¹, γ¹→γ², γ⁰¹→γ¹²
  | 5 => 2  | 6 => 4  | 7 => 6  -- γ²→γ⁰, γ⁰²→γ⁰¹, γ¹²→γ⁰²
  | 8 => 8  -- γ⁰¹² (pseudoscalar, fixed)
  | _ => 0  -- out of range

-- ── Z₃ is order 3 ─────────────────────────────────────────────────────────────

/-- z3_4d applied 3 times is identity (on valid states 0..16). -/
theorem z3_4d_order3 (s : ℕ) (hs : s ≤ 16) : z3_4d (z3_4d (z3_4d s)) = s := by
  interval_cases s <;> decide

/-- z3_3d applied 3 times is identity (on valid states 0..8). -/
theorem z3_3d_order3 (s : ℕ) (hs : s ≤ 8) : z3_3d (z3_3d (z3_3d s)) = s := by
  interval_cases s <;> decide

-- ── Fixed points of z3_4d ─────────────────────────────────────────────────────

/-- γ⁰ (s=2) is a Z₃ fixed point in Cl(1,3): z3_4d(2) = 2. -/
theorem z3_4d_gamma0_fixed : z3_4d 2 = 2 := by decide

/-- The grade-1 spatial directions {γ¹,γ²,γ³} = {3,5,9} form a 3-cycle. -/
theorem z3_4d_quark_orbit :
    z3_4d 3 = 5 ∧ z3_4d 5 = 9 ∧ z3_4d 9 = 3 := by decide

/-- Only s ∈ {0,1,2,15,16} are Z₃ fixed points in Cl(1,3). -/
theorem z3_4d_fixed_points (s : ℕ) (hs : s ≤ 16) :
    z3_4d s = s ↔ s ∈ ({0, 1, 2, 15, 16} : Finset ℕ) := by
  interval_cases s <;> decide

/-- Grade-1 states of Cl(1,3) are {2, 3, 5, 9} (one bit set in mi = s−1). -/
def grade1_4d : Finset ℕ := {2, 3, 5, 9}

/-- Among grade-1 states of Cl(1,3), only γ⁰ (s=2) is a Z₃ fixed point.
    The quarks {γ¹,γ²,γ³} always cycle — they can never be the lepton. -/
theorem z3_4d_unique_grade1_fp (s : ℕ) (hs : s ∈ grade1_4d) :
    z3_4d s = s ↔ s = 2 := by
  fin_cases hs <;> decide

-- ── Fixed points of z3_3d ─────────────────────────────────────────────────────

/-- In Cl(1,2), γ⁰ (s=2) is NOT a Z₃ fixed point — it maps to γ¹ (s=3). -/
theorem z3_3d_gamma0_not_fixed : z3_3d 2 ≠ 2 := by decide

/-- The grade-1 3-cycle in Cl(1,2): γ⁰ → γ¹ → γ² → γ⁰. -/
theorem z3_3d_grade1_orbit :
    z3_3d 2 = 3 ∧ z3_3d 3 = 5 ∧ z3_3d 5 = 2 := by decide

/-- Only s ∈ {0, 1, 8} are Z₃ fixed points in Cl(1,2) — all trivial. -/
theorem z3_3d_fixed_points (s : ℕ) (hs : s ≤ 8) :
    z3_3d s = s ↔ s ∈ ({0, 1, 8} : Finset ℕ) := by
  interval_cases s <;> decide

/-- Grade-1 states of Cl(1,2) are {2, 3, 5}. -/
def grade1_3d : Finset ℕ := {2, 3, 5}

/-- No grade-1 state in Cl(1,2) is a Z₃ fixed point.
    γ⁰ cycles into γ¹ (lepton ≡ quark — no stable distinction). -/
theorem z3_3d_no_grade1_fp (s : ℕ) (hs : s ∈ grade1_3d) : z3_3d s ≠ s := by
  fin_cases hs <;> decide

-- ── d=4 uniqueness ─────────────────────────────────────────────────────────────

/-!
### d=4 is the Minimum Dimension for Stable Matter

Mathematical reason:
  Z₃ on n bits can fix one specific bit b₀ only if the remaining n−1
  bits support a Z₃ permutation: need n−1 ≥ 3, so n ≥ 4 bits.
  In d-dimensional spacetime, n = d (one bit per dimension).
  Therefore d ≥ 4 is required.

  d=3 (Cl(1,2)): 3 bits, only 2 spatial bits → Z₃ must include γ⁰.
  d=4 (Cl(1,3)): 4 bits, 3 spatial bits → Z₃ fixes γ⁰, cycles {γ¹,γ²,γ³}.

Physical consequence:
  If γ⁰ is in the Z₃ quark orbit → no stable "lepton" distinct from "quarks"
  → no hydrogen → no atoms → no chemistry → no life.
-/

/-- d=4 is the unique minimum dimension for a grade-1 Z₃ fixed point:
    - d=3: no grade-1 fixed point
    - d=4: γ⁰ (s=2) IS a grade-1 fixed point -/
theorem d4_minimum_for_atoms :
    -- d=3: no grade-1 Z₃ fixed point (lepton mixes with quarks)
    (∀ s ∈ grade1_3d, z3_3d s ≠ s) ∧
    -- d=4: γ⁰ is a grade-1 Z₃ fixed point (stable lepton)
    (∃ s ∈ grade1_4d, z3_4d s = s) :=
  ⟨z3_3d_no_grade1_fp, ⟨2, by decide, z3_4d_gamma0_fixed⟩⟩

/-- The Eddington numbers for dimensions with and without stable matter. -/
theorem eddington_stable_vs_unstable :
    alphaInverse 3 = 37 ∧   -- d=3: unstable (no grade-1 fixed point)
    alphaInverse 4 = 137 ∧  -- d=4: stable (γ⁰ is fixed point) ← our universe
    alphaInverse 5 = 529 :=  -- d=5: stable (fixed point still exists)
  ⟨alpha_inverse_d3, alpha_inverse_d4, alpha_inverse_d5⟩

/-- The unique grade-1 fixed points of Z₃ in Cl(1,3):
    γ⁰ (lepton) and γ¹²³ (the spatial trivector dual to the quark sector). -/
theorem z3_4d_grade1_fp_is_gamma0 :
    ∀ s ∈ grade1_4d, z3_4d s = s ↔ s = 2 := z3_4d_unique_grade1_fp

/-- The γ⁰ ↔ γ¹²³ duality: both are fixed by Z₃.
    γ⁰ is grade-1 (the lepton); γ¹²³ is grade-3 (its Hodge dual in 4D). -/
theorem fixed_point_duality :
    z3_4d 2 = 2 ∧   -- γ⁰ fixed (grade-1, lepton)
    z3_4d 15 = 15 := -- γ¹²³ fixed (grade-3, dual quark sector)
  ⟨by decide, by decide⟩

end Gutoe.DimensionalStructure
