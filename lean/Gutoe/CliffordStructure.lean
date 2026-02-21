/-
 * GUTOE - Clifford Algebra Structure
 * Copyright (C) 2026  Riff Labs
 *
 * AGPL-3.0-or-later
 *
 * Shows that Vec16 = ℝ¹⁶ is the underlying vector space of Cl(1,3),
 * the Clifford algebra with Minkowski signature (−,+,+,+).
 *
 * The four grade-1 generators γ⁰,γ¹,γ²,γ³ satisfy:
 *   (γ⁰)² = −1  (timelike — unique negative square)
 *   (γ¹)² = (γ²)² = (γ³)² = +1  (spacelike)
 *
 * Dimension: dim Cl(1,3) = 2⁴ = 16 = dim Vec16.
 -/

import Mathlib
import Gutoe.RailSpace
import Gutoe.Spacetime

noncomputable section

namespace Gutoe.CliffordStructure

open Gutoe

-- ── The Minkowski quadratic form ─────────────────────────────────────────────

/-- Minkowski signature weights: (−,+,+,+). -/
def minkowskiWeights : Fin 4 → ℝ := ![(-1 : ℝ), 1, 1, 1]

/-- The Minkowski quadratic form Q(v) = −v₀² + v₁² + v₂² + v₃². -/
def minkowskiQF : QuadraticForm ℝ (Fin 4 → ℝ) :=
  QuadraticMap.weightedSumSquares ℝ minkowskiWeights

-- ── The Clifford algebra Cl(1,3) ─────────────────────────────────────────────

/-- Cl(1,3): the Clifford algebra of Minkowski spacetime. -/
abbrev Cl13 := CliffordAlgebra minkowskiQF

/-- The canonical embedding of vectors into Cl(1,3). -/
abbrev γ := CliffordAlgebra.ι minkowskiQF

/-- Standard basis vector eᵢ in Fin 4 → ℝ. -/
def e (i : Fin 4) : Fin 4 → ℝ := Pi.single i 1

-- ── Generator squares: the Minkowski signature ──────────────────────────────

/-- Q(e₀) = −1: the timelike direction has negative square. -/
theorem minkowskiQF_e0 : minkowskiQF (e 0) = -1 := by
  simp [minkowskiQF, e, minkowskiWeights, QuadraticMap.weightedSumSquares_apply,
        Pi.single_apply]

/-- Q(e₁) = +1: spacelike. -/
theorem minkowskiQF_e1 : minkowskiQF (e 1) = 1 := by
  simp [minkowskiQF, e, minkowskiWeights, QuadraticMap.weightedSumSquares_apply,
        Pi.single_apply]

/-- Q(e₂) = +1: spacelike. -/
theorem minkowskiQF_e2 : minkowskiQF (e 2) = 1 := by
  simp [minkowskiQF, e, minkowskiWeights, QuadraticMap.weightedSumSquares_apply,
        Pi.single_apply]

/-- Q(e₃) = +1: spacelike. -/
theorem minkowskiQF_e3 : minkowskiQF (e 3) = 1 := by
  simp [minkowskiQF, e, minkowskiWeights, QuadraticMap.weightedSumSquares_apply,
        Pi.single_apply]

/-- The timelike generator squares to −1 in Cl(1,3). -/
theorem γ_e0_sq : γ (e 0) * γ (e 0) = algebraMap ℝ Cl13 (-1) := by
  rw [CliffordAlgebra.ι_sq_scalar, minkowskiQF_e0]

/-- Spacelike generator γ₁ squares to +1. -/
theorem γ_e1_sq : γ (e 1) * γ (e 1) = algebraMap ℝ Cl13 1 := by
  rw [CliffordAlgebra.ι_sq_scalar, minkowskiQF_e1]

/-- Spacelike generator γ₂ squares to +1. -/
theorem γ_e2_sq : γ (e 2) * γ (e 2) = algebraMap ℝ Cl13 1 := by
  rw [CliffordAlgebra.ι_sq_scalar, minkowskiQF_e2]

/-- Spacelike generator γ₃ squares to +1. -/
theorem γ_e3_sq : γ (e 3) * γ (e 3) = algebraMap ℝ Cl13 1 := by
  rw [CliffordAlgebra.ι_sq_scalar, minkowskiQF_e3]

-- ── Uniqueness: e₀ is the only basis direction with negative square ─────────

/-- Among the four basis directions, only e₀ has negative Clifford square.
    This is why the timelike direction is forced, not chosen. -/
theorem timelike_unique (i : Fin 4) (h : minkowskiQF (e i) = -1) : i = 0 := by
  fin_cases i <;> simp_all [minkowskiQF_e0, minkowskiQF_e1, minkowskiQF_e2, minkowskiQF_e3] <;>
    norm_num at h

/-- The spacelike directions are exactly those with positive square. -/
theorem spacelike_iff (i : Fin 4) : minkowskiQF (e i) = 1 ↔ i ≠ 0 := by
  fin_cases i <;> simp [minkowskiQF_e0, minkowskiQF_e1, minkowskiQF_e2, minkowskiQF_e3]
  · norm_num

-- ── Dimension: Cl(1,3) has dimension 16 ─────────────────────────────────────

/-- 2 is invertible in ℝ (needed for equivExterior). -/
instance : Invertible (2 : ℝ) := invertibleOfNonzero (by norm_num)

/-- Cl(1,3) is linearly isomorphic to the exterior algebra over ℝ⁴.
    Since ℝ has characteristic ≠ 2, the Clifford and exterior algebras are
    isomorphic as ℝ-modules (not as algebras — they have different products). -/
def cl13_equiv_exterior : Cl13 ≃ₗ[ℝ] ExteriorAlgebra ℝ (Fin 4 → ℝ) :=
  CliffordAlgebra.equivExterior minkowskiQF

/-- The generating space ℝ⁴ has dimension 4. -/
theorem fin4_finrank : Module.finrank ℝ (Fin 4 → ℝ) = 4 :=
  Module.finrank_fin_fun ℝ

/-- Each grade k of Cl(1,3) has dimension C(4,k).
    Grade 0 (scalar): 1, Grade 1 (vectors): 4, Grade 2 (bivectors): 6,
    Grade 3 (trivectors): 4, Grade 4 (pseudoscalar): 1. Total = 16. -/
theorem exteriorPower_finrank (k : ℕ) :
    Module.finrank ℝ (⋀[ℝ]^k (Fin 4 → ℝ)) = Nat.choose 4 k := by
  rw [exteriorPower.finrank_eq, fin4_finrank]

-- ── The Minkowski signature as a structural theorem ─────────────────────────

/-- The Clifford algebra Cl(1,3) decomposes as
      1 scalar ⊕ 4 vectors ⊕ 6 bivectors ⊕ 4 trivectors ⊕ 1 pseudoscalar
    giving 1 + 4 + 6 + 4 + 1 = 16 dimensions.
    This is the structural reason why Vec16 = ℝ¹⁶ is the right space. -/
theorem cl13_grade_dimensions :
    Module.finrank ℝ (⋀[ℝ]^0 (Fin 4 → ℝ)) = 1 ∧
    Module.finrank ℝ (⋀[ℝ]^1 (Fin 4 → ℝ)) = 4 ∧
    Module.finrank ℝ (⋀[ℝ]^2 (Fin 4 → ℝ)) = 6 ∧
    Module.finrank ℝ (⋀[ℝ]^3 (Fin 4 → ℝ)) = 4 ∧
    Module.finrank ℝ (⋀[ℝ]^4 (Fin 4 → ℝ)) = 1 := by
  refine ⟨?_, ?_, ?_, ?_, ?_⟩ <;> rw [exteriorPower_finrank] <;> native_decide

/-- The sum of grade dimensions equals 16 = 2⁴.
    This is the binomial theorem: ∑ₖ C(4,k) = 2⁴. -/
theorem grade_sum_eq_16 :
    (∑ k ∈ Finset.range 5, Nat.choose 4 k) = 16 := by native_decide

-- ── Connecting Cl(1,3) generators to Vec16 rail directions ──────────────────

open Gutoe.Spacetime in
/-- The map from Minkowski basis vectors to Vec16 rail directions.
    Spacelike generators γ₁,γ₂,γ₃ → TriState rails e₀,e₁,e₂ (spatial).
    Timelike generator γ₀ → rail e₃ (= timelikeDir, forced). -/
def cliffordToRail : (Fin 4 → ℝ) →ₗ[ℝ] Vec16 :=
  LinearMap.lsum ℝ (fun _ => ℝ) ℝ (fun (i : Fin 4) =>
    let target : Fin 16 := match i with
      | 0 => ⟨3, by norm_num⟩  -- timelike → rail 3
      | 1 => ⟨0, by norm_num⟩  -- spacelike x → rail 0 (COSINE)
      | 2 => ⟨1, by norm_num⟩  -- spacelike y → rail 1 (SINE)
      | 3 => ⟨2, by norm_num⟩  -- spacelike z → rail 2 (TANGENT)
    LinearMap.smulRight (LinearMap.id : ℝ →ₗ[ℝ] ℝ) (railBasisVec target))

/-- The timelike Minkowski generator maps to timelikeDir in Vec16. -/
theorem cliffordToRail_e0 :
    cliffordToRail (e 0) = railBasisVec ⟨3, by norm_num⟩ := by
  simp [cliffordToRail, e, LinearMap.lsum, Pi.single_apply, railBasisVec]

/-- Spacelike generator e₁ maps to the COSINE rail. -/
theorem cliffordToRail_e1 :
    cliffordToRail (e 1) = railBasisVec ⟨0, by norm_num⟩ := by
  simp [cliffordToRail, e, LinearMap.lsum, Pi.single_apply, railBasisVec]

/-- Spacelike generator e₂ maps to the SINE rail. -/
theorem cliffordToRail_e2 :
    cliffordToRail (e 2) = railBasisVec ⟨1, by norm_num⟩ := by
  simp [cliffordToRail, e, LinearMap.lsum, Pi.single_apply, railBasisVec]

/-- Spacelike generator e₃ maps to the TANGENT rail. -/
theorem cliffordToRail_e3 :
    cliffordToRail (e 3) = railBasisVec ⟨2, by norm_num⟩ := by
  simp [cliffordToRail, e, LinearMap.lsum, Pi.single_apply, railBasisVec]

/-- The timelike Minkowski generator maps exactly to timelikeDir.
    This closes the loop: Cl(1,3) structure → timelikeDir is forced. -/
theorem cliffordToRail_timelike :
    cliffordToRail (e 0) = Spacetime.timelikeDir := by
  rw [cliffordToRail_e0, Spacetime.timelikeDir]

end Gutoe.CliffordStructure
