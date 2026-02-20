/-
 * GUTOE - 16D Vector Rail Space
 * Copyright (C) 2026  Riff Labs
 *
 * AGPL-3.0-or-later
 *
 * Ported and refined from VoidEmergence.RailHilbert
 * (fuck-it-lets-do-a-universe-2026).
 * Original had two sorrys (Cauchy-Schwarz, triangle inequality).
 * Both are now proven using Mathlib's EuclideanSpace instances.
 -/

import Mathlib
import Gutoe.Basic

open scoped InnerProductSpace

noncomputable section

namespace Gutoe

/-!
# 16D Vector Rail Space

Formalises the Hilbert space structure of the GUTOE vector rail space.

Each of the 16 HexStates (arising from four binary dimensions) corresponds
to a basis direction in ℝ¹⁶. The rail space is therefore `EuclideanSpace ℝ (Fin 16)`,
which inherits its full Hilbert space structure from Mathlib.

The four TriState components embed as the first four rail directions:
- VOID    → zero vector (no rail)
- COSINE  → rail direction e₀ (cosine wave axis)
- SINE    → rail direction e₁ (sine wave axis)
- TANGENT → rail direction e₂ (slope/ratio axis)

Key theorems (all PROVEN):
1. Vec16 is a real Hilbert space (complete inner product space of dimension 16)
2. Cauchy-Schwarz: |⟨v, w⟩| ≤ ‖v‖ · ‖w‖
3. Triangle inequality: ‖v + w‖ ≤ ‖v‖ + ‖w‖
4. Phase evolution is unitary (preserves inner product)
5. Constructive interference: ‖av + bv‖ = |a + b| · ‖v‖
6. TriState states embed as orthonormal rail directions
7. VOID maps to the zero rail (no direction, no coherence)
-/

-- ── Vec16: the 16-dimensional rail space ───────────────────────────────────────

/-- 16-dimensional vector rail space, the inner state of the full GUTOE geometry.
    Each of the 16 HexState directions maps to one basis direction. -/
abbrev Vec16 := EuclideanSpace ℝ (Fin 16)

/-!
All Hilbert space structure is inherited from Mathlib's `EuclideanSpace`:
- `AddCommGroup Vec16`                                   ← free
- `Module ℝ Vec16`                                       ← free
- `InnerProductSpace ℝ Vec16`                            ← free
- `FiniteDimensional ℝ Vec16`  (dim = 16)               ← free
- `CompleteSpace Vec16`         (Cauchy sequences converge) ← free
-/

/-- Vec16 has dimension 16 over ℝ — REAL -/
theorem vec16_dim : Module.finrank ℝ Vec16 = 16 :=
  finrank_euclideanSpace_fin

/-- Cauchy-Schwarz for Vec16: |⟪v, w⟫| ≤ ‖v‖ · ‖w‖ — REAL -/
theorem vec16_cauchy_schwarz (v w : Vec16) :
    abs ⟪v, w⟫_ℝ ≤ ‖v‖ * ‖w‖ :=
  abs_real_inner_le_norm v w

/-- Triangle inequality for Vec16: ‖v + w‖ ≤ ‖v‖ + ‖w‖ — REAL -/
theorem vec16_triangle (v w : Vec16) :
    ‖v + w‖ ≤ ‖v‖ + ‖w‖ :=
  norm_add_le v w

/-- Inner product is symmetric — REAL -/
theorem vec16_inner_symm (v w : Vec16) : ⟪v, w⟫_ℝ = ⟪w, v⟫_ℝ :=
  real_inner_comm w v

/-- Inner product is positive-definite — REAL -/
theorem vec16_inner_pos (v : Vec16) (h : v ≠ 0) : 0 < ⟪v, v⟫_ℝ :=
  real_inner_self_pos.mpr h

-- ── Standard basis vectors ─────────────────────────────────────────────────────

/-- The i-th standard basis vector of Vec16 -/
def railBasisVec (i : Fin 16) : Vec16 :=
  EuclideanSpace.basisFun (Fin 16) ℝ i

/-- Basis vectors are unit vectors — REAL -/
theorem railBasisVec_norm (i : Fin 16) : ‖railBasisVec i‖ = 1 :=
  (EuclideanSpace.basisFun (Fin 16) ℝ).norm_eq_one i

/-- Distinct basis vectors are orthogonal (Kronecker delta) — REAL -/
theorem railBasisVec_inner (i j : Fin 16) :
    ⟪railBasisVec i, railBasisVec j⟫_ℝ = if i = j then 1 else 0 :=
  (EuclideanSpace.basisFun (Fin 16) ℝ).inner_eq_ite i j

-- ── VectorRail structure ───────────────────────────────────────────────────────

/-- A vector rail: a direction in 16D space with wave properties.
    The direction lives in Vec16; amplitude and phase are wave scalars. -/
structure VectorRail where
  direction : Vec16   -- which of the 16 basis directions (or superposition)
  amplitude : ℝ       -- wave amplitude; energy ∝ amplitude²
  phase     : ℝ       -- wave phase ∈ [0, 2π)

/-- Extract the direction from a rail -/
def railDir (r : VectorRail) : Vec16 := r.direction

-- ── Phase evolution ────────────────────────────────────────────────────────────

/-- Phase evolution: advance the phase by 2π·frequency·dt, keeping direction/amplitude.
    This models free wave propagation along a fixed rail. -/
def phaseEvolve (dt freq : ℝ) (r : VectorRail) : VectorRail :=
  { r with phase := (r.phase + 2 * Real.pi * freq * dt) % (2 * Real.pi) }

/-- Phase evolution does not change the rail direction — REAL -/
theorem phaseEvolve_direction (dt freq : ℝ) (r : VectorRail) :
    railDir (phaseEvolve dt freq r) = railDir r := by
  simp [phaseEvolve, railDir]

/-- Phase evolution preserves the inner product (it is unitary) — REAL -/
theorem phaseEvolve_unitary (dt freq : ℝ) (r1 r2 : VectorRail) :
    ⟪railDir (phaseEvolve dt freq r1), railDir (phaseEvolve dt freq r2)⟫_ℝ
    = ⟪railDir r1, railDir r2⟫_ℝ := by
  simp [phaseEvolve, railDir]

/-- Phase evolution preserves amplitude (energy conservation) — REAL -/
theorem phaseEvolve_amplitude (dt freq : ℝ) (r : VectorRail) :
    (phaseEvolve dt freq r).amplitude = r.amplitude := by
  simp [phaseEvolve]

-- ── Superposition ──────────────────────────────────────────────────────────────

/-- Superpose two rails: weighted sum of directions by amplitude -/
def railSuperpose (r1 r2 : VectorRail) : Vec16 :=
  r1.amplitude • r1.direction + r2.amplitude • r2.direction

/-- Constructive interference: when two rails point the same way,
    their amplitudes add linearly and the norm scales accordingly — REAL -/
theorem constructive_interference (r1 r2 : VectorRail)
    (h : r1.direction = r2.direction) :
    ‖railSuperpose r1 r2‖ = abs (r1.amplitude + r2.amplitude) * ‖r1.direction‖ := by
  unfold railSuperpose
  rw [← h, ← add_smul, norm_smul, Real.norm_eq_abs]

/-- Destructive interference: when amplitudes are equal and opposite on the
    same direction, the result is zero — REAL -/
theorem destructive_interference (r1 r2 : VectorRail)
    (hdir : r1.direction = r2.direction)
    (hamp : r1.amplitude = -r2.amplitude) :
    railSuperpose r1 r2 = 0 := by
  unfold railSuperpose
  rw [← hdir, ← add_smul, hamp, neg_add_cancel, zero_smul]

-- ── Master Hilbert space theorem ───────────────────────────────────────────────

/-- Vec16 satisfies all Hilbert space axioms — REAL -/
theorem vec16_is_hilbert_space :
    -- (1) Complete inner product space of dimension 16
    Module.finrank ℝ Vec16 = 16 ∧
    -- (2) Cauchy-Schwarz inequality
    (∀ v w : Vec16, abs ⟪v, w⟫_ℝ ≤ ‖v‖ * ‖w‖) ∧
    -- (3) Orthonormal basis of 16 vectors exists
    (∃ b : Fin 16 → Vec16, ∀ i j,
      ⟪b i, b j⟫_ℝ = if i = j then 1 else 0) ∧
    -- (4) The space is complete (Cauchy sequences converge)
    CompleteSpace Vec16 :=
  ⟨vec16_dim, fun v w => vec16_cauchy_schwarz v w,
   ⟨railBasisVec, railBasisVec_inner⟩, inferInstance⟩

-- ── TriState embedding ─────────────────────────────────────────────────────────

/-!
### TriState as the first four rail directions

The four TriState wave components map into Vec16:

| TriState | Rail | Physical role                         |
|----------|------|---------------------------------------|
| VOID     | 0    | no rail direction (zero vector)       |
| COSINE   | e₀   | cosine wave axis (amplitude extremum) |
| SINE     | e₁   | sine wave axis   (zero-crossing)      |
| TANGENT  | e₂   | slope axis       (sin/cos ratio)      |

COSINE and SINE are adjacent in the hexagonal lattice with angle 60°.
In Vec16 they map to *orthogonal* rail directions (e₀ ⊥ e₁).
The veracity(SINE, COSINE) = √3/2 ≠ 0 comes from hexagonal geometry, not
from rail-space inner product (which gives 0).  The two measures are different:
veracity is a spectral overlap in the hexagonal tiling; the rail inner product
is the standard Euclidean dot product between direction vectors.
-/

/-- Embed a TriState into the 16D rail space -/
def triStateToRail : TriState → Vec16
  | TriState.VOID    => 0
  | TriState.COSINE  => railBasisVec ⟨0, by norm_num⟩
  | TriState.SINE    => railBasisVec ⟨1, by norm_num⟩
  | TriState.TANGENT => railBasisVec ⟨2, by norm_num⟩

/-- VOID maps to the zero rail — REAL -/
theorem triState_void_is_zero : triStateToRail TriState.VOID = 0 := rfl

/-- COSINE, SINE, TANGENT are unit rail vectors — REAL -/
theorem triState_basis_rail_norm (s : TriState) (hs : s ≠ TriState.VOID) :
    ‖triStateToRail s‖ = 1 := by
  cases s with
  | VOID    => exact absurd rfl hs
  | COSINE  => exact railBasisVec_norm ⟨0, by norm_num⟩
  | SINE    => exact railBasisVec_norm ⟨1, by norm_num⟩
  | TANGENT => exact railBasisVec_norm ⟨2, by norm_num⟩

/-- COSINE and SINE are orthogonal as rail directions — REAL -/
theorem triState_cosine_sine_ortho :
    ⟪triStateToRail TriState.COSINE, triStateToRail TriState.SINE⟫_ℝ = 0 := by
  simp [triStateToRail, railBasisVec_inner]

/-- COSINE and TANGENT are orthogonal as rail directions — REAL -/
theorem triState_cosine_tangent_ortho :
    ⟪triStateToRail TriState.COSINE, triStateToRail TriState.TANGENT⟫_ℝ = 0 := by
  simp [triStateToRail, railBasisVec_inner]

/-- SINE and TANGENT are orthogonal as rail directions — REAL -/
theorem triState_sine_tangent_ortho :
    ⟪triStateToRail TriState.SINE, triStateToRail TriState.TANGENT⟫_ℝ = 0 := by
  simp [triStateToRail, railBasisVec_inner]

/-- The non-VOID TriState directions are mutually orthonormal in Vec16 — REAL -/
theorem triState_embed_orthonormal :
    ∀ a b : TriState, a ≠ TriState.VOID → b ≠ TriState.VOID →
    ⟪triStateToRail a, triStateToRail b⟫_ℝ = if a = b then 1 else 0 := by
  intro a b ha hb
  cases a
  · exact absurd rfl ha
  all_goals (cases b; · exact absurd rfl hb)
  all_goals simp only [triStateToRail, railBasisVec_inner, Fin.ext_iff, ↓reduceIte]
  all_goals norm_num
  all_goals decide

/-- The TriState embedding is injective on non-VOID states — REAL -/
theorem triState_embed_injective :
    ∀ a b : TriState, a ≠ TriState.VOID → b ≠ TriState.VOID →
    triStateToRail a = triStateToRail b → a = b := by
  intro a b ha hb heq
  by_contra hne
  have h1 := triState_embed_orthonormal a b ha hb
  rw [if_neg hne, ← heq] at h1
  have h2 := triState_embed_orthonormal a a ha ha
  rw [if_pos rfl] at h2
  linarith

end Gutoe

end
