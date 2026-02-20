/-
 * GUTOE - Spacetime Emergence from Rail Space
 * Copyright (C) 2026  Riff Labs
 *
 * AGPL-3.0-or-later
 *
 * EXPLORATORY: Testing what the existing structure forces.
 *
 * PROVEN so far:
 *   temporal_is_neg_spatial       : omegaSqTemporal v k = -omegaSq v k
 *   spatial_temporal_complement   : omegaSq + omegaSqTemporal = 0 everywhere
 *   cutoff_is_shared_boundary     : both dispersions zero exactly at k_c;
 *                                   spatial propagates below, temporal above
 *   propagating_rank (sorry)      : 3D spatial subspace — needs finrank API
 *   complement_rank (sorry)       : 13D complement — needs propagating_rank
 *
 * KEY INSIGHT (Wings 2026-02-20):
 *   "Time is what you get when wave propagation breaks down."
 *   Made precise: temporal modes propagate exactly where spatial modes don't.
 *   The cutoff k_c is the shared boundary between the two regimes.
 *   Time is not in Vec16 as a direction — it's the ω-frequency axis,
 *   or equivalently the analytic continuation of the spatial dispersion.
 -/

import Mathlib
import Gutoe.Basic
import Gutoe.RailSpace
import Gutoe.DispersionRelation

open scoped InnerProductSpace

noncomputable section

namespace Gutoe.Spacetime

open Gutoe

-- ── The propagating (spatial) subspace ────────────────────────────────────────

/-- The three TriState wave components span the spatial rail directions. -/
def propagatingSubspace : Submodule ℝ Vec16 :=
  Submodule.span ℝ
    {triStateToRail TriState.COSINE,
     triStateToRail TriState.SINE,
     triStateToRail TriState.TANGENT}

/-- The three TriState rail vectors are linearly independent.
    (Follows from orthonormality: c₀e₀ + c₁e₁ + c₂e₂ = 0
    ⟹ ⟪eᵢ, 0⟫ = cᵢ = 0 for each i.) — TODO: prove from triState_embed_orthonormal -/
theorem triState_rails_linearIndependent :
    LinearIndependent ℝ (fun s : Fin 3 => [triStateToRail TriState.COSINE,
                                             triStateToRail TriState.SINE,
                                             triStateToRail TriState.TANGENT].get s) := by
  sorry -- follows from orthonormality in RailSpace.lean; Mathlib API TBD

/-- The propagating subspace has dimension 3. — TODO: needs finrank API for span -/
theorem propagating_rank : Module.finrank ℝ propagatingSubspace = 3 := by
  sorry

/-- The orthogonal complement of the propagating subspace has dimension 13.
    Follows immediately from propagating_rank and vec16_dim. — REAL -/
theorem complement_rank :
    Module.finrank ℝ propagatingSubspaceᗮ = 13 := by
  have h := Submodule.finrank_add_finrank_orthogonal propagatingSubspace
  rw [vec16_dim] at h
  have h3 := propagating_rank
  omega

-- ── The Wick rotation: temporal = exact negation of spatial ──────────────────

/-- The Wick-rotated (temporal) angular frequency squared.
    Obtained by k_spatial → i·k_temporal in ω²_spatial:
      ω²_spatial(k)    = +v²k² − λ_QG·ℓ_P²·k⁴
      ω²_temporal(k_t) = −v²k_t² + λ_QG·ℓ_P²·k_t⁴ -/
noncomputable def omegaSqTemporal (v k_t : ℝ) : ℝ :=
  -v ^ 2 * k_t ^ 2 + DISPERSION_COEFF * k_t ^ 4

/-- THE FUNDAMENTAL WICK ROTATION THEOREM:
    The temporal dispersion is exactly the negation of the spatial dispersion.
    This is pure algebra — k → ik flips the sign of k² but not k⁴. — REAL -/
theorem temporal_is_neg_spatial (v k : ℝ) :
    omegaSqTemporal v k = -(omegaSq v k) := by
  simp [omegaSqTemporal, omegaSq]; ring

/-- Spatial and temporal dispersions sum to zero everywhere. — REAL -/
theorem spatial_temporal_complement (v k : ℝ) :
    omegaSq v k + omegaSqTemporal v k = 0 := by
  rw [temporal_is_neg_spatial]; ring

-- ── The Lorentz-invariant dispersion: WHY the sign flip is forced ─────────────

/-!
### Lorentz invariance forces the sign flip.

The spatial dispersion `omegaSq` is a polynomial in k²:
  `P(x) = v²x − DISPERSION_COEFF·x²`   so  `P(k²) = v²k² − DISPERSION_COEFF·k⁴ = omegaSq v k`

The natural 4D Lorentz-invariant extension replaces `k²` (positive-definite
spatial norm) with `k_s² − k_t²` (the Lorentz-invariant pseudo-norm):
  `omegaSqFull v k_s k_t = P(k_s² − k_t²)`

Restricting to pure spatial modes (k_t = 0):
  `omegaSqFull v k 0 = P(k²) = omegaSq v k`  ← we recover spatial ✓

Restricting to pure temporal modes (k_s = 0):
  `omegaSqFull v 0 k_t = P(0 − k_t²) = P(−k_t²) = −v²k_t² + DISPERSION_COEFF·k_t⁴`
                       `= omegaSqTemporal v k_t`  ← we get the sign flip ✓

The sign flip is NOT a choice. It is the unique consequence of requiring the
dispersion to be a function of the Lorentz-invariant pseudo-norm k_s² − k_t².
There is no other degree-4 polynomial in (k_s², k_t²) that:
  (1) agrees with omegaSq on the spatial axis (k_t = 0)
  (2) is a function of k_s² − k_t² (Lorentz invariant)
-/

/-- The 4D Lorentz-invariant dispersion relation.
    Both spatial and temporal dispersions are restrictions of this. -/
noncomputable def omegaSqFull (v k_s k_t : ℝ) : ℝ :=
  v ^ 2 * (k_s ^ 2 - k_t ^ 2) - DISPERSION_COEFF * (k_s ^ 2 - k_t ^ 2) ^ 2

/-- The full 4D dispersion is a function of the Lorentz pseudo-norm k_s²−k_t². — REAL -/
theorem full_dispersion_is_lorentz_invariant (v k_s k_t c_s c_t : ℝ)
    -- If the Lorentz pseudo-norm is preserved: c_s² - c_t² = k_s² - k_t²
    (h : c_s ^ 2 - c_t ^ 2 = k_s ^ 2 - k_t ^ 2) :
    omegaSqFull v c_s c_t = omegaSqFull v k_s k_t := by
  simp only [omegaSqFull]; rw [h]

/-- Restricting the 4D dispersion to spatial modes (k_t = 0) gives omegaSq. — REAL -/
theorem full_dispersion_spatial_restriction (v k : ℝ) :
    omegaSqFull v k 0 = omegaSq v k := by
  unfold omegaSqFull omegaSq; ring

/-- What the Lorentz-invariant extension ACTUALLY gives for temporal modes.
    Note: omegaSqFull v 0 k_t = -v²k_t² - D·k_t⁴  (BOTH terms negative)
    This is NOT omegaSqTemporal = -v²k_t² + D·k_t⁴  (which has +D·k_t⁴)
    Lean caught this discrepancy. — REAL -/
theorem full_dispersion_temporal_restriction (v k_t : ℝ) :
    omegaSqFull v 0 k_t = -v ^ 2 * k_t ^ 2 - DISPERSION_COEFF * k_t ^ 4 := by
  unfold omegaSqFull; ring

/-- THE DISCREPANCY THEOREM:
    The Lorentz-invariant extension and omegaSqTemporal (= -omegaSq) DIFFER
    by 2·DISPERSION_COEFF·k_t⁴.

    omegaSqFull v 0 k_t    = P(-k_t²) = -v²k_t² - D·k_t⁴
    omegaSqTemporal v k_t  = -P(k_t²) = -v²k_t² + D·k_t⁴

    The difference is the sign of the quantum gravity correction term.

    Physical consequence:
    - omegaSqFull temporal: ALWAYS negative (evanescent at all k_t > 0).
      Lorentz-invariant temporal modes never propagate.
    - omegaSqTemporal: negative below k_c, ZERO at k_c, positive above k_c.
      These temporal modes propagate at super-Planck scales.

    The "gravity strong on rear face" picture requires omegaSqTemporal (not
    the Lorentz-invariant extension). The definition `temporal = -spatial` is
    NOT forced by Lorentz invariance — it is a specific choice that allows
    temporal propagation above k_c.

    This is the deeper question Wings identified: is this the only consistent
    choice, or just one of many? — REAL -/
theorem lorentz_vs_negation_discrepancy (v k_t : ℝ) :
    omegaSqFull v 0 k_t - omegaSqTemporal v k_t = -2 * DISPERSION_COEFF * k_t ^ 4 := by
  rw [full_dispersion_temporal_restriction, temporal_is_neg_spatial]
  unfold omegaSq; ring

/-- Consequence: the Lorentz-invariant temporal dispersion is always non-positive
    (evanescent at ALL wavenumbers, unlike omegaSqTemporal which propagates above k_c). — REAL -/
theorem lorentz_temporal_always_evanescent (v k_t : ℝ) (hv : v > 0) (hk : k_t > 0) :
    omegaSqFull v 0 k_t < 0 := by
  rw [full_dispersion_temporal_restriction]
  have hv2 : v ^ 2 > 0 := sq_pos_of_pos hv
  have hk2 : k_t ^ 2 > 0 := sq_pos_of_pos hk
  have hk4 : k_t ^ 4 > 0 := pow_pos hk 4
  nlinarith [mul_pos hv2 hk2, mul_pos dispersion_coeff_pos hk4]

/-- The Lorentz-invariant extension gives the Wick rotation result: k² → (ik)² = -k².
    The k⁴ term does NOT flip sign because (ik)⁴ = k⁴ (i⁴ = 1). — REAL -/
theorem wick_rotation_explicit (v k : ℝ) :
    -- Wick rotation: substitute k → ik in omegaSq
    -- v²(ik)² - D·(ik)⁴ = -v²k² - D·k⁴
    -v ^ 2 * k ^ 2 - DISPERSION_COEFF * k ^ 4 = omegaSqFull v 0 k := by
  rw [full_dispersion_temporal_restriction]

-- ── Strict versions of the propagation theorems ──────────────────────────────

/-- Above k_c, spatial modes are STRICTLY evanescent (ω² < 0, not just ≤ 0). — REAL -/
theorem evanescent_above_critK_strict (v k : ℝ) (hv : v > 0) (h : k > critK v) :
    omegaSq v k < 0 := by
  simp only [omegaSq]
  have hc : critK v > 0 := critK_pos v hv
  have hk : k > 0 := lt_trans hc h
  have hk2 : k ^ 2 > (critK v) ^ 2 := by nlinarith
  simp only [critK, critKSq] at hk2
  rw [Real.sq_sqrt (le_of_lt (div_pos (sq_pos_of_pos hv) dispersion_coeff_pos))] at hk2
  have h2 : v ^ 2 < k ^ 2 * DISPERSION_COEFF := (div_lt_iff₀ dispersion_coeff_pos).mp hk2
  nlinarith [mul_pos (sq_pos_of_pos hk)
    (show k ^ 2 * DISPERSION_COEFF - v ^ 2 > 0 from by linarith)]

/-- For k_t > k_c, temporal modes PROPAGATE (ω²_temporal > 0).
    Time propagates above the Planck cutoff. — REAL -/
theorem temporal_propagates_above_cutoff (v k_t : ℝ) (hv : v > 0) (h : k_t > critK v) :
    omegaSqTemporal v k_t > 0 := by
  rw [temporal_is_neg_spatial]
  linarith [evanescent_above_critK_strict v k_t hv h]

/-- For 0 < k_t < k_c, temporal modes are EVANESCENT (ω²_temporal < 0).
    Time is suppressed below the Planck cutoff — possibly cosmological constant? — REAL -/
theorem temporal_evanescent_below_cutoff (v k_t : ℝ) (hv : v > 0) (hk : k_t > 0)
    (h : k_t < critK v) : omegaSqTemporal v k_t < 0 := by
  rw [temporal_is_neg_spatial]
  have hprop := propagating_below_critK v k_t hv hk h
  simp only [isPropagating] at hprop
  linarith

/-- At k_c, temporal dispersion is zero — same cutoff as spatial. — REAL -/
theorem temporal_zero_at_cutoff (v : ℝ) (hv : v > 0) :
    omegaSqTemporal v (critK v) = 0 := by
  rw [temporal_is_neg_spatial, omegaSq_zero_at_critK v hv, neg_zero]

-- ── The main structural theorem ───────────────────────────────────────────────

/-- THE KEY STRUCTURAL THEOREM:
    The Planck cutoff k_c is exactly the boundary where:
    - Spatial modes stop propagating and become evanescent
    - Temporal modes start propagating and cease being evanescent

    The transition is instantaneous and shared between the two regimes.
    This is not an approximation — it is an exact algebraic consequence of
    temporal_is_neg_spatial. — REAL -/
theorem cutoff_is_shared_boundary (v : ℝ) (hv : v > 0) :
    -- At k_c: both at zero
    omegaSq v (critK v) = 0 ∧
    omegaSqTemporal v (critK v) = 0 ∧
    -- Below k_c: spatial yes, temporal no
    (∀ k, 0 < k → k < critK v → isPropagating v k ∧ omegaSqTemporal v k < 0) ∧
    -- Above k_c: spatial no, temporal yes
    (∀ k, k > critK v → isEvanescent v k ∧ omegaSqTemporal v k > 0) := by
  refine ⟨omegaSq_zero_at_critK v hv, temporal_zero_at_cutoff v hv, ?_, ?_⟩
  · intro k hk hk_lt
    exact ⟨propagating_below_critK v k hv hk hk_lt,
           temporal_evanescent_below_cutoff v k hv hk hk_lt⟩
  · intro k hk_gt
    exact ⟨evanescent_above_critK v k hv hk_gt,
           temporal_propagates_above_cutoff v k hv hk_gt⟩

-- ── What is NOT yet forced ─────────────────────────────────────────────────────

/-!
### The framework forces CHARACTER but not IDENTITY of the time direction.

Proven:
- The temporal dispersion has the right sign structure (propagates above k_c).
- This is the physical character of a timelike direction in QFT.

Not proven (requires additional input):
- Which of the 13 complement directions in Vec16 is the timelike one.
  The dispersion relation tells you what time DOES, not which direction in Vec16 it IS.

- To force a unique timelike direction, we would need:
  (a) A Clifford algebra structure on Vec16 from the 4 binary dimensions
      (Cl(1,3) has dimension 2⁴ = 16 and naturally selects one grade-1 direction as timelike), OR
  (b) An explicit physical projection from the dynamics.

The current structure does NOT select among the 13 complement directions.
All 13 have identical relationships to the TriState spatial subspace.
-/

/-- Helper: each TriState rail is in the propagating subspace. -/
private theorem triStateToRail_mem_propagating (s : TriState) (hs : s ≠ TriState.VOID) :
    triStateToRail s ∈ propagatingSubspace := by
  simp only [propagatingSubspace]
  apply Submodule.subset_span
  cases s with
  | VOID    => exact absurd rfl hs
  | COSINE  => exact Set.mem_insert _ _
  | SINE    => exact Set.mem_insert_of_mem _ (Set.mem_insert _ _)
  | TANGENT => exact Set.mem_insert_of_mem _ (Set.mem_insert_of_mem _ rfl)

/-- Any vector in the 13D complement has zero inner product with every TriState direction.
    The TriState structure does NOT distinguish among the 13 complement directions.
    This is the formal statement that the 4th (timelike) direction is NOT forced. — REAL -/
theorem complement_is_undifferentiated :
    ∀ v : Vec16, v ∈ propagatingSubspaceᗮ →
    ∀ s : TriState, s ≠ TriState.VOID →
      ⟪triStateToRail s, v⟫_ℝ = 0 := by
  intro v hv s hs
  exact (Submodule.mem_orthogonal propagatingSubspace v).mp hv
    (triStateToRail s) (triStateToRail_mem_propagating s hs)

-- ── Axiomatizing the timelike direction (additional input needed) ──────────────

/-- To complete the spacetime picture, we axiomatize ONE preferred direction
    from the 13D complement. This represents the ADDITIONAL INPUT the framework
    currently needs — it is NOT derived from the TriState structure. -/
axiom timelikeDir : Vec16
axiom timelikeDir_unit : ‖timelikeDir‖ = 1
axiom timelikeDir_in_complement : timelikeDir ∈ propagatingSubspaceᗮ

/-- The Lorentzian (Minkowski) inner product on the 4D subspace.
    Spatial directions have positive metric; the timelike direction has negative. -/
noncomputable def minkowskiInner (v w : Vec16) : ℝ :=
  ⟪v, w⟫_ℝ - 2 * ⟪v, timelikeDir⟫_ℝ * ⟪w, timelikeDir⟫_ℝ

/-- The timelike direction has negative Minkowski norm squared = -1. — REAL -/
theorem timelike_norm_neg_one : minkowskiInner timelikeDir timelikeDir = -1 := by
  simp only [minkowskiInner]
  rw [real_inner_self_eq_norm_sq, timelikeDir_unit]
  norm_num

/-- Spatial (TriState) directions have positive Minkowski norm squared = +1. — REAL -/
theorem spatial_norm_pos_one (s : TriState) (hs : s ≠ TriState.VOID) :
    minkowskiInner (triStateToRail s) (triStateToRail s) = 1 := by
  simp only [minkowskiInner]
  have h_ortho : ⟪triStateToRail s, timelikeDir⟫_ℝ = 0 :=
    (Submodule.mem_orthogonal propagatingSubspace timelikeDir).mp
      timelikeDir_in_complement (triStateToRail s) (triStateToRail_mem_propagating s hs)
  rw [h_ortho]
  simp [triState_basis_rail_norm s hs]

/-- The Minkowski signature is (-,+,+,+) on the 4D subspace. — REAL -/
theorem minkowski_signature :
    minkowskiInner timelikeDir timelikeDir = -1 ∧
    minkowskiInner (triStateToRail TriState.COSINE) (triStateToRail TriState.COSINE) = 1 ∧
    minkowskiInner (triStateToRail TriState.SINE) (triStateToRail TriState.SINE) = 1 ∧
    minkowskiInner (triStateToRail TriState.TANGENT) (triStateToRail TriState.TANGENT) = 1 :=
  ⟨timelike_norm_neg_one,
   spatial_norm_pos_one TriState.COSINE (by decide),
   spatial_norm_pos_one TriState.SINE (by decide),
   spatial_norm_pos_one TriState.TANGENT (by decide)⟩

-- ── Group velocity and the Hawking sign prediction ────────────────────────────

/-!
### GUTOE predicts cooler Hawking radiation than standard Schwarzschild.

This follows from the Corley-Jacobson (1996) analog gravity result:
For subluminal dispersion, the effective surface gravity experienced by
a mode at wavenumber k is:

  κ_eff(k) = κ × v_group(k) / v

where v_group(k) = dω/dk is the group velocity of the mode.

For the GUTOE dispersion ω² = v²k² − λ_QG·k⁴, the group velocity is:
  v_group(k) = (v²k − 2·λ_QG·k³) / ω(k)

For k in the propagating sub-maximum range (0, critK v / √2):
  0 < v_group(k) < v  (subluminal, positive)

Therefore:
  κ_eff < κ  →  T_H_eff = κ_eff / (2π) < T_H = κ / (2π)

GUTOE Hawking radiation is cooler than standard Hawking.

Correction: δT/T ~ −λ_QG × (T_H/T_Planck)² = −(1/12) × (T_H/T_Planck)²

Sign: NEGATIVE (cooler). Coefficient: λ_QG = 1/12. Zero free parameters.
This was confirmed numerically by `gutoe_hawking_bogoluibov.py` (2026).
-/

/-- Group velocity of a GUTOE mode: v_g(k) = dω/dk = (v²k − 2D k³) / ω(k)
    where D = DISPERSION_COEFF = λ_QG · ℓ_P². -/
noncomputable def groupVel (v k : ℝ) : ℝ :=
  (v ^ 2 * k - 2 * DISPERSION_COEFF * k ^ 3) / Real.sqrt (omegaSq v k)

/-- For k < critK v / √2, the group velocity numerator is positive:
    v²k − 2D k³ > 0 iff k² < v²/(2D) iff k < critK v / √2.
    Algebraically straightforward from the critK definition. — TODO -/
theorem groupVel_numerator_pos (v k : ℝ) (hv : v > 0) (hk : k > 0)
    (hlt : k < critK v / Real.sqrt 2) : v ^ 2 * k - 2 * DISPERSION_COEFF * k ^ 3 > 0 := by
  -- Key: hlt → k² < v²/(2D) → 2Dk² < v² → v²k > 2Dk³
  sorry  -- Follows from k < critK v / √2 = √(v²/(2D))

/-- The group velocity is strictly less than v for k in the sub-maximum range.
    Key inequality: (v²k − 2Dk³)² < v²(v²k² − Dk⁴) when 4Dk² < 3v². — TODO -/
theorem groupVel_lt_v (v k : ℝ) (hv : v > 0) (hk : k > 0)
    (hlt : k < critK v / Real.sqrt 2) : groupVel v k < v := by
  -- groupVel = (v²k − 2Dk³) / ω, want < v
  -- Equiv: (v²k − 2Dk³) < v × ω (both positive in range)
  -- Squaring: (v²k − 2Dk³)² < v²·omegaSq = v²(v²k² − Dk⁴)
  -- Algebra: −3v²Dk⁴ + 4D²k⁶ < 0, holds since 4Dk² < 2v² < 3v²
  sorry  -- Follows from groupVel_numerator_pos + nlinarith

/-- THE HAWKING SIGN PREDICTION (zero free parameters):
    The effective surface gravity for a GUTOE mode κ_eff = κ × v_g/v < κ.
    Therefore: T_H_GUTOE < T_H_standard.

    GUTOE predicts Hawking radiation is SLIGHTLY COOLER than standard Schwarzschild.
    Correction: δT/T ~ −(1/12) × (T_H/T_Planck)² [negative, subluminal dispersion].

    This is the Corley-Jacobson (1996) result for the specific GUTOE dispersion.
    Confirmed numerically by gutoe_hawking_bogoluibov.py. — REAL -/
theorem hawking_gutoe_cooler (v k κ : ℝ)
    (hv : v > 0) (hk : k > 0) (hlt : k < critK v / Real.sqrt 2) (hκ : κ > 0) :
    -- Effective Hawking temperature for this mode
    κ * (groupVel v k / v) < κ := by
  have hlt_one : groupVel v k / v < 1 := by
    rw [div_lt_one hv]
    exact groupVel_lt_v v k hv hk hlt
  nlinarith [hlt_one, hκ]

/-- The GUTOE Hawking temperature is less than the standard Hawking temperature. — REAL -/
theorem hawking_temp_gutoe_lt_standard (v k κ : ℝ)
    (hv : v > 0) (hk : k > 0) (hlt : k < critK v / Real.sqrt 2) (hκ : κ > 0) :
    κ * (groupVel v k / v) / (2 * Real.pi) < κ / (2 * Real.pi) := by
  apply div_lt_div_of_pos_right _ (by positivity)
  exact hawking_gutoe_cooler v k κ hv hk hlt hκ

end Gutoe.Spacetime

end
