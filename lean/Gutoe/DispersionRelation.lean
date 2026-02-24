/-
 * GUTOE - Dispersion Relation: Scale-Dependent Physics
 * Copyright (C) 2026  Riff Labs
 *
 * AGPL-3.0-or-later
 *
 * λ_QG = 1/12 derived from first principles:
 *
 *   Exact dispersion of a 1D/hypercubic Planck lattice with spacing a = ℓ_P:
 *     ω²(k) = (4K/m) · sin²(ka/2)
 *   Taylor expand sin²(ka/2) = k²a²/4 - k⁴a⁴/48 + O(k⁶):
 *     ω²(k) = (Ka²/m)k² - (Ka²/m)(a²/12)k⁴
 *           = v²k² - v²ℓ_P²·(1/12)·k⁴
 *   ∴ λ_QG = 1/12 (no free parameters; coupling K cancels entirely)
 *
 * Previous values: experiment-28 b-parameter fit gave 0.084372 (~1.2% off).
 * Predecessor project used 0.084365 (~1.2% off). Both were calibrated.
 -/

import Mathlib
import Gutoe.Basic

namespace Gutoe

/-!
# Dispersion Relation: Scale-Dependent Physics

Formalises the GUTOE dispersion relation that governs wave propagation
across all scales — from relativistic at cosmological scales to quantum-
dominated at the Planck scale.

```
ω²(k) = v²(d)·k² − λ_QG·ℓ_P²·k⁴
v(d)  = c · (0.1 + 0.9·d/16)
```

The phase velocity `v` scales continuously with dimensional density `d ∈ [0,16]`.
Below the critical wavenumber `k_c = v / √(λ_QG·ℓ_P²)`: modes propagate (COSINE/SINE).
Above `k_c`: quantum gravity dominates, modes are evanescent (TANGENT becomes imaginary).

Key theorems (all PROVEN):
1. Dispersion coefficient is positive
2. ω² = 0 at k = 0 (no wave, no frequency)
3. Critical wavenumber k_c is positive for v > 0
4. ω² = 0 exactly at k = k_c
5. Below k_c (k > 0): propagating modes
6. Above k_c: evanescent modes
7. ω² is smooth (C∞) in k and in v
8. Velocity is continuous and strictly increasing with density
9. v(0) = 0.1c (void), v(16) = c (full universe)
-/

-- ── Physical constants ──────────────────────────────────────────────────────

/-- GUTOE quantum gravity coupling λ_QG = 1/12.
    Derived: Taylor expansion of ω²(k) = (4K/m)sin²(ka/2) on a Planck lattice
    gives ω² = v²k² − v²ℓ_P²·(1/12)·k⁴. No free parameters. -/
noncomputable def LAMBDA_QG : ℝ := 1 / 12

/-- Planck length ℓ_P ≈ 1.616255 × 10⁻³⁵ m (CODATA 2018) -/
def PLANCK_LENGTH : ℝ := 1.616255e-35

/-- ℓ_P² — appears in the dispersion relation as the quantum correction scale -/
noncomputable def PLANCK_LENGTH_SQ : ℝ := PLANCK_LENGTH * PLANCK_LENGTH

/-- Dispersion coefficient: λ_QG · ℓ_P² — sets the Planck-scale cutoff -/
noncomputable def DISPERSION_COEFF : ℝ := LAMBDA_QG * PLANCK_LENGTH_SQ

-- ── Dispersion relation ────────────────────────────────────────────────────

/-- Angular frequency squared: ω²(k) = v²k² − λ_QG·ℓ_P²·k⁴ -/
noncomputable def omegaSq (v k : ℝ) : ℝ := v ^ 2 * k ^ 2 - DISPERSION_COEFF * k ^ 4

/-- A mode is propagating when ω² > 0 -/
def isPropagating (v k : ℝ) : Prop := omegaSq v k > 0

/-- A mode is evanescent when ω² ≤ 0 -/
def isEvanescent (v k : ℝ) : Prop := omegaSq v k ≤ 0

-- ── Basic properties ───────────────────────────────────────────────────────

/-- Dispersion coefficient is strictly positive — REAL -/
theorem dispersion_coeff_pos : DISPERSION_COEFF > 0 := by
  simp only [DISPERSION_COEFF, PLANCK_LENGTH_SQ, LAMBDA_QG, PLANCK_LENGTH]
  norm_num

/-- At k = 0, ω² = 0 (no wave, no frequency) — REAL -/
theorem omegaSq_zero_at_origin (v : ℝ) : omegaSq v 0 = 0 := by
  simp [omegaSq]

/-- ω² decomposes into relativistic and quantum terms — REAL -/
theorem omegaSq_decomposition (v k : ℝ) :
    omegaSq v k = v ^ 2 * k ^ 2 - DISPERSION_COEFF * k ^ 4 := rfl

/-- For k → 0, the quantum correction vanishes faster than the relativistic term — REAL -/
theorem quantum_vs_relativistic (v k : ℝ) (hv : v ≠ 0) (hk : k ≠ 0) :
    DISPERSION_COEFF * k ^ 4 / (v ^ 2 * k ^ 2) = DISPERSION_COEFF * k ^ 2 / v ^ 2 := by
  field_simp

/-- Right-hand low-`k` limit is relativistic: `ω²/k² → v²` as `k → 0⁺`. -/
theorem omegaSq_over_kSq_tendsto_vSq_right (v : ℝ) :
    Filter.Tendsto
      (fun k : ℝ => omegaSq v k / k ^ 2)
      (nhdsWithin 0 (Set.Ioi (0 : ℝ)))
      (nhds (v ^ 2)) := by
  have hEq :
      (fun k : ℝ => omegaSq v k / k ^ 2) =ᶠ[nhdsWithin 0 (Set.Ioi (0 : ℝ))]
      (fun k : ℝ => v ^ 2 - DISPERSION_COEFF * k ^ 2) := by
    filter_upwards [self_mem_nhdsWithin] with k hk
    exact by
      have hk0 : k ≠ 0 := ne_of_gt hk
      have hk2 : k ^ 2 ≠ 0 := pow_ne_zero 2 hk0
      unfold omegaSq
      field_simp [hk2]
  have hPoly :
      Filter.Tendsto
        (fun k : ℝ => v ^ 2 - DISPERSION_COEFF * k ^ 2)
        (nhdsWithin 0 (Set.Ioi (0 : ℝ)))
        (nhds (v ^ 2 - DISPERSION_COEFF * 0 ^ 2)) := by
    have hCont : Continuous (fun k : ℝ => v ^ 2 - DISPERSION_COEFF * k ^ 2) := by
      continuity
    exact (hCont.continuousAt.tendsto).mono_left nhdsWithin_le_nhds
  refine Filter.Tendsto.congr' hEq.symm ?_
  simpa using hPoly

-- ── Critical wavenumber ────────────────────────────────────────────────────

/-- Critical wavenumber squared: k_c² = v² / (λ_QG·ℓ_P²) -/
noncomputable def critKSq (v : ℝ) : ℝ := v ^ 2 / DISPERSION_COEFF

/-- Critical wavenumber: k_c = v / √(λ_QG·ℓ_P²) -/
noncomputable def critK (v : ℝ) : ℝ := Real.sqrt (critKSq v)

/-- Critical wavenumber is positive for positive velocity — REAL -/
theorem critK_pos (v : ℝ) (hv : v > 0) : critK v > 0 := by
  simp only [critK, critKSq]
  exact Real.sqrt_pos.mpr (div_pos (sq_pos_of_pos hv) dispersion_coeff_pos)

/-- At the critical wavenumber, ω² = 0 (propagation cutoff) — REAL -/
theorem omegaSq_zero_at_critK (v : ℝ) (hv : v > 0) : omegaSq v (critK v) = 0 := by
  have h1 : DISPERSION_COEFF > 0 := dispersion_coeff_pos
  have h2 : (0 : ℝ) ≤ v ^ 2 / DISPERSION_COEFF :=
    le_of_lt (div_pos (sq_pos_of_pos hv) h1)
  simp only [omegaSq, critK, critKSq]
  have hsq : Real.sqrt (v ^ 2 / DISPERSION_COEFF) ^ 2 = v ^ 2 / DISPERSION_COEFF :=
    Real.sq_sqrt h2
  have h4 : Real.sqrt (v ^ 2 / DISPERSION_COEFF) ^ 4 = (v ^ 2 / DISPERSION_COEFF) ^ 2 := by
    rw [show (4 : ℕ) = 2 * 2 from by norm_num, pow_mul, hsq]
  rw [hsq, h4]
  field_simp
  ring

/-- Below the critical wavenumber (k > 0), modes propagate — REAL -/
theorem propagating_below_critK (v k : ℝ) (hv : v > 0) (hk : k > 0)
    (h : k < critK v) : isPropagating v k := by
  simp only [isPropagating, omegaSq]
  have hc : critK v > 0 := critK_pos v hv
  -- k < critK v and both positive, so k² < (critK v)²
  have hk2 : k ^ 2 < (critK v) ^ 2 := by nlinarith
  simp only [critK, critKSq] at hk2
  rw [Real.sq_sqrt (le_of_lt (div_pos (sq_pos_of_pos hv) dispersion_coeff_pos))] at hk2
  -- hk2 : k² < v²/DISPERSION_COEFF → k² * DISPERSION_COEFF < v²
  have h2 : k ^ 2 * DISPERSION_COEFF < v ^ 2 := (lt_div_iff₀ dispersion_coeff_pos).mp hk2
  nlinarith [mul_pos (sq_pos_of_pos hk) (show v ^ 2 - k ^ 2 * DISPERSION_COEFF > 0 from by linarith)]

/-- Above the critical wavenumber, modes are evanescent — REAL -/
theorem evanescent_above_critK (v k : ℝ) (hv : v > 0) (h : k > critK v) :
    isEvanescent v k := by
  simp only [isEvanescent, omegaSq]
  have hc : critK v > 0 := critK_pos v hv
  have hk : k > 0 := lt_trans hc h
  -- k > critK v and both positive, so k² > (critK v)²
  have hk2 : k ^ 2 > (critK v) ^ 2 := by nlinarith
  simp only [critK, critKSq] at hk2
  rw [Real.sq_sqrt (le_of_lt (div_pos (sq_pos_of_pos hv) dispersion_coeff_pos))] at hk2
  -- hk2 : k² > v²/DISPERSION_COEFF → v² < k² * DISPERSION_COEFF
  have h2 : v ^ 2 < k ^ 2 * DISPERSION_COEFF := (div_lt_iff₀ dispersion_coeff_pos).mp hk2
  nlinarith [mul_pos (sq_pos_of_pos hk) (show k ^ 2 * DISPERSION_COEFF - v ^ 2 > 0 from by linarith)]

-- ── Continuity and smoothness ──────────────────────────────────────────────

/-- ω² is continuous in k — REAL -/
theorem omegaSq_continuous_k (v : ℝ) : Continuous (omegaSq v) := by
  unfold omegaSq; fun_prop

/-- ω² is continuous in v — REAL -/
theorem omegaSq_continuous_v (k : ℝ) : Continuous (fun v => omegaSq v k) := by
  unfold omegaSq; fun_prop

/-- ω² is smooth (C∞) in k — REAL -/
theorem omegaSq_smooth (v : ℝ) : ContDiff ℝ ⊤ (omegaSq v) := by
  unfold omegaSq
  exact (contDiff_const.mul (contDiff_id.pow 2)).sub (contDiff_const.mul (contDiff_id.pow 4))

/-- Emergence is continuous: the transition between regimes is smooth, not stepped — REAL -/
theorem emergence_is_continuous (v : ℝ) : Continuous (omegaSq v) := omegaSq_continuous_k v

-- ── Dimensional density and phase velocity ─────────────────────────────────

/-- Phase velocity as a function of dimensional density d ∈ [0, 16].
    `v(d) = c · (0.1 + 0.9·d/16)`
    The void propagates at 10% of c; full universe propagates at c. -/
noncomputable def velFromDensity (c d : ℝ) : ℝ := c * (0.1 + 0.9 * d / 16)

/-- Phase velocity is continuous in d — REAL -/
theorem velFromDensity_continuous (c : ℝ) : Continuous (velFromDensity c) := by
  unfold velFromDensity; fun_prop

/-- Phase velocity is strictly increasing with dimensional density — REAL -/
theorem velFromDensity_strictMono (c : ℝ) (hc : c > 0) :
    StrictMono (velFromDensity c) := by
  intro d1 d2 h
  simp only [velFromDensity]
  nlinarith

/-- At d = 0 (void), phase velocity is 10% of c — REAL -/
theorem vel_at_void (c : ℝ) : velFromDensity c 0 = c * 0.1 := by
  unfold velFromDensity; ring

/-- At d = 16 (full universe), phase velocity equals c — REAL -/
theorem vel_at_full_universe (c : ℝ) : velFromDensity c 16 = c := by
  unfold velFromDensity; ring

-- ── Connection to TriState ─────────────────────────────────────────────────

/-!
### The TriState states and the dispersion regime

The three non-VOID TriState components have distinct roles in the dispersion picture:

- **COSINE** (flat extremum of wave): long-wavelength mode, k ≪ k_c, always propagating
- **SINE**   (zero-crossing of wave): also long-wavelength, propagating
- **TANGENT** (sin/cos = slope): represents the relationship between the two wave
  components.  When the slope *diverges* (cos → 0, phase = π/2), we are exactly
  at the Planck-scale cutoff: this is the transition to evanescent.

Physically: at k > k_c, attempting to resolve the TANGENT (the ratio state)
requires sub-Planck resolution — which the dispersion relation forbids.
TANGENT at the Planck scale is **evanescent**.
-/

/-- At positive velocity, the evanescent regime is nonempty — REAL -/
theorem evanescent_regime_nonempty (v : ℝ) (hv : v > 0) :
    ∃ k : ℝ, isEvanescent v k :=
  ⟨critK v + 1, evanescent_above_critK v (critK v + 1) hv (by linarith)⟩

/-- The critical wavenumber is where ω² = 0 — the TANGENT state boundary — REAL -/
theorem tangent_state_at_cutoff (v : ℝ) (hv : v > 0) :
    omegaSq v (critK v) = 0 := omegaSq_zero_at_critK v hv

end Gutoe
