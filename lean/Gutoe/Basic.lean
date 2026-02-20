/-
 * GUTOE Core - Tripartite Quantum System
 * Copyright (C) 2026  Riff Labs
 *
 * AGPL-3.0-or-later
 -/

import Mathlib

/-!
# GUTOE Tripartite Quantum System

Formalization of the 4-state quantum system from the Grand Unifying Theory of Everything.

States (spectral Fourier components):
- VOID:    Zero mode — absolute nothingness, no wave content
- COSINE:  cos wave component — |0⟩, the "flat" (amplitude-extremum) part of the wave
- SINE:    sin wave component — |1⟩, the "steep" (zero-crossing) part of the wave
- TANGENT: tan = sin/cos — the slope/ratio state; diverges when the cosine component
           vanishes (at phase π/2, 3π/2).  This is the "relationship" state that
           connects the two basis modes and manifests as vector rails in the theory.

SINE and COSINE are 60° adjacent in the hexagonal lattice — NOT antipodal.
veracity(SINE, COSINE) = √3/2, not 0.

## What this file checks

The goal is NOT to blindly transcribe the theory — it is to *stress-test* it.
Theorems marked `-- REAL` are provable and hold.
Theorems marked `-- BROKEN` are false and have been corrected or removed.
-/

namespace Gutoe

/-- The four fundamental GUTOE states -/
inductive TriState
| VOID
| COSINE  -- |0⟩ cosine wave component
| SINE    -- |1⟩ sine wave component
| TANGENT -- tan = sin/cos, the relationship/slope state
deriving DecidableEq, Repr

namespace TriState

/-- Convert state to complex amplitude (real part) -/
noncomputable def toAmplitudeReal : TriState → ℝ
| VOID    => 0
| COSINE  => 1
| SINE    => 0
| TANGENT => 1 / Real.sqrt 2

/-- Convert state to complex amplitude (imaginary part) -/
noncomputable def toAmplitudeImag : TriState → ℝ
| VOID    => 0
| COSINE  => 0
| SINE    => 1
| TANGENT => 1 / Real.sqrt 2

/-- Cycle operation: SINE → COSINE → TANGENT → SINE, VOID is fixed -/
def cycle : TriState → TriState
| SINE    => COSINE
| COSINE  => TANGENT
| TANGENT => SINE
| VOID    => VOID

/-- Check if state is a basis wave component (SINE or COSINE) -/
def isBasis : TriState → Bool
| SINE   => true
| COSINE => true
| _      => false

/-- Phase factor for each state -/
noncomputable def phase : TriState → ℝ
| SINE    => 0
| COSINE  => Real.pi
| TANGENT => Real.pi / 2
| VOID    => 0

/-! ## Cycle theorems — REAL, all proven by `rfl` -/

theorem cycle_sine    : cycle SINE    = COSINE  := rfl
theorem cycle_cosine  : cycle COSINE  = TANGENT := rfl
theorem cycle_tangent : cycle TANGENT = SINE    := rfl
theorem cycle_void    : cycle VOID    = VOID    := rfl

/-! ## Basis theorems — REAL -/

theorem isBasis_sine    : isBasis SINE    = true  := rfl
theorem isBasis_cosine  : isBasis COSINE  = true  := rfl
theorem isBasis_tangent : isBasis TANGENT = false := rfl
theorem isBasis_void    : isBasis VOID    = false := rfl

/-! ## Cycle is order-3 on the non-VOID states — REAL -/

-- Applying cycle 3 times on each state returns the original.
theorem cycle_order3_sine    : cycle (cycle (cycle SINE))    = SINE    := rfl
theorem cycle_order3_cosine  : cycle (cycle (cycle COSINE))  = COSINE  := rfl
theorem cycle_order3_tangent : cycle (cycle (cycle TANGENT)) = TANGENT := rfl

-- VOID is cycle's only fixed point.
theorem cycle_fixed_iff_void (s : TriState) : cycle s = s ↔ s = VOID := by
  cases s <;> simp [cycle]

/-! ## The cycle group is Z₃ on {SINE, COSINE, TANGENT} — REAL -/

-- cycle is injective on the non-VOID subset.
theorem cycle_injective : Function.Injective cycle := by
  intro a b h
  cases a <;> cases b <;> simp_all [cycle]

/-! ## The cycle is NOT involutive — REAL -/

-- cycle ∘ cycle ≠ id (it would take 3 applications to return)
theorem cycle_not_involutive : ∃ s : TriState, cycle (cycle s) ≠ s := by
  exact ⟨SINE, by simp [cycle]⟩

end TriState

/-!
## Veracity

Veracity measures the "strength" of a relationship between two states.

Physical basis (hexagonal lattice):
- Same state = full coherence = 1
- SINE ↔ COSINE: 60° adjacent hex modes → √3/2 ≈ 0.866 (NOT orthogonal)
- TANGENT with SINE or COSINE: 1/2 (partial — ratio state)
- VOID with anything: 0 (void has no connections)
-/

/-- Veracity between two states -/
noncomputable def veracity (a b : TriState) : ℝ :=
  match a, b with
  | TriState.VOID,    _                  => 0
  | _,                TriState.VOID      => 0
  | TriState.COSINE,  TriState.COSINE    => 1
  | TriState.SINE,    TriState.SINE      => 1
  | TriState.TANGENT, TriState.TANGENT   => 1
  | TriState.COSINE,  TriState.SINE      => Real.sqrt 3 / 2
  | TriState.SINE,    TriState.COSINE    => Real.sqrt 3 / 2
  | _,                _                  => 1 / 2  -- TANGENT-SINE, TANGENT-COSINE

/-! ## Veracity theorems -/

-- Symmetry — REAL
theorem veracity_symm (a b : TriState) : veracity a b = veracity b a := by
  cases a <;> cases b <;> simp [veracity]

-- BROKEN in original: `veracity_refl` claimed `veracity s s = 1` for ALL s,
-- but veracity VOID VOID = 0.
-- Corrected: reflexivity holds for all non-VOID states.
theorem veracity_refl_basis (s : TriState) (h : s.isBasis = true) :
    veracity s s = 1 := by
  cases s <;> simp [TriState.isBasis] at h <;> simp [veracity]

-- TANGENT has full veracity with itself — REAL
theorem veracity_tangent_self : veracity TriState.TANGENT TriState.TANGENT = 1 := by
  simp [veracity]

-- VOID has zero veracity with itself — REAL
theorem veracity_void_self : veracity TriState.VOID TriState.VOID = 0 := by
  simp [veracity]

-- SINE and COSINE are 60° adjacent hex modes: veracity = √3/2 — REAL
theorem veracity_sine_cosine_hex :
    veracity TriState.SINE TriState.COSINE = Real.sqrt 3 / 2 := by
  simp [veracity]

-- Anti-symmetry for basis: veracity = 1 ↔ same basis state — REAL
-- (if two basis states have veracity 1, they are equal)
theorem veracity_antisymm_basis (a b : TriState)
    (ha : a.isBasis = true) (hb : b.isBasis = true)
    (h  : veracity a b = 1) : a = b := by
  cases a <;> cases b <;> simp_all [TriState.isBasis, veracity]
  all_goals (
    have hne : Real.sqrt 3 / 2 ≠ 1 := by
      have h3 : Real.sqrt 3 < 2 := by
        nlinarith [Real.sq_sqrt (show (3 : ℝ) ≥ 0 by norm_num),
                   sq_nonneg (Real.sqrt 3)]
      linarith
    exact absurd ‹Real.sqrt 3 / 2 = 1› hne)

end Gutoe
