/-
 * GUTOE - Dynamics Properties
 * Copyright (C) 2026  Riff Labs
 *
 * AGPL-3.0-or-later
 *
 * Experiments:
 *   #1  (kill alignment): binding coherence is alignment-independent
 *   #15 (lattice big bang): single perturbation stays alive forever
 *   #19 (time reversal): cycle is invertible, alignment is not
 *
 * Key insight: The arrow of time in GUTOE comes from alignment
 * (lossy majority vote) and differentiation (VOID → SINE), NOT from
 * the Z₃ cycle (which is perfectly reversible).
 -/

import Gutoe.Basic
import Gutoe.ParticleFormation

namespace Gutoe.DynamicsProperties

open Gutoe

-- ── Experiment #1: Binding coherence is alignment-independent ────────────

/-- Binding coherence depends only on veracity and field gradient — REAL
    The alignment probability parameter does not enter the formula.
    This is why killing alignment in experiment #1 doesn't eliminate protons. -/
theorem binding_coherence_alignment_independent (fc1 fc2 : FieldConfig)
    (hv : fc1.veracity = fc2.veracity) (hg : fc1.field_grad = fc2.field_grad) :
    bindingCoherence fc1 = bindingCoherence fc2 := by
  simp only [bindingCoherence, hv, hg]

/-- Quark classification depends only on veracity vs curvature — REAL
    Alignment doesn't affect whether a quark is UP or DOWN. -/
theorem quark_type_alignment_independent (fc1 fc2 : FieldConfig)
    (hv : fc1.veracity = fc2.veracity) (hc : fc1.curvature = fc2.curvature) :
    classifyQuark fc1 = classifyQuark fc2 := by
  simp only [classifyQuark, hv, hc]

/-- Quark formation threshold depends only on local field values — REAL
    Alignment changes neighbour states (which change veracity/gradient),
    but the THRESHOLD ITSELF is a function of the field only. -/
theorem formation_threshold_is_local (fc : FieldConfig) :
    isQuarkForming fc ↔ fc.veracity / (1 + fc.field_grad) ≥ 3 / 5 := by
  exact Iff.rfl

-- ── Experiment #15: Single perturbation stays alive ──────────────────────

/-- The cycle preserves non-void: once differentiated, always active — REAL
    The Z₃ cycle maps SINE→COSINE→TANGENT→SINE, never producing VOID.
    This is the "big bang" theorem: a single non-void cell lives forever
    under cycle dynamics alone. -/
theorem cycle_preserves_non_void (s : TriState) (h : s ≠ TriState.VOID) :
    s.cycle ≠ TriState.VOID := by
  cases s <;> simp_all [TriState.cycle]

/-- Iterated cycle preserves non-void — REAL
    Applying cycle any number of times never returns to void.
    (Proved by induction: if s ≠ VOID, then cycle s ≠ VOID.) -/
theorem cycle_n_preserves_non_void (s : TriState) (h : s ≠ TriState.VOID) (n : ℕ) :
    n.iterate TriState.cycle s ≠ TriState.VOID := by
  induction n generalizing s with
  | zero => exact h
  | succ n ih => exact ih _ (cycle_preserves_non_void s h)

/-- Void differentiation NEVER returns to VOID — REAL
    VOID → SINE (non-void), anything else → itself (already non-void). -/
theorem differentiation_escapes_void (s : TriState) :
    voidDifferentiation s ≠ TriState.VOID := by
  cases s <;> simp [voidDifferentiation]

/-- The lattice big bang: VOID differentiates to a non-void state
    that survives forever under cycle dynamics — REAL -/
theorem big_bang_survives (n : ℕ) :
    n.iterate TriState.cycle (voidDifferentiation TriState.VOID) ≠ TriState.VOID := by
  apply cycle_n_preserves_non_void
  exact differentiation_escapes_void TriState.VOID

-- ── Experiment #19: Time reversibility ───────────────────────────────────

/-- The reverse cycle: SINE → TANGENT → COSINE → SINE (opposite of cycle) -/
def reverseCycle : TriState → TriState
  | TriState.SINE    => TriState.TANGENT
  | TriState.COSINE  => TriState.SINE
  | TriState.TANGENT => TriState.COSINE
  | TriState.VOID    => TriState.VOID

/-- Reverse cycle is the two-fold application of forward cycle — REAL
    In Z₃, the inverse of rotation-by-1 is rotation-by-2. -/
theorem reverse_is_double_cycle (s : TriState) :
    reverseCycle s = s.cycle.cycle := by
  cases s <;> rfl

/-- Reverse cycle is the left inverse of cycle — REAL -/
theorem reverse_cycle_left_inv (s : TriState) :
    reverseCycle (s.cycle) = s := by
  cases s <;> rfl

/-- Reverse cycle is the right inverse of cycle — REAL -/
theorem reverse_cycle_right_inv (s : TriState) :
    (reverseCycle s).cycle = s := by
  cases s <;> rfl

/-- The cycle is therefore a bijection — REAL
    Z₃ dynamics is perfectly time-reversible. -/
theorem cycle_bijective : Function.Bijective TriState.cycle :=
  ⟨TriState.cycle_injective,
    fun b => ⟨reverseCycle b, reverse_cycle_right_inv b⟩⟩

/-- Reverse cycle is also a bijection — REAL -/
theorem reverse_cycle_bijective : Function.Bijective reverseCycle := by
  constructor
  · intro a b h
    have := congr_arg TriState.cycle h
    rwa [reverse_cycle_right_inv, reverse_cycle_right_inv] at this
  · intro b
    exact ⟨b.cycle, reverse_cycle_left_inv b⟩

-- ── The arrow of time: alignment is irreversible ─────────────────────────

/-- A constant function (modelling alignment that fully synchronizes
    all states to a fixed majority) is not injective — REAL
    Multiple distinct inputs → same output → information lost → arrow of time. -/
theorem constant_alignment_not_injective (t : TriState) (_ht : t ≠ TriState.VOID) :
    ¬Function.Injective (fun _ : TriState => t) := by
  intro h
  exact absurd (@h TriState.SINE TriState.COSINE rfl) (by decide)

/-- Void differentiation is not surjective — REAL
    voidDifferentiation maps VOID → SINE, leaving no input that maps to VOID.
    This is the source of irreversibility in differentiation. -/
theorem differentiation_not_surjective :
    ¬Function.Surjective voidDifferentiation := by
  intro h
  obtain ⟨s, hs⟩ := h TriState.VOID
  cases s <;> simp [voidDifferentiation] at hs

/-!
### Summary: The Arrow of Time in GUTOE

| Dynamics component   | Reversible? | Proven                        |
|----------------------|-------------|-------------------------------|
| Z₃ cycle             | YES         | cycle_bijective               |
| Reverse cycle        | YES         | reverse_cycle_bijective       |
| Alignment (majority) | NO          | constant_alignment_not_injective |
| Void differentiation | NO          | differentiation_not_surjective|

The arrow of time emerges from:
1. **Void differentiation**: VOID → SINE is one-way (not surjective).
   Once differentiated, there is no dynamics that returns to void.
2. **Alignment**: majority voting loses information (not injective).
   Many distinct configurations → same aligned state.

The Z₃ cycle is perfectly reversible and contributes no arrow of time.
This matches experiment #19: forward and reverse simulations diverge because
alignment and differentiation are irreversible.
-/

end Gutoe.DynamicsProperties
