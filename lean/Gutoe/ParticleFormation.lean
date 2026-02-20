/-
 * GUTOE - Particle Formation from Void
 * Copyright (C) 2026  Riff Labs
 *
 * AGPL-3.0-or-later
 -/

import Mathlib
import Gutoe.Basic
import Gutoe.RealGates

namespace Gutoe

/-!
# Particle Formation from Void

Formalises the GUTOE claim that matter (specifically quarks) can
emerge from a pure void state through tripartite field dynamics.

The formation sequence (from `void_to_hydrogen_proof.py`):
1. All lattice cells start as VOID (instability = 1, no waves)
2. Centre cell differentiates: VOID → SINE
3. Veracity rails form between the differentiated cell and its six void neighbours
4. Spacetime curvature emerges from the veracity field
5. When `binding_coherence = v / (1 + g) ≥ 0.6`, a quark forms
6. Quark type: UP if `veracity > curvature`, DOWN otherwise
7. 2 UP + 1 DOWN quarks → proton (charge +1); add an electron → hydrogen

All theorems marked `-- REAL` are fully proven.
-/

-- ── Void differentiation ───────────────────────────────────────────────────

/-- Void differentiation: the unstable void collapses to SINE (the first
    asymmetry — nothingness spontaneously becomes existence). -/
def voidDifferentiation : TriState → TriState
  | TriState.VOID => TriState.SINE
  | s             => s

/-- Void differentiates to SINE — REAL -/
theorem void_differentiates_to_sine :
    voidDifferentiation TriState.VOID = TriState.SINE := rfl

/-- Differentiation escapes the void: the output is not VOID — REAL -/
theorem void_differentiation_escapes_void :
    voidDifferentiation TriState.VOID ≠ TriState.VOID := by decide

/-- The differentiated state is a basis wave component (SINE) — REAL -/
theorem void_differentiation_is_basis :
    (voidDifferentiation TriState.VOID).isBasis = true := by decide

/-- Non-VOID states are unchanged by differentiation — REAL -/
theorem differentiation_fixes_non_void (s : TriState) (h : s ≠ TriState.VOID) :
    voidDifferentiation s = s := by
  cases s <;> simp_all [voidDifferentiation]

-- ── Field configuration ────────────────────────────────────────────────────

/-- Local field configuration at a hexagonal lattice site.
    All field values are non-negative (physically required). -/
structure FieldConfig where
  veracity   : ℝ   -- veracity field strength at this site
  curvature  : ℝ   -- spacetime curvature at this site
  field_grad : ℝ   -- average gradient from hex neighbours
  hv : veracity   ≥ 0
  hc : curvature  ≥ 0
  hg : field_grad ≥ 0

/-- Binding coherence: `v / (1 + g)`.
    Measures stability of the local field configuration.
    Threshold for quark formation: ≥ 3/5 (0.6 in the Python simulation). -/
noncomputable def bindingCoherence (fc : FieldConfig) : ℝ :=
  fc.veracity / (1 + fc.field_grad)

/-- Binding coherence is non-negative — REAL -/
theorem binding_coherence_nonneg (fc : FieldConfig) : 0 ≤ bindingCoherence fc :=
  div_nonneg fc.hv (by linarith [fc.hg])

/-- Binding coherence ≤ veracity (dividing by ≥ 1 only decreases) — REAL -/
theorem binding_coherence_le_veracity (fc : FieldConfig) :
    bindingCoherence fc ≤ fc.veracity := by
  unfold bindingCoherence
  apply div_le_self fc.hv
  linarith [fc.hg]

-- ── Quark type and formation ───────────────────────────────────────────────

/-- First-generation quark types relevant to hydrogen -/
inductive QuarkType
  | UP    -- veracity > curvature; charge +2/3
  | DOWN  -- curvature ≥ veracity; charge −1/3
  deriving DecidableEq, Repr

/-- A lattice site can form a quark when its binding coherence ≥ 3/5 -/
def isQuarkForming (fc : FieldConfig) : Prop :=
  bindingCoherence fc ≥ 3 / 5

/-- Classify quark type from the local field balance -/
noncomputable def classifyQuark (fc : FieldConfig) : QuarkType :=
  if fc.veracity > fc.curvature then QuarkType.UP else QuarkType.DOWN

/-- A site with veracity strictly exceeding curvature is UP — REAL -/
theorem classify_up (fc : FieldConfig) (h : fc.veracity > fc.curvature) :
    classifyQuark fc = QuarkType.UP := by
  unfold classifyQuark; simp only [if_pos h]

/-- A site where curvature is not exceeded by veracity is DOWN — REAL -/
theorem classify_down (fc : FieldConfig) (h : ¬fc.veracity > fc.curvature) :
    classifyQuark fc = QuarkType.DOWN := by
  unfold classifyQuark; simp only [if_neg h]

/-- UP and DOWN are the only possibilities — REAL -/
theorem quark_type_exhaustive (fc : FieldConfig) :
    classifyQuark fc = QuarkType.UP ∨ classifyQuark fc = QuarkType.DOWN := by
  unfold classifyQuark
  split_ifs with h
  · exact Or.inl rfl
  · exact Or.inr rfl

-- ── Electric charge ────────────────────────────────────────────────────────

/-- Electric charge of each quark type (rational, in units of e) -/
def quarkCharge : QuarkType → ℚ
  | QuarkType.UP   =>  2 / 3
  | QuarkType.DOWN => -1 / 3

/-- UP quark charge is +2/3 — REAL -/
theorem up_quark_charge : quarkCharge QuarkType.UP = 2 / 3 := rfl

/-- DOWN quark charge is −1/3 — REAL -/
theorem down_quark_charge : quarkCharge QuarkType.DOWN = -1 / 3 := rfl

-- ── Proton (hydrogen nucleus) ─────────────────────────────────────────────

/-!
### Proton charge from quark content

A proton is (uud): two UP quarks and one DOWN quark.
Total charge = 2 × (+2/3) + 1 × (−1/3) = 4/3 − 1/3 = 1.
-/

/-- Proton charge: 2 UP + 1 DOWN = +1 — REAL -/
theorem proton_charge :
    2 * quarkCharge QuarkType.UP + quarkCharge QuarkType.DOWN = 1 := by
  simp only [quarkCharge]; norm_num

/-- The proton charge is an integer — REAL -/
theorem proton_charge_is_integer :
    ∃ n : ℤ, (n : ℚ) = 2 * quarkCharge QuarkType.UP + quarkCharge QuarkType.DOWN :=
  ⟨1, by simp only [quarkCharge]; norm_num⟩

-- ── Mass from field configuration ─────────────────────────────────────────

/-!
### Quark mass formula

`quarkMass v c g l lam = v * c * g / l * lam²`
(defined in `RealGates.lean`, derived from the Python `particle_formation.py`).

When all parameters are strictly positive the mass is strictly positive.
-/

/-- A quark-forming configuration yields a strictly positive mass — REAL -/
theorem quark_forming_has_positive_mass
    (fc : FieldConfig) (l lam : ℝ)
    (hl  : l   > 0) (hlam : lam > 0)
    (hv  : fc.veracity   > 0)
    (hc  : fc.curvature  > 0)
    (hg  : fc.field_grad > 0) :
    quarkMass fc.veracity fc.curvature fc.field_grad l lam > 0 := by
  unfold quarkMass
  exact mul_pos (div_pos (mul_pos (mul_pos hv hc) hg) hl) (pow_pos hlam 2)

-- ── Concrete example: VOID → UP quark ────────────────────────────────────

/-!
### Existence proof: VOID → UP quark

The Python `void_to_hydrogen_proof.py` runs the centre cell with field values
matching those produced by a freshly differentiated SINE cell:

  veracity = 1.0  (full coherence with itself)
  curvature = 0.8 (spacetime curvature from neighbouring veracity rails)
  field_grad = 0.5 (gradient from 6 hex neighbours still in VOID)

Binding coherence = 1 / (1 + 0.5) = 2/3 > 0.6 → quark forms
Quark type        = UP  (veracity 1.0 > curvature 0.8)
-/

/-- Example field configuration from the Python void-to-hydrogen simulation -/
private noncomputable def exampleFC : FieldConfig :=
  { veracity   := 1
    curvature  := 4 / 5
    field_grad := 1 / 2
    hv := by norm_num
    hc := by norm_num
    hg := by norm_num }

/-- Example configuration clears the quark-formation threshold (2/3 ≥ 3/5) — REAL -/
theorem example_forms_quark : isQuarkForming exampleFC := by
  unfold isQuarkForming bindingCoherence exampleFC
  norm_num

/-- Example configuration is classified UP (veracity 1 > curvature 4/5) — REAL -/
theorem example_is_up_quark : classifyQuark exampleFC = QuarkType.UP :=
  classify_up exampleFC (by unfold exampleFC; norm_num)

/-!
### The complete VOID → UP quark path

VOID differentiates to SINE.  A SINE-sourced field configuration with
veracity > curvature meets the quark threshold and is classified UP.
Two UP quarks and one DOWN quark yield a proton with integer charge +1.
-/

/-- VOID → SINE → UP quark: the three formation steps proven together — REAL -/
theorem void_to_up_quark_path :
    -- Step 1: VOID differentiates to a non-void basis state
    voidDifferentiation TriState.VOID = TriState.SINE ∧
    -- Step 2: a SINE-sourced field configuration meets the quark threshold
    isQuarkForming exampleFC ∧
    -- Step 3: that configuration classifies as UP (veracity > curvature)
    classifyQuark exampleFC = QuarkType.UP :=
  ⟨void_differentiates_to_sine, example_forms_quark, example_is_up_quark⟩

end Gutoe
