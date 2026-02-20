/-
 * GUTOE - Gate Properties (Formal Verification)
 * Copyright (C) 2026  Riff Labs
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 -/

import Mathlib
import Gutoe.Basic

namespace Gutoe

/-!
## Gate Properties — Formal Verification

Stress-tests the gate algebra claimed by GUTOE.

Key findings, all PROVEN:
1. `TripartiteHadamard` is **not injective** — it maps both COSINE and SINE to
   UNDEFINED, destroying information.  It is therefore **not a valid quantum gate**
   (unitary operators must be bijective).

2. `cycle` has **order 3** (Z₃), not order 2 (Z₂).  Applying a cycle-based gate
   twice does NOT restore the input — three applications are required.

3. `CTripartiteNot` (controlled cycle) applied **twice does not return the target**
   to its original state.  It is not self-inverse.

Theorems marked `-- REAL` are provable.  No broken theorems are included here.
-/

-- ── TripartiteHadamard (modelled as a function) ────────────────────────────

/-- Simulation of the TripartiteHadamard gate applied to a TriState -/
def hadamard : TriState → TriState
  | TriState.COSINE  => TriState.TANGENT
  | TriState.SINE    => TriState.TANGENT
  | TriState.TANGENT => TriState.TANGENT
  | TriState.VOID    => TriState.VOID

/-- Both COSINE and SINE map to TANGENT — REAL -/
theorem hadamard_cosine : hadamard TriState.COSINE = TriState.TANGENT := rfl

/-- Both COSINE and SINE map to TANGENT — REAL -/
theorem hadamard_sine : hadamard TriState.SINE = TriState.TANGENT := rfl

/-- TripartiteHadamard is NOT injective: two distinct inputs share an output — REAL -/
theorem hadamard_not_injective :
    hadamard TriState.COSINE = hadamard TriState.SINE := rfl

/-- TripartiteHadamard is not a bijection — REAL -/
theorem hadamard_not_bijective : ¬ Function.Injective hadamard := by
  intro h
  exact absurd (h hadamard_not_injective) (by decide)

/-- Applying TripartiteHadamard twice does NOT return COSINE — REAL -/
theorem hadamard_not_involutive_cosine :
    hadamard (hadamard TriState.COSINE) ≠ TriState.COSINE := by decide

/-- Applying TripartiteHadamard twice does NOT return SINE — REAL -/
theorem hadamard_not_involutive_sine :
    hadamard (hadamard TriState.SINE) ≠ TriState.SINE := by decide

-- ── cycle is Z₃ on the non-VOID states ────────────────────────────────────

/-- Cycle does NOT have order 2 on COSINE: cycle² COSINE ≠ COSINE — REAL -/
theorem cycle_not_order_2_cosine :
    TriState.cycle (TriState.cycle TriState.COSINE) ≠ TriState.COSINE := by decide

/-- Cycle does NOT have order 2 on SINE: cycle² SINE ≠ SINE — REAL -/
theorem cycle_not_order_2_sine :
    TriState.cycle (TriState.cycle TriState.SINE) ≠ TriState.SINE := by decide

/-- Cycle does NOT have order 2 on TANGENT — REAL -/
theorem cycle_not_order_2_tangent :
    TriState.cycle (TriState.cycle TriState.TANGENT) ≠ TriState.TANGENT := by decide

-- ── Controlled cycle (analogue of CTripartiteNot) ─────────────────────────

/-- Apply cycle to target iff control = SINE, otherwise identity -/
def cCycle (control target : TriState) : TriState :=
  if control = TriState.SINE then target.cycle else target

/-- Applying cCycle with SINE control twice does NOT restore COSINE — REAL -/
theorem c_cycle_not_self_inverse_cosine :
    cCycle TriState.SINE (cCycle TriState.SINE TriState.COSINE) ≠ TriState.COSINE := by
  decide

/-- Applying cCycle with SINE control twice does NOT restore SINE — REAL -/
theorem c_cycle_not_self_inverse_sine :
    cCycle TriState.SINE (cCycle TriState.SINE TriState.SINE) ≠ TriState.SINE := by
  decide

/-- Applying cCycle with SINE control three times DOES restore the target — REAL -/
theorem c_cycle_order_3 (target : TriState) (h : target ≠ TriState.VOID) :
    cCycle TriState.SINE (cCycle TriState.SINE (cCycle TriState.SINE target)) = target := by
  cases target <;> simp_all [cCycle, TriState.cycle]

/-- Applying cCycle with SINE control twice does NOT restore TANGENT — REAL -/
theorem c_cycle_not_self_inverse_tangent :
    cCycle TriState.SINE (cCycle TriState.SINE TriState.TANGENT) ≠ TriState.TANGENT := by
  decide

/-- cCycle with any non-SINE control is the identity — REAL -/
theorem c_cycle_identity_when_control_not_sine
    (control target : TriState) (h : control ≠ TriState.SINE) :
    cCycle control target = target := by
  simp [cCycle, h]

end Gutoe
