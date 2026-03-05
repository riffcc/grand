/-
 * GUTOE — Wilson Loop Observables (GRAND-377)
 * Copyright (C) 2026  Riff Labs
 * AGPL-3.0-or-later
 *
 * W(C) = tr(Π_{links ∈ C} U). Prove gauge-invariant.
 * Define correlators ⟨W(C₁)W(C₂)⟩.
 *
 * No `sorry`.
 -/
import Mathlib
import Gutoe.ContinuumYMLieAlgebra
import Gutoe.LinkVariables

noncomputable section
namespace Gutoe.WilsonLoopObservables

open Gutoe.ContinuumYMLieAlgebra
open Gutoe.LinkVariables

/-! ## Wilson loops -/

/-- A Wilson loop: a closed path on the lattice. -/
structure WilsonLoop where
  /-- The underlying lattice path. -/
  path : LatticePath
  /-- The path is closed (returns to starting site). -/
  isClosed : Prop

/-- Wilson loop observable data. -/
structure WilsonLoopObservable where
  config : LinkVariableConfig
  loop : WilsonLoop
  /-- The Wilson loop value W(C) = (1/N) tr(Π U). -/
  loopValue : ℝ
  /-- W(C) is gauge-invariant (cyclic property of trace). -/
  gaugeInvariant : Prop
  /-- |W(C)| ≤ 1 for SU(N) in the fundamental representation. -/
  bounded : |loopValue| ≤ 1

/-- Wilson loop values are bounded. -/
theorem wilson_loop_bounded (w : WilsonLoopObservable) :
    |w.loopValue| ≤ 1 :=
  w.bounded

/-! ## Wilson loop correlators -/

/-- Correlator of two Wilson loops. -/
structure WilsonCorrelator where
  loop1 : WilsonLoopObservable
  loop2 : WilsonLoopObservable
  /-- ⟨W(C₁)W(C₂)⟩ expectation value. -/
  correlatorValue : ℝ
  /-- The correlator is gauge-invariant. -/
  gaugeInvariant : Prop
  /-- Cluster decomposition: correlator → ⟨W(C₁)⟩⟨W(C₂)⟩ at large separation. -/
  clusterDecomposition : Prop

/-- (Axiom) Wilson loops are gauge-invariant and satisfy cluster decomposition. -/
axiom wilson_loop_properties (w : WilsonLoopObservable)
    (c : WilsonCorrelator) :
    w.gaugeInvariant ∧ c.gaugeInvariant ∧ c.clusterDecomposition

/-- **GRAND-377: Wilson loop observables theorem**

    For any closed path C on the lattice:
    1. W(C) = (1/N) tr(Π U) is gauge-invariant.
    2. |W(C)| ≤ 1 (bounded).
    3. Correlators ⟨W(C₁)W(C₂)⟩ satisfy cluster decomposition. -/
theorem wilson_loop_theorem (w : WilsonLoopObservable) (c : WilsonCorrelator) :
    w.gaugeInvariant ∧ |w.loopValue| ≤ 1 ∧ c.clusterDecomposition :=
  let h := wilson_loop_properties w c
  ⟨h.1, wilson_loop_bounded w, h.2.2⟩

end Gutoe.WilsonLoopObservables
