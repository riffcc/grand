/-
 * GUTOE — GRAND-331 End-to-End OS Reconstruction Bridge
 * Copyright (C) 2026  Riff Labs
 * AGPL-3.0-or-later
 *
 * GRAND-331:
 *   Build explicit OS/Hilbert/Hamiltonian objects directly from the continuum
 *   Schwinger lane (GRAND-330) and completion/generator lane (GRAND-321).
 *
 * No `sorry`.
 -/

import Mathlib
import Gutoe.YangMillsContinuumLimit
import Gutoe.YangMillsOSTextbook
import Gutoe.YangMillsOSCompletion

noncomputable section

namespace Gutoe.YangMillsOSEndToEnd

open Gutoe.YangMillsConstructiveQFT
open Gutoe.YangMillsContinuumLimit
open Gutoe.YangMillsOSTextbook
open Gutoe.YangMillsOSCompletion
open Gutoe.YangMillsWilsonBridge
open Gutoe.YangMillsWilsonEquivalence

/-- Explicit GRAND-331 step package:
for each refinement step, expose the concrete kernel used to build the
Schwinger family, the quotient/completion Hilbert objects, and the associated
self-adjoint Hamiltonian/Wightman data. -/
structure OSEndToEndStepPackage
    (W : WilsonZ3Action)
    (a_t : ℕ → ℝ)
    (alpha : ℝ)
    (n : ℕ) where
  K : Matrix (Fin 3) (Fin 3) ℝ
  K_eq : K = wilsonKernelAt W alpha n
  schwinger_eq :
    schwingerFamilyFromKernel (fun i j => K i j) = wilsonSchwingerFamily W alpha n
  quotientNonempty : Nonempty (OSHilbertQuot K)
  completionDense : DenseRange ((↑) : OSHilbertQuot K → OSCauchyCompletion K)
  generatorSelfAdjoint : IsSelfAdjoint (osGeneratorAt W a_t alpha n)
  hamiltonianPos : 0 < osHamiltonianAt W a_t alpha n
  wightman_eq :
    ∀ t : ℕ,
      wightmanAt W a_t alpha n t =
        Real.exp (-(osHamiltonianAt W a_t alpha n) * (t : ℝ))

/-- GRAND-331 end-to-end closure:
1. An explicit continuum Schwinger family is fixed.
2. It is normalized at every refinement and n-point level.
3. Every refinement step admits an explicit OS/Hilbert/Hamiltonian package.
4. A uniform positive Hamiltonian floor survives across the full schedule. -/
theorem grand331_end_to_end_os_reconstruction_of_domain
    (W : WilsonZ3Action)
    (a_t : ℕ → ℝ)
    (alpha : ℝ)
    (hDom : WilsonEquivalenceDomain a_t alpha) :
    ∃ SF : ℕ → CorrelatorFamily,
      (∀ n, SF n = wilsonSchwingerFamily W alpha n) ∧
      (∀ n m, SF n m (fun _ => 1) = 1) ∧
      (∀ n, Nonempty (OSEndToEndStepPackage W a_t alpha n)) ∧
      (∃ c : ℝ, 0 < c ∧ ∀ n, c ≤ osHamiltonianAt W a_t alpha n) := by
  let SF : ℕ → CorrelatorFamily := wilsonSchwingerFamily W alpha
  have hSch : (∀ n m, SF n m (fun _ => 1) = 1) := by
    intro n m
    -- Pull normalization directly from the continuum Schwinger-family theorem.
    simpa [SF] using (constructive_schwinger_family_exists W a_t alpha hDom).1 n m
  refine ⟨SF, ?_, hSch, ?_, ?_⟩
  · intro n
    rfl
  · intro n
    refine ⟨{ K := wilsonKernelAt W alpha n
              K_eq := rfl
              schwinger_eq := ?_
              quotientNonempty := ?_
              completionDense := ?_
              generatorSelfAdjoint := ?_
              hamiltonianPos := ?_
              wightman_eq := ?_ }⟩
    · simpa [wilsonKernelAt, wilsonSchwingerFamily]
    · exact ⟨Quotient.mk _ (fun _ => 0)⟩
    · simpa using osQuot_dense_in_completion (wilsonKernelAt W alpha n)
    · exact osGeneratorAt_selfAdjoint W a_t alpha n
    · exact osGeneratorAt_gap_positive_of_domain W a_t alpha hDom n
    · intro t
      rfl
  · exact osGenerator_uniform_gap_floor_of_domain W a_t alpha hDom

end Gutoe.YangMillsOSEndToEnd

