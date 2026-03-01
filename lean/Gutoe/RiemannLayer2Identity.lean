import Mathlib
import Gutoe.RiemannCore
import Gutoe.RiemannBridge
import Gutoe.RiemannCounting

namespace Gutoe.RiemannLayer2Identity

open Gutoe.RiemannCore
open Gutoe.RiemannBridge
open Gutoe.RiemannCounting

noncomputable section

/-- Union-spectrum extracted from a finite truncation ladder. -/
def ladderSpec (specN : ℕ → Finset ℝ) : Set ℝ :=
  { t : ℝ | ∃ N : ℕ, t ∈ specN N }

/-- Forward (Xi-zero → finite spectral witness) side of the analytic identity. -/
def ZeroToFiniteWitness (Xi : ℂ → ℂ) (specN : ℕ → Finset ℝ) : Prop :=
  ∀ s : ℂ, Xi s = 0 → ∃ N : ℕ, ∃ t : ℝ, t ∈ specN N ∧ s = criticalLinePoint t

/-- Backward (finite spectral witness → Xi-zero) side of the analytic identity. -/
def FiniteWitnessToZero (Xi : ℂ → ℂ) (specN : ℕ → Finset ℝ) : Prop :=
  ∀ N : ℕ, ∀ t : ℝ, t ∈ specN N → Xi (criticalLinePoint t) = 0

/-- Layer-2 analytic identity assumptions:
    finite truncation ladder + two-way Xi/spectral witness equivalence. -/
structure Layer2AnalyticIdentity (Xi : ℂ → ℂ) where
  specN : ℕ → Finset ℝ
  nested : ∀ N : ℕ, specN N ⊆ specN (N + 1)
  zero_to_finite : ZeroToFiniteWitness Xi specN
  finite_to_zero : FiniteWitnessToZero Xi specN

theorem mem_ladderSpec_of_mem_level
    {specN : ℕ → Finset ℝ} {N : ℕ} {t : ℝ}
    (ht : t ∈ specN N) :
    t ∈ ladderSpec specN := by
  exact ⟨N, ht⟩

theorem countUpTo_mono_spec
    {specA specB : Finset ℝ}
    (hAB : specA ⊆ specB)
    (T : ℝ) :
    countUpTo specA T ≤ countUpTo specB T := by
  classical
  unfold countUpTo
  refine Finset.card_le_card ?_
  intro x hx
  simp at hx ⊢
  exact ⟨hAB hx.1, hx.2⟩

theorem countUpTo_mono_ladder_level
    {Xi : ℂ → ℂ}
    (hL2 : Layer2AnalyticIdentity Xi)
    (N : ℕ)
    (T : ℝ) :
    countUpTo (hL2.specN N) T ≤ countUpTo (hL2.specN (N + 1)) T := by
  exact countUpTo_mono_spec (hL2.nested N) T

theorem zero_to_ladder_point
    {Xi : ℂ → ℂ}
    (hL2 : Layer2AnalyticIdentity Xi) :
    ∀ s : ℂ, Xi s = 0 → ∃ t : ℝ, t ∈ ladderSpec hL2.specN ∧ s = criticalLinePoint t := by
  intro s hs
  rcases hL2.zero_to_finite s hs with ⟨N, t, ht, hsEq⟩
  exact ⟨t, mem_ladderSpec_of_mem_level ht, hsEq⟩

theorem ladder_point_to_zero
    {Xi : ℂ → ℂ}
    (hL2 : Layer2AnalyticIdentity Xi) :
    ∀ t : ℝ, t ∈ ladderSpec hL2.specN → Xi (criticalLinePoint t) = 0 := by
  intro t ht
  rcases ht with ⟨N, htN⟩
  exact hL2.finite_to_zero N t htN

/-- Minimal Layer-2 RH trigger:
    only the forward witness direction is required for RH itself. -/
theorem rh_of_zero_to_finite_witness
    (Xi : ℂ → ℂ)
    (specN : ℕ → Finset ℝ)
    (hzero : ZeroToFiniteWitness Xi specN) :
    RiemannHypothesisXi Xi := by
  apply rh_of_zero_parameterization Xi
  intro s hs
  rcases hzero s hs with ⟨_N, t, _ht, hsEq⟩
  exact ⟨t, hsEq⟩

/-- Layer-2 closure: finite-ladder analytic identity induces exact Xi/spectrum bridge. -/
theorem spectralBridge_of_layer2
    (Xi : ℂ → ℂ)
    (hL2 : Layer2AnalyticIdentity Xi) :
    SpectralBridge Xi (ladderSpec hL2.specN) := by
  intro s
  constructor
  · intro hs
    exact zero_to_ladder_point hL2 s hs
  · intro hs
    rcases hs with ⟨t, ht, hsEq⟩
    subst hsEq
    exact ladder_point_to_zero hL2 t ht

/-- RH reduction from the Layer-2 analytic identity assumptions. -/
theorem rh_of_layer2_identity
    (Xi : ℂ → ℂ)
    (hL2 : Layer2AnalyticIdentity Xi) :
    RiemannHypothesisXi Xi := by
  exact bridge_implies_rh Xi (ladderSpec hL2.specN) (spectralBridge_of_layer2 Xi hL2)

end

end Gutoe.RiemannLayer2Identity
