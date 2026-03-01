import Mathlib
import Gutoe.RiemannCore
import Gutoe.RiemannBridge
import Gutoe.RiemannLayer2Identity

namespace Gutoe.RiemannLimitBridge

open Gutoe.RiemannCore
open Gutoe.RiemannBridge
open Gutoe.RiemannLayer2Identity

noncomputable section

/-- Set-view of a finite spectral level. -/
def levelSpecSet (specN : ℕ → Finset ℝ) (N : ℕ) : Set ℝ :=
  fun t => t ∈ specN N

theorem mem_ladder_of_level_mem
    {specN : ℕ → Finset ℝ}
    {N : ℕ} {t : ℝ}
    (ht : t ∈ specN N) :
    t ∈ ladderSpec specN := by
  exact ⟨N, ht⟩

/-- Layer-3 limit-transfer closure theorem:
    if every finite level has an exact bridge and every `Xi`-zero appears
    as a zero of some finite level function, RH follows for `Xi`. -/
theorem rh_of_limit_transfer
    (Xi : ℂ → ℂ)
    (XiN : ℕ → (ℂ → ℂ))
    (specN : ℕ → Finset ℝ)
    (hfiniteBridge : ∀ N : ℕ, SpectralBridge (XiN N) (levelSpecSet specN N))
    (hzeroForward : ∀ s : ℂ, Xi s = 0 → ∃ N : ℕ, XiN N s = 0) :
    RiemannHypothesisXi Xi := by
  apply rh_of_zero_parameterization Xi
  intro s hs
  rcases hzeroForward s hs with ⟨N, hsN⟩
  rcases (hfiniteBridge N s).1 hsN with ⟨t, ht, hsEq⟩
  exact ⟨t, hsEq⟩

/-- Two-way limit transfer yields an exact infinite bridge for `Xi`
    against the ladder-union spectrum. -/
theorem spectralBridge_of_limit_transfer
    (Xi : ℂ → ℂ)
    (XiN : ℕ → (ℂ → ℂ))
    (specN : ℕ → Finset ℝ)
    (hfiniteBridge : ∀ N : ℕ, SpectralBridge (XiN N) (levelSpecSet specN N))
    (hzeroForward : ∀ s : ℂ, Xi s = 0 → ∃ N : ℕ, XiN N s = 0)
    (hzeroBackward : ∀ N : ℕ, ∀ s : ℂ, XiN N s = 0 → Xi s = 0) :
    SpectralBridge Xi (ladderSpec specN) := by
  intro s
  constructor
  · intro hs
    rcases hzeroForward s hs with ⟨N, hsN⟩
    rcases (hfiniteBridge N s).1 hsN with ⟨t, ht, hsEq⟩
    exact ⟨t, mem_ladder_of_level_mem ht, hsEq⟩
  · intro hs
    rcases hs with ⟨t, ht, hsEq⟩
    rcases ht with ⟨N, htN⟩
    subst hsEq
    have hzN : XiN N (criticalLinePoint t) = 0 := by
      exact (hfiniteBridge N (criticalLinePoint t)).2 ⟨t, htN, rfl⟩
    exact hzeroBackward N (criticalLinePoint t) hzN

/-- RH closure via exact bridge obtained from two-way limit transfer. -/
theorem rh_of_exact_limit_bridge
    (Xi : ℂ → ℂ)
    (XiN : ℕ → (ℂ → ℂ))
    (specN : ℕ → Finset ℝ)
    (hfiniteBridge : ∀ N : ℕ, SpectralBridge (XiN N) (levelSpecSet specN N))
    (hzeroForward : ∀ s : ℂ, Xi s = 0 → ∃ N : ℕ, XiN N s = 0)
    (hzeroBackward : ∀ N : ℕ, ∀ s : ℂ, XiN N s = 0 → Xi s = 0) :
    RiemannHypothesisXi Xi := by
  exact bridge_implies_rh Xi (ladderSpec specN)
    (spectralBridge_of_limit_transfer Xi XiN specN hfiniteBridge hzeroForward hzeroBackward)

end

end Gutoe.RiemannLimitBridge

