import Mathlib
import Gutoe.RiemannCore
import Gutoe.RiemannSelfAdjoint

namespace Gutoe.RiemannBridge

open Gutoe.RiemannCore

/-- Spectral bridge: zeros of `Xi` are exactly critical-line points with parameter in `spec`. -/
def SpectralBridge (Xi : ℂ → ℂ) (spec : Set ℝ) : Prop :=
  ∀ s : ℂ, Xi s = 0 ↔ ∃ t : ℝ, t ∈ spec ∧ s = criticalLinePoint t

theorem bridge_implies_zero_parameterization
    (Xi : ℂ → ℂ) (spec : Set ℝ)
    (hbridge : SpectralBridge Xi spec) :
    ∀ s : ℂ, Xi s = 0 → ∃ t : ℝ, s = criticalLinePoint t := by
  intro s hs
  rcases (hbridge s).1 hs with ⟨t, _ht, hsEq⟩
  exact ⟨t, hsEq⟩

/-- Core reduction theorem: exact spectral bridge implies RH for `Xi`. -/
theorem bridge_implies_rh
    (Xi : ℂ → ℂ) (spec : Set ℝ)
    (hbridge : SpectralBridge Xi spec) :
    RiemannHypothesisXi Xi := by
  apply rh_of_zero_parameterization Xi
  exact bridge_implies_zero_parameterization Xi spec hbridge

end Gutoe.RiemannBridge
