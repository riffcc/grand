import Mathlib
import Gutoe.RiemannCore
import Gutoe.RiemannBridge

namespace Gutoe.RiemannFiniteXiModel

open Gutoe.RiemannCore
open Gutoe.RiemannBridge

noncomputable section

/-- Finite spectral set as a `Set`. -/
def finiteSpecSet (spec : Finset ℝ) : Set ℝ := fun t => t ∈ spec

/-- Finite spectral Xi model:
    one zero factor per spectral ordinate. -/
def XiFinite (spec : Finset ℝ) : ℂ → ℂ :=
  fun s => Finset.prod spec (fun t => (s - criticalLinePoint t : ℂ))

theorem XiFinite_zero_of_mem
    (spec : Finset ℝ) {t : ℝ}
    (ht : t ∈ spec) :
    XiFinite spec (criticalLinePoint t) = 0 := by
  classical
  unfold XiFinite
  refine Finset.prod_eq_zero_iff.mpr ?_
  exact ⟨t, ht, by simp [criticalLinePoint]⟩

theorem XiFinite_zero_iff_exists
    (spec : Finset ℝ) (s : ℂ) :
    XiFinite spec s = 0 ↔ ∃ t : ℝ, t ∈ spec ∧ s = criticalLinePoint t := by
  classical
  constructor
  · intro hs
    unfold XiFinite at hs
    rcases (Finset.prod_eq_zero_iff.mp hs) with ⟨t, ht, hfac⟩
    refine ⟨t, ht, ?_⟩
    exact sub_eq_zero.mp hfac
  · rintro ⟨t, ht, rfl⟩
    exact XiFinite_zero_of_mem spec ht

/-- Exact finite bridge theorem for the explicit finite Xi model. -/
theorem finiteXi_spectralBridge
    (spec : Finset ℝ) :
    SpectralBridge (XiFinite spec) (finiteSpecSet spec) := by
  intro s
  constructor
  · intro hs
    rcases (XiFinite_zero_iff_exists spec s).1 hs with ⟨t, ht, hsEq⟩
    exact ⟨t, ht, hsEq⟩
  · intro hs
    rcases hs with ⟨t, ht, hsEq⟩
    rw [hsEq]
    exact XiFinite_zero_of_mem spec ht

/-- RH holds for the explicit finite Xi model by exact bridge. -/
theorem rh_XiFinite
    (spec : Finset ℝ) :
    RiemannHypothesisXi (XiFinite spec) := by
  exact bridge_implies_rh (XiFinite spec) (finiteSpecSet spec) (finiteXi_spectralBridge spec)

end

end Gutoe.RiemannFiniteXiModel
