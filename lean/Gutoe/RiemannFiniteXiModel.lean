import Mathlib
import Gutoe.RiemannCore
import Gutoe.RiemannBridge
import Gutoe.RiemannLimitBridge
import Gutoe.RiemannConvergenceTransfer

namespace Gutoe.RiemannFiniteXiModel

open Gutoe.RiemannCore
open Gutoe.RiemannBridge
open Gutoe.RiemannLimitBridge
open Gutoe.RiemannConvergenceTransfer

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

/-- Constant finite-level family used to instantiate transfer contracts for `XiFinite`. -/
def XiFiniteConst (spec : Finset ℝ) : ℕ → (ℂ → ℂ) := fun _ => XiFinite spec

/-- Constant spectral ladder used with `XiFiniteConst`. -/
def specConst (spec : Finset ℝ) : ℕ → Finset ℝ := fun _ => spec

/-- Zero tolerance profile. -/
def tolZero : ℕ → ℝ := fun _ => 0

theorem finiteBridgeFamily_XiFiniteConst
    (spec : Finset ℝ) :
    FiniteBridgeFamily (XiFiniteConst spec) (specConst spec) := by
  intro N
  simpa [XiFiniteConst, specConst, levelSpecSet, finiteSpecSet]
    using finiteXi_spectralBridge spec

theorem zeroTol_tolZero : zeroTol tolZero := by
  intro N
  simp [tolZero]

theorem approxZero_XiFiniteConst
    (spec : Finset ℝ) :
    ApproxZeroConvergence (XiFinite spec) (XiFiniteConst spec) tolZero := by
  intro s hs
  refine ⟨0, ?_⟩
  simp [XiFiniteConst, tolZero, hs]

theorem rigidity_XiFiniteConst
    (spec : Finset ℝ) :
    SpectralRigidity (XiFiniteConst spec) tolZero := by
  intro N s hs
  have hnorm0 : ‖XiFiniteConst spec N s‖ = 0 := by
    exact le_antisymm hs (norm_nonneg _)
  exact norm_eq_zero.mp hnorm0

/-- Explicit convergence-transfer contract for `XiFinite` (all obligations discharged). -/
def XiFiniteConvergenceContract
    (spec : Finset ℝ) :
    RHConvergenceTransferContract (XiFinite spec) where
  XiN := XiFiniteConst spec
  specN := specConst spec
  finiteBridge := finiteBridgeFamily_XiFiniteConst spec
  tol := tolZero
  tolNonneg := zeroTol_tolZero
  approxZero := approxZero_XiFiniteConst spec
  rigidity := rigidity_XiFiniteConst spec

/-- RH for `XiFinite` via convergence-transfer contract path. -/
theorem rh_XiFinite_via_convergence_contract
    (spec : Finset ℝ) :
    RiemannHypothesisXi (XiFinite spec) := by
  exact rh_of_convergence_transfer_contract (XiFinite spec) (XiFiniteConvergenceContract spec)

end

end Gutoe.RiemannFiniteXiModel
