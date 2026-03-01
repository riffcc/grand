import Mathlib
import Gutoe.RiemannCore
import Gutoe.RiemannLimitBridge
import Gutoe.RiemannConvergenceTransfer
import Gutoe.RiemannFiniteXiModel
import Gutoe.RiemannFinalTarget

namespace Gutoe.RiemannTargetFiniteLadder

open Gutoe.RiemannCore
open Gutoe.RiemannLimitBridge
open Gutoe.RiemannConvergenceTransfer
open Gutoe.RiemannFiniteXiModel
open Gutoe.RiemannFinalTarget

noncomputable section

/-- Finite spectral product ladder used as approximation family for `XiTarget`. -/
def XiFiniteLadder (specN : ℕ → Finset ℝ) : ℕ → (ℂ → ℂ) :=
  fun N => XiFinite (specN N)

/-- Zero tolerance profile. -/
def tolZero : ℕ → ℝ := fun _ => 0

theorem zeroTol_tolZero : zeroTol tolZero := by
  intro N
  simp [tolZero]

/-- Finite bridge family is automatic for the `XiFinite` ladder. -/
theorem finiteBridgeFamily_XiFiniteLadder
    (specN : ℕ → Finset ℝ) :
    FiniteBridgeFamily (XiFiniteLadder specN) specN := by
  intro N
  simpa [XiFiniteLadder, levelSpecSet, finiteSpecSet]
    using finiteXi_spectralBridge (specN N)

/-- Rigidity is automatic at zero tolerance for any ladder. -/
theorem rigidity_tolZero
    (XiN : ℕ → (ℂ → ℂ)) :
    SpectralRigidity XiN tolZero := by
  intro N s hs
  have hnorm0 : ‖XiN N s‖ = 0 := le_antisymm hs (norm_nonneg _)
  exact norm_eq_zero.mp hnorm0

/-- Reduced final contract:
    all structural obligations are discharged, leaving only zero-capture of `XiTarget`
    by the finite spectral ladder. -/
structure XiTargetFiniteLadderContract where
  specN : ℕ → Finset ℝ
  approxZero : ApproxZeroConvergence XiTarget (XiFiniteLadder specN) tolZero

/-- Build the full convergence-transfer contract from the reduced finite-ladder contract. -/
def toConvergenceTransferContract
    (hC : XiTargetFiniteLadderContract) :
    RHConvergenceTransferContract XiTarget where
  XiN := XiFiniteLadder hC.specN
  specN := hC.specN
  finiteBridge := finiteBridgeFamily_XiFiniteLadder hC.specN
  tol := tolZero
  tolNonneg := zeroTol_tolZero
  approxZero := hC.approxZero
  rigidity := rigidity_tolZero (XiFiniteLadder hC.specN)

/-- Final RH closure from the reduced finite-ladder contract. -/
theorem mathlibRH_of_target_finite_ladder_contract
    (hC : XiTargetFiniteLadderContract) :
    RiemannHypothesis := by
  exact mathlibRH_of_contract (toConvergenceTransferContract hC)

/-- Explicit zero-capture form of the remaining gap:
    every `XiTarget` zero appears at some finite ladder level as a listed ordinate. -/
def XiTargetLadderZeroCapture (specN : ℕ → Finset ℝ) : Prop :=
  ∀ s : ℂ, XiTarget s = 0 → ∃ N : ℕ, ∃ t : ℝ, t ∈ specN N ∧ s = criticalLinePoint t

/-- For the `XiFinite` ladder at zero tolerance, `approxZero` is equivalent to
    exact ladder zero capture. -/
theorem approxZero_tolZero_iff_zeroCapture
    (specN : ℕ → Finset ℝ) :
    ApproxZeroConvergence XiTarget (XiFiniteLadder specN) tolZero
      ↔ XiTargetLadderZeroCapture specN := by
  constructor
  · intro hApprox
    intro s hs
    rcases hApprox s hs with ⟨N, hN⟩
    have hXiN0 : XiFinite (specN N) s = 0 := by
      have hnorm0 : ‖XiFinite (specN N) s‖ = 0 := le_antisymm hN (by simp [tolZero])
      exact norm_eq_zero.mp hnorm0
    rcases (XiFinite_zero_iff_exists (specN N) s).1 hXiN0 with ⟨t, ht, hsEq⟩
    exact ⟨N, t, ht, hsEq⟩
  · intro hCap
    intro s hs
    rcases hCap s hs with ⟨N, t, ht, hsEq⟩
    refine ⟨N, ?_⟩
    subst hsEq
    have hXiN0 : XiFinite (specN N) (criticalLinePoint t) = 0 :=
      XiFinite_zero_of_mem (specN N) ht
    simp [XiFiniteLadder, tolZero, hXiN0]

/-- Direct final closure from explicit ladder zero-capture.
    This is equivalent to the zero-tolerance convergence contract in this lane. -/
theorem mathlibRH_of_target_ladder_zero_capture
    (specN : ℕ → Finset ℝ)
    (hCap : XiTargetLadderZeroCapture specN) :
    RiemannHypothesis := by
  intro s hs htriv h1
  have hsXi : XiTarget s = 0 := nontrivialZeroTransferToXiTarget s hs htriv h1
  rcases hCap s hsXi with ⟨N, t, ht, hsEq⟩
  simpa [hsEq, onCriticalLine, criticalLinePoint_re] using (criticalLinePoint_re t)

end

end Gutoe.RiemannTargetFiniteLadder
