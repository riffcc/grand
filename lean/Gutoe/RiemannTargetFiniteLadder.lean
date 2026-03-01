import Mathlib
import Gutoe.RiemannLimitBridge
import Gutoe.RiemannConvergenceTransfer
import Gutoe.RiemannFiniteXiModel
import Gutoe.RiemannFinalTarget

namespace Gutoe.RiemannTargetFiniteLadder

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

end

end Gutoe.RiemannTargetFiniteLadder

