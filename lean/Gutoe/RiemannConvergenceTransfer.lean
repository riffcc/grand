import Mathlib
import Gutoe.RiemannCore
import Gutoe.RiemannLimitBridge

namespace Gutoe.RiemannConvergenceTransfer

open Gutoe.RiemannCore
open Gutoe.RiemannLimitBridge

noncomputable section

/-- Numeric tolerance profile used to turn approximate zeros into exact finite-level zeros. -/
def zeroTol (tol : ℕ → ℝ) : Prop := ∀ N : ℕ, 0 ≤ tol N

/-- Analytic convergence obligation:
    every target-`Xi` zero is eventually (or at least once) within tolerance
    for some finite-level function. -/
def ApproxZeroConvergence
    (Xi : ℂ → ℂ)
    (XiN : ℕ → (ℂ → ℂ))
    (tol : ℕ → ℝ) : Prop :=
  ∀ s : ℂ, Xi s = 0 → ∃ N : ℕ, ‖XiN N s‖ ≤ tol N

/-- Spectral rigidity obligation:
    at each finite level, tolerance-small values are exact zeros. -/
def SpectralRigidity
    (XiN : ℕ → (ℂ → ℂ))
    (tol : ℕ → ℝ) : Prop :=
  ∀ N : ℕ, ∀ s : ℂ, ‖XiN N s‖ ≤ tol N → XiN N s = 0

/-- Core transfer theorem:
    convergence-to-tolerance + spectral rigidity implies ZeroForwardTransfer. -/
theorem zeroForward_of_convergence_and_rigidity
    (Xi : ℂ → ℂ)
    (XiN : ℕ → (ℂ → ℂ))
    (tol : ℕ → ℝ)
    (hApprox : ApproxZeroConvergence Xi XiN tol)
    (hRigid : SpectralRigidity XiN tol) :
    ZeroForwardTransfer Xi XiN := by
  intro s hs
  rcases hApprox s hs with ⟨N, hN⟩
  exact ⟨N, hRigid N s hN⟩

/-- Contract replacing direct zero-forward assumption by analytic convergence + rigidity. -/
structure RHConvergenceTransferContract (Xi : ℂ → ℂ) where
  XiN : ℕ → (ℂ → ℂ)
  specN : ℕ → Finset ℝ
  finiteBridge : FiniteBridgeFamily XiN specN
  tol : ℕ → ℝ
  tolNonneg : zeroTol tol
  approxZero : ApproxZeroConvergence Xi XiN tol
  rigidity : SpectralRigidity XiN tol

/-- Derived ZeroForwardTransfer from convergence contract. -/
theorem zeroForward_of_contract
    (Xi : ℂ → ℂ)
    (hC : RHConvergenceTransferContract Xi) :
    ZeroForwardTransfer Xi hC.XiN := by
  exact zeroForward_of_convergence_and_rigidity
    Xi hC.XiN hC.tol hC.approxZero hC.rigidity

/-- RH closure from convergence contract (no direct zero-forward axiom needed). -/
theorem rh_of_convergence_transfer_contract
    (Xi : ℂ → ℂ)
    (hC : RHConvergenceTransferContract Xi) :
    RiemannHypothesisXi Xi := by
  exact rh_of_limit_transfer
    Xi hC.XiN hC.specN hC.finiteBridge (zeroForward_of_contract Xi hC)

end

end Gutoe.RiemannConvergenceTransfer

