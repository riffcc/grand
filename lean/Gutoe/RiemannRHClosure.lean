import Mathlib
import Gutoe.RiemannCore
import Gutoe.RiemannSelfAdjoint
import Gutoe.RiemannBridge
import Gutoe.RiemannCounting
import Gutoe.RiemannLayer2Identity
import Gutoe.RiemannLimitBridge

namespace Gutoe.RiemannRHClosure

open Gutoe.RiemannCore
open Gutoe.RiemannSelfAdjoint
open Gutoe.RiemannBridge
open Gutoe.RiemannCounting
open Gutoe.RiemannLayer2Identity
open Gutoe.RiemannLimitBridge

/-- Program-level assumptions required for the RH closure reduction lane.
    This is intentionally explicit so no hidden assumptions are smuggled in. -/
structure RHProgramAssumptions (Xi : ℂ → ℂ) where
  n : ℕ
  selfAdjointFinite : finiteSelfAdjoint (structuralRiemannMatrix n)
  spec : Set ℝ
  bridge : SpectralBridge Xi spec

/-- Main closure reduction theorem:
    if the exact spectral bridge is established, RH follows for `Xi`. -/
theorem rh_of_program_assumptions
    (Xi : ℂ → ℂ)
    (hProg : RHProgramAssumptions Xi) :
    RiemannHypothesisXi Xi := by
  exact bridge_implies_rh Xi hProg.spec hProg.bridge

/-- Constructive witness: finite self-adjointness of the structural matrix
    is always available for any matrix size. -/
theorem finite_selfAdjoint_witness (n : ℕ) :
    finiteSelfAdjoint (structuralRiemannMatrix n) := by
  exact structuralRiemannMatrix_finiteSelfAdjoint n

/-- Packaging theorem: if an exact bridge is supplied, one obtains RH from
    an explicit program record carrying all named obligations. -/
theorem rh_from_explicit_bridge
    (Xi : ℂ → ℂ)
    (n : ℕ)
    (spec : Set ℝ)
    (hbridge : SpectralBridge Xi spec) :
    RiemannHypothesisXi Xi := by
  let hProg : RHProgramAssumptions Xi :=
    { n := n
      selfAdjointFinite := finite_selfAdjoint_witness n
      spec := spec
      bridge := hbridge }
  exact rh_of_program_assumptions Xi hProg

/-- Counting-side compatibility transport theorem for finite proxy sets. -/
theorem counting_transport
    (specA specB specC : Finset ℝ)
    (hAB : FiniteCountingMatch specA specB)
    (hBC : FiniteCountingMatch specB specC) :
    FiniteCountingMatch specA specC := by
  exact finiteCountingMatch_trans hAB hBC

/-- Layer-2 bridge packaging theorem:
    an explicit finite-ladder analytic identity yields RH directly. -/
theorem rh_from_layer2_analytic_identity
    (Xi : ℂ → ℂ)
    (hL2 : Layer2AnalyticIdentity Xi) :
    RiemannHypothesisXi Xi := by
  exact rh_of_layer2_identity Xi hL2

/-- Minimal Layer-2 closure form:
    proving only zero-to-finite-witness already forces RH. -/
theorem rh_from_layer2_forward_only
    (Xi : ℂ → ℂ)
    (specN : ℕ → Finset ℝ)
    (hzero : ZeroToFiniteWitness Xi specN) :
    RiemannHypothesisXi Xi := by
  exact rh_of_zero_to_finite_witness Xi specN hzero

/-- Layer-3 limit transfer packaging:
    RH follows once finite exact bridges + zero-forward transfer are supplied. -/
theorem rh_from_limit_transfer
    (Xi : ℂ → ℂ)
    (XiN : ℕ → (ℂ → ℂ))
    (specN : ℕ → Finset ℝ)
    (hfiniteBridge : FiniteBridgeFamily XiN specN)
    (hzeroForward : ZeroForwardTransfer Xi XiN) :
    RiemannHypothesisXi Xi := by
  exact rh_of_limit_transfer Xi XiN specN hfiniteBridge hzeroForward

/-- Contract API: one record carries the named remaining obligations. -/
theorem rh_from_limit_transfer_contract
    (Xi : ℂ → ℂ)
    (hC : RHLimitTransferContract Xi) :
    RiemannHypothesisXi Xi := by
  exact rh_of_limit_transfer_contract Xi hC

end Gutoe.RiemannRHClosure
