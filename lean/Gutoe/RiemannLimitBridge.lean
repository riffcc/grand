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

/-- Obligation A: every finite level has an exact XiN↔Spec bridge. -/
def FiniteBridgeFamily
    (XiN : ℕ → (ℂ → ℂ))
    (specN : ℕ → Finset ℝ) : Prop :=
  ∀ N : ℕ, SpectralBridge (XiN N) (levelSpecSet specN N)

/-- Obligation B(fwd): every target-Xi zero appears at some finite level. -/
def ZeroForwardTransfer
    (Xi : ℂ → ℂ)
    (XiN : ℕ → (ℂ → ℂ)) : Prop :=
  ∀ s : ℂ, Xi s = 0 → ∃ N : ℕ, XiN N s = 0

/-- Optional strengthening B(back): finite-level zeros are genuine Xi zeros. -/
def ZeroBackwardTransfer
    (Xi : ℂ → ℂ)
    (XiN : ℕ → (ℂ → ℂ)) : Prop :=
  ∀ N : ℕ, ∀ s : ℂ, XiN N s = 0 → Xi s = 0

/-- Canonical endgame contract:
    fill this record to close RH-for-Xi in this lane. -/
structure RHLimitTransferContract (Xi : ℂ → ℂ) where
  XiN : ℕ → (ℂ → ℂ)
  specN : ℕ → Finset ℝ
  finiteBridge : FiniteBridgeFamily XiN specN
  zeroForward : ZeroForwardTransfer Xi XiN

/-- Strong endgame contract: adds backward transfer for exact infinite bridge. -/
structure RHExactLimitTransferContract (Xi : ℂ → ℂ) extends RHLimitTransferContract Xi where
  zeroBackward : ZeroBackwardTransfer Xi toRHLimitTransferContract.XiN

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
    (hfiniteBridge : FiniteBridgeFamily XiN specN)
    (hzeroForward : ZeroForwardTransfer Xi XiN) :
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
    (hfiniteBridge : FiniteBridgeFamily XiN specN)
    (hzeroForward : ZeroForwardTransfer Xi XiN)
    (hzeroBackward : ZeroBackwardTransfer Xi XiN) :
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
    (hfiniteBridge : FiniteBridgeFamily XiN specN)
    (hzeroForward : ZeroForwardTransfer Xi XiN)
    (hzeroBackward : ZeroBackwardTransfer Xi XiN) :
    RiemannHypothesisXi Xi := by
  exact bridge_implies_rh Xi (ladderSpec specN)
    (spectralBridge_of_limit_transfer Xi XiN specN hfiniteBridge hzeroForward hzeroBackward)

/-- Contract-packaged RH closure theorem (minimal endgame). -/
theorem rh_of_limit_transfer_contract
    (Xi : ℂ → ℂ)
    (hC : RHLimitTransferContract Xi) :
    RiemannHypothesisXi Xi := by
  exact rh_of_limit_transfer Xi hC.XiN hC.specN hC.finiteBridge hC.zeroForward

/-- Contract-packaged exact bridge theorem (strong endgame). -/
theorem spectralBridge_of_exact_limit_transfer_contract
    (Xi : ℂ → ℂ)
    (hC : RHExactLimitTransferContract Xi) :
    SpectralBridge Xi (ladderSpec hC.toRHLimitTransferContract.specN) := by
  exact spectralBridge_of_limit_transfer Xi
    hC.toRHLimitTransferContract.XiN
    hC.toRHLimitTransferContract.specN
    hC.toRHLimitTransferContract.finiteBridge
    hC.toRHLimitTransferContract.zeroForward
    hC.zeroBackward

end

end Gutoe.RiemannLimitBridge
