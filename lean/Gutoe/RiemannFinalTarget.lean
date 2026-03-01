import Mathlib
import Mathlib.NumberTheory.LSeries.RiemannZeta
import Gutoe.RiemannCore
import Gutoe.RiemannConvergenceTransfer

namespace Gutoe.RiemannFinalTarget

open Gutoe.RiemannCore
open Gutoe.RiemannConvergenceTransfer

noncomputable section

/-- Final analytic target in this lane: entire completed zeta with poles removed. -/
def XiTarget : ℂ → ℂ := completedRiemannZeta₀

/-- RH-for-`XiTarget` from convergence-transfer contract. -/
theorem rhXiTarget_of_contract
    (hC : RHConvergenceTransferContract XiTarget) :
    RiemannHypothesisXi XiTarget := by
  exact rh_of_convergence_transfer_contract XiTarget hC

/-- Transfer obligation from nontrivial `ζ`-zeros to `XiTarget`-zeros. -/
def NontrivialZeroTransferToXiTarget : Prop :=
  ∀ s : ℂ, riemannZeta s = 0 → (¬ ∃ n : ℕ, s = -2 * (n + 1)) → s ≠ 1 → XiTarget s = 0

/-- Final closure theorem surface:
    contract + nontrivial-zero transfer implies Mathlib's RH statement. -/
theorem mathlibRH_of_contract_and_transfer
    (hC : RHConvergenceTransferContract XiTarget)
    (hTransfer : NontrivialZeroTransferToXiTarget) :
    RiemannHypothesis := by
  intro s hs htriv h1
  have hsXi : XiTarget s = 0 := hTransfer s hs htriv h1
  exact rhXiTarget_of_contract hC s hsXi

end

end Gutoe.RiemannFinalTarget

