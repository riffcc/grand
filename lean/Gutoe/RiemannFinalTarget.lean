import Mathlib
import Mathlib.NumberTheory.LSeries.RiemannZeta
import Mathlib.Analysis.SpecialFunctions.Gamma.Deligne
import Gutoe.RiemannCore
import Gutoe.RiemannConvergenceTransfer

namespace Gutoe.RiemannFinalTarget

open Gutoe.RiemannCore
open Gutoe.RiemannConvergenceTransfer

noncomputable section

/-- Final analytic target in this lane: entire completed zeta with poles removed. -/
def XiTarget : ℂ → ℂ := completedRiemannZeta

/-- RH-for-`XiTarget` from convergence-transfer contract. -/
theorem rhXiTarget_of_contract
    (hC : RHConvergenceTransferContract XiTarget) :
    RiemannHypothesisXi XiTarget := by
  exact rh_of_convergence_transfer_contract XiTarget hC

/-- Transfer obligation from nontrivial `ζ`-zeros to `XiTarget`-zeros. -/
def NontrivialZeroTransferToXiTarget : Prop :=
  ∀ s : ℂ, riemannZeta s = 0 → (¬ ∃ n : ℕ, s = -2 * (n + 1)) → s ≠ 1 → XiTarget s = 0

/-- Nontrivial-zero transfer is derivable for `XiTarget = completedRiemannZeta`. -/
theorem nontrivialZeroTransferToXiTarget :
    NontrivialZeroTransferToXiTarget := by
  intro s hs htriv _h1
  have hs0 : s ≠ 0 := by
    intro hs0Eq
    have hz0 : riemannZeta 0 = 0 := by simpa [hs0Eq] using hs
    have hz0_ne : riemannZeta 0 ≠ 0 := by
      norm_num [riemannZeta_zero]
    exact hz0_ne hz0
  have hGammaR_ne : Complex.Gammaℝ s ≠ 0 := by
    intro hGammaR0
    rcases (Complex.Gammaℝ_eq_zero_iff).1 hGammaR0 with ⟨n, hn⟩
    cases n with
    | zero =>
        have hsEq0 : s = 0 := by simpa using hn
        exact hs0 hsEq0
    | succ k =>
        have hsEven : s = -2 * (k + 1) := by
          simpa [Nat.succ_eq_add_one, mul_comm, mul_left_comm, mul_assoc] using hn
        exact htriv ⟨k, hsEven⟩
  have hzDiv : completedRiemannZeta s / Complex.Gammaℝ s = 0 := by
    calc
      completedRiemannZeta s / Complex.Gammaℝ s = riemannZeta s := by
        simpa using (riemannZeta_def_of_ne_zero (s := s) hs0).symm
      _ = 0 := hs
  rcases (div_eq_zero_iff).1 hzDiv with hnum | hden
  · exact hnum
  · exact (hGammaR_ne hden).elim

/-- Final closure theorem surface:
    contract + nontrivial-zero transfer implies Mathlib's RH statement. -/
theorem mathlibRH_of_contract_and_transfer
    (hC : RHConvergenceTransferContract XiTarget)
    (hTransfer : NontrivialZeroTransferToXiTarget) :
    RiemannHypothesis := by
  intro s hs htriv h1
  have hsXi : XiTarget s = 0 := hTransfer s hs htriv h1
  exact rhXiTarget_of_contract hC s hsXi

/-- Final closure with transfer discharged. Remaining assumption is only the
    convergence-transfer contract at the target function. -/
theorem mathlibRH_of_contract
    (hC : RHConvergenceTransferContract XiTarget) :
    RiemannHypothesis := by
  exact mathlibRH_of_contract_and_transfer hC nontrivialZeroTransferToXiTarget

end

end Gutoe.RiemannFinalTarget
