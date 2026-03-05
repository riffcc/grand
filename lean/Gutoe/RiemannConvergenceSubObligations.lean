import Mathlib
import Gutoe.RiemannFinalTarget
import Gutoe.RiemannFiniteXiModel
import Gutoe.RiemannHurwitzKernel
import Gutoe.RiemannOperatorLadder

namespace Gutoe.RiemannConvergenceSubObligations

open Gutoe.RiemannFinalTarget
open Gutoe.RiemannFiniteXiModel
open Gutoe.RiemannHurwitzKernel
open Gutoe.RiemannOperatorLadder

noncomputable section
open scoped Topology

/-- Sub-obligation package for the Hadamard partial-product lane:
`XiN` is the finite Hadamard family, converging locally uniformly to `XiTarget`
(`completedRiemannZeta`) on compact subsets of `ℂ`. -/
structure HadamardXiNConvergence where
  XiN : ℕ → (ℂ → ℂ)
  isHadamardFamily : Prop
  convergesLocallyUniformly :
    TendstoLocallyUniformly XiN XiTarget (Filter.atTop : Filter ℕ)
  differentiable_XiN : ∀ N : ℕ, Differentiable ℂ (XiN N)
  nontrivialAtZero : ∀ N : ℕ, XiN N 0 ≠ 0

/-- Sub-obligation package for the Hurwitz transfer lane:
1. a zero-order witness interface for `XiTarget` zeros,
2. nearby-zero production for approximants under local-uniform convergence, and
3. a rigidity-upgrade slot specialized to the operator ladder, turning nearby
analytic transfer into exact finite-level operator zero capture. -/
structure HurwitzZeroTransfer where
  zeroOrder : ℂ → ℕ → Prop
  zeroOrder_of_zero : ∀ s : ℂ, XiTarget s = 0 → ∃ k : ℕ, zeroOrder s k
  eventually_zero_near :
    ∀ {XiN : ℕ → (ℂ → ℂ)},
      TendstoLocallyUniformly XiN XiTarget (Filter.atTop : Filter ℕ) →
      (∀ N : ℕ, Differentiable ℂ (XiN N)) →
      (∀ s : ℂ, XiTarget s = 0 → DifferentiableAt ℂ XiTarget s) →
      (∀ N : ℕ, XiN N 0 ≠ 0) →
      ∀ s : ℂ, ∀ k : ℕ, XiTarget s = 0 → zeroOrder s k →
        ∀ ε : ℝ, 0 < ε →
          ∃ N : ℕ, ∃ z : ℂ, XiN N z = 0 ∧ ‖s - z‖ < ε
  operatorExactZeroUpgrade :
    ∀ s : ℂ, XiTarget s = 0 → ∃ N : ℕ, operatorXiFiniteLadder N s = 0

/-- Main reducer for the RH convergence gap:
if the Hadamard convergence package and Hurwitz transfer package are available,
then `OperatorApproxZeroConvergence` follows. -/
theorem operatorApproxZeroConvergence_of_subObligations
    (_hHad : HadamardXiNConvergence)
    (hHurwitz : HurwitzZeroTransfer) :
    OperatorApproxZeroConvergence := by
  refine (operatorApproxZeroConvergence_iff_eventual_exact_zero).2 ?_
  intro s hs
  exact hHurwitz.operatorExactZeroUpgrade s hs

end

end Gutoe.RiemannConvergenceSubObligations
