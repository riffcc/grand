import Mathlib
import Mathlib.NumberTheory.LSeries.RiemannZeta
import Gutoe.RiemannCore
import Gutoe.RiemannFinalTarget

namespace Gutoe.RiemannHurwitzKernel

open Gutoe.RiemannCore
open Gutoe.RiemannFinalTarget

noncomputable section
open scoped Topology

/-- Abstract Hurwitz output surface for a ladder `XiN` converging to a target
function `Xi`: every target zero is approximated by finite-level zeros. -/
def HurwitzZeroApproxTransfer
    (Xi : ℂ → ℂ)
    (XiN : ℕ → (ℂ → ℂ)) : Prop :=
  ∀ s : ℂ, Xi s = 0 → ∀ ε : ℝ, 0 < ε →
    ∃ N : ℕ, ∃ z : ℂ, XiN N z = 0 ∧ ‖s - z‖ < ε

/-- Finite-level zero witness interface: each zero of `XiN N` is represented by
a critical-line point indexed by `specN N`. -/
def FiniteZeroWitness
    (XiN : ℕ → (ℂ → ℂ))
    (specN : ℕ → Finset ℝ) : Prop :=
  ∀ N : ℕ, ∀ z : ℂ, XiN N z = 0 →
    ∃ t : ℝ, t ∈ specN N ∧ z = criticalLinePoint t

/-- Approximate critical-line capture surface induced by finite spectral sets. -/
def ApproximateCriticalCapture
    (Xi : ℂ → ℂ)
    (specN : ℕ → Finset ℝ) : Prop :=
  ∀ s : ℂ, Xi s = 0 → ∀ ε : ℝ, 0 < ε →
    ∃ N : ℕ, ∃ t : ℝ, t ∈ specN N ∧ ‖s - criticalLinePoint t‖ < ε

/-- Engineering kernel for Hurwitz-style transfer.
This isolates Step-2 as a reusable module obligation independent of RH plumbing. -/
def HurwitzZeroApproxKernel
    (Xi : ℂ → ℂ)
    (XiN : ℕ → (ℂ → ℂ)) : Prop :=
  TendstoLocallyUniformly XiN Xi (Filter.atTop : Filter ℕ) →
    (∀ N : ℕ, Differentiable ℂ (XiN N)) →
    (∀ s : ℂ, Xi s = 0 → DifferentiableAt ℂ Xi s) →
    (∀ N : ℕ, XiN N 0 ≠ 0) →
    HurwitzZeroApproxTransfer Xi XiN

/-- Kernel instantiation theorem: once a Hurwitz kernel is provided for `(Xi, XiN)`
and convergence/regularity hypotheses are supplied, one gets concrete transfer. -/
theorem hurwitzTransfer_of_kernel
    {Xi : ℂ → ℂ}
    {XiN : ℕ → (ℂ → ℂ)}
    (hKernel : HurwitzZeroApproxKernel Xi XiN)
    (hconv : TendstoLocallyUniformly XiN Xi (Filter.atTop : Filter ℕ))
    (hDiffN : ∀ N : ℕ, Differentiable ℂ (XiN N))
    (hDiffXi : ∀ s : ℂ, Xi s = 0 → DifferentiableAt ℂ Xi s)
    (hNontriv : ∀ N : ℕ, XiN N 0 ≠ 0) :
    HurwitzZeroApproxTransfer Xi XiN := by
  exact hKernel hconv hDiffN hDiffXi hNontriv

/-- Hurwitz transfer + finite zero witnesses implies approximate critical-line
capture through the ladder `specN`. -/
theorem approximateCriticalCapture_of_hurwitzTransfer_and_witness
    {Xi : ℂ → ℂ}
    {XiN : ℕ → (ℂ → ℂ)}
    {specN : ℕ → Finset ℝ}
    (hHurwitz : HurwitzZeroApproxTransfer Xi XiN)
    (hWitness : FiniteZeroWitness XiN specN) :
    ApproximateCriticalCapture Xi specN := by
  intro s hsXi ε hε
  rcases hHurwitz s hsXi ε hε with ⟨N, z, hz0, hdist⟩
  rcases hWitness N z hz0 with ⟨t, htN, hzEq⟩
  refine ⟨N, t, htN, ?_⟩
  simpa [hzEq] using hdist

/-- `XiTarget = completedRiemannZeta` is differentiable at each of its zeros. -/
theorem differentiableAt_XiTarget_of_zero
    {s : ℂ} (hs : XiTarget s = 0) :
    DifferentiableAt ℂ XiTarget s := by
  have h0 : s ≠ 0 := by
    intro hs0
    exact xiTargetNonvanishingObligations.1 (by simpa [hs0] using hs)
  have h1 : s ≠ 1 := by
    intro hs1
    exact xiTargetNonvanishingObligations.2.1 (by simpa [hs1] using hs)
  simpa [XiTarget] using differentiableAt_completedZeta h0 h1

end

end Gutoe.RiemannHurwitzKernel
