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

/-- Reverse transfer target for the endgame:
    every `XiTarget`-zero is a nontrivial `riemannZeta`-zero. -/
def XiTargetZeroTransferToNontrivialZeta : Prop :=
  ∀ s : ℂ, XiTarget s = 0 →
    riemannZeta s = 0 ∧ (¬ ∃ n : ℕ, s = -2 * (n + 1)) ∧ s ≠ 1

/-- Atomic nonvanishing obligations that imply reverse transfer:
    `XiTarget` is nonzero at `0`, at `1`, and at every trivial-zero location. -/
def XiTargetNonvanishingObligations : Prop :=
  XiTarget 0 ≠ 0 ∧ XiTarget 1 ≠ 0 ∧ ∀ n : ℕ, XiTarget (-2 * (n + 1)) ≠ 0

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

/-- Reverse transfer derived from explicit nonvanishing obligations. -/
theorem xiTargetZeroTransferToNontrivialZeta_of_nonvanishing
    (hNv : XiTargetNonvanishingObligations) :
    XiTargetZeroTransferToNontrivialZeta := by
  intro s hsXi
  have hs0 : s ≠ 0 := by
    intro hs0Eq
    exact hNv.1 (by simpa [hs0Eq] using hsXi)
  have hsZeta : riemannZeta s = 0 := by
    have hnum : completedRiemannZeta s = 0 := by
      simpa [XiTarget] using hsXi
    have hzDiv : completedRiemannZeta s / Complex.Gammaℝ s = 0 := by
      simp [hnum]
    calc
      riemannZeta s = completedRiemannZeta s / Complex.Gammaℝ s := by
        simpa using (riemannZeta_def_of_ne_zero (s := s) hs0)
      _ = 0 := hzDiv
  have htriv : ¬ ∃ n : ℕ, s = -2 * (n + 1) := by
    intro h
    rcases h with ⟨n, hsEq⟩
    exact (hNv.2.2 n) (by simpa [hsEq] using hsXi)
  have h1 : s ≠ 1 := by
    intro hs1
    exact hNv.2.1 (by simpa [hs1] using hsXi)
  exact ⟨hsZeta, htriv, h1⟩

/-- `XiTarget = completedRiemannZeta` is nonvanishing at `0`, `1`, and all trivial-zero locations. -/
theorem xiTargetNonvanishingObligations :
    XiTargetNonvanishingObligations := by
  have h1 : XiTarget 1 ≠ 0 := by
    intro hXi1
    have hz1 : riemannZeta (1 : ℂ) = 0 := by
      rw [riemannZeta_def_of_ne_zero (s := (1 : ℂ)) one_ne_zero]
      simpa [XiTarget, hXi1]
    exact riemannZeta_one_ne_zero hz1
  have h0 : XiTarget 0 ≠ 0 := by
    -- Functional equation: `Λ(1) = Λ(0)`.
    have hfe0 : XiTarget 1 = XiTarget 0 := by
      simpa [XiTarget] using (completedRiemannZeta_one_sub (0 : ℂ))
    simpa [hfe0] using h1
  have htriv : ∀ n : ℕ, XiTarget (-2 * (n + 1)) ≠ 0 := by
    intro n
    let s : ℂ := -2 * (n + 1)
    let u : ℂ := (1 + 2 * (n + 1) : ℂ)
    have hu_re_ge1 : (1 : ℝ) ≤ u.re := by
      have hu_re : u.re = (1 + 2 * (n + 1) : ℝ) := by
        norm_num [u]
      linarith
    have hz_u_ne : riemannZeta u ≠ 0 := riemannZeta_ne_zero_of_one_le_re hu_re_ge1
    have hu_ne_zero : u ≠ 0 := by
      have hu_re_pos : (0 : ℝ) < u.re := by
        have hu_re : u.re = (1 + 2 * (n + 1) : ℝ) := by
          norm_num [u]
        linarith
      intro hu0
      have hu_re_zero : u.re = 0 := by simpa [hu0]
      linarith
    have hXi_u_ne : XiTarget u ≠ 0 := by
      intro hXi_u
      have hz_u_zero : riemannZeta u = 0 := by
        rw [riemannZeta_def_of_ne_zero (s := u) hu_ne_zero]
        have hnum : completedRiemannZeta u = 0 := by
          simpa [XiTarget] using hXi_u
        simp [hnum]
      exact hz_u_ne hz_u_zero
    have hu_eq : u = 1 - s := by
      apply Complex.ext <;> norm_num [u, s]
    have hs_eq : XiTarget (1 - s) = XiTarget s := by
      simpa [XiTarget] using (completedRiemannZeta_one_sub s)
    have hs_ne : XiTarget s ≠ 0 := by
      have hOneMinus_ne : XiTarget (1 - s) ≠ 0 := by
        simpa [hu_eq] using hXi_u_ne
      simpa [hs_eq] using hOneMinus_ne
    simpa [s] using hs_ne
  exact ⟨h0, h1, htriv⟩

/-- Reverse transfer discharged by the explicit nonvanishing theorem for `XiTarget`. -/
theorem xiTargetZeroTransferToNontrivialZeta :
    XiTargetZeroTransferToNontrivialZeta :=
  xiTargetZeroTransferToNontrivialZeta_of_nonvanishing xiTargetNonvanishingObligations

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
