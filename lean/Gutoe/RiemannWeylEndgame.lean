import Mathlib
import Gutoe.RiemannCore
import Gutoe.RiemannTargetFiniteLadder
import Gutoe.RiemannConvergenceTransfer
import Gutoe.RiemannFinalTarget

namespace Gutoe.RiemannWeylEndgame

open Gutoe.RiemannCore
open Gutoe.RiemannTargetFiniteLadder
open Gutoe.RiemannConvergenceTransfer
open Gutoe.RiemannFinalTarget

noncomputable section

/-- `2π` helper used in the Riemann–von Mangoldt main term. -/
def twoPi : ℝ := 2 * Real.pi

/-- Main asymptotic term in the Riemann–von Mangoldt zero-counting formula. -/
def rvmMain (T : ℝ) : ℝ :=
  (T / twoPi) * Real.log (T / (twoPi * Real.exp 1))

/-- Big-Oh style envelope around the Riemann–von Mangoldt main term. -/
def RiemannVonMangoldtEnvelope (N : ℝ → ℝ) : Prop :=
  ∃ C T0 : ℝ, 0 ≤ C ∧ 0 ≤ T0 ∧
    ∀ T : ℝ, T0 ≤ T → |N T - rvmMain T| ≤ C * Real.log (T + 2)

/-- Herglotz-type positivity condition for an m-function candidate. -/
def HerglotzLike (m : ℂ → ℂ) : Prop :=
  ∀ z : ℂ, 0 < z.im → 0 ≤ (m z).im

/-- Abstract m-function identity placeholder: same analytic object up to equality. -/
def MFunctionIdentity (mH mXi : ℂ → ℂ) : Prop := mH = mXi

/-- Canonical critical-line pole-like support extracted from an m-function.
    This is now definitionally tied to the function, not carried as abstract data. -/
def poleSet (m : ℂ → ℂ) : Set ℝ := { t : ℝ | m (criticalLinePoint t) = 0 }

/-- Counting semantics for a spectrum-like point set:
    `N` is explicitly realized as cardinality of finite cuts. -/
structure StepCountingSemantics (N : ℝ → ℝ) (S : Set ℝ) where
  cut : ℝ → Finset ℝ
  cut_spec :
    ∀ T : ℝ, ∀ t : ℝ, t ∈ cut T ↔ t ∈ S ∧ 0 ≤ t ∧ t ≤ T
  count_def :
    ∀ T : ℝ, N T = ((cut T).card : ℝ)

/-- Endgame contract encoding the Weyl + m-function attack path.
    The remaining work is now decomposed into explicit semantic bridge obligations:
    zero→pole (`mXi`), while pole transfer and ordinate enumeration are now
    canonical consequences of `m_identity` and countability. -/
structure RiemannWeylIdentityContract where
  N_H : ℝ → ℝ
  N_xi : ℝ → ℝ
  mH : ℂ → ℂ
  mXi : ℂ → ℂ
  counting_H : StepCountingSemantics N_H (poleSet mH)
  counting_xi : StepCountingSemantics N_xi (poleSet mXi)
  weyl_H : RiemannVonMangoldtEnvelope N_H
  weyl_xi : RiemannVonMangoldtEnvelope N_xi
  count_exact : ∀ T : ℝ, 0 ≤ T → N_H T = N_xi T
  mH_herglotz : HerglotzLike mH
  mXi_herglotz : HerglotzLike mXi
  m_identity : MFunctionIdentity mH mXi
  poleH_countable : (poleSet mH).Countable
  zero_to_poleXi :
    ∀ s : ℂ, riemannZeta s = 0 →
      (¬ ∃ n : ℕ, s = -2 * (n + 1)) →
      s ≠ 1 →
      ∃ t : ℝ, s = criticalLinePoint t ∧ t ∈ poleSet mXi

/-- Explicit semantic bridge theorem:
    from zero→pole and canonical pole/countability semantics,
    obtain an ordinate enumeration. -/
theorem exists_ordinateEnumeration_of_semantic_bridge
    (hC : RiemannWeylIdentityContract) :
    ∃ ρ : ℕ → ℝ, RiemannNontrivialZeroOrdinateEnumeration ρ := by
  by_cases hne : (poleSet hC.mH).Nonempty
  · rcases hC.poleH_countable.exists_eq_range hne with ⟨ρ, hρ⟩
    refine ⟨ρ, ?_⟩
    intro s hs htriv h1
    rcases hC.zero_to_poleXi s hs htriv h1 with ⟨t, hsEq, htXi⟩
    have htH : t ∈ poleSet hC.mH := by
      change hC.mH (criticalLinePoint t) = 0
      rw [hC.m_identity]
      exact htXi
    have htRange : t ∈ Set.range ρ := by simpa [hρ] using htH
    rcases htRange with ⟨n, hn⟩
    refine ⟨n, ?_⟩
    simpa [hn] using hsEq
  · refine ⟨fun _ => 0, ?_⟩
    intro s hs htriv h1
    exfalso
    rcases hC.zero_to_poleXi s hs htriv h1 with ⟨t, _hsEq, htXi⟩
    have htH : t ∈ poleSet hC.mH := by
      change hC.mH (criticalLinePoint t) = 0
      rw [hC.m_identity]
      exact htXi
    exact hne ⟨t, htH⟩

/-- RH closure from the Weyl/m-function endgame contract. -/
theorem mathlibRH_of_weyl_identity_contract
    (hC : RiemannWeylIdentityContract) :
    RiemannHypothesis := by
  rcases exists_ordinateEnumeration_of_semantic_bridge hC with ⟨ρ, hEnum⟩
  exact mathlibRH_of_ordinate_enumeration ρ hEnum

/-- Instantiation helper for the remaining bridge obligation:
    if `mXi` is identified with `XiTarget` and the target convergence-transfer
    contract is available, then `zero_to_poleXi` is derivable automatically. -/
theorem zero_to_poleXi_of_target_contract
    (mXi : ℂ → ℂ)
    (hmXi : mXi = XiTarget)
    (hTarget : RHConvergenceTransferContract XiTarget) :
    ∀ s : ℂ, riemannZeta s = 0 →
      (¬ ∃ n : ℕ, s = -2 * (n + 1)) →
      s ≠ 1 →
      ∃ t : ℝ, s = criticalLinePoint t ∧ t ∈ poleSet mXi := by
  intro s hs htriv h1
  have hsXi : XiTarget s = 0 := nontrivialZeroTransferToXiTarget s hs htriv h1
  have hsLine : onCriticalLine s := rhXiTarget_of_contract hTarget s hsXi
  have hsEq : s = criticalLinePoint s.im := by
    apply Complex.ext
    · simpa [onCriticalLine, criticalLinePoint] using hsLine
    · simp [criticalLinePoint]
  refine ⟨s.im, hsEq, ?_⟩
  · change mXi (criticalLinePoint s.im) = 0
    rw [hmXi, ← hsEq]
    exact hsXi

/-- Canonical constructor: when `mXi` is identified with `XiTarget`, the
    `zero_to_poleXi` bridge field is auto-derived from the target contract. -/
def mkWeylIdentityContract_of_target
    (N_H N_xi : ℝ → ℝ)
    (mH mXi : ℂ → ℂ)
    (counting_H : StepCountingSemantics N_H (poleSet mH))
    (counting_xi : StepCountingSemantics N_xi (poleSet mXi))
    (weyl_H : RiemannVonMangoldtEnvelope N_H)
    (weyl_xi : RiemannVonMangoldtEnvelope N_xi)
    (count_exact : ∀ T : ℝ, 0 ≤ T → N_H T = N_xi T)
    (mH_herglotz : HerglotzLike mH)
    (mXi_herglotz : HerglotzLike mXi)
    (m_identity : MFunctionIdentity mH mXi)
    (poleH_countable : (poleSet mH).Countable)
    (hmXi : mXi = XiTarget)
    (hTarget : RHConvergenceTransferContract XiTarget) :
    RiemannWeylIdentityContract where
  N_H := N_H
  N_xi := N_xi
  mH := mH
  mXi := mXi
  counting_H := counting_H
  counting_xi := counting_xi
  weyl_H := weyl_H
  weyl_xi := weyl_xi
  count_exact := count_exact
  mH_herglotz := mH_herglotz
  mXi_herglotz := mXi_herglotz
  m_identity := m_identity
  poleH_countable := poleH_countable
  zero_to_poleXi := zero_to_poleXi_of_target_contract mXi hmXi hTarget

/-- One-shot RH closure surface for the target-identified Weyl contract:
    all Weyl-side obligations plus target convergence-transfer contract. -/
theorem mathlibRH_of_weyl_identity_contract_of_target
    (N_H N_xi : ℝ → ℝ)
    (mH mXi : ℂ → ℂ)
    (counting_H : StepCountingSemantics N_H (poleSet mH))
    (counting_xi : StepCountingSemantics N_xi (poleSet mXi))
    (weyl_H : RiemannVonMangoldtEnvelope N_H)
    (weyl_xi : RiemannVonMangoldtEnvelope N_xi)
    (count_exact : ∀ T : ℝ, 0 ≤ T → N_H T = N_xi T)
    (mH_herglotz : HerglotzLike mH)
    (mXi_herglotz : HerglotzLike mXi)
    (m_identity : MFunctionIdentity mH mXi)
    (poleH_countable : (poleSet mH).Countable)
    (hmXi : mXi = XiTarget)
    (hTarget : RHConvergenceTransferContract XiTarget) :
    RiemannHypothesis := by
  exact mathlibRH_of_weyl_identity_contract
    (mkWeylIdentityContract_of_target
      N_H N_xi mH mXi
      counting_H counting_xi
      weyl_H weyl_xi
      count_exact
      mH_herglotz mXi_herglotz
      m_identity poleH_countable
      hmXi hTarget)

end

end Gutoe.RiemannWeylEndgame
