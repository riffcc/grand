import Mathlib
import Gutoe.RiemannCore
import Gutoe.GaugeGroupSU3

namespace Gutoe.RiemannCounting

open Gutoe.GaugeGroupSU3

noncomputable section

/-- Finite counting proxy: number of spectral parameters with `|t| ≤ T`. -/
def countUpTo (spec : Finset ℝ) (T : ℝ) : ℕ :=
  (spec.filter (fun t => |t| ≤ T)).card

/-- Finite counting-match predicate used as the theorem-level analog
    of empirical `N(T)` agreement checks. -/
def FiniteCountingMatch (specA specB : Finset ℝ) : Prop :=
  ∀ T : ℝ, countUpTo specA T = countUpTo specB T

theorem countUpTo_mono (spec : Finset ℝ) {T₁ T₂ : ℝ} (hT : T₁ ≤ T₂) :
    countUpTo spec T₁ ≤ countUpTo spec T₂ := by
  classical
  unfold countUpTo
  refine Finset.card_le_card ?_
  intro x hx
  simp at hx ⊢
  exact ⟨hx.1, le_trans hx.2 hT⟩

theorem finiteCountingMatch_refl (spec : Finset ℝ) :
    FiniteCountingMatch spec spec := by
  intro T
  rfl

theorem finiteCountingMatch_symm {specA specB : Finset ℝ}
    (h : FiniteCountingMatch specA specB) :
    FiniteCountingMatch specB specA := by
  intro T
  symm
  exact h T

theorem finiteCountingMatch_trans {specA specB specC : Finset ℝ}
    (hAB : FiniteCountingMatch specA specB)
    (hBC : FiniteCountingMatch specB specC) :
    FiniteCountingMatch specA specC := by
  intro T
  calc
    countUpTo specA T = countUpTo specB T := hAB T
    _ = countUpTo specC T := hBC T

theorem finiteCountingMatch_of_eq {specA specB : Finset ℝ}
    (hEq : specA = specB) :
    FiniteCountingMatch specA specB := by
  intro T
  simpa [hEq]

/-- Structural sanity anchor for RH lanes: quark orbit count remains 3. -/
theorem z3_structural_count_anchor : quarkOrbit.card = 3 := by
  simpa using quarkOrbit_card

end

end Gutoe.RiemannCounting
