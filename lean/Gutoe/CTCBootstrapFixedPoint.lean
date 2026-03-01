import Mathlib
import Gutoe.DynamicTopologyCreation

/-!
GUTOE — CTC bootstrap fixed-point lane

We model the "send energy back to create the same door" idea as a recurrence:

`E_past = η * E_future - loss`

A bootstrap self-consistency point is a fixed point `E*` with
`E* = η*E* - loss`.

This module proves the sharp closed-cycle condition:
- If `η ≤ 1` and `loss ≥ 0`, any strictly positive fixed point forces
  `η = 1` and `loss = 0`.
- So a positive bootstrap point under non-amplifying closed-cycle dynamics is
  possible only in the lossless ideal limit.
-/

namespace Gutoe.CTCBootstrapFixedPoint

open Gutoe.DynamicTopologyCreation

/-- Closed-cycle return map for bootstrap energy circulation. -/
def bootstrapReturn (η loss Efuture : ℝ) : ℝ := η * Efuture - loss

/-- Fixed-point condition for bootstrap circulation. -/
def bootstrapFixed (η loss E : ℝ) : Prop := bootstrapReturn η loss E = E

/-- If `η < 1`, a positive fixed point forces strictly negative loss. -/
theorem eta_lt_one_positive_fixed_implies_negative_loss
    (η loss E : ℝ)
    (hEta : η < 1)
    (hE : 0 < E)
    (hFix : bootstrapFixed η loss E) :
    loss < 0 := by
  have hLossLinear : loss = (η - 1) * E := by
    unfold bootstrapFixed bootstrapReturn at hFix
    nlinarith
  have hMulNeg : (η - 1) * E < 0 := by
    exact mul_neg_of_neg_of_pos (sub_lt_zero.mpr hEta) hE
  simpa [hLossLinear] using hMulNeg

/-- Under `η ≤ 1` and `loss ≥ 0`, any positive fixed point is possible only if
`η = 1` and `loss = 0`. -/
theorem positive_fixed_under_eta_le_one_loss_nonneg_is_lossless
    (η loss E : ℝ)
    (hEtaLe : η ≤ 1)
    (hLoss : 0 ≤ loss)
    (hE : 0 < E)
    (hFix : bootstrapFixed η loss E) :
    η = 1 ∧ loss = 0 := by
  by_cases hEtaLt : η < 1
  · have hLossNeg : loss < 0 :=
      eta_lt_one_positive_fixed_implies_negative_loss η loss E hEtaLt hE hFix
    have hFalse : False := (not_lt_of_ge hLoss) hLossNeg
    exact False.elim hFalse
  · have hOneLe : 1 ≤ η := le_of_not_gt hEtaLt
    have hEtaEq : η = 1 := le_antisymm hEtaLe hOneLe
    have hLossEq : loss = 0 := by
      unfold bootstrapFixed bootstrapReturn at hFix
      rw [hEtaEq] at hFix
      linarith
    exact ⟨hEtaEq, hLossEq⟩

/-- Threshold-targeted bootstrap: if a structural threshold energy is a positive
fixed point under `η ≤ 1` and `loss ≥ 0`, closed-cycle dynamics must be exactly
lossless (`η=1`, `loss=0`). -/
theorem threshold_bootstrap_requires_lossless_closed_cycle
    (η loss radius period : ℝ)
    (hEtaLe : η ≤ 1)
    (hLoss : 0 ≤ loss)
    (hThrPos : 0 < structuralCreationThreshold radius period)
    (hFix : bootstrapFixed η loss (structuralCreationThreshold radius period)) :
    η = 1 ∧ loss = 0 := by
  exact positive_fixed_under_eta_le_one_loss_nonneg_is_lossless
    η loss (structuralCreationThreshold radius period) hEtaLe hLoss hThrPos hFix

end Gutoe.CTCBootstrapFixedPoint
