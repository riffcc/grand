import Mathlib
import Gutoe.ProjectionFibers
import Gutoe.ContainmentScope

/-!
GUTOE — Fiber Channel Capacity / Compression Lane

This module separates three statements:
1. Structural hidden dimension per event is `12` (kernel rank).
2. Passive decode from projection alone is impossible (non-injective).
3. Active decode with a supplied kernel key is exact.
-/

namespace Gutoe.FiberChannelCapacity

open Gutoe
open Gutoe.ProjectionFibers
open Gutoe.ContainmentScope

noncomputable section

/-- Visible projection dimension (`Fin 4 -> ℝ`). -/
def visibleDim : ℕ := Module.finrank ℝ (Fin 4 → ℝ)

/-- Hidden fiber dimension (`ker grade1Projection`). -/
def hiddenDim : ℕ := Module.finrank ℝ (LinearMap.ker grade1Projection)

theorem visible_dim_eq_4 : visibleDim = 4 := by
  unfold visibleDim
  simpa using (finrank_euclideanSpace_fin : Module.finrank ℝ (Fin 4 → ℝ) = 4)

theorem hidden_dim_eq_12 : hiddenDim = 12 := by
  unfold hiddenDim
  exact grade1Projection_kernel_finrank

/-- Hidden/visible dimensional ratio. -/
def hiddenVisibleRatioQ : ℚ := (hiddenDim : ℚ) / (visibleDim : ℚ)

theorem hidden_visible_ratio_eq_3 : hiddenVisibleRatioQ = 3 := by
  unfold hiddenVisibleRatioQ
  rw [hidden_dim_eq_12, visible_dim_eq_4]
  norm_num

/-- Quantized visible payload bits per event if each visible axis has `b` bits. -/
def visibleBitsPerEvent (b : ℕ) : ℕ := visibleDim * b

/-- Quantized hidden payload bits per event if each hidden axis has `b` bits. -/
def hiddenBitsPerEvent (b : ℕ) : ℕ := hiddenDim * b

/-- Quantized total payload bits per event. -/
def totalBitsPerEvent (b : ℕ) : ℕ := visibleBitsPerEvent b + hiddenBitsPerEvent b

theorem visible_bits_per_event_closed_form (b : ℕ) :
    visibleBitsPerEvent b = 4 * b := by
  unfold visibleBitsPerEvent
  rw [visible_dim_eq_4]

theorem hidden_bits_per_event_closed_form (b : ℕ) :
    hiddenBitsPerEvent b = 12 * b := by
  unfold hiddenBitsPerEvent
  rw [hidden_dim_eq_12]

theorem total_bits_per_event_closed_form (b : ℕ) :
    totalBitsPerEvent b = 16 * b := by
  unfold totalBitsPerEvent
  rw [visible_bits_per_event_closed_form, hidden_bits_per_event_closed_form]
  omega

theorem hidden_bits_are_three_times_visible (b : ℕ) :
    hiddenBitsPerEvent b = 3 * visibleBitsPerEvent b := by
  rw [hidden_bits_per_event_closed_form, visible_bits_per_event_closed_form]
  ring

/-- Passive decode impossibility from projection alone. -/
theorem passive_decode_impossible :
    ¬ ∃ rec : (Fin 4 → ℝ) → Vec16, Function.LeftInverse rec grade1Projection := by
  exact no_global_state_reconstructor

/-- Active decode with explicit kernel key. -/
def decodeWithKey (x : Fin 4 → ℝ) (k : LinearMap.ker grade1Projection) : Vec16 :=
  fiberBase x + k.1

/-- Decoding with key always projects back to the requested visible state. -/
theorem decode_with_key_projects_back
    (x : Fin 4 → ℝ) (k : LinearMap.ker grade1Projection) :
    grade1Projection (decodeWithKey x k) = x := by
  unfold decodeWithKey
  have hk : grade1Projection k.1 = 0 := k.2
  calc
    grade1Projection (fiberBase x + k.1)
        = grade1Projection (fiberBase x) + grade1Projection k.1 := by simp
    _ = x + 0 := by rw [grade1Projection_fiberBase, hk]
    _ = x := by simp

end
end Gutoe.FiberChannelCapacity
