import Mathlib
import Gutoe.FineStructure

namespace Gutoe.EverydayExtremes

open Gutoe.FineStructure

/-- Leading-order electromagnetic coupling used by the everyday lanes. -/
def alphaLeadingOrderQ : ℚ := (1 : ℚ) / (alphaInverse 4 : ℚ)

theorem alpha_leading_order_q :
    alphaLeadingOrderQ = 1 / 137 := by
  unfold alphaLeadingOrderQ
  rw [alpha_inverse_d4]
  norm_num

/-- Warm-ice friction coefficient snapshot from the Rust lane at `T=-2°C`. -/
def iceMuWarmQ : ℚ := 11506984942670542 / 100000000000000000

/-- Cold-ice friction coefficient snapshot from the Rust lane at `T=-20°C`. -/
def iceMuColdQ : ℚ := 2564129732919839 / 10000000000000000

/-- Friction drop from `-20°C` to `-2°C` at fixed pressure/speed. -/
def iceMuDropQ : ℚ := iceMuColdQ - iceMuWarmQ

theorem ice_slipperiness_drop_positive :
    iceMuDropQ > 1 / 10 := by
  native_decide

/-- Pressure-ready temperature snapshot (°C) from the dynamic popcorn lane. -/
def popcornReadyTempCQ : ℚ := 152519 / 1000

/-- Burst temperature snapshot (°C) from the dynamic popcorn lane. -/
def popcornBurstTempCQ : ℚ := 184938 / 1000

/-- Thermal hysteresis (°C) between pressure-ready and actual burst. -/
def popcornHysteresisCQ : ℚ := popcornBurstTempCQ - popcornReadyTempCQ

/-- Expansion ratio snapshot from the popcorn lane. -/
def popcornExpansionRatioQ : ℚ := 46521 / 1000

theorem popcorn_ready_temperature_gate :
    145 ≤ popcornReadyTempCQ ∧ popcornReadyTempCQ ≤ 175 := by
  native_decide

theorem popcorn_burst_temperature_gate :
    170 ≤ popcornBurstTempCQ ∧ popcornBurstTempCQ ≤ 195 := by
  native_decide

theorem popcorn_ready_before_burst :
    popcornReadyTempCQ < popcornBurstTempCQ := by
  native_decide

theorem popcorn_hysteresis_gate :
    popcornHysteresisCQ ≥ 5 := by
  native_decide

theorem popcorn_expansion_gate :
    popcornExpansionRatioQ ≥ 10 := by
  native_decide

/-- Optimal raindrop diameter snapshot (mm) from the shape sweep. -/
def raindropOptimalDiameterMmQ : ℚ := 59 / 10

/-- Optimal raindrop aspect ratio snapshot (minor/major axis). -/
def raindropOptimalAspectQ : ℚ := 6731128988915706 / 10000000000000000

theorem raindrop_diameter_window :
    5 / 2 ≤ raindropOptimalDiameterMmQ ∧ raindropOptimalDiameterMmQ ≤ 6 := by
  native_decide

theorem raindrop_aspect_window :
    3 / 5 ≤ raindropOptimalAspectQ ∧ raindropOptimalAspectQ ≤ 4 / 5 := by
  native_decide

/-- Default Mpemba-case freeze time snapshot for the initially hot sample (minutes). -/
def mpembaHotFreezeTimeMinQ : ℚ := 199047 / 1000

/-- Default Mpemba-case freeze time snapshot for the initially cold sample (minutes). -/
def mpembaColdFreezeTimeMinQ : ℚ := 274865 / 1000

/-- Fraction of hot-faster outcomes in the fixed small regime sweep. -/
def mpembaSweepHotFasterFractionQ : ℚ := 1

theorem mpemba_default_ordering :
    mpembaHotFreezeTimeMinQ < mpembaColdFreezeTimeMinQ := by
  native_decide

theorem mpemba_sweep_fraction_gate :
    mpembaSweepHotFasterFractionQ ≥ 1 / 5 := by
  native_decide

/-- GRAND constraint spine for the everyday-extremes lane.
    This packages the operational gate inequalities into Lean. -/
theorem everyday_extremes_constraint_spine :
    alphaLeadingOrderQ = 1 / 137 ∧
    iceMuDropQ > 1 / 10 ∧
    (145 ≤ popcornReadyTempCQ ∧ popcornReadyTempCQ ≤ 175) ∧
    (170 ≤ popcornBurstTempCQ ∧ popcornBurstTempCQ ≤ 195) ∧
    popcornReadyTempCQ < popcornBurstTempCQ ∧
    popcornHysteresisCQ ≥ 5 ∧
    popcornExpansionRatioQ ≥ 10 ∧
    (5 / 2 ≤ raindropOptimalDiameterMmQ ∧ raindropOptimalDiameterMmQ ≤ 6) ∧
    (3 / 5 ≤ raindropOptimalAspectQ ∧ raindropOptimalAspectQ ≤ 4 / 5) ∧
    mpembaHotFreezeTimeMinQ < mpembaColdFreezeTimeMinQ ∧
    mpembaSweepHotFasterFractionQ ≥ 1 / 5 := by
  exact ⟨alpha_leading_order_q,
    ice_slipperiness_drop_positive,
    popcorn_ready_temperature_gate,
    popcorn_burst_temperature_gate,
    popcorn_ready_before_burst,
    popcorn_hysteresis_gate,
    popcorn_expansion_gate,
    raindrop_diameter_window,
    raindrop_aspect_window,
    mpemba_default_ordering,
    mpemba_sweep_fraction_gate⟩

end Gutoe.EverydayExtremes
