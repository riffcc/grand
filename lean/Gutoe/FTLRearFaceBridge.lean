import Mathlib
import Gutoe.FTLGeometry
import Gutoe.VacuumEnergyBounds

/-!
GUTOE — FTL Rear-Face Bridge

Bridges two previously separate lanes:

1. `FTLGeometry`: coordinate-effective superluminal travel is logically compatible
   with local signal bound `≤ c` when shortcut factor `0 < s < 1`.
2. `VoidRearFace`/`VacuumEnergyBounds`: rear-channel suppression is exactly `1/10`.

This module identifies the shortcut factor with the rear-face suppression factor
and derives the resulting coordinate-effective `10c` scaling in that bridged lane.
-/

namespace Gutoe.FTLRearFaceBridge

open Gutoe.FTLGeometry
open Gutoe.VacuumEnergyBounds

/-- Bridged shortcut factor: use rear-face suppression as the geometry shortcut factor. -/
noncomputable def rearShortcutFactor : ℝ := rearFaceSuppressionR

theorem rear_shortcut_factor_eq_one_tenth :
    rearShortcutFactor = (1 : ℝ) / 10 := by
  simpa [rearShortcutFactor] using rear_face_suppression_eq_one_tenth

theorem rear_shortcut_factor_pos : 0 < rearShortcutFactor := by
  rw [rear_shortcut_factor_eq_one_tenth]
  norm_num

theorem rear_shortcut_factor_lt_one : rearShortcutFactor < 1 := by
  rw [rear_shortcut_factor_eq_one_tenth]
  norm_num

/-- Closed-form coordinate-effective speed in the bridged rear-face shortcut lane. -/
theorem coordinate_speed_rear_shortcut_eq_ten_c
    {d c : ℝ} (hd : d ≠ 0) (hc : c ≠ 0) :
    coordinateSpeed d (shortcutTravelTime d c rearShortcutFactor) = 10 * c := by
  have hs : rearShortcutFactor ≠ 0 := ne_of_gt rear_shortcut_factor_pos
  rw [coordinate_speed_shortcut_closed_form hd hc hs]
  rw [rear_shortcut_factor_eq_one_tenth]
  field_simp

/-- Local causal bound with rear-face shortcut still yields coordinate-effective
superluminal speed; this specializes the generic geometry-loophole theorem. -/
theorem local_bound_compatible_with_rear_coordinate_superluminal
    {d c : ℝ} (hd : d > 0) (hc : c > 0) :
    LocalSignalBound c c ∧
      coordinateSpeed d (shortcutTravelTime d c rearShortcutFactor) > c := by
  exact local_bound_compatible_with_coordinate_superluminal
    hd hc rear_shortcut_factor_pos rear_shortcut_factor_lt_one

end Gutoe.FTLRearFaceBridge
