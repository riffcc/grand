import Mathlib

/-!
GUTOE — Geometry Shortcut Lane for FTL Questions

This module formalizes the distinction between:
1) local signal-speed bounds (`v_local <= c`), and
2) coordinate-effective travel speed through geometry/topology shortcuts.

It does **not** claim a physical engine exists; it states the logical structure:
shortcut geometry can yield coordinate-effective speed `> c` without requiring
local propagation `> c`.
-/

namespace Gutoe.FTLGeometry

/-- Local signal-speed bound (causal lane): local propagation does not exceed `c`. -/
def LocalSignalBound (c vLocal : ℝ) : Prop := vLocal ≤ c

/-- Coordinate-effective speed over coordinate distance `d` and elapsed coordinate time `dt`. -/
noncomputable def coordinateSpeed (d dt : ℝ) : ℝ := d / dt

/-- Shortcut factor `s` rescales coordinate distance to a shorter effective path:
effective path length = `s * d`, with `0 < s < 1` for a true shortcut. -/
noncomputable def shortcutTravelTime (d c s : ℝ) : ℝ := (s * d) / c

/-- With positive distance, speed scale, and shortcut factor, shortcut travel time is positive. -/
theorem shortcut_travel_time_pos {d c s : ℝ}
    (hd : d > 0) (hc : c > 0) (hs : s > 0) :
    shortcutTravelTime d c s > 0 := by
  unfold shortcutTravelTime
  positivity

/-- Coordinate-effective speed induced by shortcut time is exactly `c / s`. -/
theorem coordinate_speed_shortcut_closed_form {d c s : ℝ}
    (hd : d ≠ 0) (hc : c ≠ 0) (hs : s ≠ 0) :
    coordinateSpeed d (shortcutTravelTime d c s) = c / s := by
  unfold coordinateSpeed shortcutTravelTime
  field_simp [hd, hc, hs]

/-- If `0 < s < 1`, coordinate-effective speed is strictly greater than `c`. -/
theorem coordinate_speed_gt_c_of_shortcut {d c s : ℝ}
    (hd : d > 0) (hc : c > 0) (hs0 : s > 0) (hs1 : s < 1) :
    coordinateSpeed d (shortcutTravelTime d c s) > c := by
  have hd0 : d ≠ 0 := ne_of_gt hd
  have hc0 : c ≠ 0 := ne_of_gt hc
  have hs0' : s ≠ 0 := ne_of_gt hs0
  rw [coordinate_speed_shortcut_closed_form hd0 hc0 hs0']
  have hdiv : c / s > c := by
    have hcs : c * s < c := by
      nlinarith
    exact (lt_div_iff₀ hs0).2 hcs
  exact hdiv

/-- Core geometry-loophole theorem:
there exist regimes with local bound `vLocal <= c` and coordinate-effective
speed `> c` simultaneously (for `0 < s < 1`). -/
theorem local_bound_compatible_with_coordinate_superluminal {d c s : ℝ}
    (hd : d > 0) (hc : c > 0) (hs0 : s > 0) (hs1 : s < 1) :
    LocalSignalBound c c ∧ coordinateSpeed d (shortcutTravelTime d c s) > c := by
  refine ⟨?_, ?_⟩
  · unfold LocalSignalBound
    exact le_rfl
  · exact coordinate_speed_gt_c_of_shortcut hd hc hs0 hs1

end Gutoe.FTLGeometry
