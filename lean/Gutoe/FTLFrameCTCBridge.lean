import Mathlib
import Gutoe.FTLRearFaceBridge
import Gutoe.CTCLegality

/-!
GUTOE — Frame/CTC bridge for coordinate-effective superluminal segments

This module formalizes a standard tachyon-antitelephone-style ingredient:

- If a segment has coordinate speed `u > c`, then there exists a subluminal
  boost parameter `0 < v < c` such that the boosted time numerator
  `Δt - vΔx/c²` is negative.

It then specializes this to the rear-face shortcut lane (`s = 1/10`) from
`FTLRearFaceBridge`, and pairs it with the existing CTC legality witness.
-/

namespace Gutoe.FTLFrameCTCBridge

open Gutoe.FTLGeometry
open Gutoe.FTLRearFaceBridge
open Gutoe.CTCLegality

/-- Boosted-frame time numerator (up to positive Lorentz factor). -/
noncomputable def boostedTimeNumerator (dt dx c v : ℝ) : ℝ :=
  dt - (v * dx) / (c ^ 2)

/-- For a coordinate-effective superluminal segment `dx = u*dt` with `u > c > 0`,
there exists a subluminal boost `0 < v < c` with negative boosted time numerator. -/
theorem exists_subluminal_boost_with_negative_time_numerator
    {c u dt : ℝ}
    (hc : 0 < c)
    (hu : u > c)
    (hdt : 0 < dt) :
    ∃ v : ℝ,
      0 < v ∧ v < c ∧
      boostedTimeNumerator dt (u * dt) c v < 0 := by
  have hu_pos : 0 < u := lt_trans hc hu
  have hc2_pos : 0 < c ^ 2 := by positivity
  have hcut_lt_c : c ^ 2 / u < c := by
    have hmul : c ^ 2 < c * u := by nlinarith
    exact (div_lt_iff₀ hu_pos).2 hmul
  let v : ℝ := (c + c ^ 2 / u) / 2
  have hv_pos : 0 < v := by
    have hcut_nonneg : 0 ≤ c ^ 2 / u := by
      exact div_nonneg hc2_pos.le hu_pos.le
    unfold v
    nlinarith [hc, hcut_nonneg]
  have hv_lt_c : v < c := by
    unfold v
    nlinarith [hcut_lt_c]
  have hv_gt_cut : c ^ 2 / u < v := by
    unfold v
    nlinarith [hcut_lt_c]
  refine ⟨v, hv_pos, hv_lt_c, ?_⟩
  unfold boostedTimeNumerator
  have huv : c ^ 2 < v * u := by
    have hmul : (c ^ 2 / u) * u < v * u := by
      nlinarith [hv_gt_cut, hu_pos]
    have hcut_mul : (c ^ 2 / u) * u = c ^ 2 := by
      field_simp [hu_pos.ne']
    nlinarith [hmul, hcut_mul]
  have hratio : 1 < (v * u) / (c ^ 2) := by
    have hscaled : 1 * (c ^ 2) < v * u := by nlinarith [huv]
    exact (lt_div_iff₀ hc2_pos).2 hscaled
  have hfac : 1 - (v * u) / (c ^ 2) < 0 := by nlinarith
  have hmul_neg : dt * (1 - (v * u) / (c ^ 2)) < 0 := by
    exact mul_neg_of_pos_of_neg hdt hfac
  have hrew : dt - (v * (u * dt)) / (c ^ 2) = dt * (1 - (v * u) / (c ^ 2)) := by
    ring
  rw [hrew]
  exact hmul_neg

/-- Rear-face shortcut specialization:
`s = 1/10` implies a frame with negative boosted-time numerator. -/
theorem rear_shortcut_has_predeparture_frame
    {d c : ℝ}
    (hd : d > 0)
    (hc : c > 0) :
    ∃ v : ℝ,
      0 < v ∧ v < c ∧
      boostedTimeNumerator
        (shortcutTravelTime d c rearShortcutFactor)
        d c v < 0 := by
  let dt : ℝ := shortcutTravelTime d c rearShortcutFactor
  let u : ℝ := coordinateSpeed d dt
  have hdt_pos : 0 < dt := by
    simpa [dt] using shortcut_travel_time_pos hd hc rear_shortcut_factor_pos
  have hu_gt_c : u > c := by
    simpa [u, dt] using
      coordinate_speed_gt_c_of_shortcut hd hc
        rear_shortcut_factor_pos rear_shortcut_factor_lt_one
  rcases exists_subluminal_boost_with_negative_time_numerator
      (c := c) (u := u) (dt := dt) hc hu_gt_c hdt_pos with
    ⟨v, hv0, hvc, hvneg⟩
  refine ⟨v, hv0, hvc, ?_⟩
  have hdt_ne : dt ≠ 0 := ne_of_gt hdt_pos
  have hdx : u * dt = d := by
    unfold u coordinateSpeed
    field_simp [hdt_ne]
  have hvneg_d : boostedTimeNumerator dt d c v < 0 := by
    simpa [hdx] using hvneg
  simpa [dt] using hvneg_d

/-- Combined statement:
the rear-shortcut lane admits a predeparture frame witness, and the current
periodic timelike-identification lane admits a CTC witness. -/
theorem rear_shortcut_predeparture_and_ctc_witness
    {d c T : ℝ}
    (hd : d > 0)
    (hc : c > 0)
    (hT : T > 0) :
    (∃ v : ℝ,
      0 < v ∧ v < c ∧
      boostedTimeNumerator
        (shortcutTravelTime d c rearShortcutFactor)
        d c v < 0) ∧
    (∃ a b : Event, Timelike a b ∧ sameOnTimeCylinder T a b) := by
  refine ⟨rear_shortcut_has_predeparture_frame hd hc, ?_⟩
  exact ctc_exists_on_time_cylinder T hT

/-- End-to-end possibility statement for the current bridged lane:
local causal bound is preserved, coordinate-effective superluminal transport is
present, a predeparture boosted-frame witness exists, and the CTC legality lane
provides a timelike identified-loop witness. -/
theorem rear_lane_time_travel_possible
    {d c T : ℝ}
    (hd : d > 0)
    (hc : c > 0)
    (hT : T > 0) :
    (LocalSignalBound c c ∧
      coordinateSpeed d (shortcutTravelTime d c rearShortcutFactor) > c) ∧
    (∃ v : ℝ,
      0 < v ∧ v < c ∧
      boostedTimeNumerator
        (shortcutTravelTime d c rearShortcutFactor)
        d c v < 0) ∧
    (∃ a b : Event, Timelike a b ∧ sameOnTimeCylinder T a b) := by
  have hLocal :
      LocalSignalBound c c ∧
        coordinateSpeed d (shortcutTravelTime d c rearShortcutFactor) > c :=
    local_bound_compatible_with_rear_coordinate_superluminal hd hc
  have hBridge :
      (∃ v : ℝ,
        0 < v ∧ v < c ∧
        boostedTimeNumerator
          (shortcutTravelTime d c rearShortcutFactor)
          d c v < 0) ∧
      (∃ a b : Event, Timelike a b ∧ sameOnTimeCylinder T a b) :=
    rear_shortcut_predeparture_and_ctc_witness hd hc hT
  exact ⟨hLocal, hBridge.1, hBridge.2⟩

end Gutoe.FTLFrameCTCBridge
