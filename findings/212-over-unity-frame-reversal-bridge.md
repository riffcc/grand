# 212 — Over-unity coordinate segment implies predeparture frame witness

## New theorem module
- `lean/Gutoe/FTLFrameCTCBridge.lean`

## Wiring
- Added `Gutoe.FTLFrameCTCBridge` to `lean/lakefile.lean` roots.

## New formal results
1. `exists_subluminal_boost_with_negative_time_numerator`
   - For `u > c > 0` and `dt > 0`, there exists `v` with `0 < v < c`
     such that `dt - v*dx/c^2 < 0` (with `dx = u*dt`).
   - This is the frame-reversal witness (tachyon-antitelephone ingredient).

2. `rear_shortcut_has_predeparture_frame`
   - Specializes the above to the rear-shortcut lane:
     `dt = shortcutTravelTime d c rearShortcutFactor` and `dx = d`.
   - Uses previously proven `rearShortcutFactor = 1/10` and
     `coordinateSpeed > c` from `FTLRearFaceBridge`.

3. `rear_shortcut_predeparture_and_ctc_witness`
   - Produces conjunction:
     - predeparture-frame witness for rear shortcut segment, and
     - CTC witness from periodic timelike identification (`CTCLegality`).

## Interpretation
- This module closes the “if over-unity in one frame...” step at theorem level:
  coordinate-effective superluminal segment implies existence of a subluminal
  boost frame with reversed coordinate time ordering witness.
- It remains a logical/kinematic bridge, not a standalone engineering proof.

## Build status
- `lake build Gutoe.FTLFrameCTCBridge` passed.
- `lake build Gutoe` passed (warnings only).
