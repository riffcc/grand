# 211 — FTLGeometry × VoidRearFace bridge (explicit `s = 1/10` lane)

## New theorem module
- `lean/Gutoe/FTLRearFaceBridge.lean`

## Wiring
- Added `Gutoe.FTLRearFaceBridge` to `lean/lakefile.lean` roots.

## What is now formally connected

Previously separate:
1. `Gutoe.FTLGeometry`:
   - local bound and coordinate-effective superluminal compatibility for `0 < s < 1`.
2. `Gutoe.VoidRearFace` / `Gutoe.VacuumEnergyBounds`:
   - rear-channel suppression/cost factor exactly `1/10`.

Now bridged:
- `rearShortcutFactor := rearFaceSuppressionR`
- `rearShortcutFactor = 1/10`

## New theorems
1. `rear_shortcut_factor_eq_one_tenth`
2. `rear_shortcut_factor_pos`
3. `rear_shortcut_factor_lt_one`
4. `coordinate_speed_rear_shortcut_eq_ten_c`
   - `coordinateSpeed d (shortcutTravelTime d c rearShortcutFactor) = 10 * c`
     (under `d ≠ 0`, `c ≠ 0`)
5. `local_bound_compatible_with_rear_coordinate_superluminal`
   - local causal bound (`v_local ≤ c`) and coordinate-effective `> c`
     hold simultaneously in the rear-shortcut-specialized lane.

## Interpretation boundary
- This proves logical/mathematical compatibility and the explicit `10c` coordinate
  scaling under the bridged shortcut factor.
- It does **not** by itself prove a physical engine implementation.

## Build status
- `lake build Gutoe.FTLRearFaceBridge` passed.
- `lake build Gutoe` passed (warnings only).
