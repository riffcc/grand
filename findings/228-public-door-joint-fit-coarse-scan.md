# Finding 228 — Public Door Joint Fit (Coarse Scan)

## Summary

Ran a direct "look for doors" synthesis over current public lanes with available
quantitative residuals.

Binary:

- `crates/gutoe-physics/src/bin/ctc_public_door_joint_fit.rs`

Outputs:

- `/tmp/bh_renders/ctc_public_door_joint_fit/ctc_public_door_joint_fit.txt`
- `/tmp/bh_renders/ctc_public_door_joint_fit/ctc_public_door_joint_fit.json`

## Run

```bash
cargo run -p gutoe-physics --bin ctc_public_door_joint_fit
```

## Result

- `n_numeric_lanes_included = 3`
- `z_stouffer = 0.796188`
- `max_abs_z = 1.135135`
- `chi2_null = 1.348020`
- `door_detected = false`

## Lanes

Included numeric lanes:

1. `cosmic_birefringence_anisotropic_A_CB`
2. `eht_shadow_fractional_deviation_from_kerr`
3. `electroweak_running_lowQ_weak_charge_residual`

Tracked but not included in joint scalar yet:

4. `desi_void_topology_statistics` (`pending_quantitative_scalar`)

## Notes

- This is intentionally a coarse pass and not a publication-grade global fit.
- It gives a direct immediate answer to "look for doors now" with current
  accessible lane-level residuals.
- Current coarse public-lane synthesis shows **no robust spontaneous-door
  detection**.
