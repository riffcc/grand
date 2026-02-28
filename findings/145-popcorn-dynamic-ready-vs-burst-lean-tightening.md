# 145 - Popcorn Dynamic Lane: Ready vs Burst + Lean Constraint Tightening

Date: 2026-02-28  
Scope: Replace static popcorn threshold with dynamic thermal/kinetic/fracture model and tighten Lean gates accordingly.

## Why This Change

The earlier popcorn lane used a single static criterion (`P_int >= P_thresh`) and returned a pop temperature near `155°C`.
That was physically plausible as a **pressure-ready** point, but it conflated:

- readiness to fail (threshold crossed),
- and observed burst (delayed rupture after viscoelastic/damage effects).

This update introduces explicit `T_ready` and `T_burst`.

## Model Upgrades

In `crates/gutoe-physics/src/everyday_extremes.rs`:

1. Thermal lag:
- Core kernel temperature follows external heating with first-order lag (`tau`).

2. Vapor-rate kinetics:
- Vapor fraction evolves with Arrhenius dynamics (`A * exp(-Ea/RT)`), not instant equilibrium.

3. Hull softening + damage integral:
- Temperature-dependent strength reduction.
- Overstress-driven damage accumulation.
- Burst occurs when cumulative damage crosses threshold.

This yields:
- `T_ready`: first time `P_int >= P_thresh(T)`,
- `T_burst`: first time damage criterion is met under overstress.

## Returned Popcorn Numbers (Current Default Run)

From `/tmp/bh_renders/everyday_extremes/everyday_extremes_report.txt`:

- `ready_temperature_c = 152.519`
- `burst_temperature_c = 184.938`
- `hysteresis_delta_c = 32.419`
- `ready_time_s = 119.600`
- `burst_time_s = 144.600`
- `internal_pressure_ready_mpa = 1.020484`
- `internal_pressure_burst_mpa = 2.300463`
- `rupture_threshold_ready_mpa = 1.013830`
- `rupture_threshold_burst_mpa = 0.951586`
- `estimated_expansion_ratio = 46.521`

Interpretation:
- The old `~155°C` result now maps cleanly to **pressure-ready**.
- Observed burst aligns with the `~180°C` regime via delayed fracture dynamics.

## CI Gate Tightening

Updated `crates/gutoe-physics/src/bin/everyday_extremes_ci_gate.rs`:

- old single window:
  - `150 ≤ T_pop ≤ 220`
- new mechanistic windows:
  - `145 ≤ T_ready ≤ 175`
  - `170 ≤ T_burst ≤ 195`
  - `T_burst > T_ready`
  - `hysteresis_delta_c ≥ 5`
  - expansion still constrained (`≥ 10`)

Current gate output:
- `overall_pass = true`

## Lean Tightening

Updated `lean/Gutoe/EverydayExtremes.lean` with new theorem-level constraints:

- `popcorn_ready_temperature_gate`
- `popcorn_burst_temperature_gate`
- `popcorn_ready_before_burst`
- `popcorn_hysteresis_gate`
- retained `popcorn_expansion_gate`

Spine theorem `everyday_extremes_constraint_spine` now encodes two-temperature popcorn behavior, not a loose one-temperature band.

## Verification

- `cargo run -q -p gutoe-physics --bin everyday_extremes_report` passed
- `cargo run -q -p gutoe-physics --bin everyday_extremes_ci_gate` passed
- `cargo test -q -p gutoe-physics --lib everyday_extremes` passed (4/4)
- `cd lean && lake build Gutoe.EverydayExtremes` passed
- `cd lean && lake build Gutoe` passed (warnings only)
