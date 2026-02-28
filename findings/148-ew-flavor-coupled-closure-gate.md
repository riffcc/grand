# 148 - EW + Flavor Coupled Closure Gate

## Goal
Track the coupled closure target in one lane:
1. `sin²(theta_W)` at `M_Z`
2. neutrino splitting ratio `R = Δm²32 / Δm²21`

Implemented as:
- `crates/gutoe-physics/src/bin/ew_flavor_coupled_ci_gate.rs`
- output: `/tmp/bh_renders/ew_flavor_coupled_ci_gate.json`

## Current snapshot
Run:
```bash
cargo run -q -p gutoe-physics --bin ew_flavor_coupled_ci_gate
```

Result:
- `overall_pass = false`
- `ew_ok = true`
- `ratio_ok = false`
- `hierarchy_ok = true`
- `ordering_ok = true`

### Electroweak target
- target `sin²(theta_W)(M_Z) = 0.23122`
- bridge value `0.231195465518`
- abs error `2.4534481981e-5`
- within tolerance (`5e-4`)

### Flavor-ratio target
- target ratio `R_target = 32.576361221780`
- observed ratio `R_model = 0.297103211006`
- relative error `-0.990879791362` (about `-99.09%`)
- outside tolerance (`5%`)

## Interpretation
This confirms the current structural state:
- EW uses the clean structural bridge `3/13 + α²·d/2` and lands at 0.01% error.
- flavor mass splitting remains compressed.

The closure problem is therefore concentrated in flavor-to-mass separation,
not in the overall neutrino mass scale cap or hierarchy sign.

## Next closure criterion
This lane closes when both hold simultaneously:
- `|sin²(theta_W)(M_Z) - 0.23122| <= 5e-4`
- `|R_model / R_target - 1| <= 5%`

No separate success claims until both pass in the same artifact.
