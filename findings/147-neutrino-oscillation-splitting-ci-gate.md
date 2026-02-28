# 147 - Neutrino Oscillation Splitting CI Gate (Closure Target)

## Why this gate exists
The absolute-mass transduction lane must be consistent with oscillation
measurements, not only hierarchy/type and upper caps.

This gate explicitly checks the two measured mass-squared splittings:
- Solar: `Δm²21 ≈ 7.53e-5 eV²`
- Atmospheric (normal ordering): `Δm²32 ≈ 2.453e-3 eV²`

## Implementation
New binary:
- `crates/gutoe-em/src/bin/neutrino_oscillation_ci_gate.rs`

Output artifact:
- `/tmp/bh_renders/neutrino_oscillation_ci_gate.json`

Current gate window:
- Relative tolerance: `±5%`

## Current result (expected fail)
From current absolute masses:
- `m1 = 8.497119214462e-4 eV`
- `m2 = 6.952371162470e-3 eV`
- `m3 = 7.904528763902e-3 eV`

Derived splittings:
- `Δm²21 = 4.761345443130e-5 eV²` (rel err `-36.768%`)
- `Δm²32 = 1.414611019861e-5 eV²` (rel err `-99.423%`)

Gate status:
- `hierarchy_ok = true`
- `ordering_ok = true`
- `solar_ok = false`
- `atmospheric_ok = false`
- `overall_pass = false`

## Interpretation
This is now a hard quantitative closure target for the neutrino mass-scale lane.
The binary blind lock remains hierarchy/type only, while this gate tracks the
remaining oscillation mismatch explicitly.

## Repro
```bash
cargo run -q -p gutoe-em --bin neutrino_oscillation_ci_gate
```
