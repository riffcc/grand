# 146 - Blind Neutrino Register (Zero-Parameter Lock)

## Why this ticket exists
A serious critique is that many strong results can still be interpreted as structured postdiction unless at least one high-value prediction is frozen *before* external resolution and CI-locked against drift.

## Candidate selection
Chosen blind candidate: `BLIND-NEUTRINO-001`

- Observable A: neutrino hierarchy (`normal` vs `inverted`)
- Observable B: mass character (`dirac` vs `majorana`-like)

Why this lane:
- low knob count (binary outcomes)
- no seeded Ki/drug-anchor dependence
- already connected to core Cl(1,3)->flavor structure
- clean falsification path via external measurements

Lock scope (updated): hierarchy + mass-character only.
Absolute mass outputs are retained as advisory diagnostics and are not part
of the blind lock pass/fail.

## Frozen outputs
From `/tmp/bh_renders/blind_prediction_register/blind_prediction_register.json`:

- `hierarchy_prediction = normal`
- `mass_character_prediction = dirac`
- `majorana_symmetry_residual = 0.9948906419732648`
- `m1_ev = 8.497119214461528e-4`
- `m2_ev = 6.952371162469944e-3`
- `m3_ev = 7.904528763902007e-3`
- `sum_ev = 1.5706611847818103e-2`

## Falsification criteria
1. Measured hierarchy is not `normal` -> lane falsified.
2. Robust neutrinoless double-beta signal consistent with Majorana neutrinos -> this `dirac` lock falsified.

## CI enforcement
New bins:
- `blind_prediction_register`
- `blind_prediction_register_ci_gate`

Gate JSON:
- `/tmp/bh_renders/blind_prediction_register_ci_gate.json`
- `overall_pass` must stay true for hierarchy/type checks only.

Global integration:
- `global_gate_report` now runs both blind bins and requires `blind_prediction_pass=true` in `overall_pass`.

## Repro commands
```bash
cargo run -q -p gutoe-physics --bin blind_prediction_register
cargo run -q -p gutoe-physics --bin blind_prediction_register_ci_gate
cargo run -q -p gutoe-physics --bin global_gate_report
```
