# Finding 258 — Remaining12 Unified CI Gate

Date: 2026-03-01  
Runner: `cargo run -q -p gutoe-em --bin remaining12_gate`

## What was added

- Unified gate binary:
  - `crates/gutoe-em/src/bin/remaining12_gate.rs`

This gate combines both previously split endgame lanes:

1. **Neutrino endgame checks**
   - structural no-fit lane (hierarchy + tiny + ratio gate)
   - triangulated oscillation closure lane (ratio + absolute splitting machine-precision gates)
2. **Absolute-scale endgame checks**
   - electron anchor
   - lattice-vs-Fermi VEV agreement
   - lattice branch W/Z/H mass envelope
   - Fermi branch W/Z/H mass envelope

The gate emits unified artifacts and returns nonzero on failure.

## Artifacts

- `/tmp/bh_renders/remaining12_gate.txt`
- `/tmp/bh_renders/remaining12_gate.json`

## Gate result

- `overall_pass = true`
- Neutrino: `no_fit_pass=true`, `triangulated_pass=true`
- Absolute-scale: `overall_pass=true`

