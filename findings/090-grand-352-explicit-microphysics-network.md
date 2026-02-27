# Finding 090 — GRAND-352 Explicit BBN/Recombination Microphysics

Date: 2026-02-27
Status: COMPLETE

## Summary
Added an explicit microphysics lane (time-evolving BBN reaction network + recombination/opacity evolution) and wired it into universe assembly and CI.

New files:
- `crates/gutoe-physics/src/microphysics.rs`
- `crates/gutoe-physics/src/bin/microphysics_report.rs`
- `crates/gutoe-physics/src/bin/microphysics_ci_gate.rs`

Universe integration:
- `UniverseScorecard` now includes:
  - `microphysics_ok`
  - `microphysics: MicrophysicsScorecard`
- `passes_early_universe()` now requires `microphysics_ok`

## Network behavior
BBN lane evolves baryon fractions in:
- free neutrons/protons
- deuterium
- helium-4

Recombination lane evolves:
- ionization fraction `x_e(z)`
- opacity/optical-depth proxy
- visibility peak redshift

## Verified outputs (default run)
From `cargo run -q -p gutoe-physics --bin microphysics_report`:
- `Y_p_network = 0.24841`
- `D/H_network = 2.035e-5`
- `z_visibility_peak = 1063.0`
- Gate pass: `true`

Artifacts:
- `/tmp/bh_renders/microphysics_report.txt`
- `/tmp/bh_renders/microphysics_report.json`
- `/tmp/bh_renders/microphysics_ci_gate.json`

## Notes
- This is an explicit dynamical lane (not a checkpoint-only anchor).
- Default windows remain falsifiable (`Y_p`, `D/H`, recombination visibility region).
