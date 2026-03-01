# 192 — Escher Stair CTC: Finite Local, Unbounded Global

## What was added

### Lean (formal)
File: `lean/Gutoe/CTCLegality.lean`

Added theorems:
- `fixed_space_n_steps_timelike`
- `fixed_space_n_steps_interval_exact`
- `same_on_time_cylinder_after_n_periods`
- `escher_stair_n_laps`
- `escher_stair_cover_time_unbounded`
- `ctc_step_forward_but_identified`

Interpretation:
- Local steps remain timelike (`ds² < 0`) and finite.
- Time-cylinder identification closes the loop after integer laps.
- Covering-space coordinate time is unbounded with lap count.

### Rust (quant)
File: `crates/gutoe-physics/src/bin/ctc_time_cylinder_sim.rs`

Added:
- `GUTOE_CTC_LAPS` input
- multi-lap worldline generation
- `escher_stair` report block in txt/json:
  - `proper_time_per_lap`
  - `cover_time_total`
  - `finite_local_lap`
  - `unbounded_global_cover_witness`

## Verification

- Lean module build: `lake build Gutoe.CTCLegality` ✅
- Full Lean build: `lake build Gutoe` ✅ (`8139 jobs`)
- Rust run:
  - `GUTOE_CTC_LAPS=10000 GUTOE_CTC_PERIOD_T=1 GUTOE_CTC_STEPS=200 cargo run -q -p gutoe-physics --bin ctc_time_cylinder_sim`

Observed report:
- `timelike_all_segments = true`
- `identified_closed = true`
- `proper_time_per_lap = 1.0`
- `cover_time_total = 10000.0`
- `unbounded_global_cover_witness = true`

## Honest boundary

This lane proves/measures **kinematic topology behavior** (finite local timelike laps with unbounded covering progression under identification). It does **not** yet prove a physical mechanism for dynamic topology creation.
