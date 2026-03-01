# 213 — Faithful time-travel simulation (bridge lane)

## New executable
- `crates/gutoe-physics/src/bin/ctc_faithful_bridge_sim.rs`

## Purpose
Single-run composite simulator for the already-proven theorem lanes:
1. Rear shortcut factor `s = 1/10`.
2. Local-causal bound (`v_local <= c`) with coordinate-effective `u = c/s`.
3. Subluminal boosted-frame predeparture witness (`dt' < 0`).
4. Timelike + identified closure witness on a time-cylinder step.
5. Path-A style effective-arrival loop summary.

## Default run result
- Output text:
  - `/tmp/bh_renders/ctc_faithful_bridge_sim/ctc_faithful_bridge_sim.txt`
- Output JSON:
  - `/tmp/bh_renders/ctc_faithful_bridge_sim/ctc_faithful_bridge_sim.json`

Key values (default):
- `s = 0.1`
- `u/c = 10.0`
- `local_bound_ok = true`
- `dynamic_gate_pass = true`
- `dt'_boosted < 0` (predeparture frame witness true)
- `timelike_step = true`, `identified_closed = true`
- `predeparture_effective = true`
- `faithful_sim_possible = true`

## Scope boundary
- This is a theorem-faithful consistency simulation.
- It is not a physical-engine claim.
