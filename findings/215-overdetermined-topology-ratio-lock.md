# 215 — Overdetermined topology ratio lock (checkmate constraints lane)

## New Lean module
- `lean/Gutoe/OverdeterminedTopologyRatios.lean`

## New Rust probe
- `crates/gutoe-physics/src/bin/topology_overdetermination_probe.rs`

## What is now proven in Lean
1. Structural ratios are recovered from shared counting primitives:
   - `branching = 3`
   - `void = 3/16`
   - `eta = 4/6 = 2/3`
   - `infra = 16/6 = 8/3`
2. Topology closure identity:
   - `G = branching * void * eta * infra = 1`
   - equivalent `144/144 = 1`
3. Compatibility with existing lanes:
   - `topologyGainQ = geffZ3VoidSplitQ`
   - `topologyGainQ = 1` and `geffCanonicalQ = 1`
4. Overdetermination theorems:
   - If `G=1` and three factors are fixed, the fourth is uniquely forced.
   - Implemented for all four factors:
     - `infra_forced_by_unit_gain`
     - `eta_forced_by_unit_gain`
     - `void_forced_by_unit_gain`
     - `branching_forced_by_unit_gain`

Interpretation:
- Ratios are no longer tunable knobs; they are constrained variables in a
  solved system.

## Runtime confirmation
- Output text:
  - `/tmp/bh_renders/topology_overdetermination_probe/topology_overdetermination_probe.txt`
- Output JSON:
  - `/tmp/bh_renders/topology_overdetermination_probe/topology_overdetermination_probe.json`

Default probe result:
- `G = 1.0` with residual `0.0`
- inferred values from three-factor inversion exactly match structural values:
  - inferred `infra = 8/3`
  - inferred `eta = 2/3`
  - inferred `void = 3/16`
  - inferred `branching = 3`

## Build status
- `lake build Gutoe.OverdeterminedTopologyRatios` passed.
- `lake build Gutoe` passed (warnings only).
