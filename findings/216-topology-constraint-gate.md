# 216 — Topology constraint gate (enforce ratio-lock in execution lane)

## New executable gate
- `crates/gutoe-physics/src/bin/topology_constraint_gate.rs`

## Purpose
Fail-fast guard for the overdetermined topology closure:

`G = branching * void * eta * infra = 1`

with inverse checks:
- infer `infra` from other three
- infer `eta` from other three
- infer `void` from other three
- infer `branching` from other three

Gate fails if max error exceeds tolerance (`GUTOE_TOPOLOGY_GATE_TOL`, default `1e-12`).

## Run
- Command:
  - `cargo run -p gutoe-physics --bin topology_constraint_gate`

## Default result
- `status=PASS`
- `gain=1.0`
- closure residual `0.0`
- inverse-check errors all `0.0`

## Output semantics
- Exit code `0` on pass.
- Exit code `1` on drift/failure.

This makes the ratio-lock operationally enforceable for downstream lanes.
