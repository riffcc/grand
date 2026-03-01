# 191 — Lentz Soliton Hard Zoom: Compact Energy Floor + Near-Core Dead End

Date: 2026-02-28

## Scope

Focused follow-up on the only prior matrix-passing engine (`Lentz-style warp soliton`) to answer:

1. Is it still viable under a quantitative compact-shell energy gate?
2. Does near-core exotic curvature help or hurt?

## What changed

### 1) Upgraded one-way Lentz probe

Extended:
- `crates/gutoe-physics/src/bin/ftl_one_way_lentz_probe.rs`

New outputs:
- Classical vs GUTOE-`G_eff` (scale-radius and Planck upper bound) requirement lanes.
- Tube geometry integrated-energy zoom.
- Compact self-contained shell sweep:
  - radius model (`R_curv = bubble radius`),
  - thickness model (`R_curv = wall thickness`),
  - core model (`R_curv = r_core = sqrt(C_inf) * l_P`).

Artifacts:
- `/tmp/bh_renders/ftl_one_way_lentz_probe/ftl_one_way_lentz_probe.txt`
- `/tmp/bh_renders/ftl_one_way_lentz_probe/ftl_one_way_lentz_probe.json`

### 2) Updated exhaustion matrix with quantitative Lentz gate

Extended:
- `crates/gutoe-physics/src/bin/ftl_engine_report.rs`

New framework gate fields:
- `lentz_compact_budget_met`
- `lentz_compact_floor_j`
- `lentz_compact_shortfall_ratio`
- `near_core_dead_end_multiplier`

And Lentz now includes:
- `requires_compact_positive_budget = true`

So matrix pass now requires not only topology/energy-condition assumptions, but also clearing a compact-shell positive-energy budget floor.

Artifact:
- `/tmp/bh_renders/ftl_engine_report/ftl_engine_report.{txt,json}`

## Key numbers

From compact-shell sweep (`rear10`, `v_eff = 2c`):

- Best radius-model compact shell:
  - `E ≈ 6.057e41 J`
  - `~5.014e7` solar-years
  - `(R=100 m, thickness=0.1 m)`

- Best thickness-model compact shell:
  - `E ≈ 1.916e44 J`
  - `~1.586e10` solar-years
  - `(R=10 m, thickness=5 m)`

- Near-core model (`R_curv = r_core ≈ 1.195e-35 m`):
  - `E ≈ 4.280e113 J`
  - `~3.543e79` solar-years
  - catastrophic increase relative to compact radius-model.

From matrix gate:

- `lentz_compact_floor_j ≈ 6.057e41 J`
- `lentz_compact_shortfall_ratio ≈ 1.483e41`
- `near_core_dead_end_multiplier ≈ 7.067e71`
- `lentz_compact_budget_met = false`

## Interpretation

1. The previous Lentz matrix pass was assumption-level (`allows_ctc=true`, no negative-energy/ANEC blockers).
2. Under quantitative compact-shell budgeting, Lentz fails by ~`1e41` against the optimistic positive-energy source lane.
3. Near-core exotic curvature is a dead end in this lane:
   the same UV/singularity-resolution structure that regularizes geometry makes extreme-curvature propulsion prohibitively expensive.

This is a consistency result, not a bug:
- singularity resolution protects spacetime from exploitative infinite-curvature channels.

## Build/run verification

Executed:
- `cargo run -q -p gutoe-physics --bin ftl_one_way_lentz_probe`
- `cargo run -q -p gutoe-physics --bin ftl_engine_report`
- `GUTOE_ALLOW_CTC=1 cargo run -q -p gutoe-physics --bin ftl_engine_report`

Observed:
- `ftl matrix: pass=0 fail=9 (total=9)` after quantitative Lentz gate integration.

