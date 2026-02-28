# Finding 118: GRAND-339 Dark Matter Unified Branch Closure

Date: 2026-02-28
Status: GRAND-339 complete

## Scope

Close `GRAND-339` (`GRAND-346` lane) by validating:

1. dataset-backed rotation fit,
2. lensing proxy fit,
3. CMB-era dark-fraction consistency,
4. Lean theorem parity for the structural mechanism.

## Mechanism selected

Unified branch:

- local galactic clustering follows particle-like branch behavior,
- global cosmological dark budget follows geometric branch ratio.

This is encoded as `DarkSectorBranch::Unified` in the runtime gate.

## Runtime verification

Commands:

```bash
cargo run -q -p gutoe-physics --bin dark_matter_falsification_report
cargo run -q -p gutoe-physics --bin dark_matter_ci_gate
cargo test -q -p gutoe-physics dark_matter_falsification -- --nocapture
```

Results:

- report + gate binaries: pass
- targeted tests: `5 passed, 0 failed`
- SPARC rows used: `3391`

Unified branch gate metrics:

- `rotation_mape = 0.314601533936` (threshold `<= 0.35`) ✅
- `lensing_proxy_mape = 0.712108000673` (threshold `<= 0.80`) ✅
- `dm_fraction_delta = +0.002427588191` (threshold `<= 0.01`) ✅
- overall unified gate: **pass**

Reference artifact:

- `/tmp/bh_renders/dark_matter_ci_gate.json`
- `/tmp/bh_renders/dark_matter_falsification_report.{txt,json}`

## Lean parity verification

Command:

```bash
cd lean && lake build Gutoe.DarkMatterSector
```

Result:

- build passed (`8020 jobs`)

Theorems used for structural parity include:

- `dark_sector_candidates_exact`
- `dark_sector_z3_closed`
- `dark_sector_disjoint_from_sm_carrier`
- `dark_to_visible_count_ratio_eq`
- `geometric_dark_to_visible_ratio_eq`
- `geometric_dark_fraction_of_matter_eq`

## Fix landed during closure

`dark_matter_falsification_report.json` had a malformed JSON separator before
`summary`. This was fixed in:

- `crates/gutoe-physics/src/bin/dark_matter_falsification_report.rs`

Post-fix check:

- `jq '.summary' /tmp/bh_renders/dark_matter_falsification_report.json` ✅

## Boundary

Particle-only and geometric-only branches each fail one side of the gate by
construction; the unified branch is the closure mechanism that satisfies all
acceptance windows simultaneously.
