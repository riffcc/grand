# 241 — Dual Scoreboard Encoded: Core Algebra vs Speculative Detection

## Scope
Implemented explicit dual scoreboards in the main falsification harness artifact:
- `core_algebraic_score`
- `speculative_detection_score`

Updated file:
- `crates/gutoe-physics/src/bin/ctc_falsification_20_harness.rs`

Artifact:
- `/tmp/bh_renders/ctc_falsification_20/ctc_falsification_20.json`
- `/tmp/bh_renders/ctc_falsification_20/ctc_falsification_20.txt`

## Current result snapshot
- `summary`: PASS=11, FAIL=4, OPEN=7, TOTAL=22
- `core_algebraic_score`: PASS=11, FAIL=0, OPEN=0, TOTAL=11
- `speculative_detection_score`: PASS=0, FAIL=4, OPEN=7, TOTAL=11

## Interpretation boundary
- Core algebraic lanes are cleanly separated and currently all pass.
- Speculative detection lanes are tracked independently and include null/blocked outcomes.
- This split prevents mixed interpretation of “theory-derived constants” vs “new-physics detection” status.

## Notes
- Weak-angle slope-lock (`T14`) is explicitly in `core_algebraic` and passes.
- Null replication (`T16`) is explicitly in `core_algebraic` and passes.
