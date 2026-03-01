# Finding 231 — Falsification-20 Live-Fire Harness

## Summary

Implemented and executed a 20-test rapid falsification harness that encodes:

- dataset lane
- test metric
- explicit pass threshold
- explicit kill threshold
- current status (`PASS` / `FAIL` / `OPEN`)

Binary:

- `crates/gutoe-physics/src/bin/ctc_falsification_20_harness.rs`

Outputs:

- `/tmp/bh_renders/ctc_falsification_20/ctc_falsification_20.txt`
- `/tmp/bh_renders/ctc_falsification_20/ctc_falsification_20.json`

## Run

```bash
cargo run -p gutoe-physics --bin ctc_falsification_20_harness
```

## Current Scoreboard

- `PASS = 3`
- `FAIL = 1`
- `OPEN = 16`
- `TOTAL = 20`

Closed tests (auto-populated from existing artifacts):

1. `T01` Kill naive weak-angle `3/13` lane: **PASS**  
   - observed `reduced_chi2(base_fixed) = 16.910485`

2. `T02` Zero-free weak-angle survival: **PASS**  
   - observed `reduced_chi2(zero_free_formula) = 1.155080`

3. `T03` Zero-free weak-angle max pull guard: **PASS**  
   - observed `max_abs_pull_sigma = 2.142615`

4. `T04` Spontaneous-door global excess: **FAIL**  
   - observed `z_stouffer = 0.796188`, `door_detected = false`

## Notes

- This harness is intentionally strict about turning ideas into kill/pass gates.
- Most of the campaign remains in `OPEN` state pending additional public-lane
  extraction (especially CB/void scalarization and EHT joint residual lanes).
- The harness is designed for fast iteration: as lane artifacts land, statuses
  move from `OPEN` to `PASS`/`FAIL` without changing gate definitions.
