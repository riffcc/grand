# Finding 217 — Accelerated 50-Year CTC Campaign Fast Sim

## Summary

Added an integrated program-level simulator that compresses a long-horizon
time-travel/transport campaign into one executable lane:

- device physics constraints (rear shortcut, local bound, threshold lane),
- infrastructure growth and budget allocation,
- fabrication bottlenecks,
- safety maturation and mission operations,
- transport of simulated beings across campaign phases.

Binary:

- `crates/gutoe-physics/src/bin/ctc_50y_campaign_fast.rs`

## Run

```bash
cargo run -p gutoe-physics --bin ctc_50y_campaign_fast
```

Outputs:

- `/tmp/bh_renders/ctc_50y_campaign_fast/ctc_50y_campaign_fast.txt`
- `/tmp/bh_renders/ctc_50y_campaign_fast/ctc_50y_campaign_fast.json`

## Default-run highlights

- `years = 50`
- `missions_total = 158802.089`
- `transported_sim_beings_total = 7502385.083`
- `predeparture_fraction = 0.835663`
- `first predeparture_enabled year = 31`
- `safety_index_final = 1.0`

## Notes

- This is a systems/campaign simulation lane, not a physical engine claim.
- Predeparture operations are gated to late phases (`year >= 31`) and high
  safety maturity.
- Fabrication is bottlenecked by an explicit capacity term to avoid naive
  budget-only runaway device counts.

