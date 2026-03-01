# Finding 218 — Retro-Retry Closure Theorem + Campaign Eventual Success Metrics

## Summary

Added a formal retry-closure lane in Lean and mirrored it in the accelerated
campaign simulator as explicit eventual/asymptotic success metrics.

## Lean

New module:

- `lean/Gutoe/CTCRetroRetryClosure.lean`

Added to roots in:

- `lean/lakefile.lean`

Key theorems:

- `retry_partition`
- `retry_success_monotone`
- `retry_success_strict_step`
- `retry_success_tendsto_one`
- `exists_retry_count_for_target`
- `structural_seed_eventual_closure`

Interpretation:

- With `0 < p < 1`, retry success `S_n = 1 - (1-p)^n` is increasing and tends to `1`.
- For any target `< 1`, there exists finite retry depth exceeding that target.

## Rust

Updated binary:

- `crates/gutoe-physics/src/bin/ctc_50y_campaign_fast.rs`

New per-year metrics:

- `retry_depth`
- `eventual_success_prob`
- `asymptotic_success_prob`
- `missions_eventual_success`
- `missions_asymptotic_success`

New campaign totals:

- `mission_success_rate_first_pass`
- `mission_success_rate_eventual`
- `mission_success_rate_asymptotic`
- `transported_sim_beings_total_first_pass`
- `transported_sim_beings_total_eventual`
- `transported_sim_beings_total_asymptotic`

## Verification

- `lake build Gutoe.CTCRetroRetryClosure` passed.
- `lake build Gutoe` passed.
- `cargo run -p gutoe-physics --bin ctc_50y_campaign_fast` passed.

Default run excerpt:

- `mission_success_rate_first_pass = 0.944872`
- `mission_success_rate_eventual = 0.959440`
- `mission_success_rate_asymptotic = 0.959440`
- `transported_sim_beings_total_first_pass = 7502385.083`
- `transported_sim_beings_total_eventual = 7618055.118`
- `transported_sim_beings_total_asymptotic = 7618055.118`

