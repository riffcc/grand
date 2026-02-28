# GRAND-346 — Alpha Two-Term Lane Closure (`137 + 5α - 9α²`)

## Goal

Promote the two-term correction lane from exploratory output to formal/runtime parity with an explicit regression gate.

## Lean Formalization

File:

- `lean/Gutoe/FineStructure.lean`

Added/verified theorem chain:

- `triangular_clifford_dim_plus_one_eq_137`
- `alpha_first_order_explicit`
- `alpha_second_order_explicit`
- `alpha_second_order_within_2e5_band`
- `alpha_second_order_closer_than_first`

Interpretation:

- leading term is fixed by Cl(1,3): `T(2^4)+1=137`
- first-order lane: `137 + 5/137`
- second-order lane: `137 + 5/137 - 9/137^2`
- second-order remains closer than first-order to the decimal reference lane.

Build proof gate:

- `cd lean && lake build Gutoe` -> success.

## Runtime Parity + Regression Gate

File:

- `crates/gutoe-physics/src/bin/alpha_web_ci_report.rs`

Added runtime gate term:

- `TWO_TERM_ABS_ERR_MAX = 1.0e-5`
- boolean `alpha_second_order_within_band`
- `passes_all` now requires:
  - identity checks,
  - lane sanity,
  - second-order improvement,
  - second-order residual within band.

Repro command:

- `cargo run -q -p gutoe-physics --bin alpha_web_ci_report`

Artifacts:

- `/tmp/bh_renders/alpha_web_ci_report/alpha_web_ci_report.txt`
- `/tmp/bh_renders/alpha_web_ci_report/alpha_web_ci_report.json`

Observed key values (fresh run):

- `delta_target = 3.5999084e-2`
- `delta_first_order_5alpha = 3.6486763e-2`
- `delta_second_order_5alpha_minus_9alpha2 = 3.6007501e-2`
- `first_abs_error = 4.8767885e-4`
- `second_abs_error = 8.4166557e-6`
- `second_abs_error_band_max = 1.0e-5`
- `second_order_improves = true`
- `second_order_within_band = true`
- `ci_gate.passes_all = true`

## Acceptance Check

- `lake build Gutoe` passes with theorem scaffold and no new `sorry`: **done**
- one command reproduces two-term artifact: **done**
- quantitative improvement vs `137 + 5α` shown: **done**
