# 002 — Alpha Consistency (Rust vs Lean)

## Claim

- Lean theorem chain (`lean/Gutoe/FineStructure.lean`) proves the **leading-order structural** result:
  - `alphaInverse 4 = 137`
  - therefore `alpha_LO = 1/137`.
- Rust runtime constants use the **measured low-energy** value:
  - `ALPHA = 1/137.036...` (encoded as `7.2973525693e-3`).

These are not contradictory; they represent different orders of approximation.

## Quantitative check

Relative difference between runtime and leading-order alpha:

`|alpha_runtime - alpha_LO| / alpha_LO ≈ 2.63e-4 ≈ 0.026%`

This matches the expected small correction scale discussed in the Lean comments
(higher-order QED/radiative effects beyond the leading structural value).

## Code parity hooks

- `crates/gutoe-physics/src/constants.rs`
  - `ALPHA` (runtime measured value)
  - `ALPHA_LEADING_ORDER` (proof-aligned value `1/137`)
- Unit test `test_alpha_runtime_vs_leading_order_offset_is_small` enforces that the offset remains small.
