# Finding 019: Weinberg M_Z Correction Budget

Date: 2026-02-24
Issue: GRAND-61

## Goal

Quantify the exact correction budget required to move from the structural
Lean/Rust value `sin²(theta_W)=3/13` to the observed `M_Z` value.

## Command

- `cargo run -p gutoe-em --bin weinberg_mz_report`

## Output

- Structural: `sin²θ_W = 0.230769231`
- Observed (M_Z): `sin²θ_W = 0.231220000`
- Delta: `+0.000450769`
- Relative deviation: `0.194953%`

Artifacts:
- `/tmp/bh_renders/weinberg_mz_report.csv`
- `/tmp/bh_renders/weinberg_mz_summary.txt`

## Interpretation

This is the quantitative target for GRAND-61: the RG/loop flow model must
account for a +0.000450769 shift from the structural tree value at M_Z.
