# Finding 018: Gauge Coupling Running / Unification Snapshot

Date: 2026-02-24
Issues: GRAND-63, GRAND-133

## Method

Executed one-loop gauge-coupling scan binary:

- `cargo run -p gutoe-physics --bin coupling_unification_report`

This scans representative MS-bar couplings from `M_Z` to high scale and reports
where the three inverse couplings come closest.

## Result

From the fresh run:

- Best unification-like scale: `mu* = 2.5119e14 GeV`
- Minimum spread in inverse couplings: `3.688117`
- Values at `mu*`:
  - `alpha1^-1 = 40.399758`
  - `alpha2^-1 = 44.074874`
  - `alpha3^-1 = 40.386757`

## Artifacts

- `/tmp/bh_renders/coupling_unification_report.csv`
- `/tmp/bh_renders/coupling_unification_summary.txt`

## Interpretation

At one-loop with current coefficients/input values, couplings approach but do not
exactly meet at a single point. The tool now gives a reproducible quantitative
answer to “if/where they approximately unify,” with concrete scale and residual spread.
