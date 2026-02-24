# 014 — Coupling Unification Scan, Locality Gate, and Decoherence Gate

Issues addressed:
- GRAND-63 (`Coupling unification at GUT scale`)
- GRAND-118 (`Locality or emergent locality from lattice`)
- GRAND-120 (`Classical emergence / decoherence`)

## GRAND-63: Coupling unification scan

Added executable:
- `crates/gutoe-physics/src/bin/coupling_unification_report.rs`

Run:
- `cargo run -p gutoe-physics --bin coupling_unification_report`

Artifacts:
- `/tmp/bh_renders/coupling_unification_report.csv`
- `/tmp/bh_renders/coupling_unification_summary.txt`

Observed summary (current one-loop scan):
- best unification-like scale `mu ~= 2.5119e14 GeV`
- minimum spread in `alpha^-1`: `~3.688`
- `alpha^-1(mu*) ~= (40.3998, 44.0749, 40.3868)`

Interpretation:
- One-loop trajectories approach but do not exactly meet.
- This provides a reproducible baseline and identifies the remaining mismatch to close in higher-order/threshold passes.

## GRAND-118: Locality gate

Code anchor:
- `crates/gutoe-em/src/geometry.rs`

Runs:
- `cargo test -p gutoe-em neighbours_are_intra_layer -- --nocapture`
- `cargo test -p gutoe-em each_site_has_six_neighbours -- --nocapture`

Result:
- Both pass.
- Confirms finite local stencil (`k=6`) and no unintended cross-layer coupling in `mesh_neighbours`.

## GRAND-120: Classical emergence / decoherence gate

Code anchor:
- `crates/gutoe-em/src/quantum_lepton.rs`

Runs:
- `cargo test -p gutoe-em arrow_of_time_free_vs_bound -- --nocapture`
- `cargo test -p gutoe-em entropy_extremes -- --nocapture`

Result:
- Both pass.
- Quantitative behavior:
  - free evolution: entropy increases (`ΔS > 0`)
  - bound Coulomb evolution: entropy decreases (`ΔS < 0`)

Interpretation:
- The runtime demonstrates an explicit entropy-direction split supporting classical emergence/decoherence narratives in the current model.

