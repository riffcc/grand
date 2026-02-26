# Finding 083: GRAND-346 Dark Sector Candidate Split from Z₃ Orbits

Date: 2026-02-26
Status: GRAND-346 in progress

## Goal

Start a concrete dark-matter lane from existing Cl(1,3)/Z₃ primitives without
adding free parameters or fake interaction assumptions.

## What landed

New Lean module:

- `lean/Gutoe/DarkMatterSector.lean`

Wired into root build:

- `lean/lakefile.lean` adds `Gutoe.DarkMatterSector`

Core definitions/theorems:

- `darkSectorCandidates := dualEmTriplet ∪ {15,16}`
- `smInteractionCarrier := leptonState ∪ quarkTriplet ∪ emTriplet ∪ magneticTriplet`
- `dark_sector_z3_closed`
- `dark_sector_disjoint_from_sm_carrier`
- `visible_dark_state_count_split` (`11` visible vs `5` dark candidates; disjoint; total `16`)
- `dark_to_visible_count_ratio_eq` (`5/11`)

Rust harness additions:

- `crates/gutoe-physics/src/dark_sector.rs`
- `crates/gutoe-physics/src/bin/dark_matter_report.rs`

Core runtime outputs:

- particle branch density map: `ρ_dark = (5/11) ρ_visible`
- geometric branch density map: `ρ_dark = (5/11) κ(r) ρ_visible`
- rotation/lensing proxies from enclosed mass
- CMB-era matter-fraction consistency report from structural ratio

Initial run artifact:

- `/tmp/bh_renders/dark_matter_report.txt`
- `/tmp/bh_renders/dark_matter_report.json`

Current headline from the report:

- Structural particle branch gives dark matter fraction `5/16 = 0.3125`
- Observed matter dark fraction baseline is `~0.8426`
- Delta is `~ -0.5301` (large mismatch; this lane is not closed yet)

## Why this is real progress

This isolates a structurally defined sector that is:

1. Z₃-stable
2. disjoint from the current SM interaction carrier orbits in the finite lane
3. counted from existing shared primitives only

No new physics constants were introduced.

## Honest boundary

This is a **candidate-sector isolation result**, not yet a cosmology fit.
Remaining GRAND-346 closure work:

1. Upgrade the report from proxy checks to dataset-backed scoring (rotation/lensing curves).
2. Tighten the geometric branch by deriving `κ(r)` from the cosmology/metric lane instead of using a report-level proxy.
3. Add structural weighting from energy-per-state dynamics (not just state counts) to test whether the 0.3125 fraction can move toward observed matter fractions without fit knobs.
4. Push branch comparison into a falsification gate with explicit pass/fail thresholds.

## Build sanity

- `lake build Gutoe.DarkMatterSector` ✅
- `lake build Gutoe` ✅
- `cargo check -p gutoe-physics --bin dark_matter_report` ✅

No `sorry` introduced.
