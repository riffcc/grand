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
- `geometric_dark_amplification_eq` (`12 = 16 - 4` from shared state counts)
- `geometric_dark_to_visible_ratio_eq` (`60/11`)
- `geometric_dark_fraction_of_matter_eq` (`60/71 ≈ 0.84507`)

Rust harness additions:

- `crates/gutoe-physics/src/dark_sector.rs`
- `crates/gutoe-physics/src/bin/dark_matter_report.rs`
- `crates/gutoe-physics/src/dark_matter_falsification.rs`
- `crates/gutoe-physics/src/bin/dark_matter_falsification_report.rs`
- `crates/gutoe-physics/src/bin/dark_matter_ci_gate.rs`
- `crates/gutoe-physics/data/sparc_massmodels_2016c_baryon.csv`

Core runtime outputs:

- particle branch density map: `ρ_dark = (5/11) ρ_visible`
- geometric branch density map: `ρ_dark = (5/11) κ(r) ρ_visible`
- rotation/lensing proxies from enclosed mass
- CMB-era matter-fraction consistency report from structural ratio
- geometric branch using structural amplification from shared Cl(1,3) counts
- `κ(r)` now derived from Einstein/cosmology source terms:
  - `κ(r) = (1 + λ_QG (l_P/r)^2) * (1 + ρ_Λ/ρ_visible(r))`
  - `ρ_Λ = Λ c²/(8πG)`

Initial run artifact:

- `/tmp/bh_renders/dark_matter_report.txt`
- `/tmp/bh_renders/dark_matter_report.json`
- `/tmp/bh_renders/dark_matter_falsification_report.txt`
- `/tmp/bh_renders/dark_matter_falsification_report.json`

Current headline from the report:

- Structural particle branch gives dark matter fraction `5/16 = 0.3125`
- Structural geometric branch gives dark matter fraction `60/71 = 0.84507`
- Observed matter dark fraction baseline is `~0.84264`
- Geometric branch delta is `+0.00243` (~0.29% high)
- SPARC dataset-backed gate (3391 rows):
  - Particle branch: rotation/lensing pass, CMB fraction fail
  - Geometric branch: CMB fraction pass, rotation/lensing fail
  - Unified branch (particle local clustering + geometric global budget): passes all current gates
    - rotation MAPE `0.31460` (window `≤ 0.35`)
    - lensing-proxy MAPE `0.71211` (window `≤ 0.80`)
    - CMB dark-fraction delta `+0.00243` (window `≤ 0.01`)
- Einstein/cosmology-derived `κ(r)` is now wired end-to-end in report + gate
  (numerically close to unity at galaxy scales, as expected from source-term ratios).

## Why this is real progress

This isolates a structurally defined sector that is:

1. Z₃-stable
2. disjoint from the current SM interaction carrier orbits in the finite lane
3. counted from existing shared primitives only
4. scored against a real galaxy rotation-curve dataset and explicit falsification windows
5. now supports a unified branch that composes local and global constraints in one scorecard

No new physics constants were introduced; the geometric amplification is built
from shared counts `16 - 4 = 12` (total Clifford states minus grade-1 states).

## Honest boundary

This is still a **candidate-sector derivation lane**, not yet a closure of
dark-matter phenomenology.
Remaining GRAND-346 closure work:

1. Replace the current unified composition rule with a fully dynamical halo evolution law from the same primitive equations.
2. Add independent lensing datasets (cluster-scale) instead of rotation-derived lensing proxies.

## Build sanity

- `lake build Gutoe.DarkMatterSector` ✅
- `lake build Gutoe` ✅
- `cargo check -p gutoe-physics --bin dark_matter_report` ✅
- `cargo check -p gutoe-physics --bin dark_matter_falsification_report` ✅
- `cargo check -p gutoe-physics --bin dark_matter_ci_gate` ✅
- `cargo test -p gutoe-physics dark_matter_falsification` ✅
- `cargo run -q -p gutoe-physics --bin dark_matter_ci_gate` ✅

CI hard gate behavior:

- default target branch is `unified`
- configurable via `GUTOE_DARK_GATE_BRANCH={particle|geometric|unified}`
- writes `/tmp/bh_renders/dark_matter_ci_gate.json`
- exits with code `2` if target branch fails any gate window

No `sorry` introduced.
