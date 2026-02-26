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

## Why this is real progress

This isolates a structurally defined sector that is:

1. Z₃-stable
2. disjoint from the current SM interaction carrier orbits in the finite lane
3. counted from existing shared primitives only

No new physics constants were introduced.

## Honest boundary

This is a **candidate-sector isolation result**, not yet a cosmology fit.
Remaining GRAND-346 closure work:

1. Add a phenomenology harness (rotation-curve/lensing/CMB-era matter fraction tests).
2. Prove/derive the gravitational-source map from this sector into the cosmology pipeline.
3. Test particle-vs-geometric branch explicitly and keep both falsifiable.

## Build sanity

- `lake build Gutoe.DarkMatterSector` ✅
- `lake build Gutoe` ✅

No `sorry` introduced.
