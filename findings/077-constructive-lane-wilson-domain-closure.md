# Finding 077: GRAND-302 Constructive Lane Closure via Wilson-Equivalence Domain

Date: 2026-02-26
Status: GRAND-302 complete

## Updated module

- `lean/Gutoe/YangMillsConstructiveQFT.lean`

## What landed

### 1) Constructive embedding tied directly to Theorem-C domain

New theorem:
- `mass_gap_embedded_of_wilson_equivalence_domain`

This threads:
- constructive targets (`ConstructiveYMModel` checklist)
- Wilson-equivalence domain (`a_t > 0`, bounded cap, `alpha > 0`)

into a single constructive-lane mass-gap embedding with
`eps n = minorizationEps (wilsonRowTotalsSchedule W n) alpha`.

### 2) Wilson-specialized closure theorem in constructive lane

New theorem:
- `constructive_lane_gap_closure_of_wilson_equivalence_domain`

This yields both:
- per-step positivity of Doeblin lower-bound gap
- non-vanishing sequence result (`¬ TendsToZeroSeq ...`)

in the same constructive framework.

## Why this matters

- GRAND-302 is no longer just an interface checklist.
- The constructive lane is now explicitly connected to GRAND-300's Wilson
  equivalence domain and GRAND-299 continuum-survival machinery.
- Mass-gap consequences are proved *inside* the constructive lane, not externally.

## Explicit boundary (honest)

This closes the constructive-lane theorem plumbing and embedding requirements.
It does **not** claim full standalone constructive YM existence has been proved
from scratch; the OS/milestone targets remain explicit assumptions in the model
interface by design.

## Build sanity

- `lake build Gutoe.YangMillsConstructiveQFT` ✅
- `lake build Gutoe` ✅

No `sorry` introduced.
