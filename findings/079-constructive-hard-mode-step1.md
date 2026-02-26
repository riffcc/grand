# Finding 079: GRAND-316 Constructive Hard-Mode Step 1

Date: 2026-02-26
Status: GRAND-316 complete

## Goal

Start replacing the broad `hTargets : constructiveTargetsSatisfied M` input
with a canonical constructive model where target checklist items are
discharged directly from Wilson-equivalence and Wilson-kernel theorems.

## Module added

- `lean/Gutoe/YangMillsConstructiveHardMode.lean`

Also added to roots in:

- `lean/lakefile.lean`

## What landed

### 1) Canonical hard-mode model + sequences

- `hardModeEpsSeq`
- `hardModeGapSeq`
- `HardModeCoreObligations`
- `hardModeInterfaces`
- `hardModeModel`

The model now encodes:

- reflection-positivity proxy as gap nonnegativity
- cluster proxy as non-vanishing gap sequence
- Schwinger-existence proxy as per-step positive gap
- Euclidean-invariance proxy from `wilson_kernel_row_offset_invariant`
- regularity proxy from `wilson_kernel_row_sum_one`

with only these residual core obligations explicit:

- `osReconstruction`
- `wightmanCompatibility`

### 2) Domain-driven discharge theorems

- `hard_mode_reflection_positivity_of_domain`
- `hard_mode_schwinger_exists_of_domain`
- `hard_mode_cluster_property_of_domain`

These derive three checklist items from `WilsonEquivalenceDomain` via
`c3_gap_correspondence_of_domain` and existing non-vanishing lemmas.

Additional unconditional discharge theorems:

- `hard_mode_euclidean_invariance`
- `hard_mode_regularity`

### 3) Target-assumption reduction theorem

- `constructive_targets_satisfied_of_hard_mode_core`

From Wilson-domain assumptions plus residual `HardModeCoreObligations`,
it constructs full `constructiveTargetsSatisfied` for the canonical hard-mode
model, without separately assuming:

- reflection positivity
- Euclidean invariance
- regularity
- cluster property
- Schwinger existence

### 4) Hard-mode closure in the constructive lane

- `mass_gap_embedded_of_hard_mode_core`
- `constructive_lane_gap_closure_of_hard_mode_core`

These instantiate the existing GRAND-302 closure path with the hard-mode model.

## Why this matters

This is a genuine assumption-surface reduction step for GRAND-302:

- before: all constructive targets entered as a single external assumption
- now: three targets are discharged from the theorem chain itself, and only
  core OS/Wightman obligations remain explicit.

## Build sanity

- `lake build Gutoe.YangMillsConstructiveHardMode` ✅
- `lake build Gutoe` ✅

No `sorry` introduced.
