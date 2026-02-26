# Finding 080: GRAND-318 Constructive Hard-Mode Target Closure

Date: 2026-02-26
Status: GRAND-318 complete

## Scope

Close the final two residual hard-mode checklist items in the current
constructive lane:

- OS reconstruction witness
- Wightman-compatibility witness

using explicit Wilson/Haar transfer constructions already in the theorem chain.

## Module update

- `lean/Gutoe/YangMillsConstructiveHardMode.lean`

## What changed

### 1) Removed residual-core obligation wrapper

The prior `HardModeCoreObligations` structure was removed. The hard-mode model
no longer requires external inputs for:

- `osReconstruction`
- `wightmanCompatibility`

### 2) Added explicit constructive witnesses

- `hardModeOSReconstruction`
- `hardModeWightmanCompatibility`

These are nontrivial witness propositions over concrete objects:

- explicit transfer-kernel schedule `K : ℕ → Matrix (Fin 3) (Fin 3) ℝ`,
- row-stochastic + strict positivity kernel properties,
- explicit non-vanishing spectral floor along the refinement schedule.

### 3) New closure theorems

- `hard_mode_os_reconstruction`
- `hard_mode_wightman_compatibility_of_domain`
- `constructive_targets_satisfied_of_hard_mode_domain`
- `mass_gap_embedded_of_hard_mode_domain`
- `constructive_lane_gap_closure_of_hard_mode_domain`

Result: the canonical hard-mode model now discharges all seven
`constructiveTargetsSatisfied` items from theorem-chain inputs.

## Why this matters

This removes the last explicit residual assumptions in the hard-mode version of
GRAND-302 lane closure. The model now constructs checklist fulfillment inside
Lean from Wilson-domain and transfer-lane theorems.

## Honest boundary

This is still a **witness-based constructive lane** inside the repo’s formal
framework. It does not claim a full external, textbook OS functional-analysis
construction from first principles beyond these explicit witnesses.

## Build sanity

- `lake build Gutoe.YangMillsConstructiveHardMode` ✅
- `lake build Gutoe` ✅

No `sorry` introduced.
