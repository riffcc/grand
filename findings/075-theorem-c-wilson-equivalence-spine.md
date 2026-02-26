# Finding 075: GRAND-300 Theorem-C Wilson Equivalence Spine

Date: 2026-02-26
Status: GRAND-300 closure candidate

## New module

- `lean/Gutoe/YangMillsWilsonEquivalence.lean`

## What landed

### 1) Precise equivalence domain/limits statement

- `WilsonEquivalenceDomain`

This encodes exactly the bridge domain:
- `a_t n > 0` (all refinement steps)
- bounded refinement cap `∃ aCap > 0, a_t n ≤ aCap`
- `alpha > 0`

### 2) Action correspondence lemma

- `action_correspondence_of_domain`

For each refinement step:

`wilsonKernel 1 (centerPlaquetteActionSchedule W alpha n)`
`= smoothedTransition (z3NearestNeighborCounts ...) (rowTotalsFromCounts ...) alpha`

This is the Wilson-action ↔ transfer-kernel action layer in theorem form.

### 3) Measure correspondence lemma

- `measure_correspondence_of_center_quotient_normalization`

This imports GRAND-312's discharge and proves normalized expectation equality
(full lane vs center quotient lane) from quotient-normalization data.

### 4) Correlator correspondence lemma

- `twoPointKernelObservable`
- `correlator_correspondence_finite_of_row_scale_orbit`

Finite two-point correlator collapse is proved on the transfer lane under the
row-scale orbit hypothesis (same structural condition as finite fiber collapse).

### 5) Consolidated Theorem-C spine

- `theorem_c_wilson_equivalence_domain_limits`

Produces in one theorem:
- action correspondence across all refinements
- non-vanishing continuum Doeblin gap lower bound on the Wilson lane

## Explicit residual assumptions (honest boundary)

Remaining assumptions are explicit and localized:

- Measure correspondence still requires quotient-normalization data:
  - `hFiberObs`, `hFiberMass`, `hc`, `hMassQ`
- Finite correlator correspondence requires row-scale orbit hypothesis:
  - `hscale`

No hidden assumptions are introduced; all are theorem arguments.

## Build sanity

- `lake build Gutoe.YangMillsWilsonEquivalence` ✅
- `lake build Gutoe` ✅

No `sorry` introduced.
