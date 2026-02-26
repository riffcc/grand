# Finding 074: GRAND-312 Common-Factor Discharge from Quotient Normalization

Date: 2026-02-26
Status: GRAND-312 complete

## Updated module

- `lean/Gutoe/HaarFiberCollapse.lean`

## What landed

### 1) Derived common-factor hypotheses (no direct `hInt/hMass` assumption)

New theorem:
- `common_factor_hypotheses_of_center_quotient_normalization`

Given center-quotient decomposition and explicit fiber normalization equalities:
- `fiberExpectation f = (fun q => c * fQ q)`
- `fiberExpectation 1 = (fun _ => c)`

it derives:
- `expectation μG f = c * expectation μQ fQ`
- `(μG Set.univ).toReal = c * (μQ Set.univ).toReal`

with `μQ = quotientFiberMeasure μG 𝓕`.

### 2) Closed normalized reduction from derived data

New theorem:
- `normalized_expectation_reduce_to_center_of_quotient_normalization`

This discharges GRAND-311's remaining explicit `hInt/hMass` seam from quotient
normalization structure, then applies the existing normalized collapse theorem.

## Why this matters

- GRAND-311 isolated the seam (`hInt/hMass`).
- GRAND-312 now proves that seam from the center-quotient normalization layer.
- The collapse path is tighter and more structural: decomposition + normalization
  data -> common factors -> normalized equality.

## Build sanity

- `lake build Gutoe.HaarFiberCollapse` ✅
- `lake build Gutoe` ✅

No `sorry` introduced.
