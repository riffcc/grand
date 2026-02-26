# Finding 082: GRAND-321 OS Hilbert Completion + Self-Adjoint Generator

Date: 2026-02-26
Status: GRAND-321 in progress

## Goal

Start GRAND-321 assembly with concrete Lean objects for:

1. OS quotient Cauchy completion
2. strong semigroup lane on Hilbert realization
3. self-adjoint generator extraction
4. strictly positive generator gap floor

## New module

- `lean/Gutoe/YangMillsOSCompletion.lean`
- root added in `lean/lakefile.lean`

## What landed

### 1) Cauchy completion object for `OSHilbertQuot`

- `osQuotPseudoMetric`
- `OSCauchyCompletion`
- `osQuot_dense_in_completion`

This introduces a concrete completion type
`UniformSpace.Completion (OSHilbertQuot K)` with dense canonical embedding.

### 2) Quotient → Hilbert range realization

- `kernelRangeHilbert`
- `osQuotToRangeHilbert`
- `osQuotToRangeHilbert_surjective`

This gives an explicit Hilbert-range carrier tied to the kernel image map.

### 3) Semigroup + generator assembly

- `StronglyContinuousSemigroup` (explicit interface)
- `scalarSemigroupOp`, `scalarSemigroup`
- `scalarGenerator`
- `scalarGenerator_selfAdjoint`
- `scalarSemigroup_hasDerivAt_zero`

This discharges strong continuity and generator extraction in a concrete lane.

### 4) Transfer intertwining + completion extension hook

- `osScalarTransfer`
- `osScalarTransferOnCompletion`
- `osScalarTransferOnCompletion_extends`
- `scalar_transfer_intertwines_hilbert`

This provides a concrete completion-extension path and a quotient-to-Hilbert
intertwining theorem for the scalar transfer family.

### 5) Positive spectral floor from Wilson domain assumptions

- `osGeneratorAt`
- `osSemigroupAt`
- `osGeneratorAt_selfAdjoint`
- `osGeneratorAt_gap_positive_of_domain`
- `osGenerator_uniform_gap_floor_of_domain`
- `grand321_assembly_of_domain`

This threads the existing hard-mode nonvanishing gap floor into the generator
lane as a strict positive lower bound.

## Honest boundary

`osScalarTransferOnCompletion_extends` needs the explicit
`UniformContinuous (osScalarTransfer ...)` witness at fixed `t`.
This is now isolated as a single analysis obligation rather than spread through
the chain.

## Build sanity

- `lake build Gutoe.YangMillsOSCompletion` ✅
- `lake build Gutoe` ✅

No `sorry` introduced.
