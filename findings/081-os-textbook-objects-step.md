# Finding 081: GRAND-320 Concrete OS/Wightman Objects (Textbook Step)

Date: 2026-02-26
Status: GRAND-320 in progress

## Goal

Move the constructive lane from structural witness-only statements toward
concrete functional-analytic objects (quotient carrier, inner product, transfer
semigroup) on the Wilson hard-mode lane.

## New module

- `lean/Gutoe/YangMillsOSTextbook.lean`
- added root in `lean/lakefile.lean`

## What landed

### 1) Concrete OS quotient carrier and inner product

- `EuclideanTestSpace := Fin 3 -> R`
- `OSRel`, `osSetoid`, `OSHilbertQuot`
- `osInnerRep`, `osInner`
- `osInner_self_nonneg`

This gives an explicit quotient-level inner-product candidate built from kernel
images, not a free witness placeholder.

### 2) Concrete transfer semigroup

- `osTransfer`
- `osSemigroup`
- `osSemigroup_add`

So the semigroup law is proved on the concrete quotient carrier.

### 3) Stepwise kernel regularity/positivity and textbook OS package

- `wilsonKernelAt`
- `osRowStochasticAt`
- `osKernelPositiveAt`
- `osReconstructionTextbookAt`
- `os_reconstruction_textbook_at`

This packages row-stochasticity, strict positivity, quotient nonemptiness,
inner nonnegativity, and semigroup composition per refinement step.

### 4) Hard-mode constructive closure with textbook object layer

- `hard_mode_os_reconstruction_from_textbook`
- `constructive_targets_and_textbook_objects_of_domain`

This threads textbook-object facts with existing hard-mode domain closure and
Wightman floor witness in one theorem output.

## Why this matters

It replaces "OS reconstruction exists" as an opaque proposition with explicit
constructed objects and proved laws in the lane’s concrete finite setting.

## Honest boundary

This is a major step, but not yet a full external textbook OS reconstruction in
the infinite-dimensional analytic form. Remaining work for GRAND-320 is the
full functional-analytic completion in that stronger sense.

## Build sanity

- `lake build Gutoe.YangMillsOSTextbook` ✅
- `lake build Gutoe` ✅

No `sorry` introduced.
