# Finding 081: GRAND-320 Concrete OS/Wightman Objects (Textbook Completion)

Date: 2026-02-26
Status: GRAND-320 done

## Goal

Complete the constructive lane milestone by replacing witness-level OS/Wightman
claims with concrete object-level constructions in Lean on the Wilson hard-mode
lane.

## New module

- `lean/Gutoe/YangMillsOSTextbook.lean`
- added root in `lean/lakefile.lean`

## What landed (final)

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

### 4) Hamiltonian, spectral condition, and Wightman threading

- `osHamiltonianAt`
- `osSpectralCondition`
- `wightmanAt`
- `wightmanAt_nonneg`
- `wightmanAt_semigroup_time`

This makes the time-evolution and spectral positivity layer explicit at each
refinement step rather than proxy-only.

### 5) Full concrete textbook package at fixed step

- `OSTextbookPackageAt`
- `os_textbook_package_at_of_domain`
- `os_spectral_condition_of_domain`

This packages kernel, quotient carrier, positivity, semigroup law, Hamiltonian,
spectral nonnegativity, and Wightman semigroup behavior as concrete objects.

### 6) Hard-mode constructive closure with textbook object layer

- `hard_mode_os_reconstruction_from_textbook`
- `constructive_targets_and_textbook_objects_of_domain`

This threads textbook-object facts with existing hard-mode domain closure and
Wightman compatibility in one theorem output, including nonempty concrete
package witnesses at every refinement step.

## Why this matters

It replaces "OS reconstruction exists" as an opaque proposition with explicit
constructed objects and proved laws in the lane’s concrete finite setting.

## Scope boundary (tracked separately)

Infinite-dimensional OS completion and full self-adjoint reconstruction polish
are now tracked as follow-on work (GRAND-321). They are not blockers for
GRAND-320 acceptance, which is the concrete object-level closure in the current
hard-mode lane.

## Build sanity

- `lake build Gutoe.YangMillsOSTextbook` ✅
- `lake build Gutoe` ✅

No `sorry` introduced.
