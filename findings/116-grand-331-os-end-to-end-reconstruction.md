# Finding 116: GRAND-331 Full OS Reconstruction End-to-End Theorem

Date: 2026-02-28
Status: GRAND-331 complete

## Scope

GRAND-331 required a full theorem-level bridge from continuum Schwinger functions
(GRAND-330 lane) to explicit OS/Hilbert/Hamiltonian objects (GRAND-321 lane),
without reverting to standalone existential interface assumptions.

## What landed

New Lean module:

- `lean/Gutoe/YangMillsOSEndToEnd.lean`

Integrated into root build:

- `lean/lakefile.lean` (`Gutoe.YangMillsOSEndToEnd` root)

### Core objects

- `OSEndToEndStepPackage`
  - per-step explicit kernel `K`
  - Schwinger-family/kernel identity
  - quotient nonemptiness
  - completion dense embedding
  - self-adjoint generator
  - strictly positive Hamiltonian
  - Wightman exponential threading

### Closure theorem

- `grand331_end_to_end_os_reconstruction_of_domain`

For `WilsonEquivalenceDomain a_t alpha`, this theorem yields:

1. explicit schedule `SF : ℕ → CorrelatorFamily`
2. `SF n = wilsonSchwingerFamily W alpha n`
3. normalization `SF n m (fun _ => 1) = 1`
4. nonempty `OSEndToEndStepPackage` at every refinement step
5. uniform positive Hamiltonian floor `∃ c > 0, ∀ n, c ≤ osHamiltonianAt ... n`

## Additional unblock applied

`lean/Gutoe/YangMillsContinuumLimitKolmogorov.lean` was failing root builds due
stale/incorrect Mathlib API usage. It was repaired to current APIs:

- replaced invalid `Kernel.ofFintype` usage with `Kernel.ofFunOfCountable` + PMF rows
- aligned trajectory index type with `Finset.Iic` (required by `trajMeasure`)
- fixed `Measure.dirac` / `IsProbabilityMeasure` usage
- repaired Markov-kernel inheritance through `Kernel.comap`
- restored GRAND-322 package theorem with normalized expectation + mass-gap floor

This removed the global Lean blocker so GRAND-331 can be verified in full root CI.

## Verification

Executed successfully:

- `cd lean && lake build Gutoe.YangMillsOSEndToEnd`
- `cd lean && lake build Gutoe.YangMillsContinuumLimitKolmogorov`
- `cd lean && lake build Gutoe`

All succeeded (warnings only; no proof errors, no `sorry` additions).
