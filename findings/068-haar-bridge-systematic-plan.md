# Finding 068: Haar Bridge Systematic Execution Plan

Date: 2026-02-26
Status: Active
Scope: Path-2 bridge from Cl(1,3) -> Z3 -> SU(3) to continuous Haar/coset expectation statements

## Why this split

The mathematical target is:

- decompose SU(3) integration into center-orbit and coset factors
- prove gauge-invariant observables are fiber-constant on SU(3)/Z3 fibers
- collapse normalized expectations to center-sector observables

In Lean, this is too large to land monolithically. We split into four dependency-ordered tickets.

## Tickets (dependency chain)

1. GRAND-308 (In Progress): Lean SU(3) + quotient scaffold
2. GRAND-309 (Todo): Haar measure existence/uniqueness hooks
3. GRAND-310 (Todo): SU(3) expectation decomposition over center/coset fibers
4. GRAND-311 (Todo): Gauge-invariant fiber collapse + parity checks

Dependency:

GRAND-308 -> GRAND-309 -> GRAND-310 -> GRAND-311

## Acceptance criteria by ticket

### GRAND-308 (structural scaffold)

- [ ] Introduce explicit SU(3) carrier/group object used by bridge layer.
- [ ] Define Z3 center embedding and prove normality/centrality statements.
- [ ] Define quotient object SU(3)/Z3 and projection maps.
- [ ] Prove descent lemma: functions constant on cosets factor through quotient.
- [ ] Keep theorem statements non-vacuous and trace bridge anchors to existing Cl(1,3)->Z3->SU(3) theorems.

### GRAND-309 (Haar hooks)

- [ ] Pin the Haar measure object on SU(3) (existence/uniqueness reference in Lean form).
- [ ] Pin quotient/coset measure object and the projection measurability assumptions.
- [ ] State exact normalization conventions used by later expectation theorems.
- [ ] If any heavy theorem is not fully formalized, isolate as explicit hypothesis with ticket reference.

### GRAND-310 (measure decomposition)

- [ ] State and prove a usable decomposition theorem for expectation/integration:
      `E_SU3[f] = E_center[ E_fiber[f] ]` (precise Lean form).
- [ ] Verify compatibility with normalized kernels already used in Yang-Mills bridge modules.
- [ ] Add a finite analog parity check that matches existing transfer-lane theorems.

### GRAND-311 (fiber collapse)

- [ ] Prove gauge-invariant observables are fiber-constant on cosets.
- [ ] Collapse the coset integral to a normalization factor that cancels in normalized expectations.
- [ ] Produce parity theorem tying continuous collapse to existing finite transfer collapse.
- [ ] Record explicit assumptions (if any) remaining after collapse.

## Pragmatic proof policy

- Use full formalization where mathlib machinery is already mature.
- Use explicit bridge hypotheses (not hidden axioms) where Lie/measure machinery is still heavy.
- Every hypothesis must be named, scoped, and attached to one of GRAND-308..311.
- No vacuous theorems (`True`, tautologies, or disguised restatements).

## Current code anchor

`lean/Gutoe/YangMillsWilsonBridge.lean` already contains finite transfer-lane analogs:

- row-scale gauge redundancy equivalence
- kernel fiber constancy
- finite-fiber expectation collapse

These are the discrete parity anchors for GRAND-310/311 and should remain source-of-truth comparators.
