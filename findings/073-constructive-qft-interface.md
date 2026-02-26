# Finding 073: GRAND-302 Constructive QFT Interface (Scaffold)

Date: 2026-02-26
Status: GRAND-302 in progress (interface + non-vanishing consequences landed)

## New module

- `lean/Gutoe/YangMillsConstructiveQFT.lean`

## What landed

### 1) Explicit axiom set and theorem targets

- `OSAxiomInterface`
- `ConstructiveMilestones`
- `ConstructiveYMModel`
- `constructiveTargetsSatisfied`
- `constructive_targets_unpacked`

This gives an explicit Lean object for OS/Wightman-compatible target conditions.

### 2) Lean-facing constructive existence lane interface

- `massGapEmbeddedInConstructiveLane`

This states that constructive targets and continuum-survival hypotheses live in one theorem lane.

### 3) Mass-gap statement embedded in same framework

- `mass_gap_embedded_of_continuum_survival`
- `embedded_mass_gap_extract`

This threads the existing non-vanishing continuum gap theorem directly into the constructive-lane package.

### 4) Non-trivial dynamical consequences from the embedding

- `embedded_gap_positive_each_step`
- `TendsToZeroSeq`
- `not_tends_to_zero_of_uniform_positive_floor`
- `embedded_gap_not_tends_to_zero`
- `constructive_lane_gap_closure`

This upgrades GRAND-302 beyond checklist scaffolding: once embedded, the gap
lane is provably positive at every refinement step and cannot collapse to zero
as a sequence.

## Build sanity

- `lake build Gutoe.YangMillsConstructiveQFT` ✅
- `lake build Gutoe` ✅

No `sorry` introduced.

## Remaining for GRAND-302 closure

- Replace milestone `Prop`s with concrete constructive objects/morphisms
  (Schwinger families, reconstruction maps, operator algebras).
- Prove a constructive compatibility theorem from those objects into the
  existing mass-gap lane without interface assumptions.
- Pin dossier-grade statement for constructive existence that can be referenced
  directly from GRAND-304.
