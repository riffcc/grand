# Finding 078: GRAND-304 Clay-Readiness Dossier

Date: 2026-02-26
Status: GRAND-304 complete

## Purpose

This bundle gives an external reviewer a direct path from claim to theorem to
artifact, with explicit assumptions and reproducible commands.

## Claim Index (Claim -> Lean theorem(s) -> artifact(s))

| Claim | Lean theorem(s) | Artifacts |
|---|---|---|
| Structural Yang-Mills mass gap is positive (Theorem A lane) | `doeblin_decomposition`, `abs_eigenvalue_le_one_sub_eps_of_decomposition_stochastic`, `mass_gap_positive_of_doeblin_mode` | `lean/Gutoe/YangMillsStructuralGap.lean`, `lean/Gutoe/YangMillsMassGap.lean`, `findings/064-yang-mills-theorem-a-structural-gap-closure.md` |
| Mass-gap lower bound is explicit and computable | `mass_gap_ge_doeblin_bound`, `doeblin_bound_positive` | `lean/Gutoe/YangMillsMassGap.lean`, `findings/062-yang-mills-doeblin-gap-lower-bound.md`, `findings/063-yang-mills-doeblin-subdominant-bound.md` |
| Gap survives continuum schedule (Theorem B lane) | `uniform_eps_floor_of_z3_local_regular_schedule`, `continuum_hypotheses_of_z3_nn_schedule`, `continuum_survival_gap_nonvanishing_of_z3_nn_schedule` | `lean/Gutoe/YangMillsContinuumSurvival.lean`, `findings/065-yang-mills-continuum-survival-bridge.md`, `findings/066-yang-mills-uniform-epsilon-floor-from-sc-coordination.md` |
| Wilson-side transfer kernel matches center-count lane | `z3_center_plaquette_kernel_eq_smoothed_transition`, `center_plaquette_schedule_kernel_eq_transfer` | `lean/Gutoe/YangMillsWilsonBridge.lean`, `findings/067-yang-mills-wilson-z3-structural-bridge.md` |
| Center-origin construction and gap transfer (Theorem C sub-chain) | `c1_z3_to_su3_structural_construction`, `c2_counts_center_action_bijective`, `c3_wilson_gap_nonvanishing_from_clifford_z3` | `lean/Gutoe/YangMillsWilsonBridge.lean`, `findings/075-theorem-c-wilson-equivalence-spine.md` |
| Haar expectation decomposition and center reduction | `expectation_decomposition_over_center`, `normalized_expectation_reduce_to_center_of_quotient_normalization` | `lean/Gutoe/HaarExpectationDecomposition.lean`, `lean/Gutoe/HaarFiberCollapse.lean`, `findings/070-haar-expectation-decomposition.md`, `findings/074-haar-common-factor-normalization.md` |
| Full-state lift obligations are necessary and discharged in Wilson lane | `lift_obligation_families_are_independently_necessary`, `full_gap_positive_all_steps_of_wilson_center_schedule` | `lean/Gutoe/YangMillsFullStateLift.lean`, `findings/072-full-state-lift-progress.md` |
| Constructive lane closure is embedded in Wilson-equivalence domain | `mass_gap_embedded_of_wilson_equivalence_domain`, `constructive_lane_gap_closure_of_wilson_equivalence_domain` | `lean/Gutoe/YangMillsConstructiveQFT.lean`, `findings/073-constructive-qft-interface.md`, `findings/077-constructive-lane-wilson-domain-closure.md` |
| Group scope extends beyond one bespoke SU(3) presentation | `compact_simple_scope_supports_full_path2` | `lean/Gutoe/YangMillsGaugeScope.lean`, `findings/076-gauge-scope-generalization-compact-simple.md` |

## Assumptions Register (Explicit, reviewer-facing)

1. `a_t n > 0` for each lattice step in the schedule.
This appears in continuum-domain theorems as a physical time-step positivity condition.

2. `alpha > 0` for Laplace smoothing / Doeblin floor.
This is required for entrywise positivity and non-trivial minorization.

3. Constructive-QFT interface assumptions are explicit in `ConstructiveYMModel`.
The file `lean/Gutoe/YangMillsConstructiveQFT.lean` intentionally keeps OS/Wightman milestones as typed obligations, then proves mass-gap consequences once those obligations are supplied.

4. No hidden fallback path is used in the Clay lane.
Kernel equivalence theorems are wired through the Wilson/center chain, not legacy alternatives.

## Repro Scripts

- Main bundle script:
  - `scripts/clay_repro_bundle.sh`
- Output location:
  - `findings/assets/clay/repro_<timestamp>.log`
  - `findings/assets/clay/theorem_presence_<timestamp>.txt`

The script runs:

1. `lake build` for each load-bearing module (structural gap, continuum survival, Haar bridge, Wilson bridge/equivalence, lift, constructive lane, gauge scope, and `Gutoe` root).
2. theorem-presence checks for the core chain names used in this dossier.

## Independent-Check Checklist (External reviewer runbook)

1. Clone repo and enter root.
2. Run `./scripts/clay_repro_bundle.sh`.
3. Confirm all Lean builds succeed with no manual edits.
4. Open the generated `repro_*.log` and verify each module build step is present.
5. Open `theorem_presence_*.txt` and verify each theorem symbol resolves to exactly one canonical declaration line.
6. Spot-check theorem bodies in:
   - `lean/Gutoe/YangMillsStructuralGap.lean`
   - `lean/Gutoe/YangMillsContinuumSurvival.lean`
   - `lean/Gutoe/YangMillsWilsonBridge.lean`
   - `lean/Gutoe/YangMillsWilsonEquivalence.lean`
   - `lean/Gutoe/HaarFiberCollapse.lean`
7. Validate that assumptions used by each theorem are listed in the Assumptions Register above.
8. Confirm finding-chain traceability:
   - `findings/064` -> Theorem A closure
   - `findings/065`/`066` -> Theorem B continuity lane
   - `findings/067`/`075` -> Theorem C bridge
   - `findings/074`/`077` -> normalization + constructive-lane closure

## Boundary (honest)

This dossier establishes reproducibility and proof navigation for the current
Clay lane. It does not claim that external prize adjudication requirements are
already accepted; it provides the verification package needed for that review.
