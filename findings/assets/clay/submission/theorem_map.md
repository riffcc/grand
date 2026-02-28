# Theorem Map — GRAND-330

Date: 2026-02-28

This map links the Clay-lane requirements used by GRAND-330 to concrete Lean
symbols, files, and verification commands.

## Requirement mapping

| Requirement | Plane issue | Lean theorem(s) | Primary file |
|---|---|---|---|
| Continuum Schwinger function construction | Sequence 323 (`GRAND-330`) | `constructive_schwinger_family_exists` | `lean/Gutoe/YangMillsContinuumLimit.lean` |
| End-to-end OS reconstruction in continuum lane | Sequence 325 (`GRAND-331`) | `grand331_end_to_end_os_reconstruction_of_domain` | `lean/Gutoe/YangMillsOSEndToEnd.lean` |
| Bridge completeness to standard Yang-Mills formulation | Sequence 326 (`GRAND-332`) | `theorem_c_wilson_equivalence_domain_limits` | `lean/Gutoe/YangMillsWilsonEquivalence.lean` |
| Positive mass gap in reconstructed continuum model | Sequence 328 (`GRAND-333`) | `osGenerator_uniform_gap_floor_of_domain` | `lean/Gutoe/YangMillsOSCompletion.lean` |

## Supporting theorem notes

- The reconstructed-model gap statement is also threaded in
  `grand331_end_to_end_os_reconstruction_of_domain` as the final conjunct:
  `∃ c > 0, ∀ n, c ≤ osHamiltonianAt ... n`.
- The continuity of the gap transfer from Wilson-equivalence assumptions is
  provided in `c3_gap_correspondence_of_domain` and consumed by the OS lane.

## Repro entrypoint

Run:

```bash
bash scripts/clay_submission_bundle.sh
```

Outputs:

- `findings/assets/clay/submission/submission_repro_<timestamp>.log`
- `findings/assets/clay/submission/submission_theorem_presence_<timestamp>.txt`
