# Reviewer Guide — GRAND-330

Date: 2026-02-28

This runbook is for an independent reviewer validating the Clay-lane package
without local manual theorem patching.

## 1) Repository baseline

1. Clone and enter repository root.
2. Confirm tools:
   - Lean/Lake available
   - `rg` available

## 2) Build + theorem-presence bundle

Run:

```bash
bash scripts/clay_submission_bundle.sh
```

This executes:

1. Lean builds for:
   - `Gutoe.YangMillsContinuumLimit`
   - `Gutoe.YangMillsOSEndToEnd`
   - `Gutoe.YangMillsWilsonEquivalence`
   - `Gutoe.YangMillsOSCompletion`
   - `Gutoe` (root)
2. Symbol checks for canonical theorem names.
3. `no-sorry` checks in the mapped theorem files.

## 3) Inspect generated artifacts

Check the newest files in:

- `findings/assets/clay/submission/submission_repro_<timestamp>.log`
- `findings/assets/clay/submission/submission_theorem_presence_<timestamp>.txt`

Pass criteria:

1. All listed `lake build` steps succeed.
2. Every mapped theorem symbol resolves to one declaration line.
3. No `FOUND_SORRY` lines appear in the log.

## 4) Spot-check theorem bodies

Open and inspect:

1. `lean/Gutoe/YangMillsContinuumLimit.lean`
   - `constructive_schwinger_family_exists`
2. `lean/Gutoe/YangMillsOSEndToEnd.lean`
   - `grand331_end_to_end_os_reconstruction_of_domain`
3. `lean/Gutoe/YangMillsWilsonEquivalence.lean`
   - `theorem_c_wilson_equivalence_domain_limits`
4. `lean/Gutoe/YangMillsOSCompletion.lean`
   - `osGenerator_uniform_gap_floor_of_domain`

## 5) Boundary

This package is a reproducibility and traceability artifact for independent
review. It is not a substitute for external adjudication.
