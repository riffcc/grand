# Finding 117: GRAND-330 Clay Submission Package

Date: 2026-02-28
Status: GRAND-330 complete

## Scope

Produce the reviewer-facing package requested by GRAND-330:

1. cover letter,
2. theorem map,
3. step-by-step reviewer guide,
4. reproducibility artifacts for the mapped Clay requirements.

## What landed

Submission package directory:

- `findings/assets/clay/submission/README.md`
- `findings/assets/clay/submission/cover_letter.md`
- `findings/assets/clay/submission/theorem_map.md`
- `findings/assets/clay/submission/reviewer_guide.md`

New reproducibility script:

- `scripts/clay_submission_bundle.sh`

It verifies the four requirement lanes through:

- module builds,
- theorem-presence checks,
- `no-sorry` checks in the mapped Lean files.

## Requirement map used

1. Sequence 323 (`GRAND-330`): `constructive_schwinger_family_exists`
2. Sequence 325 (`GRAND-331`): `grand331_end_to_end_os_reconstruction_of_domain`
3. Sequence 326 (`GRAND-332`): `theorem_c_wilson_equivalence_domain_limits`
4. Sequence 328 (`GRAND-333`): `osGenerator_uniform_gap_floor_of_domain`

## Verification command

Run:

```bash
bash scripts/clay_submission_bundle.sh
```

Outputs:

- `findings/assets/clay/submission/submission_repro_<timestamp>.log`
- `findings/assets/clay/submission/submission_theorem_presence_<timestamp>.txt`

## Boundary

This package is a reproducibility and reviewer-navigation artifact. It does not
claim adjudication outcome.
