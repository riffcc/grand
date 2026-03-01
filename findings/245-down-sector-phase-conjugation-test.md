# Finding 245 — Down-Sector Phase-Conjugation Test

## Question
Does down-sector closure improve if we use phase-conjugated Z3 harmonics instead of CKM amplitude-space rotation?

Tested transforms:
- `delta -> -delta`
- `delta -> pi - delta`
- both with free fitted `s` and fixed `s = sqrt(2)`

Closure metric:
- RMS log error against structural down targets:
  - `ms/md = 19`
  - `mb/ms = (8/3)*19*(67/66) = 51.434343...`

## Inputs
- `m_d = 4.67 MeV`
- `m_s = 93.0 MeV`
- `m_b = 4180.0 MeV`

Fitted down parameters from baseline extraction:
- `M = 25.485862847`
- `s = 1.546135152`
- `delta = 126.294242914 deg`

## Results
- Baseline observed down masses:
  - `ms/md = 19.914346895`
  - `mb/ms = 44.946236559`
  - `rms_log = 0.100972056`  (best)

- Free-s, `delta -> -delta`:
  - same closure values (permutation-equivalent)
  - `rms_log = 0.100972056`

- Free-s, `delta -> pi - delta`:
  - `ms/md = 9.125199379`
  - `mb/ms = 1.394889649`
  - `rms_log = 2.603062074`

- Fixed `s = sqrt(2)`, `delta -> -delta`:
  - `ms/md = 7.015799213`
  - `mb/ms = 31.092695822`
  - `rms_log = 0.789273988`

- Fixed `s = sqrt(2)`, `delta -> pi - delta`:
  - `ms/md = 14.949344680`
  - `mb/ms = 1.371723218`
  - `rms_log = 2.568325709`

## Verdict
Phase conjugation does **not** close the down sector. Baseline remains best, and all nontrivial conjugation variants are significantly worse.

Artifacts:
- `/tmp/bh_renders/yukawa_down_phase_conjugation_test.txt`
- `/tmp/bh_renders/yukawa_down_phase_conjugation_test.json`
