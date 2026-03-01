# 234 — Correction Hunts: α⁻¹, mp/me, mW/mZ (2026-03-01)

## Scope
Priority correction hunts requested after strict snapshot:
1. `alpha^-1`
2. `mp/me`
3. `mW/mZ`

Runner:
- `crates/gutoe-physics/src/bin/ctc_correction_hunts_alpha_mpme_mwmz.rs`

Artifacts:
- `/tmp/bh_renders/ctc_correction_hunts_alpha_mpme_mwmz/ctc_correction_hunts_alpha_mpme_mwmz.txt`
- `/tmp/bh_renders/ctc_correction_hunts_alpha_mpme_mwmz/ctc_correction_hunts_alpha_mpme_mwmz.json`

## Results

### 1) α⁻¹ correction hunt
Target offset from LO:
- `ALPHA_INV_PHYS - 137 = 0.03599917700001`

Candidates:
- `137 + 5α = 137.036486762822` (3.558 ppm)
- `137 + 5α - 9α² = 137.036007500632` (0.061 ppm)

Best current structural lane in this hunt:
- `137 + 5α - 9α²`

### 2) mp/me correction hunt
LO base:
- `6π^5 = 1836.118108711688`

Target offset from LO:
- `1836.15267343 - 6π^5 = 0.03456471831169`

Candidates:
- `6π^5 + (5α - 9α²) = 1836.154116212320` (0.786 ppm)
- integer scan in `6π^5 + (5α - c α²)`, `c in [0,128]`:
  - best `c = 36`
  - `1836.152678425750` (0.003 ppm)

Observation:
- `c=36` is structurally suggestive (`36 = 3*12`).

### 3) mW/mZ correction hunt
Experimental ratio:
- `mW/mZ = 0.881361062250`

Tree-level from weak-angle base (`sin^2 = 508/2197`):
- `sqrt(1 - 508/2197) = 0.876798496289` (5176.727 ppm)

Using running coefficient:
- `delta1 = α ln(10)/(4π) = 1.337122367993e-3`

Scheme-candidate lane:
- `sin^2_on-shell ~= sin^2_MSbar - 6*delta1`
- ratio prediction: `0.881361638207` (0.653 ppm)
- inferred needed shift ratio: `delta_scheme_needed / delta1 = 5.999240718`

Observation:
- Needed shift is numerically ~`6*delta1`.

## Status
- No new free-fit parameter introduced in these reported candidate lanes.
- These are correction-hunt candidates, not final derivations.
