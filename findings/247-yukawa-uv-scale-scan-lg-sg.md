# Finding 247 — Yukawa UV Scale Scan with `L_g / S_g` Decomposition

## Goal
Implement the requested pipeline:
1. Run all quark masses to common `μ` using one-loop, threshold-matched QCD.
2. Decompose into generation modes:
   - `L_g = sqrt(m_up,g * m_down,g)`
   - `S_g = 0.5 * ln(m_up,g / m_down,g)`
3. Evaluate:
   - Z3 closure on `L_g`
   - functional form of `S_g` vs generation index.
4. Scan `μ = m_t, 1e4, 1e8, 1e12, 1e16 GeV`.

## Implementation
- New binary: `crates/gutoe-em/src/bin/yukawa_uv_scale_scan.rs`
- Artifacts:
  - `/tmp/bh_renders/yukawa_uv_scale_scan.txt`
  - `/tmp/bh_renders/yukawa_uv_scale_scan.csv`
  - `/tmp/bh_renders/yukawa_uv_scale_scan.json`
  - `/tmp/bh_renders/yukawa_uv_scale_scan_plot.png`

## Core numeric results
At `μ = m_t`:
- `L_g = [1.997742, 207.238607, 22535.574194] MeV`
- `S_g = [-0.385525, 1.264887, 2.036808]`
- `S` sign pattern: `[-, +, +]` (first-generation sign flip preserved)

`L_g` Z3 fit:
- Free `(M,s,δ)` fit: exact by construction, `s^2 = 2.956734`, RMS `~1.8e-15`
- Constrained `s^2 = 2` fit: RMS relative error `~42.19` (very poor)

`S_g` vs generation (linear fit, g=1..3):
- slope `1.211166819`
- intercept `-1.450276881`
- RMSE `0.207062591`
- `R^2 = 0.957999737`

## Scale scan behavior
Across the full scan (`m_t` to `1e16 GeV`):
- `S_g` values are effectively invariant in this one-loop setup.
- `L_g` rescales downward with `μ`, but shape parameters stay invariant:
  - free-fit `s^2(L_g)` constant at `2.956734...`
  - constrained `s^2=2` RMS remains `~42.19`.

Monotonic flags:
- `lg_fixed_rms_nonincreasing = false`
- `sg_linear_rmse_nonincreasing = true` (trivially stable because `S_g` invariant here)

## Interpretation
This one-loop QCD pipeline does **not** produce UV tightening of down-sector closure.

Important caveat: at one loop, quark-mass anomalous dimension is flavor-universal for fixed `n_f`. That makes many mass ratios nearly scale-invariant (after fixed threshold offsets), so this setup cannot generate the kind of scale-dependent compression needed to collapse the down-sector mismatch.

## Next technical step
To test the UV-native hypothesis nontrivially, upgrade the runner to:
- two-loop (or higher) QCD mass running,
- explicit matching conditions beyond simple one-loop continuity,
- consistent MSbar input set at one reference scale for all six quarks.
