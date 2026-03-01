# Finding 255 — Down-Sector Casimir Discriminator via PDG Uncertainty Propagation

## Question
Is the `0.85%`-level mismatch between the down-sector Casimir candidate and extracted `s²_down` a real structural miss, or consistent with quark-mass uncertainty?

Target tested:
- `s²_target = 2 + 4/9`.

## Method
Added a reproducible runner:
- `crates/gutoe-em/src/bin/yukawa_down_uncertainty_propagation.rs`

Default inputs (broad PDG-style light-quark uncertainties discussed in-lane):
- `m_d = 4.67 ± 0.48 MeV`
- `m_s = 93.0 ± 11.0 MeV`
- `m_b = 4180 ± 30 MeV`

Procedure:
1. Compute center value of `s²_down = 6K(d,s,b)-2`.
2. Gaussian Monte Carlo propagation (`N=600,000`).
3. Uniform ±1σ box scan (corners + dense interior).

Artifacts:
- `/tmp/bh_renders/yukawa_down_uncertainty_propagation.txt`
- `/tmp/bh_renders/yukawa_down_uncertainty_propagation.json`

## Result

From the run:

- Center:
  - `s²_center = 2.390533910`
  - `s²_target = 2.444444444`
  - center offset: `-2.205%`

- Monte Carlo:
  - mean `2.393439929`
  - sd `0.056409882`
  - 1σ interval `[2.337361766, 2.449363393]`
  - 95% interval `[2.290014172, 2.511060157]`
  - target pull `z = 0.904`
  - `target_in_1sigma = true`
  - `target_in_95 = true`

- ±1σ box:
  - range `[2.322731462, 2.464346379]`
  - target inside box: `true`.

## Discriminator verdict

Under broad PDG-style quark-mass uncertainties, `2 + 4/9` is **consistent** with propagated `s²_down` at better than 1σ.

So for this discriminator:
- the Casimir candidate is **supported / not falsified** by current uncertainty-dominated inputs.
