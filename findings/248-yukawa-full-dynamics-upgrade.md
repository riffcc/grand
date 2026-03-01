# Finding 248 — Yukawa Tooling Upgrade to Full Dynamics

## What changed
Added a new runner:
- `crates/gutoe-em/src/bin/yukawa_full_dynamics_scan.rs`

This upgrades the earlier one-loop universal lane to a coupled dynamics lane with:
- two-loop QCD running for `alpha_s`,
- running `alpha_em`,
- flavor-sensitive quark mass flow with QCD + QED + Yukawa self term,
- threshold-aware active flavor counting,
- preserved `L_g / S_g` diagnostics and constrained-Z3 closure checks.

## Outputs
- `/tmp/bh_renders/yukawa_full_dynamics_scan.txt`
- `/tmp/bh_renders/yukawa_full_dynamics_scan.csv`
- `/tmp/bh_renders/yukawa_full_dynamics_scan.json`
- `/tmp/bh_renders/yukawa_full_dynamics_scan_plot.png`

## Key differences vs prior one-loop lane
The previous one-loop lane showed almost exact scale-invariance of shape diagnostics.
The full-dynamics lane now shows nontrivial scale dependence:

- `s^2(L_g)` evolves with scale:
  - at `mu = m_t`: `2.917522`
  - at `mu = 1e16 GeV`: `2.961018`

- `S_g` evolves with scale (not flat):
  - `S1`: `-0.3869 -> -0.3998`
  - `S2`: ` 1.2640 ->  1.2511`
  - `S3`: ` 1.9693 ->  2.0570`

- Constrained closure on `L_g`:
  - fixed `s^2=2` remains poor (`~36.9 -> 42.8` RMS rel)
  - fixed `s^2=3` is dramatically better (`0.445 -> 0.260` RMS rel), but not yet precision closure.

## Interpretation
This upgrade does what was requested: it breaks the universal-ratio artifact and restores genuine dynamical behavior in the Yukawa scan tooling.

It does **not** yet close the quark sector, but it provides the correct instrumentation layer for the next phase:
- direct lattice-side extraction of effective Yukawas (`dm_f/dv`) and
- RG transport to observable scales.

## Next actionable step
Add an in-sim Yukawa-response probe:
- perturb lattice order parameter `v` in controlled windows,
- measure `dm_f/dv` per flavor sector,
- feed those effective UV Yukawas into this full-dynamics transport lane.
