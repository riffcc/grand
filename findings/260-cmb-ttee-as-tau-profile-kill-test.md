# Finding 260 — CMB TT+TE+EE `A_s`-`tau` Profile Kill Test (Current Scorer)

Date: 2026-03-01

## Scope
Built and executed a minimum kill-test harness that profiles `tau_reio` while re-fitting `A_s` at each point, using the current in-repo CLASS-to-Planck scorer over full TT/TE/EE spectra.

New runner:
- `crates/gutoe-physics/src/bin/cmb_ttee_as_tau_profile.rs`

Artifacts:
- `/tmp/bh_renders/cmb_ttee_as_tau_profile.csv`
- `/tmp/bh_renders/cmb_ttee_as_tau_profile.json`

## Method (explicit)
Grid scan:
- `tau_reio ∈ [0.040, 0.080]` (17 points)
- `A_s` scale factor `∈ [0.90, 1.10]` (17 points)
- fixed remaining cosmology from structural base inputs

Scoring:
- channel χ² from `compare_class_to_planck` per spectrum
- combined χ² = `χ²_TT + χ²_TE + χ²_EE`
- combined reduced χ² uses total points and `(N_total - 1)` dof

## Result
Best combined profile point:
- `tau_reio = 0.0500`
- `A_s = 2.108320e-9`
- `chi2_total = 7406.100`
- `reduced_total = 1.141507`

Structural tau check (`tau_struct = 0.067637719634`):
- profiled/interpolated combined `chi2_total = 7418.403`
- `Δχ² = +12.303` versus profile best

Channel contribution to `Δχ²` at structural `tau`:
- TT: `+1.094`
- TE: `+9.819`
- EE: `+1.390`

## Interpretation boundary (critical)
This is a valid kill test **within the current scorer**, but it is **not** the official Planck likelihood comparison.

Current scorer limitations:
- diagonal-per-point pull χ² (no full covariance matrix)
- no Planck nuisance/foreground calibration parameter block
- no official `-2 ln L` (clik/CamSpec/Plik) backend integration

Therefore this finding supports:
- structural `tau≈0.0676` is disfavored relative to the profile minimum in the present approximate scorer,

and does **not** support:
- any claim of “beating Planck/ΛCDM likelihood” on official apples-to-apples likelihood machinery.

## Status
Kill Test 4 (minimum viable: `tau` profile with `A_s` re-fit) is now implemented and executed for combined TT+TE+EE in current infrastructure.
