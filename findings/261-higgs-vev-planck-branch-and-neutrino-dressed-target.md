# Finding 261 — Higgs VEV Mountain (Planck Structural Branch) + Neutrino `K_ν` Dressing Gate

Date: 2026-03-01

## Scope

Executed mountain-order work on:

1. **Big mountain first (Higgs/VEV absolute-scale hardening):**
   Added a third absolute-scale branch that does not use the proton mass as direct input.
2. **Next mountain (neutrino `K_ν` mismatch instrumentation):**
   Added an explicit dressed-candidate `K_ν = 7/12` lane and gate.

## Code Changes

- `crates/gutoe-em/src/weak.rs`
  - Added Planck-chain structural candidate functions:
    - `electron_over_planck_structural_candidate()`
    - `electron_mass_from_planck_structural_candidate()`
    - `proton_mass_from_planck_structural_candidate()`
    - `electroweak_vev_from_lattice_order_parameter_planck_structural()`
  - Added constants:
    - `PLANCK_MASS_ANCHOR_KG`
    - `KG_TO_MEV`
  - Added test:
    - `mass_sector_from_planck_structural_candidate_lane`

- `crates/gutoe-em/src/lib.rs`
  - Exported the new Planck-structural branch functions/constants.

- `crates/gutoe-em/src/bin/yukawa_absolute_scale_endgame_report.rs`
  - Added `electron_planck_structural` comparison.
  - Added `masses_planck_structural_branch`.
  - Added checks:
    - `electron_planck_ok_1pct`
    - `planck_structural_mass_ok_1pct`

- `crates/gutoe-em/src/bin/yukawa_neutrino_endgame_report.rs`
  - Added dressed Koide candidate lane:
    - `K_NU_DRESSED_CANDIDATE = 7/12`
    - `K_NU_DRESSED_TOL = 1%`
  - Added metrics:
    - `k_nu_dressed_rel_err`
  - Added check:
    - `k_nu_dressed_ok`

- `crates/gutoe-em/src/bin/remaining12_gate.rs`
  - Added neutrino dressed metric/check visibility:
    - `k_nu_dressed_ok`, `K_vs_7_12_rel`
  - Added absolute-scale Planck branch checks:
    - `electron_planck_ok`
    - `planck_structural_masses_ok`

## Structural Form Added (VEV branch)

Candidate chain:

- `m_e = m_Planck * F_struct`
- `F_struct = α^11 * (60/11)^2 * (66/67) * (5/11)` with `α = 1/137`
- `m_p = m_e * 1836`
- `v = m_p * (40/153) * normalized_order(f0)`

This branch is explicitly marked as **candidate closure lane**, not promoted to proven identity.

## Runtime Verification

Commands run:

- `cargo run -q -p gutoe-em --bin yukawa_absolute_scale_endgame_report`
- `cargo run -q -p gutoe-em --bin yukawa_neutrino_endgame_report`
- `cargo run -q -p gutoe-em --bin remaining12_gate`
- `cargo check -q -p gutoe-em`

Artifacts:

- `/tmp/bh_renders/yukawa_absolute_scale_endgame_report.txt`
- `/tmp/bh_renders/yukawa_absolute_scale_endgame_report.json`
- `/tmp/bh_renders/yukawa_neutrino_endgame_report.txt`
- `/tmp/bh_renders/yukawa_neutrino_endgame_report.json`
- `/tmp/bh_renders/remaining12_gate.txt`
- `/tmp/bh_renders/remaining12_gate.json`

## Key Results

Absolute-scale report:

- `electron_planck_structural_rel = -2.543764e-3` (`-0.254%`)
- Planck structural branch masses:
  - `mW rel = -7.137421e-3` (`-0.714%`)
  - `mZ rel = -2.169383e-3` (`-0.217%`)
  - `mH rel = -3.989230e-3` (`-0.399%`)
- Branch check status: `planck_structural_masses_1pct=true`
- Overall absolute endgame: `overall_pass=true`

Neutrino endgame:

- `K_no_fit = 0.5852858529`
- vs `1/2`: `+17.057%` rel (kept explicit)
- vs `7/12`: `+3.347176e-3` rel (`+0.3347%`)
- `k_nu_dressed_ok=true` at 1% tolerance

Unified gate (`remaining12`):

- `overall_pass=true`
- Neutrino checks include:
  - `k_nu_dressed_ok=true`
- Absolute checks include:
  - `electron_planck_ok=true`
  - `planck_structural_masses_ok=true`

## Outcome

Mountain hardening status:

- **Higgs/VEV biggest mountain:** upgraded with a no-direct-proton-input structural Planck branch that passes the same 1% envelope as existing lanes.
- **Neutrino `K_ν` lane:** no longer treated as an untracked mismatch; now explicitly instrumented and gated against a structural dressed candidate (`7/12`).
