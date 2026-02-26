# Finding 052 — EWSB/Higgs Origin and Mass-Sector Closure Slice

## Scope
Advance GRAND-80 and GRAND-131 from conceptual notes to an executable + formally verified derivation slice.

Chain implemented:
`Cl(1,3) shared counts -> Higgs quartic λ_H and critical fraction f_c -> quartic potential + broken-phase branch -> electroweak mass outputs (v, m_W, m_Z, m_H)`.

## Structural derivation (Lean)
New module:
- `lean/Gutoe/EWSBHiggs.lean`

Added to build roots:
- `lean/lakefile.lean` (`Gutoe.EWSBHiggs`)

Key theorems:
- `ew_grade_sum_eq_10`: `|grade1| + |grade2| = 10`
- `higgs_quartic_eq_13_100`: `λ_H = (16-3)/(10^2) = 13/100`
- `critical_void_fraction_eq_3_16`: `f_c = 3/16`
- `nontrivial_vev_sq_pos`: broken phase (`f₀ > f_c`) implies non-trivial branch `φ² = μ²/(2λ) > 0`
- `higgs_deriv_zero_at_nontrivial_stationary`: that branch is stationary for `V(φ;f₀) = -μ² φ² + λ φ⁴`
- `higgs_mass_over_vev_sq_eq_13_50`: `(m_H/v)^2 = 2λ_H = 13/50`
- `ewsb_structural_closure`: combined closure theorem for quartic, broken branch, mass-ratio coefficient

## Runtime parity (Rust)
Updated:
- `crates/gutoe-em/src/weak.rs`

New structural constants/functions:
- `HIGGS_QUARTIC_LAMBDA = 13/100`
- `HIGGS_CRITICAL_VOID_FRACTION = 3/16`
- `higgs_mu_sq`, `higgs_potential`, `higgs_potential_derivative`, `higgs_nontrivial_vev`
- `higgs_mass_from_vev`
- `electroweak_vev_from_fermi`
- `weak_coupling_from_alpha`
- `w_mass_from_vev_and_alpha`, `z_mass_from_vev_and_alpha`

Re-exported via:
- `crates/gutoe-em/src/lib.rs`

New report binary:
- `crates/gutoe-em/src/bin/ewsb_mass_report.rs`
- outputs:
  - `/tmp/bh_renders/ewsb_mass_report.txt`
  - `/tmp/bh_renders/ewsb_mass_report.json`

## Numerical closure slice (current run)
From `/tmp/bh_renders/ewsb_mass_report.json`:
- `λ_H = 0.130000`
- `f_c = 0.187500`
- `sin²θ_W = 3/13 = 0.230769`
- `m_W/m_Z = 0.877058`

Using `G_F = 1.1663787e-5` and `α(m_Z)^-1 = 127.95`:
- `v = 246.2197 GeV`
- `m_W = 80.3135 GeV` (`Δ = -0.0635` vs 80.377)
- `m_Z = 91.5715 GeV` (`Δ = +0.3839` vs 91.1876)
- `m_H = 125.5479 GeV` (`Δ = +0.2979` vs 125.25)

## Verification
- `cd lean && lake build Gutoe.EWSBHiggs` ✅
- `cd lean && lake build Gutoe` ✅
- `cargo test -p gutoe-em weak -- --nocapture` ✅
- `cargo run -q -p gutoe-em --bin ewsb_mass_report` ✅

## Boundary (what this does NOT claim yet)
- This closes the **origin/mechanism slice** (GRAND-80/131): structural quartic + phase split + stationary branch + mass-ratio coefficient.
- Full absolute-mass closure tickets remain separate:
  - `GRAND-78` (`m_Z` absolute)
  - `GRAND-79` (`m_H` absolute)

Those now have a stronger structural base and a tighter runtime calibration harness.
