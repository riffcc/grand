# Finding 086: GRAND-347 Inflation Amplitude + Reheating + CMB Likelihood

Date: 2026-02-27  
Status: GRAND-347 upgraded beyond slow-roll-only closure

## Goal

Take inflation from a `(N, n_s, r)` skeleton to a fuller quantitative lane:

1. scalar amplitude `A_s`
2. reheating map (`N_reh`, `w_reh`, `T_reh`)
3. combined CMB proxy likelihood over `(n_s, r, A_s)`

without introducing free fit knobs.

## What landed

### Rust lane upgrades

Updated:

- `crates/gutoe-physics/src/inflation.rs`
- `crates/gutoe-physics/src/bin/inflation_report.rs`
- `crates/gutoe-physics/src/bin/inflation_ci_gate.rs`

New structural pieces:

- `inflation_hubble_ratio_structural()`
- `scalar_amplitude(...)`
- `reheating_w_structural()`
- `reheating_efolds_structural()`
- `rho_end_planck_units(...)`
- `rho_reheat_planck_units(...)`
- `reheating_temperature_gev(...)`
- `cmb_proxy_chi2(...)`
- `cmb_proxy_loglike(...)`

Current structural composition:

- `H/M_pl = α_LO² * (60/11) * (1 - λ_QG) * (3/6) * 1/sqrt(486)`
- `A_s = (H/M_pl)² / (8π²ε)`
- `w_reh = 5/16`, `N_reh = N/12 = 5`
- `ρ_reh = ρ_end * exp(-3(1+w_reh)N_reh)`

Inflation gate now checks:

- `N` window
- `n_s` window
- `r` upper bound
- `A_s` window
- reheating floor (`T_reh > 1 MeV`)
- CMB proxy likelihood bound (`χ² <= 9`)
- graceful exit condition

### Lean parity upgrades

Updated:

- `lean/Gutoe/Inflation.lean`

New parity theorems:

- `inflation_hubble_ratio_eq`
- `inflation_hubble_ratio_pos`
- `scalar_amplitude_pos`
- `reheating_w_eq`
- `reheating_efolds_eq`

No `sorry`.

## Quantitative result

From `/tmp/bh_renders/inflation_report.txt`:

- `N = 60`
- `n_s = 0.965417`
- `r = 0.003333`
- `A_s = 2.219e-9` (target `2.10e-9`)
- `χ²_proxy = 0.649`
- `T_reh = 2.357e13 GeV`
- `passes_all = true`

## Why this matters

This moves GRAND-347 from a shape-only inflation lane to one that includes:

- perturbation amplitude scale (`A_s`)
- thermodynamic handoff (`T_reh`)
- one-number combined CMB consistency score (`χ²_proxy`)

It is still a proxy lane, but it is now meaningfully closer to reviewer-grade
cosmology diagnostics.

## Honest boundary

Still open for hard closure:

1. Replace proxy χ² with direct Planck likelihood inputs (`A_s`, `n_s`, `r`, covariance).
2. Tie reheating to explicit Boltzmann evolution (not only closed-form structural map).
3. Add tensor tilt / running checks (`n_t`, `dn_s/dlnk`) if this lane is promoted.

## Build sanity

- `cargo check -p gutoe-physics --bin inflation_report --bin inflation_ci_gate` ✅
- `cargo test -p gutoe-physics inflation -- --nocapture` ✅
- `cargo run -q -p gutoe-physics --bin inflation_report` ✅
- `cargo run -q -p gutoe-physics --bin inflation_ci_gate` ✅
- `cd lean && lake build Gutoe.Inflation` ✅
- `cd lean && lake build Gutoe` ✅

No new `sorry`.
