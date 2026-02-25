# Finding 040 — Strict Positive pp Weak Rate (SU(2) + Gamow)

Date: 2026-02-25  
Scope: GRAND-280

## What Landed

We formalized a strict-positive pp weak reaction-rate kernel in Lean, then mirrored it in Rust runtime code.

### Lean (`lean/Gutoe/StellarFusion.lean`)

Added:
- `weakFermiPrefactor (f0) = 1 / (2 f0^2)`
- bridge theorem from SU(2) mass relation:
  `weak_fermi_prefactor_from_su2_relation`
- positivity theorem:
  `weak_fermi_prefactor_positive`
- kernel:
  `ppWeakRateFromSU2 g f0 rho_p m_r E`

Main theorem:
- `pp_weak_rate_positive_from_su2_and_gamow`

This theorem returns both:
1. weak charged-current vertex existence (from the existing SU(2) chain), and
2. strict positivity of the pp weak-rate kernel under finite/positive physical inputs.

### Rust (`crates/gutoe-physics/src/equations.rs`)

Added parity functions:
- `weak_fermi_prefactor`
- `weak_prefactor_from_su2`
- `sommerfeld_parameter`
- `gamow_factor`
- `pp_weak_rate_from_su2_and_gamow`

Added tests:
- SU(2) prefactor identity and positivity
- Gamow factor bounded in (0,1)
- strict positive pp weak-rate kernel under physical inputs

## Verification

- `lake build Gutoe.StellarFusion` ✅
- `lake build Gutoe` ✅
- `cargo test -p gutoe-physics su2_fermi_prefactors_match_and_are_positive -- --nocapture` ✅
- `cargo test -p gutoe-physics gamow_factor_is_between_zero_and_one_for_positive_inputs -- --nocapture` ✅
- `cargo test -p gutoe-physics pp_weak_rate_kernel_is_strictly_positive_under_physical_inputs -- --nocapture` ✅

## Follow-up

- GRAND-282: Maxwell-Boltzmann averaged pp reaction-rate integral (from point-kernel positivity to thermal-rate envelopes).
