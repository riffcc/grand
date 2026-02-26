# Finding 056 — FRW / H(z) Λ Phenomenology Harness (GRAND-294)

## Status
Implemented and exercised.

## Scope
Built a runtime harness that propagates the derived cosmological-term branches into
FRW expansion observables:

- `H(z)` and normalized `E(z) = H(z)/H0`
- inferred `H0` from `(Λ, Ω_Λ)` under flat closure assumptions
- inferred `Ω_Λ` from `(Λ, H0)` at reference `H0` values
- curvature residual `Ω_k = 1 - Ω_m - Ω_r - Ω_Λ`

## Implementation
- `crates/gutoe-physics/src/bin/frw_hz_report.rs`
- output artifacts:
  - `/tmp/bh_renders/frw_hz_report.txt`
  - `/tmp/bh_renders/frw_hz_report.json`

Assumptions (env-overridable):
- `Ω_m0 = 0.315`
- `Ω_r0 = 9.0e-5`
- `Ω_k0 = 0.0`
- flat target `Ω_Λ = 1 - Ω_m0 - Ω_r0 - Ω_k0 = 0.68491`

## Key outputs
Using the GRAND-295 full candidate `Λ_full = 1.105602561045e-52 m^-2`:

- `H0(full, flat Ω_Λ)` = `67.8568 km/s/Mpc`
- At `H0 = 67.4`, inferred `Ω_Λ(full)` = `0.694225`
- At `H0 = 73.0`, inferred `Ω_Λ(full)` = `0.591799`

For comparison:
- `Λ_observed` branch gives essentially the same `H0(flat)` (`67.8567`)
- pure structural `Λ_struct` branch gives `H0(flat) = 80.6127` (inconsistent without extra correction)

## Interpretation
- The corrected Λ branches (`signature`, `full`) produce an FRW-normalized expansion
  scale in the high-60s km/s/Mpc for standard matter/radiation assumptions.
- The harness cleanly exposes tension geometry through inferred `Ω_Λ` and `Ω_k` at
  different `H0` reference values.

## Notes
- This is a phenomenology harness, not yet a full derivation of matter-sector
  cosmological abundances from first principles.
- Next closure step is deriving `Ω_m(z)` inputs from GUTOE matter content rather than
  fixed external assumptions.
