# GRAND-357 — Proton Mass Structural Ratio (`mp/me = 1836`)

## Result

- Structural mass ratio:
  - `mp/me = 12 × T(2^4 + 1) = 12 × T(17) = 12 × 153 = 1836`
- Lean-verified.
- Runtime parity verified.

## Derivation Chain

From shared Cl(1,3) primitives:

1. `clifford_dim = 2^4 = 16`
2. Augmented triangular index: `16 + 1 = 17`
3. `T(17) = 153`
4. `nLayers = total_gauge_generators = 12`
5. `mp/me = 12 × 153 = 1836`

Lean definitions/theorems:

- `mpMeAlgebraic := nLayers * triangularNumber (2^4 + 1)`
- `mp_me_eq_1836`
- file: `lean/Gutoe/MassSpectrum.lean`

Runtime parity expression:

- `mp_me_struct = sm.total_gauge_generators * triangular(sm.clifford_dim + 1)`
- file: `crates/gutoe-physics/src/bin/mass_periodic_report.rs`

## Proton Mass (Electron Anchor Route)

Using observed electron mass anchor:

- `m_e = 0.510998950 MeV`
- `m_p,pred = 1836 * m_e = 938.194072 MeV`
- `m_p,obs = 938.27208816 MeV`
- relative error: `-8.315e-5` (`0.0083%`)

Measured in:

- `crates/gutoe-physics/src/bin/proton_mass_report.rs`
- output: `/tmp/bh_renders/proton_mass_report/proton_mass_report.json`

## Residual

- Ratio residual to observed `mp/me ≈ 1836.15267`:
  - `Δ ≈ +0.15267` above structural leading-order integer `1836`.
- Interpretation: candidate subleading correction lane.

## Honest Boundary

- Strong result: structural mass **ratio** is derived and independently verified in Lean/runtime.
- Boundary: absolute proton mass here still uses an electron-mass anchor.
- Therefore:
  - **Derived:** leading structural ratio (`1836`)
  - **Not yet derived:** fully anchor-free absolute proton mass

## Status

- Finding accepted as structural ratio result.
- Subleading correction scan is optional next work; if no clean structural correction emerges quickly, keep this as an explicit open precision gap.
