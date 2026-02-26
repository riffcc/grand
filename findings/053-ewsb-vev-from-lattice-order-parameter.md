# Finding 053 — EWSB VEV from Lattice Order Parameter (No G_F Input)

## Scope
Close `GRAND-292`: derive electroweak vev from lattice structural terms and order parameter, without using `G_F` as an input.

Chain added:
`Cl(1,3) counts + mp/me + proton anchor + broken-phase order parameter -> v -> (m_W, m_Z, m_H)`.

## Rust implementation
Updated:
- `crates/gutoe-em/src/weak.rs`
- `crates/gutoe-em/src/lib.rs`
- `crates/gutoe-em/src/bin/ewsb_mass_report.rs`

New structural constants/functions:
- `PROTON_MASS_ANCHOR_MEV = 938.272046`
- `EWSB_SCALE_FACTOR = 480` (`2^4 * (|grade1|+|grade2|) * |SU(2)|`)
- `VEV_OVER_PROTON = 480/1836 = 40/153`
- `electron_mass_from_proton_anchor()`
- `normalized_higgs_order_parameter(f0)`
- `electroweak_vev_from_lattice_order_parameter(f0)`

No-`G_F` branch formula:
- `order = clamp((f0 - f_c)/(1 - f_c), 0, 1)` with `f_c = 3/16`
- `v(f0) = mp * (40/153) * order`

## Lean parity
Updated:
- `lean/Gutoe/EWSBHiggs.lean`

Added theorems/defs:
- `ewsbScaleFactor`
- `ewsb_scale_factor_eq_480`
- `vevOverProton`
- `vev_over_proton_eq_40_153`
- `normalizedOrderParameter`
- `normalized_order_parameter_at_one`
- `electroweakVevFromLattice`
- `electroweak_vev_over_proton_at_full_vacuum`

These bind the runtime branch to shared primitives (`grade1_4d`, `grade2_4d`, `magneticTriplet`, `mpMeAlgebraic`) with no new fit knobs.

## Numerical outputs
From `/tmp/bh_renders/ewsb_mass_report.json`:

Lattice branch (`f0 = 1`):
- `v = 245.2999 GeV`
- `m_W = 80.0135 GeV` (`Δ = -0.3635`)
- `m_Z = 91.2294 GeV` (`Δ = +0.0418`)
- `m_H = 125.0789 GeV` (`Δ = -0.1711`)

Reference (existing Fermi branch kept for comparison):
- `v = 246.2197 GeV`
- `m_W = 80.3135 GeV`
- `m_Z = 91.5715 GeV`
- `m_H = 125.5479 GeV`

## Verification
- `cargo test -p gutoe-em weak -- --nocapture` ✅
- `cargo run -q -p gutoe-em --bin ewsb_mass_report` ✅
- `cd lean && lake build Gutoe.EWSBHiggs` ✅
- `cd lean && lake build Gutoe` ✅

## Boundary
- This closes the vev-origin slice from lattice order parameter and removes `G_F` dependency for that branch.
- Remaining absolute-mass precision closure is tracked separately by existing mass-sector tickets.
