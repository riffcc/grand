# Finding 050 — GRAND-85 PMNS Derivation Closure

## Scope
Close GRAND-85 (PMNS matrix: lepton/neutrino mixing observables) with a traceable derivation chain:

`Cl(1,3) primitives -> grade/Z3 texture entries -> Hermitian lepton/neutrino textures -> diagonalization (U_l, U_nu) -> PMNS = U_l^† U_nu -> observables + J`.

## Acceptance checklist
- [x] Lepton/neutrino mass textures built from Cl(1,3) primitives only (no observable seeding)
- [x] PMNS computed as `U_l^† U_nu` from texture diagonalization
- [x] PMNS observables + J pass PDG envelope gate
- [x] Lean parity for PMNS structural values preserved
- [x] Findings artifact includes explicit deltas and chain trace

## Code paths
- Rust texture derivation + diagonalization:
  - `crates/gutoe-em/src/flavor.rs`
- CI falsification gate:
  - `crates/gutoe-em/src/bin/flavor_ci_gate.rs`
- Lean structural parity:
  - `lean/Gutoe/FlavorMixing.lean`

## Why PMNS angles are large (vs CKM)
The PMNS closure follows the same diagonalization machinery as CKM, but the structural inputs differ:
- Quark sector uses product suppressions (`1/(4*6)`, `1/(16*17)`) that force small mixing.
- Lepton sector uses ratio structures (`4/13`, `4/7`, `1/7`) that keep angles O(1) in sine space.

This separation is captured in Lean theorem:
- `quark_lepton_mixing_gap` in `lean/Gutoe/FlavorMixing.lean`

## PMNS explicit deltas (from `/tmp/bh_renders/flavor_mix_report.json`)
Target used by harness: `theta12=33.4°`, `theta23=49.0°`, `theta13=8.5°`, `delta=197.0°`, `J=-1.0e-2`

### Direct algebraic branch
- `theta12 = 33.690068°` (`+0.290068°`)
- `theta23 = 49.106605°` (`+0.106605°`)
- `theta13 = 8.213211°` (`-0.286789°`)
- `delta   = 198.434949°` (`+1.434949°`)
- `J       = -1.010759e-2` (`-1.075896e-4`)

### Texture-diagonalization branch (load-bearing for GRAND-85)
- `theta12 = 32.500499°` (`-0.899501°`)
- `theta23 = 50.274732°` (`+1.274732°`)
- `theta13 = 9.104071°` (`+0.604071°`)
- `delta   = 192.707949°` (`-4.292051°`)
- `J       = -7.559176e-3` (`+2.440824e-3`)

## Verification
- `cd lean && lake build Gutoe.FlavorMixing` ✅
- `cargo test -p gutoe-em pmns -- --nocapture` ✅
- `cargo run -q -p gutoe-em --bin flavor_ci_gate` ✅

Gate output snapshot:
- `pmns_direct`: pass
- `pmns_texture`: pass

## Residual risk
PMNS closure is now robust at the envelope level, but tighter phase/J constraints and convention handling still feed directly into GRAND-86 (CP-violation closure).
