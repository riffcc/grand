# Finding 049 — GRAND-84 CKM Derivation Closure

## Scope
Close GRAND-84 (CKM matrix: 3 angles + 1 CP phase) with an explicit traceable chain:

`Cl(1,3) primitives -> grade/Z3 texture entries -> Hermitian mass textures -> diagonalization (U_u, U_d) -> CKM = U_u^† U_d -> observables + J`.

## Acceptance checklist
- [x] Mass textures built from Cl(1,3) primitives only (no observable seeding)
- [x] CKM computed as `U_u^† U_d` from texture diagonalization
- [x] All four CKM observables + J within PDG envelope gate
- [x] Lean parity for structural texture coefficients
- [x] Findings artifact with explicit deltas

## Code paths
- Rust texture derivation + diagonalization:
  - `crates/gutoe-em/src/flavor.rs`
- CI falsification gate:
  - `crates/gutoe-em/src/bin/flavor_ci_gate.rs`
- Lean structural parity:
  - `lean/Gutoe/FlavorMixing.lean`

## CKM explicit deltas (from `/tmp/bh_renders/flavor_mix_report.json`)
Target used by harness: `theta12=13.0°`, `theta23=2.4°`, `theta13=0.2°`, `delta=68.0°`, `J=3.0e-5`

### Direct algebraic branch
- `theta12 = 13.262676°` (`+0.262676°`)
- `theta23 = 2.388015°` (`-0.011985°`)
- `theta13 = 0.210647°` (`+0.010647°`)
- `delta   = 68.130102°` (`+0.130102°`)
- `J       = 3.171628e-5` (`+1.716285e-6`)

### Texture-diagonalization branch (load-bearing for GRAND-84)
- `theta12 = 13.009816°` (`+0.009816°`)
- `theta23 = 2.380149°` (`-0.019851°`)
- `theta13 = 0.197027°` (`-0.002973°`)
- `delta   = 64.838985°` (`-3.161015°`)
- `J       = 2.832702e-5` (`-1.672976e-6`)

## Verification
- `cd lean && lake build Gutoe.FlavorMixing` ✅
- `cargo test -p gutoe-em ckm -- --nocapture` ✅
- `cargo run -q -p gutoe-em --bin flavor_ci_gate` ✅

Gate output snapshot:
- `ckm_direct`: pass
- `ckm_texture`: pass

## Residual risk
CKM closure is now robust for GRAND-84, but PMNS/CP-wide ticket closure still depends on maintaining envelope thresholds across future texture refactors (`GRAND-85`, `GRAND-86`, and trend persistence in `GRAND-290`).
