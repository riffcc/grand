# Finding 051 — GRAND-86 CP-Violation Closure (CKM + PMNS)

## Scope
Close GRAND-86 by hardening CP-violation checks from structural phase sources through runtime observables.

Chain:
`Cl(1,3) + Z3 complex phases -> texture entries carry complex arguments -> diagonalization -> non-CP-conserving delta + nonzero Jarlskog`.

## Acceptance checklist
- [x] Complex phase sources are explicit in CKM/PMNS texture construction
- [x] CKM and PMNS branches both exhibit CP-violation witness (`|J| > threshold` and `delta` away from CP-conserving branches)
- [x] CI gate now enforces CP witness in addition to PDG envelope checks
- [x] Lean structural support remains green (`ckm_delta_in_open_pi`, `ckm_jarlskog_positive`)
- [x] Findings artifact records explicit witness values

## Code changes
- CP witness helpers and thresholds:
  - `crates/gutoe-em/src/flavor.rs`
    - `cp_violation_witness(...)`
    - `CKM_CP_J_MIN = 1e-6`
    - `PMNS_CP_J_MIN = 1e-3`
    - `CP_PHASE_TOL_DEG = 5.0`
- CP witness tests:
  - `crates/gutoe-em/src/flavor.rs`
    - `ckm_direct_has_cpv_witness`
    - `ckm_texture_has_cpv_witness`
    - `pmns_direct_has_cpv_witness`
    - `pmns_texture_has_cpv_witness`
- CI gate upgrade:
  - `crates/gutoe-em/src/bin/flavor_ci_gate.rs`
    - now enforces envelope + CP witness per branch

## Witness values (current run)
From `/tmp/bh_renders/flavor_mix_report.json` + `flavor_ci_gate` output:

### CKM
- direct: `delta=68.130°`, `J=+3.171628e-5`
- texture: `delta=64.839°`, `J=+2.832702e-5`

### PMNS
- direct: `delta=198.435°`, `J=-1.010759e-2`
- texture: `delta=192.708°`, `J=-7.559176e-3`

All branches pass CP witness thresholds and are outside CP-conserving phase neighborhoods (`0°, 180°` within `±5°`).

## Honesty note on derivation quality
The texture-diagonalization branch (the load-bearing derivation path) is less precise than the direct closed-form branch:
- CKM texture keeps angle precision tight but shows larger phase drift (`delta ~ -3.16°` vs target center).
- PMNS texture carries degree-scale angle/phase drifts (still envelope-safe).

This is expected from the extra approximation layer (matrix texture construction + numerical diagonalization) and identifies the next refinement hotspot: phase fidelity in texture construction, especially `delta`.

## Lean support
- `lean/Gutoe/FlavorMixing.lean`
  - `ckm_delta_in_open_pi`
  - `ckm_jarlskog_positive`
  - structural PMNS definitions remain unchanged and verified.

## Verification
- `cd lean && lake build Gutoe.FlavorMixing` ✅
- `cargo test -p gutoe-em flavor -- --nocapture` ✅ (12/12 flavor tests)
- `cargo run -q -p gutoe-em --bin flavor_ci_gate` ✅

## Residual risk
For PMNS, phase convention and global-fit uncertainty remain broader than CKM. The current closure is a robust CP witness gate, not yet a narrow-uncertainty precision claim.
