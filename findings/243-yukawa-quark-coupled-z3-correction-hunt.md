# 243 — Yukawa Quark Coupled-Z3 Correction Hunt

## Scope
Execute the requested quark correction hunt in four parts:
1. Separate up/down Z3 harmonic fits `(s_u, δ_u)` and `(s_d, δ_d)`.
2. Compare `δ_u - δ_d` against structural Cabibbo prediction.
3. Split mismatch into sector-only vs cross-sector ratio lanes.
4. Test QCD-like scalar correction `s -> s * (1 + α_s / π)`.

Binary:
- `crates/gutoe-em/src/bin/yukawa_quark_coupled_z3_hunt.rs`

Artifacts:
- `/tmp/bh_renders/yukawa_quark_coupled_z3_hunt.txt`
- `/tmp/bh_renders/yukawa_quark_coupled_z3_hunt.json`

## Core outputs

### Phase-lock lane
- Structural Cabibbo: `θ_C(struct) = 13.262676°`
- Best wrapped phase gap: `δ_u - δ_d = 10.555847°`
- Gap to Cabibbo: `2.706829°`

Best permutation pair found:
- up fit perm: `[1, 2, 0]`
- down fit perm: `[0, 2, 1]`

Extracted sector parameters:
- up: `s_u^2 = 3.094038785`
- down: `s_d^2 = 2.390533910`

Interpretation: there is a nontrivial phase relation close to (but not equal to) Cabibbo.

### Sector mismatch decomposition (vs current structural ratio targets)

Up-only ratios:
- `m_c/m_u`: `0.275%`
- `m_t/m_c`: `0.023%`

Down-only ratios:
- `m_s/m_d`: `4.812%`
- `m_b/m_s`: `12.614%`

Cross-only ratios:
- `m_u/m_d`: `1.713%`
- `m_c/m_s`: `14.370%`
- `m_t/m_b`: `1.101%`

Interpretation: up sector is already tight, while down and especially `m_c/m_s` are the dominant residual lanes.

### QCD-like scalar correction test
Scan:
- `s -> s * (1 + α_s / π)`, `α_s ∈ [0, 0.6]`
- objective: RMS log-mismatch vs structural ratio targets

Result:
- baseline at `α_s=0`: `rms_log=0.074496070`
- best scan point: `α_s=0`, same score

Interpretation: this scalar `s` rescaling does not improve closure in this lane.

## What this means
- The coupled-sector story is partially right: a meaningful up/down phase relation exists.
- The simple QCD scalar correction is insufficient.
- The primary correction burden is not uniform; it is concentrated in down-sector + cross-sector couplings.

## Immediate next lane
Promote from scalar `s` correction to a coupled phase-correction model:
- keep separate `(s_u, s_d)`,
- introduce structured phase bridge `δ_u - δ_d = f(θ_C, Clifford counts)` with no free ad hoc terms,
- target `m_c/m_s` and `m_b/m_s` first, then re-run the same hunt.
