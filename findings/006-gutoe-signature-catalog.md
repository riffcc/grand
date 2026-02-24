# 006 — GUTOE Signature Catalogue (GR Comparison Targets)

This is a practical checklist of where GUTOE differs from pure GR in the current code/proof stack.

## 1) UV lattice correction in dispersion / metric coefficients

- Shared constant: `lambda_qg = 1/12`
- Lean source: `lean/Gutoe/LambdaQG.lean`
- Runtime sources:
  - `crates/gutoe-physics/src/constants.rs`
  - `crates/gutoe-gpu/src/metric.rs`

Observable impact: high-curvature geodesic and transfer terms receive a small positive UV correction versus pure continuum GR.

## 2) Singularity-free core (lattice floor)

- GUTOE path includes a core floor instead of letting trajectories run to `r=0` singular behavior.
- Runtime rendering supports interior/core modes and toggles in bh_render/bh_viewer.

Observable impact: inside-horizon trajectories can terminate/interact with finite core structure rather than singular crash behavior.

## 3) GUTOE-vs-GR render toggle

- Runtime rendering supports GR comparison mode (`gr=true`) and GUTOE mode.
- Existing reference issue/output path: `GRAND-160` difference-map work.

Observable impact: direct image-space deltas in ring shape/brightness and interior behavior under matched camera setup.

## 4) Transfer-model signatures (thin vs RIAF, tau scaling)

- Transfer parity report path: `/tmp/bh_renders/transfer_parity_m87star.csv`
- Current parity report now includes per-backend mean-luminance and delta-from-baseline columns.

Observable impact: transfer activation produces measurable intensity deltas in both backends; this helps isolate transfer physics from base geodesic mismatch.

## 5) Theorem/runtime coefficient parity

- Harness binary: `crates/gutoe-physics/src/bin/theorem_parity.rs`
- Output: `/tmp/bh_renders/theorem_runtime_parity.csv`

Observable impact: catches drift between theorem-level constants and runtime coefficients before visual regressions are interpreted as physics.

## Suggested reporting set for papers/demo notes

1. Always publish GR and GUTOE matched-camera pairs.
2. Include transfer parity CSV + theorem parity CSV as provenance artifacts.
3. Quote which signatures are structural (Lean-proven) vs numerical (runtime/model choices).
