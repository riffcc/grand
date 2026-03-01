# 208 — Origin closure probe wired with `even_subalgebra` structural mode

## Code change
- Updated: `crates/gutoe-physics/src/bin/ctc_origin_energy_closure_probe.rs`

New mode:
- `GUTOE_CTC_ORIGIN_MODE=even_subalgebra`
  - `branching = 3`
  - `merge_fraction = 3/16`
  - `eta = 2/3`
  - `infra_gain = 8/3`
  - so `b_eff = 3 * (3/16) * (2/3) * (8/3) = 1` (exact arithmetic in code path)

Also added:
- `mode` and `mode_note` to output payload/text for auditable runs.

## Lean parity
- Updated `lean/Gutoe/EvenSubalgebraSuppression.lean` with:
  - split-product theorem: `(2/3)*(8/3) = 16/9`
  - full closure theorem: `3*(3/16)*(2/3)*(8/3) = 1`
- `lake build Gutoe.EvenSubalgebraSuppression` passes.
- `lake build Gutoe` passes (warnings only).

## Probe outputs

`even_subalgebra` mode run:
- output dir: `/tmp/bh_renders/ctc_origin_energy_closure_probe_even_sub`
- `b_eff = 1.0`
- `bmin_for_finite_horizon = 1.0000000000000102`
- `finite_horizon_reaches_target = false`

`legacy` mode regression check:
- output dir: `/tmp/bh_renders/ctc_origin_energy_closure_probe_legacy_after_mode`
- unchanged baseline: `b_eff = 0.375`

## Interpretation
- The structural expression is now first-class executable and theorem-backed.
- At current finite-horizon defaults, exact unit gain is still just below the
  measured finite-depth threshold (`~1 + 1.02e-14`), matching prior criticality.
