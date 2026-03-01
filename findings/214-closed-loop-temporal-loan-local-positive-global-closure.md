# 214 — Closed-loop temporal loan: local positive window with global closure

## New Lean module
- `lean/Gutoe/CTCClosedLoopTemporalLoan.lean`

## New Rust probe
- `crates/gutoe-physics/src/bin/ctc_closed_loop_temporal_loan_probe.rs`

## Proven in Lean
1. `closed_cycle_export_eq_drawdown_minus_loss`
   - In a closed packet cycle, export is exactly drawdown minus losses.
2. `closed_cycle_positive_export_iff_drawdown_gt_loss`
   - Positive export in a closed packet cycle is equivalent to drawdown exceeding losses.
3. `two_phase_temporal_loan_exists`
   - There exists a two-phase closed cycle with:
     - phase A positive local export and drawdown,
     - phase B negative export (repayment),
     - zero net export and restored door state.
4. `persistent_positive_export_blocked_under_nodraw`
   - Restates hard guard: no persistent per-step positive export under
     no-drawdown and nonnegative loss.

Interpretation:
- Closed-loop "works" as **temporal liquidity shifting** (local positive window),
  not as persistent net creation.

## Runtime check
- Output text:
  - `/tmp/bh_renders/ctc_closed_loop_temporal_loan_probe/ctc_closed_loop_temporal_loan_probe.txt`
- Output JSON:
  - `/tmp/bh_renders/ctc_closed_loop_temporal_loan_probe/ctc_closed_loop_temporal_loan_probe.json`

Default probe values:
- Phase A export: `+1 J` (closed packet, exact conservation residual `0`)
- Phase B export: `-1 J` (closed packet, exact conservation residual `0`)
- Door restored: `true`
- Net export over cycle: `0 J`
- Pattern flag: `closed_loop_temporal_loan_pattern = true`

## Build status
- `lake build Gutoe.CTCClosedLoopTemporalLoan` passed.
- `lake build Gutoe` passed (warnings only).
