# GRAND — Triangulation Solver Promotion (Shared API + CI)

Date: 2026-02-28

## Scope

Moved multiplicative triangulation from a standalone report into shared `gutoe-em` APIs and wired a dedicated CI gate.

This turns the residual analysis into reusable infrastructure for downstream gates and closure work.

## Files

- `crates/gutoe-em/src/flavor.rs`
  - Added:
    - `NeutrinoTriangulatedSolution`
    - `EwShiftTriangulatedSolution`
    - `neutrino_splitting_ratio_from_exponent`
    - `solve_neutrino_exponent_for_ratio`
    - `triangulate_neutrino_from_splittings`
    - `triangulate_ew_shift_for_target`
  - Added tests for machine-precision ratio closure and absolute-splitting reconstruction.
- `crates/gutoe-em/src/lib.rs`
  - Re-exported triangulation types/functions.
- `crates/gutoe-physics/src/bin/triangulate_params.rs`
  - Refactored to consume shared `gutoe-em` triangulation APIs.
- `crates/gutoe-physics/src/bin/triangulation_ci_gate.rs`
  - New CI gate to lock solver integrity.
- `crates/gutoe-physics/src/bin/global_gate_report.rs`
  - Added invocation of `triangulation_ci_gate`.

## Validation

Commands run:

- `cargo test -q -p gutoe-em flavor`
- `cargo run -q -p gutoe-physics --bin triangulate_params`
- `cargo run -q -p gutoe-physics --bin triangulation_ci_gate`
- `cd lean && lake build Gutoe`

All passed.

## Key Numbers (current lane)

From `/tmp/bh_renders/triangulate_params_report.txt`:

- `p_ratio = 13.688110433760`
- `kappa_geo = 34.697396055505`
- `ew_coeff_required = 8.460487692308`
- `ratio_fit_rel_err ≈ 1.77e-11`

Structural references:

- `p_structural = 13.7`
- `kappa_structural = 60/11 = 5.454545...`
- `ew_coeff_structural = 8.0`

Interpretation:

- Ratio pattern closure is effectively exact under triangulation.
- The dominant unresolved gap remains absolute neutrino normalization (`kappa` uplift vs structural baseline).
- EW bridge requires a modest coefficient uplift over `d/2`.

## Honesty

This is a forcing/diagnostic lane, not a final zero-parameter closure claim.

It quantifies the exact residual targets that the next structural derivation must explain.
