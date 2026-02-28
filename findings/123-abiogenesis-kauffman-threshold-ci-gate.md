# Abiogenesis Kauffman Threshold CI Gate

## Scope

Implemented a theorem-style abiogenesis closure lane that answers a binary gate:

- Does structurally derived prebiotic chemistry exceed the autocatalytic closure threshold?

The lane uses:

- Structural monomer count: `N = 20` (`16 + 4`)
- Conservative catalytic lower bound from derived kinetics
- Kauffman threshold: `N * p > 1`
- Uncertainty-aware robustness check using a 3σ lower bound

## Physics Outputs

From `/tmp/bh_renders/abiogenesis_ci_gate.json`:

- `p_min = 0.0678523696926082`
- `N * p = 1.357047393852164`
- `closure_excess = 0.3570473938521641`
- `N * p lower_3sigma = 1.3417679226399386`
- `robust_margin = 0.3417679226399386`
- `pved_delta_e_ev = 1.798586691734443e-17`
- `overall_pass = true`

Interpretation:

- The closure inequality is not marginal. It clears threshold by ~35.7%.
- The 3σ-lowered control still clears threshold with ~34.2% margin.

## Rust Lane Changes

- Added abiogenesis lane implementation:
  - `crates/gutoe-physics/src/abiogenesis.rs`
- Added report + CI binaries:
  - `crates/gutoe-physics/src/bin/abiogenesis_report.rs`
  - `crates/gutoe-physics/src/bin/abiogenesis_ci_gate.rs`
- Exported lane in crate API:
  - `crates/gutoe-physics/src/lib.rs`
- Wired gate into global CI spine:
  - `crates/gutoe-physics/src/bin/global_gate_report.rs`

Global gate now includes:

- Execution of `abiogenesis_ci_gate`
- Parsing of abiogenesis gate JSON outputs
- Inclusion in `overall_pass`
- Emission in both text and JSON report payloads

## Lean Parity

Added formal closure module:

- `lean/Gutoe/AbiogenesisThreshold.lean`

Key theorem:

- `abiogenesis_kauffman_closure_exceeds_threshold :
   abiogenesisClosureControlQ > kauffmanClosureThresholdQ`

Derived rational chain includes:

- `p_min = 660/9727`
- `N*p = 13200/9727`
- `margin = 3473/9727 > 1/4`

Registered in Lean build roots:

- `lean/lakefile.lean` includes `Gutoe.AbiogenesisThreshold`

## Verification

Executed:

- `cargo run -q -p gutoe-physics --bin abiogenesis_report`
- `cargo run -q -p gutoe-physics --bin abiogenesis_ci_gate`
- `cargo run -q -p gutoe-physics --bin global_gate_report`
- `cargo test -q -p gutoe-physics abiogenesis`
- `cd lean && lake build Gutoe`

Status:

- Abiogenesis lane and CI gate pass
- Global gate remains green with abiogenesis included
- Lean root build succeeds with new theorem module
