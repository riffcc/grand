# Finding 264 — Riemann Phase-3 All-Lanes Report (Train/Hold/Freeze)

Date: 2026-03-01  
Scope: Execute the full phase-3 RH exploratory suite in one reproducible lane.

## Added

- Shared 120-zero reference dataset:
  - `crates/gutoe-physics/data/zeta_zeros_first_120.txt`
- `riemann_lane` now loads reference zeros from that data file (via `include_str!` + `OnceLock`):
  - `crates/gutoe-physics/src/riemann_lane.rs`
- New full-lane report binary:
  - `crates/gutoe-physics/src/bin/riemann_phase3_all_lanes_report.rs`

## Lane Coverage

The report implements all requested lanes:

1. Determinant lock proxy
2. de Branges proxy
3. Weil positivity proxy
4. Prime trace probe
5. Inverse spectral reconstruction
6. Functional-equation symmetry embedding
7. Hilbert–Pólya bridge status lane
8. Train/holdout/freeze protocol (`1..40`, `41..80`, `81..120`)

## Run

```bash
cargo run -q -p gutoe-physics --bin riemann_phase3_all_lanes_report
```

Artifacts:

- `/tmp/bh_renders/riemann_phase3_all_lanes_report.txt`
- `/tmp/bh_renders/riemann_phase3_all_lanes_report.json`

## Key Results

Baseline structural map (n=512):
- Holdout MAPE: `5.8191e-2` (5.82%)
- Freeze MAPE: `1.2953e-1` (12.95%)

Lane 1 — Determinant-lock proxy:
- Hold mean nearest-root distance: `3.430e-1`
- Freeze mean nearest-root distance: `1.4527e1`

Lane 2 — de Branges proxy:
- Minimum separation: `9.943e-1`
- HB pass fraction: `0.0`

Lane 3 — Weil positivity proxy:
- Target min eig: `-7.12e-16`
- Pred min eig: `-7.07e-16`
- Practically matched at numerical epsilon.

Lane 4 — Prime trace probe:
- Target mean prime-peak z: `-1.652`
- Pred mean prime-peak z: `+0.084`

Lane 5 — Inverse spectral reconstruction:
- Fitted affine map (train):
  - slope `≈ 6.116706655796e-1`
  - shift `≈ 3.124709244024e2`
- Holdout MAPE: `5.9109e-2`
- Freeze MAPE: `1.2995e-1`
- Very close to structural-map baseline; no out-of-sample gain yet.

Lane 6 — Symmetry embedding:
- Pair-symmetry residual: `7.79e-12` (excellent spectral ± symmetry)
- Holdout MAPE: `2.050e-1`
- Freeze MAPE: `4.832e-1`
- Symmetry is enforced but predictive fit degrades materially.

Lane 7 — Hilbert–Pólya bridge:
- Status: `exploratory_report_lane`
- Numerical/operator prerequisites are tracked; formal proof lane remains open.

## Interpretation

- The “do all lanes” phase is now real, runnable, and logged.
- Baseline operator + structural map remains strongest predictive lane in this sweep.
- Symmetry embedding gives clean functional-equation-style spectral symmetry but currently at a large fit cost.
- The freeze slice (`81..120`) is now wired as a true out-of-sample discipline.

