# Finding 263 — Riemann Phase-2: Family Scan + Spacing Objective

Date: 2026-03-01  
Scope: Execute phase-2 RH exploratory lane:
1) operator-family expansion,  
2) spacing-stat objective,  
3) high-fidelity scan and ranking.

## What Changed

- Added shared lane module:
  - `crates/gutoe-physics/src/riemann_lane.rs`
- Exported it via:
  - `crates/gutoe-physics/src/lib.rs`
- Upgraded report binary to include spacing metrics/objective:
  - `crates/gutoe-physics/src/bin/riemann_straight_shot_report.rs`
- Added family scan binary:
  - `crates/gutoe-physics/src/bin/riemann_operator_family_scan.rs`

## Objective Function

For each candidate:

```text
objective_total =
  mape
  + 0.50 * spacing_ks
  + 0.25 * spacing_mape
  + 0.25 * spacing_var_abs_err_to_gue
```

Where:
- `mape` = ordinate MAPE vs first 80 zeta-zero targets
- `spacing_ks` = KS distance between normalized spacing samples (pred vs target)
- `spacing_mape` = spacing-shape MAPE (pred vs target)
- `spacing_var_abs_err_to_gue` = |Var(spacing_pred) − (3π/8 − 1)|

## Operator Family

Shared operator form:

```text
H_ii     = ln(i + 1 + 13/16) + potential_scale / (i + 1 + 13/16)
H_i,i+1  = hop1_scale * sqrt((i+1)(i+2))
H_i,i+2  = hop2_scale * sqrt((i+1)(i+3))
```

Mapped ordinates:

```text
γ_pred = shift + slope * λ_raw
slope = 11/18
shift = 13*24 + 8/17
```

## High-Fidelity Scan Result

Run:

```bash
cargo run -q -p gutoe-physics --bin riemann_operator_family_scan
```

Configuration (default):
- `dimension = 512`
- `k = 80`
- candidates = 18
- neighborhood: `hop1 ∈ {0.50, 0.52}`, `hop2 ∈ {0, -0.02}`, `pot ∈ {0, -0.10}`

Best candidate:
- `family = baseline`
- `hop1 = 0.50`
- `hop2 = 0.00`
- `pot = 0.00`
- `mape = 4.0166%`
- `spacing_ks = 0.1772`
- `objective_total = 0.2237`

Top-ranked table confirms baseline remains the optimum in the tested family neighborhood.

## Interpretation

- Phase-2 objective successfully added and operational.
- Expanded family did **not** beat baseline under spacing-aware scoring at `n=512`.
- This is useful: it narrows the next move to changing operator structure (not parameter jitter).

## Artifacts

- `/tmp/bh_renders/riemann_straight_shot_report.txt`
- `/tmp/bh_renders/riemann_straight_shot_report.json`
- `/tmp/bh_renders/riemann_operator_family_scan.txt`
- `/tmp/bh_renders/riemann_operator_family_scan.json`

