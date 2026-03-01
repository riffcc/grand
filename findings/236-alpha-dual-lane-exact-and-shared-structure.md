# 236 — Alpha Dual-Lane Run: Exact Coefficients and Shared Higher-Order Structure

## Scope
Executed exactly the two requested lanes.

- **Lane A**: exact two-term coefficient extraction
  - `alpha^-1 = 137 + 5alpha - b alpha^2`
  - `mp/me = 6pi^5 + 5alpha - B alpha^2`
- **Lane B**: shared higher-order structure
  - `alpha^-1 = 137 + 5alpha - b alpha^2 + c alpha^3`
  - `mp/me    = 6pi^5 + 5alpha - B alpha^2 + C alpha^3`
  - with tie: `B = g b`, `C = g c`

Runner:
- `crates/gutoe-physics/src/bin/ctc_alpha_dual_lane_search.rs`

Artifacts:
- `/tmp/bh_renders/ctc_alpha_dual_lane_search/ctc_alpha_dual_lane_search.txt`
- `/tmp/bh_renders/ctc_alpha_dual_lane_search/ctc_alpha_dual_lane_search.json`

## Lane A (exact two-term)
Exact implied coefficients:
- `b_alpha_exact = 9.156308354932`
- `B_mp_exact    = 36.093814510993`
- `B/b ratio     = 3.941961444707`

Top compact approximants near `b_alpha_exact` from constrained structural-rational family:
- `9 + 5/32 = 9.15625` (abs err `5.835e-05`)
- nearby alternatives listed in artifact.

## Lane B (shared higher-order)
Exact tied scale from residual ratio:
- `g_exact = 3.941961444707`

Representative scenarios:
- `b=9, g=4`: solve `c` from alpha equation;
  - alpha lane closes exactly by construction,
  - mp lane residual is small but nonzero.
- `b=9, g=g_exact`: solve `c` from alpha equation;
  - both alpha and mp lanes close exactly by construction.
- `b=9+5/32` variants included similarly.

## Strict interpretation
- The non-rounded coefficients (`9.156...`, `36.093...`, ratio `3.94196...`) are now explicitly surfaced and treated as primary signal.
- Integer-near structure (`9`, `36`, `4`) is strongly suggested but not yet a uniqueness proof.
- The dual-lane tool isolates exactly where higher-order structure must enter (through `c` and/or non-integer effective scale `g`).
