# Finding 249 — Full-Dynamics `s²(L_g)` Planck Extension

## Scope
Extended the full-dynamics Yukawa UV scan to include Planck-near anchors:
- `mu = 1e18 GeV`
- `mu = 1e19 GeV`

Code update:
- `crates/gutoe-em/src/bin/yukawa_full_dynamics_scan.rs`
  - `scan_mu` now includes `[m_t, 1e4, 1e8, 1e12, 1e16, 1e18, 1e19]`
  - summary string now reports the final scan scale dynamically.

## Outputs
- `/tmp/bh_renders/yukawa_full_dynamics_scan.txt`
- `/tmp/bh_renders/yukawa_full_dynamics_scan.csv`
- `/tmp/bh_renders/yukawa_full_dynamics_scan.json`
- `/tmp/bh_renders/yukawa_full_dynamics_scan_s2_planck_plot.png`

## Result: `s²(L_g)` trajectory

From the extended run:

- `mu = m_t`: `s²(L_g) = 2.917521825`
- `mu = 1e4`: `2.925051167`
- `mu = 1e8`: `2.938855312`
- `mu = 1e12`: `2.950503142`
- `mu = 1e16`: `2.961017988`
- `mu = 1e18`: `2.965986213`
- `mu = 1e19`: `2.968411312`

and fixed-`s²=3` closure RMS improves monotonically:

- `0.445229 -> 0.422888 -> 0.364955 -> 0.315458 -> 0.260399 -> 0.225460 -> 0.212560`

## Asymptotic read

The slope per decade in `log10(mu)` decreases monotonically:

- `[m_t,1e4]`: `0.004271829`
- `[1e4,1e8]`: `0.003451036`
- `[1e8,1e12]`: `0.002911958`
- `[1e12,1e16]`: `0.002628711`
- `[1e16,1e18]`: `0.002484112`
- `[1e18,1e19]`: `0.002425099`

Interpretation at current depth:
- trend is still upward at `1e19`,
- growth rate is decelerating,
- no hard flattening by `1e19`,
- data is consistent with an approach toward `s²=3` but does not yet prove a strict fixed-point lock at this scan depth.
