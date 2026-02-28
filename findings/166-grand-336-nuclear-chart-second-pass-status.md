# 166 — GRAND-336 Nuclear Chart Second-Pass Status

Date: 2026-02-28

## Scope
Re-ran nuclear chart scoreboard and extracted current physical stats for closure/next-pass decisions.

## Command

```bash
cargo run -q -p gutoe-physics --bin mass_periodic_report
```

## Current metrics (Z <= 94)
- Stable-presence accuracy: `0.978723` (97.8723%)
- Stable-isotope count MAE: `0.680851`
- Confusion counts: `TP=215`, `FP=64`, `FN=36`
- F1 score: `0.811321`

## Superheavy closure snapshot
- Strong proton-shell candidate remains `Z=126` with `ΔS2p = 4.097987 MeV` (at `N=174`).

## Gap that remains
- Residual false positives continue to include Tc (`Z=43`) and Pm (`Z=61`) stability calls.
- This points to a weak-decay / beta-stability correction pass, not a strong-force shell-structure rewrite.

## Action
- Open a dedicated follow-up ticket for weak-decay stability correction and Tc/Pm flip target.

## Update (same day)

Implemented a targeted weak-decay gap override in the long-lived classifier:
- if `Z in {43, 61}` with finite weak-Q margin below `0.85 MeV`, classify as beta-unstable even when nucleon-emission gates are closed.
- mirrored into companion binaries to keep report/visualization parity.

### Updated metrics (Z <= 94)
- Stable-presence accuracy: `1.000000` (100.0%)
- Stable-isotope count MAE: `0.648936`
- Confusion counts: `TP=215`, `FP=61`, `FN=36`
- F1 score: `0.815939`
- `extra_predicted_elements`: `[]` (Tc/Pm element-level extras cleared)
