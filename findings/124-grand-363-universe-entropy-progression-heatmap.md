# GRAND-363 — Universe Entropy Progression Heatmap

Date: 2026-02-28
Status: Closed (Rust lane + Lean parity + CI gate integrated)

## Goal
Emit a heatmap over the full simulated universe history that visualizes entropy-production channels from bare-rock baseline through intelligence-era amplification, and gate the progression so silent disconnections cannot regress unnoticed.

## What shipped
- New Rust lane: `entropy_progression`
  - Uses `evaluate_universe_gate_with_depth` across full cosmic history (`z_max = 1e9`, `history_points = 768`).
  - Computes per-area and horizon-surface integrated entropy-production proxies.
  - Splits channels by stage:
    - bare rock
    - prebiotic chemistry
    - autocatalytic life
    - photosynthetic biosphere
    - multicellular ecosystem
    - technological intelligence
- New report binary:
  - `cargo run -q -p gutoe-physics --bin entropy_progression_report`
  - Writes:
    - `/tmp/bh_renders/entropy_progression/entropy_progression_report.txt`
    - `/tmp/bh_renders/entropy_progression/entropy_progression_report.json`
    - `/tmp/bh_renders/entropy_progression/entropy_progression_heatmap.png`
- New CI gate binary:
  - `cargo run -q -p gutoe-physics --bin entropy_progression_ci_gate`
  - Writes `/tmp/bh_renders/entropy_progression_ci_gate.json`
- Global gate integration:
  - `global_gate_report` now executes `entropy_progression_ci_gate` and includes its metrics in `overall_pass`.
- Lean parity:
  - `Gutoe.LifeProgressionEntropy` added and registered in `lean/lakefile.lean`.
  - Proves strict stage ordering and intelligence-step dominance using shared structural definitions.

## Key outputs (current run)
From `/tmp/bh_renders/entropy_progression/entropy_progression_report.json`:

- `universe_age_gyr = 13.626893663495203`
- `h0_km_s_mpc = 68.01633117530825`
- `hubble_radius_m = 1.3600599308517463e26`
- `hubble_surface_area_m2 = 2.3244807601613887e53`
- `final_total_per_area_w_m2_k = 0.7121672062350028`
- `final_total_universe_w_k = 1.6554189689111519e53`
- `local_maxima_count = 5`
- `local_minima_count = 4`
- `max_positive_step_age_gyr = 12.861350927038746`
- `max_positive_step_w_m2_k = 0.5286465226921973`

Gate verdicts:
- `monotone_stage_plateaus = true`
- `intelligence_step_dominant = true`
- `extrema_present = true`
- `passes_all = true`

## Stage activations and gains
- prebiotic chemistry: age `0.8516808539684502` Gyr, gain `0.08029197080291971`
- autocatalytic life: age `4.25840426984225` Gyr, gain `0.3570473938521641`
- photosynthetic biosphere: age `9.368489393652952` Gyr, gain `0.8450704225352113`
- multicellular ecosystem: age `11.071851101589852` Gyr, gain `0.45454545454545453`
- technological intelligence: age `12.775212809526753` Gyr, gain `5.454545454545454`

## Verification
- `cargo run -q -p gutoe-physics --bin entropy_progression_ci_gate` ✅
- `cargo run -q -p gutoe-physics --bin entropy_progression_report` ✅
- `cargo run -q -p gutoe-physics --bin global_gate_report` ✅ (`overall_pass=true`)
- `cd lean && lake build Gutoe.LifeProgressionEntropy` ✅
- `cd lean && lake build Gutoe` ✅

## Notes
- Heatmap rows are channel-separated and log-color-scaled, with vertical overlays at stage activation epochs.
- Plateau monotonicity is evaluated on effective multipliers (total/baseline) to avoid false negatives from late-epoch baseline decline.
