# Finding 115: GRAND-349 Channel-Specific Li-7 Dark Coupling Closure

Date: 2026-02-28
Status: GRAND-349 complete (predictive Li-7 lane + Lean parity)

## Scope

GRAND-349 required replacing global Li-7 suppression with a channel-level derivation:

- derive absolute `Li7/H` from BBN reaction channels (not from `LI7H_OBSERVED` normalization),
- apply dark coupling only to the Be-7 mediated channel,
- keep D/H and He-4 lanes structurally unchanged,
- close Rust/Lean parity with no theorem statement drift.

Files touched:

- `crates/gutoe-physics/src/bbn.rs`
- `crates/gutoe-physics/src/bin/bbn_report.rs`
- `crates/gutoe-physics/src/bin/lithium7_report.rs`
- `lean/Gutoe/BBN.lean`

## Physical stats (runtime)

Source: `/tmp/bh_renders/bbn_report.json` and `/tmp/bh_renders/lithium7_report/lithium7_report.json`

Core BBN outputs at structural `eta10 = 6.312987028056`:

- `Yp = 2.455216450468e-1`
- `D/H = 2.347983877421e-5`
- `3He/H = 1.066946141622e-5`
- `Li7/H (predictive) = 1.817947130892e-10`

Residuals / gate metrics:

- `Yp delta = 5.216450467602e-4`
- `D/H rel error = 0.078137464695`
- `3He/H rel error = 0.030048962162`
- `Li tension ratio (predictive / observed) = 1.136216956807`
- windows used: `Yp<=0.010`, `D/H<=0.15`, `3He/H<=0.15`, `Li ratio in [0.8, 1.4]`
- gate status: `yp_ok=true`, `dh_ok=true`, `he3_ok=true`, `li_tension_ok=true`, `passes_primary=true`, `passes_all=true`

Li-7 channel decomposition:

- visible fraction: `11/16 = 0.687500000000`
- reaction-network gain: `33/16 = 2.062500000000`
- source before branch closure: `5.177919161251e-10`
- direct component (`1/16`): `3.236199475782e-11`
- Be-7 component raw (`15/16`): `4.854299213673e-10`
- Be-7 dark suppression factor: `165/536 = 0.307835820896`
- Be-7 component dark-coupled: `1.494327183313e-10`
- channel-coupled factor: `3011/8576 = 0.351096082090`

Observed-anchored diagnostic retained for comparison only:

- observed-anchored Li-7 ratio: `3.321150434700`
- predictive channel-derived ratio: `1.136216956807`

## Lean parity closure

`lean/Gutoe/BBN.lean` now mirrors the Rust mechanism exactly:

- `primordialLithium7RatioObservedAnchored` (diagnostic)
- `lithium7VisibleFraction = 11/16`
- `lithium7ReactionNetworkGain = 33/16`
- `lithium7ReactionNetworkSource`
- `lithium7DirectComponent`
- `lithium7Be7ComponentUnsuppressed`
- `lithium7Be7ComponentDarkCoupled`
- `primordialLithium7RatioChannelCoupled`
- predictive alias `primordialLithium7Ratio` routes to channel-coupled lane

Key closure theorems include:

- `lithium7_visible_fraction_eq`
- `lithium7_reaction_network_gain_eq`
- `lithium7_reaction_network_source_at_reference`
- `lithium7_channel_coupled_factorization`
- `lithium7_tension_ratio_at_reference`
- `lithium7_channel_coupled_tension_ratio_at_reference`
- `lithium7_channel_coupled_ratio_reference_window`

No new `sorry` in `Gutoe/BBN.lean`.

## Verification

Executed successfully:

- `cargo check -p gutoe-physics --bin bbn_report --bin lithium7_report --bin bbn_ci_gate`
- `cargo test -p gutoe-physics --lib bbn -- --nocapture`
- `cargo run -q -p gutoe-physics --bin bbn_report`
- `cargo run -q -p gutoe-physics --bin lithium7_report`
- `cargo run -q -p gutoe-physics --bin bbn_ci_gate`
- `cd lean && lake build Gutoe.BBN`
- `cd lean && lake build Gutoe`

Builds contain existing linter warnings in unrelated files, but no GRAND-349 blocking errors.
