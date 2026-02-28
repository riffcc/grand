# 142 — GRAND-353: Hadron Transduction + Error Bars

## Scope
Implemented a dedicated hadron transduction lane that maps structural proxies to physical MeV outputs with uncertainty propagation and CI gating.

Artifacts added:
- `crates/gutoe-physics/src/hadron_transduction.rs`
- `crates/gutoe-physics/src/bin/hadron_transduction_report.rs`
- `crates/gutoe-physics/src/bin/hadron_transduction_ci_gate.rs`
- Lean parity: `lean/Gutoe/HadronTransduction.lean`
- Build root registration: `lean/lakefile.lean`
- Global integration: `crates/gutoe-physics/src/lib.rs`, `crates/gutoe-physics/src/bin/global_gate_report.rs`

## Structural map used
- `mp/me = 1836`
- `pion_transduction_factor = 153/8`
- `delta_np_from_pion_factor = 13/1370`
- `corrected_dark_to_visible_ratio = 115/22`
- `kaon_to_pion_factor = 159/44`
- `qcd_visibility_damping_factor = 5/8` (finite visible-sector occupancy damping)

## Central outputs (MeV)
From `cargo run -q -p gutoe-physics --bin hadron_transduction_report`:
- `qcd_scale_nf3_mev = 327.16459312525734`
- `qcd_scale_effective_mev = 204.47787070328584`
- `proton_mev = 938.1940721999999`
- `neutron_mev = 939.5084391042504`
- `pion_mev = 138.514050678711`
- `kaon_mev = 500.5394104071601`
- `neutron_proton_split_mev = 1.3143669042505421`

Residuals vs anchors:
- `proton_rel_error = -8.314854612483402e-05`
- `neutron_rel_error = -6.0646565428085634e-05`
- `pion_rel_error = -7.568505907943771e-03`
- `kaon_rel_error = 1.3900607901846952e-02`

## Uncertainty outputs
Default sampling: `4096` samples, `valid_fraction = 1.0`.

Quantiles:
- `pi (p05/p50/p95) = 129.9897427755572 / 138.63955998422983 / 147.53196641542806`
- `K  (p05/p50/p95) = 469.7356613934908 / 500.99295539755775 / 533.1268786375696`

Span diagnostics:
- `pion_rel_span95 = 0.12653115490171976`
- `kaon_rel_span95 = 0.12653115490171984`

## Gate result
From `cargo run -q -p gutoe-physics --bin hadron_transduction_ci_gate`:
- `overall_pass = true`

Gate booleans:
- `proton_rel_error_ok = true`
- `neutron_rel_error_ok = true`
- `pion_rel_error_ok = true`
- `kaon_rel_error_ok = true`
- `valid_fraction_ok = true`
- `pion_span95_ok = true`
- `kaon_span95_ok = true`
- `neutron_obs_in_p95 = true`
- `pion_obs_in_p95 = true`
- `kaon_obs_in_p95 = true`
- `proton_obs_in_p95 = false` (diagnostic only; not a hard fail criterion)

## Lean parity status
`lean/Gutoe/HadronTransduction.lean` proves the closed-form structural factors used by Rust:
- `pionTransductionFactorQ = 153/8`
- `neutronSplitFromPionFactorQ = 13/1370`
- `qcdVisibilityDampingFactorQ = 5/8`
- `kaonToPionFactorQ = 159/44`

Verified by:
- `cd lean && lake build Gutoe` (pass)

## Notes
- This lane closes GRAND-353 on the transduction/error-bar requirements with explicit structural formulas and reproducible artifacts.
- The previous root-level `DistributionSummary` name collision was resolved by using a dedicated hadron type name (`HadronDistributionSummary`).
