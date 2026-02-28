# 141 — Phage Host-Matching Bypass Lane

## Scope
Built a computational phage-matching lane that scores phage tail-fiber binding to bacterial receptor profiles and ranks best matches per strain.

Implemented artifacts:
- `crates/gutoe-physics/src/phage_host_matching.rs`
- `crates/gutoe-physics/src/bin/phage_host_matching_report.rs`
- `crates/gutoe-physics/src/bin/phage_host_matching_ci_gate.rs`
- `crates/gutoe-physics/src/bin/phage_host_matching_ingest_report.rs`

## Model intent
This lane bypasses beta-lactamase resistance by construction:
- Inputs: receptor expression profile + phage tail-fiber interaction proxies
- Outputs: predicted KD, attachment probability, lysis potential score
- Resistance marker is carried for reporting but excluded from the binding path.

## Core readout
From `cargo run -q -p gutoe-physics --bin phage_host_matching_report`:
- `pair_count = 16` (4 strains x 4 phages)
- `mean_best_lysis_score = 0.6474`
- `resistance_independence_probe_abs_delta = 0.0`

Best matches:
- `kp_ndm1_clinical` -> `phi_kp_omp` (lysis score 0.8796)
- `kp_kpc_clinical` -> `phi_kp_omp` (lysis score 0.8796)
- `ec_tem1_clinical` -> `phi_lambda_lamb` (lysis score 0.8159)
- `pa_mdr_clinical` -> `phi_lps_broad` (lysis score 0.0147)

Interpretation:
- NDM-1 Klebsiella is strongly matchable in this lane via Omp-targeted phage.
- The resistance-independence probe confirms this path is decoupled from beta-lactamase class.

## CI gate
From `cargo run -q -p gutoe-physics --bin phage_host_matching_ci_gate`:
- `overall_pass = true`
- Guardrails include pair count, minimum mean best lysis, exact resistance-independence probe bound, and NDM top-match expectation.

Lane CI is green and ready for global-gate wiring.
