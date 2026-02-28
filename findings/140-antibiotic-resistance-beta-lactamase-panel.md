# 140 — Antibiotic Resistance Rescue: Beta-lactamase Inhibitor Panel

## Scope
Built a reduced-order ranking lane for inhibitor-enzyme binding energetics:
- Inhibitors: clavulanic acid, sulbactam, tazobactam, avibactam, vaborbactam
- Enzymes: TEM-1, KPC, NDM-1
- Outputs: pairwise predicted potency (nM), occupancy at 1 uM, and per-enzyme ranking.

Implemented artifacts:
- `crates/gutoe-physics/src/antibiotic_resistance.rs`
- `crates/gutoe-physics/src/bin/antibiotic_resistance_report.rs`
- `crates/gutoe-physics/src/bin/antibiotic_resistance_ci_gate.rs`
- global gate integration in `crates/gutoe-physics/src/bin/global_gate_report.rs`

## Data anchoring
Anchors are ChEMBL snapshot priors (Ki/IC50-derived medians under assay-description filters) for the 5x3 matrix.
One NDM-1 sulbactam cell is explicitly marked imputed due missing filtered direct record.

## Core readout
From `cargo run -q -p gutoe-physics --bin antibiotic_resistance_report`:
- `pair_count = 15`
- `mean_abs_log10_error_pred_vs_anchor = 0.6196`
- `ndm_max_predicted_occupancy_at_1uM = 0.00448`

Best-by-enzyme (predicted):
- TEM-1: `avibactam` (74.99 nM predicted)
- KPC: `vaborbactam` (114.86 nM predicted)
- NDM-1: `avibactam` (222,092.69 nM predicted)

Interpretation:
- The lane recovers strong serine-class inhibition (TEM-1/KPC) and weak NDM-1 inhibition across this inhibitor set.
- NDM-1 remains a low-occupancy regime at 1 uM in this model.

## CI gate
From `cargo run -q -p gutoe-physics --bin antibiotic_resistance_ci_gate`:
- `overall_pass = true`
- Guardrails include pair count, mean log10 error, low NDM occupancy, and expected TEM/KPC winner constraints.

Global gate now includes this lane and remains green.
