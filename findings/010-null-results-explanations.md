# 010 — Null-Result Explanations (Proton Decay, EDM, Fifth Force, LV)

Status: baseline explanation map for GRAND-122.

## Why null results matter

A TOE must explain both positive detections and durable nulls.

## Null classes and GUTOE framing

1. Proton decay non-observation
- Experimental context: very long proton lifetime bounds.
- GUTOE framing: baryon-number-violating channels are not generic low-order lattice/Cl(1,3) operators in the current chain.
- Required follow-up: enumerate leading allowed effective operators and prove/order their suppression.

2. EDM nulls (electron/neutron)
- Experimental context: very strong EDM upper limits.
- GUTOE framing: CP phases in low-energy effective sector must remain sufficiently suppressed/cancelled.
- Required follow-up: explicit EDM operator estimate from GUTOE coefficients and comparison to bounds.

3. Fifth-force nulls
- Experimental context: no robust extra long-range force in tested windows.
- GUTOE framing: no unscreened additional light mediator is currently required by the proved gauge + metric chain.
- Required follow-up: bound emergent mediator mass/coupling space from model terms.

4. Lorentz-violation nulls
- Experimental context: tight limits on low-energy LV.
- GUTOE framing: Lorentz sector emerges from grade-2 bivector structure and must recover low-energy Lorentz invariance.
- References: `lean/Gutoe/LorentzInvariance.lean`, `lean/Gutoe/ContinuumLimit.lean`.
- Required follow-up: derive explicit LV residual scaling with energy/lattice scale and compare to constraints.

## Short risk register

- High risk: EDM and LV require explicit numeric bounds, not narrative arguments.
- Medium risk: proton-decay suppression requires operator-level closure.
- Medium risk: fifth-force constraints require full parameter exclusion plots.

## Deliverable quality bar

For each null class, produce: (a) operator or mechanism, (b) scaling law, (c) numerical bound, (d) pass/fail statement.
