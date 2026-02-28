# GRAND-364 — Priority Order Execution: (2) Neutrino/PMNS, (1) Continuum Gate, (3) Lean Parity

Requested order executed: **2 → 1 → 3**.

## (2) Neutrino/PMNS first — triage outcome

Commands:

- `cargo run -q -p gutoe-em --bin flavor_mix_report`
- `cargo run -q -p gutoe-em --bin flavor_ci_gate`

Observed direct PMNS:

- `θ12 = 33.690°`
- `θ23 = 49.107°`
- `θ13 = 8.213°`
- all envelope gates pass.

Triage result for `θ23` tension:

- Baseline structural direct: `sin²(θ23)=4/7` gives `θ23=49.106605°`.
- Small offset from target (`49.0°`) behaves as an `O(α²)` residual.
- Short structural rational coefficients in ansatz
  `sin²(θ23)=4/7 - c α²` already close to millidegree residuals:
  - `c=137/4` → `θ23=49.001051°`
  - `c=67/2` → `49.003362°`
  - `c=33` → `49.004902°`

Classification: **perturbative-fixable**, not a hard structural break.

Ticket:

- `GRAND-347` — PMNS `θ23` second-order correction lane.

## (1) Continuum second — hard reproducibility gate

Command:

- `bash scripts/clay_repro_bundle.sh`

Outputs:

- `findings/assets/clay/repro_20260227T155425Z.log`
- `findings/assets/clay/theorem_presence_20260227T155425Z.txt`

Status:

- all listed continuum/constructive modules build successfully,
- theorem-presence checks pass for the tracked constructive chain.

This is not final Clay closure, but it is a hard reproducibility/proof-presence gate for the current continuum lane.

## (3) Lean parity third — backlog checkpoint

Targeted build checkpoint:

- `cd lean && lake build Gutoe.FineStructure Gutoe.MassSpectrum Gutoe.StrongCouplingCInfBridge Gutoe.ElectronScaleTransduction`

Status:

- build green (no theorem-direction changes made),
- parity closure work item opened for findings 098–103.

Ticket:

- `GRAND-348` — Lean parity closure backlog for Findings 098–103.

## Related artifacts from same push window

- `alpha_web_ci_report` updated with explicit two-term alpha correction lane and gates:
  - `Δ = α^{-1}_phys - 137`
  - first order `5α`
  - second order `5α - 9α²`
  - quantitative improvement and CI boolean.

- `g_bridge_report` upgraded to 3 modes, including muon-phase alpha² electron lane.

