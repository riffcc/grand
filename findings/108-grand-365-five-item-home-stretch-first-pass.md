# GRAND-365 — Five-Item Home Stretch (First Pass)

Requested sequence: neutrino hierarchy, continuum lane, Lean parity backlog, hadron extension run, third-order alpha precision.

## 1) Neutrino mass hierarchy prediction

Implemented:

- `crates/gutoe-em/src/flavor.rs`
  - `neutrino_texture_eigenvalues()`
  - `neutrino_hierarchy_prediction()`
- `crates/gutoe-em/src/bin/neutrino_hierarchy_report.rs`

Run:

- `cargo run -q -p gutoe-em --bin neutrino_hierarchy_report`

Artifact:

- `/tmp/bh_renders/neutrino_hierarchy_report/neutrino_hierarchy_report.txt`
- `/tmp/bh_renders/neutrino_hierarchy_report/neutrino_hierarchy_report.json`

Output:

- hierarchy prediction: **normal**
- texture eigenvalues (current lane): `m1=-6.299e-1, m2=7.699e-2, m3=7.162e-1`

Note: this is a lane-level hierarchy call from current texture normalization, not yet an absolute-neutrino-mass closure.

## 2) Continuum limit hard gate

Run:

- `bash scripts/clay_repro_bundle.sh`

Artifacts:

- `findings/assets/clay/repro_20260227T155425Z.log`
- `findings/assets/clay/theorem_presence_20260227T155425Z.txt`

Status:

- module build chain passes,
- theorem-presence checks pass for tracked constructive/continuum obligations.

## 3) Lean parity backlog checkpoint (098–107)

Build checkpoint run:

- `cd lean && lake build Gutoe.FineStructure Gutoe.MassSpectrum Gutoe.StrongCouplingCInfBridge Gutoe.ElectronScaleTransduction`

Status:

- green build for targeted parity modules.

Operational tracking tickets created earlier in this push:

- GRAND-348: Lean parity closure backlog for Findings 098–103
- GRAND-346: alpha two-term correction formalization

## 4) Hadron extension run (beyond proton)

Runs:

- `cargo run -q -p gutoe-physics --bin mass_periodic_report`
- `cargo run -q -p gutoe-physics --bin qcd_scale_report`

Artifacts:

- `/tmp/nuclear_chart/mass_periodic_report.json`
- `/tmp/bh_renders/qcd_scale_report/qcd_scale_report.{txt,json}`

Extracted neutron extension (from `mass_periodic_report.json`):

- `neutron_minus_proton_struct_mev = 1.000000000`
- `neutron_mass_mev_pred = 939.194072200`
- `neutron_mass_mev_obs = 939.565420520`

QCD/hadron context (from `qcd_scale_report`):

- threshold-matched 3-loop `lambda_nf3 = 0.327165 GeV`
- structural-corrected 3-loop `lambda_nf3 = 0.333761 GeV`
- `mp/lambda_nf3 (3-loop structural-corrected) = 2.8112`

## 5) Third-order alpha precision artifact

Generated third-order scan for:

- `Δ = α^{-1}_phys - 137`
- model: `Δ = 5α - 9α² + c3 α³`

Artifact:

- `/tmp/bh_renders/alpha_third_order_scan/alpha_third_order_scan.txt`

Results:

- exact fitted coefficient: `c3_exact = -21.659258031216655`
- best short structural-rational candidate from sweep:
  - `c3 = -153/7 = -21.8571428571`
  - absolute error in `Δ`: `7.689683798745017e-08`

## Net

All five requested items were executed to first-pass artifact stage in one session.
