# 134 — MS Cyclosporine PK Bridge

## Scope
Bridge the mechanistic lane's effect-site target (`20 nM`) to measurable blood concentration via an explicit uncertainty transduction model.

This lane is a translational proxy and does **not** prescribe dose.

## Model
- Site target concentration: `20 nM`
- Molecular weight: `1202.61 g/mol`
- Blood/site gain factor: lognormal (`median = 8.0`, `GSD = 1.6`)
- Monte Carlo samples: `50,000`
- Seed: `1337`

## Results
Concentration quantiles:
- `p05 = 74.03 nM` (`89.03 ng/mL`)
- `p25 = 117.25 nM` (`141.00 ng/mL`)
- `p50 = 160.67 nM` (`193.23 ng/mL`)
- `p75 = 220.17 nM` (`264.78 ng/mL`)
- `p95 = 346.96 nM` (`417.26 ng/mL`)
- mean: `215.46 ng/mL`

## Artifact paths
- `/tmp/bh_renders/ms_cyclosporine_pk_bridge/ms_cyclosporine_pk_bridge_report.txt`
- `/tmp/bh_renders/ms_cyclosporine_pk_bridge/ms_cyclosporine_pk_bridge_report.json`
- `/tmp/bh_renders/ms_cyclosporine_pk_bridge/ms_cyclosporine_pk_bridge_quantiles.csv`
