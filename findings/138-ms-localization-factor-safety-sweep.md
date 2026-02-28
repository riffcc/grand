# 138 — MS Localization-Factor Safety Sweep (Concept Probe)

## Scope
Conceptual probe for CNS-localization intuition:
- Hold effect-site target fixed (`20 nM`) and efficacy mapping unchanged.
- Reduce systemic exposure via localization factor `lf` applied to PK gain median.

Definition:
- baseline gain median = `8.0`
- localized gain median = `8.0 * lf`

Interpretation:
- lower `lf` means less systemic blood exposure for same modeled site-level effect.

## Results
| lf | gain median | p50 ng/mL | p95 ng/mL | P(>renal caution) | P(>renal high) | P(>neuro caution) | gate pass |
|---:|------------:|----------:|----------:|------------------:|---------------:|------------------:|:---------|
| 1.00 | 8.0 | 193.23 | 417.26 | 0.10206 | 0.02136 | 0.03576 | true |
| 0.80 | 6.4 | 154.58 | 333.81 | 0.04068 | 0.00586 | 0.01136 | true |
| 0.60 | 4.8 | 115.94 | 250.35 | 0.00930 | 0.00084 | 0.00178 | true |
| 0.50 | 4.0 | 96.61 | 208.63 | 0.00270 | 0.00020 | 0.00042 | true |
| 0.40 | 3.2 | 77.29 | 166.90 | 0.00058 | 0.00006 | 0.00010 | false |
| 0.30 | 2.4 | 57.97 | 125.18 | 0.00008 | 0.00000 | 0.00002 | false |
| 0.20 | 1.6 | 38.65 | 83.45 | 0.00000 | 0.00000 | 0.00000 | false |

## Readout
- Safety-risk probabilities drop rapidly as `lf` decreases.
- Current gate flips to `false` below `lf~0.5`, mostly because the current *target-zone lower bound* (`80 ng/mL`) is no longer met, not because toxicity rises.
- This is consistent with your conceptual point: localization can decouple efficacy from systemic exposure constraints.

## Notes
- This is a conceptual translational probe, not dose guidance.
- To model true CNS-local modulation, the gate itself should be split into:
  1) efficacy-site attainment criterion
  2) systemic safety criterion
  rather than relying on a single blood-zone proxy.
