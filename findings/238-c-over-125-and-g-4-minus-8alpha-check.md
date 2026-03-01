# 238 — Check: `c = -1/5^3` and `g = 4 - 8α`

## Setup
Locked values:
- `b = 9 + 5/32 = 9.15625`
- `alpha^-1_phys = 137.035999177`
- `alpha = 1/alpha^-1_phys`
- `mp/me_phys = 1836.15267343`

Tested hypotheses:
1. `c = -1/125`
2. `g = 4 - 8alpha`

Formulas:
- `alpha^-1 = 137 + 5alpha - b alpha^2 + c alpha^3`
- `mp/me = 6pi^5 + 5alpha - g b alpha^2 + g c alpha^3`

## Exact outputs
### Cubic candidate
- `c_exact` from strict alpha closure with locked `b`:
  - `c_exact = -0.007996757405957876626533392673425759...`
- candidate:
  - `c_125 = -1/125 = -0.008`
- difference:
  - `c_125 - c_exact = -3.242594042123373e-6`

### Alpha residual for `c = -1/125`
- `alpha^-1_pred = 137.035999176998739947709320116065796...`
- residual:
  - `Delta = -1.260052290679883934203596923907e-12`
- normalized:
  - `9.195045814584537e-9 ppm`
  - `9.195045814584537e-6 ppb`

Interpretation: `c = -1/125` is effectively exact at this precision level for alpha lane.

### mp/me residuals
#### A) `g = 4`, `c = -1/125`
- residual:
  - `Delta = -2.82987812920081685462040935508e-5`
- normalized:
  - `0.015411997978983423 ppm`
  - `15.411997978983423 ppb`

#### B) `g = 4 - 8alpha`, `c = -1/125`
- `g = 3.941621179485348599758033637882466...`
- residual:
  - `Delta = +1.65903949562334745510538154617e-7`
- normalized:
  - `0.00009035411486367315 ppm`
  - `0.09035411486367315 ppb`

#### C) `g = 4 - 8alpha`, `c = c_exact` (reference)
- residual:
  - `Delta = +1.65908916211130948369977376741e-7`
- normalized:
  - `0.00009035681978514731 ppm`
  - `0.09035681978514731 ppb`

## Verdict
- `c = -1/125` passes extremely strongly for the alpha lane.
- Applying `g = 4 - 8alpha` drops `mp/me` residual from `15.412 ppb` to `0.09035 ppb`.
- This is below `1 ppb` by more than an order of magnitude.
