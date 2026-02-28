# 137 — MS Cyclosporine Fragility Probe

## Scope
Technical fragility probe focused on the most sensitive mapping:
`Ki -> effective shift (efficiency scalar) -> ARR proxy`.

Requested analyses:
1. Efficiency sensitivity sweep at fixed values:
   - `0.15, 0.20, 0.25, 0.30, 0.35, 0.40`
2. Off-target penalty amplification with stochastic adverse-event disability noise.

## 1) Efficiency sensitivity sweep
ARR reduction (combo vs standard, 2-year horizon):
- `eff=0.15 -> 6.97%`
- `eff=0.20 -> 9.90%`
- `eff=0.25 -> 13.18%`
- `eff=0.30 -> 16.86%`
- `eff=0.35 -> 20.88%`
- `eff=0.40 -> 20.88%`

10-year lesion reduction (combo vs standard):
- `18.86%, 26.71%, 35.43%, 45.19%, 55.76%, 55.76%`

Sample-size consequence (Poisson ARR approx, 2y, 80% power):
- `n/arm` falls from `5679` at `eff=0.15` to `921` at `eff=0.30`,
  then `587` at `eff>=0.35`.

### Linearity/threshold diagnosis
- `0.15 -> 0.35`: near-linear monotone gain.
- `>=0.35`: hard saturation/ceiling (combo drive clipped to zero in this model).
- Conclusion: not chaotic-threshold fragile in the requested band, but **ceiling-fragile** due model clipping above ~0.35.

## 2) Off-target adverse-event noise amplification (stochastic)
Reference efficiency: `0.30` (the ~45% lesion-reduction point).

Noise model:
- Event rate scaled by `off_target_occupancy * amplification`.
- Monthly stochastic events inject lesion and disability perturbations.
- `10,000` Monte Carlo samples.

Results:
- `amp=1.00`: lesion reduction mean `44.80%` (p05 `43.01%`), disability mean `0.02181`
- `amp=1.25`: lesion reduction mean `44.71%` (p05 `42.72%`), disability mean `0.02192`
- `amp=1.50` (mild): lesion reduction mean `44.60%` (p05 `42.41%`), disability mean `0.02204`

Against standard disability (`0.03874`), probability combo remains better:
- `amp=1.00`: `1.0000`
- `amp=1.25`: `1.0000`
- `amp=1.50`: `0.9999`

### Robustness call
At mild off-target noise amplification, the ~45% lesion-reduction signal **survives** with only modest degradation (~0.6 percentage points from deterministic 45.19% to noisy mean 44.60%).

## Plot and artifacts
- ARR plot:
  - `/tmp/bh_renders/ms_cyclosporine_transduction_sweep/arr_reduction_vs_efficiency.png`
- Efficiency sweep CSV/JSON:
  - `/tmp/bh_renders/ms_cyclosporine_transduction_sweep/ms_cyclosporine_transduction_sweep.csv`
  - `/tmp/bh_renders/ms_cyclosporine_transduction_sweep/ms_cyclosporine_transduction_sweep.json`
- Off-target noise sweep CSV:
  - `/tmp/bh_renders/ms_cyclosporine_transduction_sweep/ms_cyclosporine_offtarget_noise_sweep.csv`
- Text summary:
  - `/tmp/bh_renders/ms_cyclosporine_transduction_sweep/ms_cyclosporine_transduction_sweep.txt`
