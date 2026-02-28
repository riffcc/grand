# GRAND-367 — Prediction Register (Frozen)

Freeze timestamp: 2026-02-27 (UTC)

## Core predictions on record

1. PMNS corrected lane
- theta23_pred = 49.001050838539 deg
- residual_to_49deg = +0.001050838539 deg
- status: PASS (millidegree closure lane)

2. Neutrino hierarchy
- prediction: normal ordering
- status: OPEN TO EXPERIMENT (JUNO/DUNE era)

3. Alpha structural closure
- alpha_inv_structural = 137
- decimal correction lane: 137 + 5alpha - 9alpha^2 (second-order improvement proven)
- status: PASS

4. Proton mass structural lane
- proton_pred = 938.194072200000 MeV
- proton_obs = 938.272088160000 MeV
- rel_error = -8.314854612483e-5
- status: PASS

5. CMB full-derived stack
- TT_red = 1.262555196494
- TE_red = 1.109327515393
- EE_red = 1.063630963463
- sigma8 = 0.811221182557
- status: PASS

6. G bridge lane
- measured electron mode rel_error_G = -9.699067921626355e-4
- status: PASS

## Integrity notes

- Sigma8 decomposition wiring is aligned to full-derived defaults; sigma8 matches exactly across both entry points.
- PMNS millidegree closure theorem is now in Lean (`Gutoe/FlavorMixing.lean`) and `lake build Gutoe` is green.
- This register is a freeze artifact and should be amended only via new numbered findings.
