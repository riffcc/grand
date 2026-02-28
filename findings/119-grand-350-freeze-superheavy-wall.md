# Finding 119 — GRAND-350 Superheavy Wall Freeze

Status: GRAND-350 closure data frozen from extended chart scan

## Goal

Freeze the prediction:
- no secondary superheavy island at `Z=164`
- practical stability wall by `Z≈145–151`

## Repro Command

```bash
GUTOE_NUCLEAR_Z_MAX=220 GUTOE_NUCLEAR_N_MAX=420 \
  cargo run -q -p gutoe-physics --bin nuclear_chart_scan
```

Primary artifacts:
- `/tmp/nuclear_chart/nuclides.csv`
- `/tmp/nuclear_chart/top_islands.csv`
- `/tmp/nuclear_chart/summary.txt`

## Extracted Statistics

Top-10 best `stability_score` by proton number for `Z>=126`:

1. `Z=126`, `N=198`, `stability=7.276975`
2. `Z=127`, `N=202`, `stability=7.213158`
3. `Z=128`, `N=204`, `stability=7.148586`
4. `Z=129`, `N=208`, `stability=7.092773`
5. `Z=130`, `N=212`, `stability=7.054990`
6. `Z=131`, `N=214`, `stability=7.026514`
7. `Z=132`, `N=218`, `stability=7.007930`
8. `Z=133`, `N=222`, `stability=6.986465`
9. `Z=134`, `N=226`, `stability=6.966395`
10. `Z=135`, `N=230`, `stability=6.940615`

`Z=164` best point:
- `N=350`, `A=514`
- `stability_score=5.976536`
- `sf_log10_half_life_s=-21.969322`
- `fission_barrier_mev=0.003305`
- rank `39/95` among `Z>=126` best-per-Z points

## Wall Signature

Best-per-Z stability declines monotonically through the proposed wall region:

- `Z=145`: `6.637655`
- `Z=150`: `6.470545`
- `Z=151`: `6.433786`
- `Z=164`: `5.976536`

Interpretation:
- no revival/secondary peak emerges near `Z=164`
- stability degrades continuously past the `145–151` band
- fission barrier at `Z=164` is effectively collapsed (`~3.3e-3 MeV`)

## Conclusion

Freeze accepted:
- **No secondary superheavy island at `Z=164`**
- **Practical wall near `Z≈145–151`**
