# 004 — Neighbor Sweep (Combinatorial vs Mesh Degree)

## Objective

Run `gutoe_nbr_sweep.py` to test whether proton emergence is tied to a specific neighbor count or to topology/structure.

## Run setup

- Lattice: `12x12x12` (1728 cells)
- Seeds per mesh: 20
- Steps: 200
- Metric: peak proton count and UP/DN ratio

## Results

- `hex-6 (planar only)`: `278.1 ± 6.2`, `UP/DN = 2.03`
- `hex-12 (HCP)`: `360.8 ± 5.4`, `UP/DN = 1.53`
- `hex-16 (hexadecimal)`: `368.9 ± 9.8`, `UP/DN = 1.30`
- `hex-20 (full prism)`: `385.8 ± 11.5`, `UP/DN = 1.31`
- `square-6 (NO triangles)`: `0.0 ± 0.0`, `UP/DN = 0.00`

## Interpretation

- Proton count is **not invariant** under neighbor count (6→20 changes yield strongly).
- The square control is decisive: removing triangle topology collapses proton formation.
- This run supports: triangle-rich hex topology is load-bearing; degree alone is not the whole story.
