# Finding 250 — “Duet” Mode vs UV Triplet Trend in Yukawa Dynamics

## Thesis
In the current full-dynamics lane, quark Yukawa structure separates cleanly into two modes:

- **Ladder mode** `L_g = [sqrt(m_u m_d), sqrt(m_c m_s), sqrt(m_t m_b)]`  
  This is the Z3-harmonic lane with free-fit `s²(L_g)` that rises toward `3` in the UV.

- **Duet mode** `S_g = [0.5 ln(m_u/m_d), 0.5 ln(m_c/m_s), 0.5 ln(m_t/m_b)]`  
  This is the **isospin split** lane (the “duet”: up/down pair in each generation).

So, in physics language:  
**duet = SU(2) up/down pairing signal (`S_g`)**  
while the UV `s²` trend in `L_g` tracks triplet-like structure.

## Data (extended to 1e19 GeV)
From `/tmp/bh_renders/yukawa_full_dynamics_scan.csv`:

- `s²(L_g)`:
  - `m_t`: `2.917521825`
  - `1e4`: `2.925051167`
  - `1e8`: `2.938855312`
  - `1e12`: `2.950503142`
  - `1e16`: `2.961017988`
  - `1e18`: `2.965986213`
  - `1e19`: `2.968411312`

- Fixed `s²=3` closure RMS:
  - `0.445229 -> 0.422888 -> 0.364955 -> 0.315458 -> 0.260399 -> 0.225460 -> 0.212560`

All are monotonic in the favorable direction.

## Shape read
- `s²(L_g)` increases monotonically.
- Gap `Δ = 3 - s²` decreases monotonically.
- Increment per decade decreases (`ds²/dlog10(mu)` decelerates), consistent with an asymptotic approach.

At present depth:
- evidence supports **approach toward 3 from below**,
- but no strict fixed-point lock at `3` is proven by `1e19`.

## Interpretation boundary
What is supported now:
- “Duet” is a useful physical mode label for `S_g` (generation-wise up/down splitting).
- UV trend in `L_g` is compatible with triplet-like limit behavior.

What is not yet proven:
- exact identity `s² = 3` at finite scale,
- unique structural law for residual gap `3 - s²`,
- full first-principles derivation of `S_g` slope/intercept from Cl(1,3) counts alone.
