# GRAND — Triangulation: Because vs Fit (Formal Split)

Date: 2026-02-28

## New Formal Provenance Module

- `lean/Gutoe/TriangulatedTermProvenance.lean`

Added theorem-level provenance for the key terms:

- `7 = |grade₂| + 1 = C(4,2)+1`
- `19/3 = (d + |SU(2)|)/|SU(2)| = (16+3)/3`
- `1/36 = 1/|grade₂|²`
- `1/(7*13*136)` from
  `1 / ((|grade₂|+1) * (d-|SU(2)|) * T(16))`

And decomposition theorems:

- `kappa` candidate is exactly geometric dark/visible ratio times the sum of
  the three provenance terms.
- EW uplift term is exactly
  `|grade₂|/(d-|SU(2)|) - 1/((|grade₂|+1)T(16))`.
- `p` candidate is structural baseline `α⁻¹/10` minus finite-lattice term
  `1/((|grade₂|+1)*N_gauge)`.

## Build Verification

- `cd lean && lake build Gutoe`
- Success: `Gutoe.TriangulatedTermProvenance` built, full library green.

## Honest Status

Because (now formal):

- each term has explicit Cl(1,3)-count provenance;
- each candidate constant is decomposed into those provenance operations.

Still fit (not yet uniqueness-closed):

- uniqueness proof that this combination is the only admissible combination
  under a constrained operator grammar is not formalized yet.
- current state is "structurally justified reconstruction" with theorem-backed
  provenance, not final uniqueness closure.
