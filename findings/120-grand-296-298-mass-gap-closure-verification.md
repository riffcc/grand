# Finding 120: GRAND-296/298 Mass-Gap Lane Verification and Closure

Date: 2026-02-28
Status: GRAND-296 complete, GRAND-298 complete

## Scope

Re-verify the finite-volume Yang-Mills mass-gap scaffold + lower-bound extraction lane and close `GRAND-296` and `GRAND-298`.

## Runtime verification

Commands:

```bash
cargo run -q -p gutoe-physics --bin ym_mass_gap_report
cargo test -q -p gutoe-physics mass_gap -- --nocapture
```

Artifacts:

- `/tmp/bh_renders/ym_mass_gap_report.txt`
- `/tmp/bh_renders/ym_mass_gap_report.json`

Key reported values (from JSON):

- `monotone_nonincreasing_in_volume = true`
- `continuum_stability_band = [0.00131704512638, 0.02952447387407]`
- `continuum_fit_a2.intercept = 0.01556022861296` (`intercept_positive = true`)
- Finite-volume selected gap estimates:
  - `L=6`: `gap = 0.01774386560597`
  - `L=8`: `gap = 0.01770454945220`
  - `L=10`: `gap = 0.01689593830653`
  - `L=12`: `gap = 0.01542075950022`

This satisfies the ticket requirements for:

- computable finite-volume glueball-channel lower-bound extraction,
- explicit error-bearing outputs,
- monotone-in-volume trend check,
- continuum-stability pathway.

## Lean parity verification

Command:

```bash
cd lean && lake build Gutoe.YangMillsMassGap
```

Result:

- Build completed successfully (`8017 jobs`).

## Notes

This verification is consistent with earlier implementation findings:

- `findings/057-yang-mills-mass-gap-scaffold.md`
- `findings/058-yang-mills-empirical-z3-transfer-matrix.md`
- `findings/059-lean-validation-yang-mills-gevp-gap.md`

No new theorem-direction/sign changes were made in this closure pass.
