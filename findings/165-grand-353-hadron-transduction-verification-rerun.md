# 165 — GRAND-353 Verification Rerun (Hadron Transduction + Error Bars)

Date: 2026-02-28

## Scope
Re-verified the GRAND-353 lane deliverables in current head after flavor-lane updates.

## Commands

```bash
cargo run -q -p gutoe-physics --bin hadron_transduction_ci_gate
cd lean && lake build Gutoe.HadronTransduction
```

## Result
- `hadron_transduction_ci_gate overall_pass=true`
- Lean module build: `Build completed successfully`

## Notes
- Structural transduction and uncertainty lane remains green.
- Existing artifact dossier remains `findings/142-grand-353-hadron-transduction-error-bars.md`.
