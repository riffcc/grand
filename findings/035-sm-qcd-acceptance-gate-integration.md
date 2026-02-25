# 035 — SM×QCD Acceptance Gate Integrated Into Falsifiability Catalog

Status: complete.

## What changed

- Updated `lean/Gutoe/FalsifiabilityCatalog.lean` to import
  `Gutoe.SMQCDUnification`.
- Added integrated gate:
  - `smQcdAcceptanceGate`
  - `sm_qcd_acceptance_gate_holds`

## Meaning

The catalog now has a single top-level acceptance theorem that combines:

1. SM×QCD unification bundle (including Strong-CP general-case split), and
2. existing SM+GR limit recovery bundle.

This closes the loop from isolated theorem islands to one catalog-level gate.

## Verification

- `lake build Gutoe` passes after integration.
