# Finding 272 — RH Named Obligation Contracts (Lean)

Date: 2026-03-01  
Scope: Make the remaining RH gap a single, explicit Lean contract object.

## Updated

- `lean/Gutoe/RiemannLimitBridge.lean`
- `lean/Gutoe/RiemannRHClosure.lean`

## New named obligations

In `RiemannLimitBridge`:

- `FiniteBridgeFamily XiN specN`
  - every finite level has an exact bridge.
- `ZeroForwardTransfer Xi XiN`
  - every `Xi` zero appears in some finite level.
- `ZeroBackwardTransfer Xi XiN`
  - (optional strengthening) finite-level zero implies `Xi` zero.

## New contract records

- `RHLimitTransferContract Xi`
  - fields:
    - `XiN`
    - `specN`
    - `finiteBridge`
    - `zeroForward`
- `RHExactLimitTransferContract Xi`
  - extends the above with `zeroBackward`.

## Closure theorems from contracts

In `RiemannLimitBridge`:

- `rh_of_limit_transfer_contract`
- `spectralBridge_of_exact_limit_transfer_contract`

In `RiemannRHClosure`:

- `rh_from_limit_transfer_contract`

## Why this matters

The unresolved RH frontier is now one importable Lean interface:

- fill `RHLimitTransferContract Xi` to get RH-for-`Xi`,
- fill `RHExactLimitTransferContract Xi` to get exact infinite bridge + RH.

This removes ambiguity about “what remains” and prevents proof drift.

## Build verification

Executed:

```bash
cd lean
lake build Gutoe
```

Result: **passes** (`8168` jobs, warnings only).

