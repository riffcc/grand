# Finding 038 — pp-Chain Q-Values Derived from Shared Mass Infrastructure

Date: 2026-02-25  
Scope: GRAND-278

## What Changed

We removed direct pp-chain Q-value literals from the Lean fusion theorem path and from Rust's reaction-network source, replacing them with derivations from a shared rounded mass table.

### Lean (`lean/Gutoe/StellarFusion.lean`)

Added a shared mass layer:
- `FusionNucleus` enum (`H1`, `H2`, `He3`, `He4`)
- `fusionNuclearRestMassMeV : FusionNucleus → ℚ`
- `protonRestMassMeV`, `electronRestMassMeV`
- `positronAnnihilationThermalMeV`

Reworked pp-chain Q definitions to be mass-derived:
- `q_pp1_mev = 2 m_p - m_d - m_e + 2 m_e`
- `q_pp2_mev = m_d + m_p - m_He3`
- `q_pp3_mev = 2 m_He3 - m_He4 - 2 m_p`

`ppChainNetQMeV` remains exactly `6683/250 = 26.732 MeV` and `pp_chain_exothermic` remains proven.

### Rust (`crates/gutoe-physics/src/stellar_reactions.rs`)

Added the same rounded mass table and helper functions:
- `q_pp1_thermalized_mev()`
- `q_pp2_mev()`
- `q_pp3_mev()`

Reaction graph now uses these derived values instead of literals.

## Why This Matters

- Removes duplicated magic-number Q literals in both Lean and runtime code.
- Keeps Lean/Rust parity on the pp-chain baseline while preserving existing network values.
- Makes the next upgrades (GRAND-279/280) cleaner by centralizing reaction energetics.

## Verification

- `lake build Gutoe.StellarFusion` ✅
- `lake build Gutoe` ✅
- `cargo test -p gutoe-physics pp_q_values_match_baseline -- --nocapture` ✅

## Remaining Follow-ups

- GRAND-279: replace linear ignition model with a polytropic/Lane-Emden style theorem chain.
- GRAND-280: promote weak+Gamow structure to a strict positive pp reaction-rate theorem.
