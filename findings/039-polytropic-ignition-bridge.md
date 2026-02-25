# Finding 039 — Polytropic Ignition Bridge (Lean + Rust Parity)

Date: 2026-02-25  
Scope: GRAND-279

## What Landed

We upgraded the stellar ignition leg of the fusion chain from a purely linear mass-threshold model to a Lane-Emden-style **polytropic compression bridge** with Lean/Rust parity.

### Lean (`lean/Gutoe/StellarFusion.lean`)

Added polytropic ignition primitives:
- `laneEmdenCompressionProxy (M rho_c) = M * rho_c`
- `coreTemperaturePolytropic (G μ ξ M rho_c) = ξ * G * μ * sqrt(M * rho_c)`
- `minimumPolytropicCompression (G μ ξ T_ign) = (T_ign / (ξ G μ))^2`

Proved:
- `polytropic_ignition_from_compression`:
  if compression exceeds threshold, then `coreTemperaturePolytropic >= T_ign`.

Integrated into a full fusion witness theorem:
- `stellar_ignition_equilibrium_exists_polytropic_from_lattice_params`
  combines pp energetics, weak-vertex existence, Gamow positivity, and the new polytropic ignition condition.

### Rust (`crates/gutoe-physics/src/equations.rs`)

Added matching runtime functions:
- `lane_emden_compression_proxy`
- `core_temperature_polytropic`
- `minimum_polytropic_compression`
- `polytropic_ignition_condition`

Added tests:
- threshold equality case (`T_core == T_ign` at exact compression threshold)
- monotonic temperature increase with stronger compression

## Verification

- `lake build Gutoe.StellarFusion` ✅
- `lake build Gutoe` ✅
- `cargo test -p gutoe-physics polytropic -- --nocapture` ✅

## Caveat (Explicit)

This is a **polytropic proxy bridge**, not yet a full Lane-Emden ODE solution proof.
To keep the gap explicit and non-hidden, follow-up ticket opened:

- GRAND-281: Upgrade proxy to full Lane-Emden ODE bridge (mass-radius + ODE-backed ignition bounds)

