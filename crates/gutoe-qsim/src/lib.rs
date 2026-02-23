// GUTOE QSIM — Content-Addressed Quantum State Simulation
// Copyright (C) 2026 Riff Labs, AGPL-3.0-or-later
//
// Empirically measures state-space compression in a Cl(1,3) lattice
// by evolving CoW clones with BLAKE3 content-addressed deduplication.
//
// Key insight: physics produces polynomially many distinguishable histories
// (area law). Identical starting states + same deterministic gate = one new
// state, not many. This crate measures that compression.

pub mod actor;
pub mod gate;
pub mod metrics;
pub mod orchestrator;
pub mod state;
pub mod store;

pub use actor::Actor;
pub use gate::{
    CompositeGate, EmHopGate, FullStepGate, Gate, XorProductGate, Z3CycleGate,
};
pub use metrics::{ExperimentResult, RoundMetrics};
pub use orchestrator::Orchestrator;
pub use state::{CowLattice, LatticeState};
pub use store::StateStore;
