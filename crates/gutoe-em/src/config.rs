// GUTOE EM — Lattice configuration and charge constants
// Copyright (C) 2026 Riff Labs, AGPL-3.0-or-later

pub const VOID: u8 = 0;
pub const LEPTON_SEED: u8 = 2; // γ⁰ — grade-1, mi=0b0001
pub const QUARK_SEED: u8 = 3;  // γ¹ — grade-1, mi=0b0010

pub const UP_CHARGE: f64 = 2.0 / 3.0;
pub const DOWN_CHARGE: f64 = -1.0 / 3.0;
pub const LEPTON_CHARGE: f64 = -1.0;

/// Quark classification from local field balance.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum QuarkType {
    /// veracity > curvature  →  charge +2/3
    Up,
    /// curvature ≥ veracity  →  charge −1/3
    Down,
}

/// Runtime parameters for the 12×12×12 Cl(1,3) hex-toroid lattice.
#[derive(Clone, Debug)]
pub struct LatticeConfig {
    pub hex_rows: usize,
    pub hex_cols: usize,
    pub layers: usize,
    pub differentiation_prob: f64,
    pub cycle_prob: f64,
    pub clifford_prob: f64,
    pub alignment_strength: f64,
    pub quark_threshold: f64,
    /// k=4: majority vote threshold for void activation / quark alignment
    pub void_votes: usize,
    pub em_prob: f64,
    /// Wave speed c=0.4: stable c < 1/√6 ≈ 0.408 for hex-6 lattice
    pub photon_c: f64,
    pub photon_coupling: f64,
    pub poisson_iters: usize,
}

impl Default for LatticeConfig {
    fn default() -> Self {
        Self {
            hex_rows: 12,
            hex_cols: 12,
            layers: 12,
            differentiation_prob: 0.02,
            cycle_prob: 0.05,
            clifford_prob: 0.03,
            alignment_strength: 0.15,
            quark_threshold: 0.6,
            void_votes: 4,
            em_prob: 0.5,
            photon_c: 0.4,
            photon_coupling: 0.05,
            poisson_iters: 80,
        }
    }
}

impl LatticeConfig {
    /// Total number of lattice sites.
    pub fn n_sites(&self) -> usize {
        self.hex_rows * self.hex_cols * self.layers
    }

    /// Sites per layer (hex_rows × hex_cols).
    pub fn layer_stride(&self) -> usize {
        self.hex_rows * self.hex_cols
    }
}
