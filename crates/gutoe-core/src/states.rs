/*!
 * GUTOE Core - Tripartite Quantum States
 * Copyright (C) 2026  Riff Labs
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

//! Core tripartite quantum states
//!
//! The GUTOE framework uses a 4-state system instead of traditional 2-state qubits:
//! - VOID: Pure undefined state (the null reference from which all emerges)
//! - COSINE: |0⟩ state (passive non-existence, quantum vacuum)
//! - SINE: |1⟩ state (active existence)
//! - TANGENT: tan = sin/cos — the slope/ratio state, diverges when cos = 0

use num_complex::Complex64;
use rust_decimal::Decimal;
use std::fmt;

/// Fundamental GUTOE constants derived from vector rail simulations
pub mod constants {
    /// Quantum gravity coupling constant - derived from simulation
    /// This is a KEY value to verify - the framework claims λ_QG ≈ 0.084372
    pub const LAMBDA_QG: f64 = 0.084372;

    /// Maximum virtualized qubits supported
    pub const MAX_QUBITS: usize = 3_000_000;

    /// Default quantum coherence (99.9999%)
    pub const DEFAULT_COHERENCE: f64 = 0.999999;

    /// Planck length approximation (in simulation units)
    pub const PLANCK_LENGTH: f64 = 1.0e-35;

    /// Wave velocity (should reduce to c in appropriate units)
    pub const WAVE_VELOCITY: f64 = 299_792_458.0;
}

/// The four fundamental states in the GUTOE tripartite system
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TriState {
    /// VOID - Pure undefined state (the null reference)
    VOID,
    /// COSINE - |0⟩ passive non-existence
    COSINE,
    /// SINE - |1⟩ active existence
    SINE,
    /// TANGENT - tan = sin/cos, the slope/ratio state; diverges when cos = 0
    TANGENT,
}

impl TriState {
    /// Get all possible states
    pub fn variants() -> [TriState; 4] {
        [TriState::VOID, TriState::COSINE, TriState::SINE, TriState::TANGENT]
    }

    /// Convert to amplitude representation
    /// Returns (real, imag) components for state vector
    pub fn to_amplitude(&self) -> (Complex64, Complex64) {
        match self {
            // |0⟩ = COSINE
            TriState::COSINE => (Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)),
            // |1⟩ = SINE
            TriState::SINE => (Complex64::new(0.0, 0.0), Complex64::new(1.0, 0.0)),
            // Superposition = (|0⟩ + |1⟩)/√2
            TriState::TANGENT => (
                Complex64::new(1.0 / 2.0_f64.sqrt(), 0.0),
                Complex64::new(1.0 / 2.0_f64.sqrt(), 0.0),
            ),
            // VOID - special case, represents null reference
            TriState::VOID => (Complex64::new(0.0, 0.0), Complex64::new(0.0, 0.0)),
        }
    }

    /// Cycle order: SINE → COSINE → TANGENT → SINE (Z₃ group)
    /// This is the TripartiteNOT operation
    pub fn cycle(&self) -> TriState {
        match self {
            TriState::SINE => TriState::COSINE,
            TriState::COSINE => TriState::TANGENT,
            TriState::TANGENT => TriState::SINE,
            TriState::VOID => TriState::VOID, // Void is fixed point
        }
    }

    /// Check if state is a basis state (not superposition)
    pub fn is_basis(&self) -> bool {
        matches!(self, TriState::SINE | TriState::COSINE)
    }

    /// Get the phase factor for this state
    pub fn phase(&self) -> f64 {
        match self {
            TriState::SINE => 0.0,
            TriState::COSINE => std::f64::consts::PI,
            TriState::TANGENT => std::f64::consts::PI / 2.0,
            TriState::VOID => 0.0,
        }
    }
}

impl fmt::Display for TriState {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            TriState::VOID => write!(f, "VOID"),
            TriState::COSINE => write!(f, "|0⟩"),
            TriState::SINE => write!(f, "|1⟩"),
            TriState::TANGENT => write!(f, "|T⟩"),
        }
    }
}

/// A quantum register of tripartite states
#[derive(Debug, Clone)]
pub struct Register {
    states: Vec<TriState>,
    coherence: f64,
}

impl Register {
    /// Create a new register with n qubits in |0⟩ (COSINE) state
    pub fn new(n: usize) -> Self {
        Self {
            states: vec![TriState::COSINE; n],
            coherence: constants::DEFAULT_COHERENCE,
        }
    }

    /// Get the number of qubits
    pub fn len(&self) -> usize {
        self.states.len()
    }

    /// Check if register is empty
    pub fn is_empty(&self) -> bool {
        self.states.is_empty()
    }

    /// Apply a state transition (gate) to a specific qubit
    pub fn apply_gate(&mut self, gate: &dyn Fn(TriState) -> TriState, index: usize) {
        if index < self.states.len() {
            self.states[index] = gate(self.states[index]);
        }
    }

    /// Get state at index
    pub fn get(&self, index: usize) -> Option<TriState> {
        self.states.get(index).copied()
    }

    /// Set coherence
    pub fn set_coherence(&mut self, coherence: f64) {
        self.coherence = coherence.clamp(0.0, 1.0);
    }

    /// Get coherence
    pub fn coherence(&self) -> f64 {
        self.coherence
    }
}

/// Quantum amplitude with complex number representation
#[derive(Debug, Clone)]
pub struct Amplitude {
    pub re: f64,
    pub im: f64,
}

impl Amplitude {
    pub fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }

    pub fn magnitude(&self) -> f64 {
        (self.re * self.re + self.im * self.im).sqrt()
    }

    pub fn phase(&self) -> f64 {
        self.im.atan2(self.re)
    }

    /// Multiply by another amplitude
    pub fn multiply(&self, other: &Amplitude) -> Amplitude {
        Amplitude {
            re: self.re * other.re - self.im * other.im,
            im: self.re * other.im + self.im * other.re,
        }
    }

    /// Add another amplitude
    pub fn add(&self, other: &Amplitude) -> Amplitude {
        Amplitude {
            re: self.re + other.re,
            im: self.im + other.im,
        }
    }

    /// Complex conjugate
    pub fn conj(&self) -> Amplitude {
        Amplitude {
            re: self.re,
            im: -self.im,
        }
    }

    /// Normalize to unit length
    pub fn normalize(&self) -> Amplitude {
        let mag = self.magnitude();
        if mag > 0.0 {
            Amplitude {
                re: self.re / mag,
                im: self.im / mag,
            }
        } else {
            self.clone()
        }
    }
}

impl Default for Amplitude {
    fn default() -> Self {
        Amplitude::new(0.0, 0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── Cycle group structure ────────────────────────────────────────────────

    #[test]
    fn cycle_map_is_correct() {
        assert_eq!(TriState::SINE.cycle(),      TriState::COSINE);
        assert_eq!(TriState::COSINE.cycle(),    TriState::TANGENT);
        assert_eq!(TriState::TANGENT.cycle(), TriState::SINE);
    }

    #[test]
    fn cycle_is_order_3_on_non_void() {
        // Applying cycle 3 times must return the original for each non-VOID state.
        for s in [TriState::SINE, TriState::COSINE, TriState::TANGENT] {
            assert_eq!(s.cycle().cycle().cycle(), s,
                "cycle³({s:?}) ≠ {s:?} — cycle is not Z₃");
        }
    }

    #[test]
    fn void_is_only_fixed_point() {
        // cycle(VOID) = VOID, and no other state is fixed.
        assert_eq!(TriState::VOID.cycle(), TriState::VOID);
        for s in [TriState::SINE, TriState::COSINE, TriState::TANGENT] {
            assert_ne!(s.cycle(), s, "cycle({s:?}) = {s:?} — unexpected fixed point");
        }
    }

    #[test]
    fn cycle_is_not_involutive() {
        // If cycle were its own inverse, cycle(cycle(s)) = s.
        // This is FALSE for this system (it's Z₃ not Z₂).
        assert_ne!(TriState::SINE.cycle().cycle(), TriState::SINE,
            "cycle is involutive — that would make it Z₂, not Z₃");
    }

    // ── Amplitude / probability interpretation ───────────────────────────────

    #[test]
    fn amplitudes_are_normalized_for_basis_states() {
        // For basis states, the probability should sum to 1.
        for s in [TriState::COSINE, TriState::SINE] {
            let (a0, a1) = s.to_amplitude();
            let prob = a0.norm_sqr() + a1.norm_sqr();
            assert!((prob - 1.0).abs() < 1e-12,
                "{s:?}: probability sum = {prob}, expected 1.0");
        }
    }

    #[test]
    fn tangent_amplitudes_are_equal_superposition() {
        // TANGENT reference amplitude: (1/√2)|0⟩ + (1/√2)|1⟩, each with prob 0.5.
        let (a0, a1) = TriState::TANGENT.to_amplitude();
        let p0 = a0.norm_sqr();
        let p1 = a1.norm_sqr();
        assert!((p0 - 0.5).abs() < 1e-12, "p(|0⟩) = {p0}, expected 0.5");
        assert!((p1 - 0.5).abs() < 1e-12, "p(|1⟩) = {p1}, expected 0.5");
    }

    #[test]
    fn void_amplitude_is_zero_probability() {
        // VOID represents the null reference — zero probability of observation.
        let (a0, a1) = TriState::VOID.to_amplitude();
        assert_eq!(a0.norm_sqr() + a1.norm_sqr(), 0.0,
            "VOID has nonzero probability — breaks physical interpretation");
    }

    // ── Phase structure ──────────────────────────────────────────────────────

    #[test]
    fn phases_are_ordered_correctly() {
        // SINE=0, TANGENT=π/2, COSINE=π
        // The ordering should reflect real→imaginary→negative-real on the unit circle.
        use std::f64::consts::{FRAC_PI_2, PI};
        assert!((TriState::SINE.phase() - 0.0).abs() < 1e-12);
        assert!((TriState::TANGENT.phase() - FRAC_PI_2).abs() < 1e-12);
        assert!((TriState::COSINE.phase() - PI).abs() < 1e-12);
    }

    // ── λ_QG constant ────────────────────────────────────────────────────────

    #[test]
    fn lambda_qg_matches_claimed_value() {
        // GUTOE claims λ_QG ≈ 0.084372.  This test pins that number down.
        // If the constant changes, it must be a deliberate decision.
        assert!((constants::LAMBDA_QG - 0.084372).abs() < 1e-6,
            "λ_QG = {} ≠ 0.084372 — the claimed value has shifted",
            constants::LAMBDA_QG);
    }

    // ── Amplitude arithmetic ─────────────────────────────────────────────────

    #[test]
    fn amplitude_normalization_is_idempotent() {
        let amp = Amplitude::new(3.0, 4.0);
        let n1  = amp.normalize();
        let n2  = n1.normalize(); // normalizing a unit vector should be identity
        assert!((n1.re - n2.re).abs() < 1e-12 && (n1.im - n2.im).abs() < 1e-12);
    }

    #[test]
    fn amplitude_magnitude_after_normalize_is_one() {
        let amp = Amplitude::new(3.0, 4.0);
        assert!((amp.normalize().magnitude() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn amplitude_multiply_preserves_magnitude() {
        // |a·b| = |a|·|b|
        let a = Amplitude::new(3.0, 4.0);
        let b = Amplitude::new(1.0, 2.0);
        let product = a.multiply(&b);
        let expected_mag = a.magnitude() * b.magnitude();
        assert!((product.magnitude() - expected_mag).abs() < 1e-10);
    }

    #[test]
    fn amplitude_conjugate_times_self_is_real() {
        // a * conj(a) = |a|² + 0i — imaginary part must be zero
        let a = Amplitude::new(3.0, 4.0);
        let product = a.multiply(&a.conj());
        assert!(product.im.abs() < 1e-12,
            "a·ā has nonzero imaginary part: {}", product.im);
        assert!((product.re - a.magnitude().powi(2)).abs() < 1e-10);
    }

    // ── Register ─────────────────────────────────────────────────────────────

    #[test]
    fn register_initializes_to_cosine() {
        let reg = Register::new(4);
        assert_eq!(reg.len(), 4);
        for i in 0..4 {
            assert_eq!(reg.get(i), Some(TriState::COSINE),
                "qubit {i} not in |0⟩ = COSINE state on init");
        }
    }

    #[test]
    fn register_coherence_clamped_to_unit_interval() {
        let mut reg = Register::new(1);
        reg.set_coherence(1.5);
        assert_eq!(reg.coherence(), 1.0);
        reg.set_coherence(-0.1);
        assert_eq!(reg.coherence(), 0.0);
    }
}
