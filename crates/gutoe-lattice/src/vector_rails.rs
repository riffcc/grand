/*!
 * GUTOE Lattice - Vector Rails
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

//! Vector rails - wave-like connections for state propagation
//!
//! From GUTOE.md:
//! "Vector rails are fundamentally wave equations that manifest as:
//!  - Oscillatory patterns of connection strength
//!  - Propagating waves of state information
//!  - Coherent field-like structures spanning entities"

use super::hex_lattice::HexCoord;
use std::collections::HashMap;

/// A vector rail connection between two lattice points
#[derive(Debug, Clone)]
pub struct VectorRail {
    pub from: HexCoord,
    pub to: HexCoord,
    pub veracity: f64,      // Strength of the connection
    pub phase: f64,        // Current phase of oscillation
    pub frequency: f64,    // Oscillation frequency
}

impl VectorRail {
    pub fn new(from: HexCoord, to: HexCoord) -> Self {
        Self {
            from,
            to,
            veracity: 1.0,
            phase: 0.0,
            frequency: 1.0,
        }
    }

    /// Update phase based on time step
    pub fn tick(&mut self, dt: f64) {
        self.phase += 2.0 * std::f64::consts::PI * self.frequency * dt;
        self.phase = self.phase % (2.0 * std::f64::consts::PI);
    }

    /// Get current oscillation amplitude
    pub fn amplitude(&self) -> f64 {
        self.veracity * self.phase.sin()
    }

    /// Get current state value (complex-like)
    pub fn state_value(&self) -> (f64, f64) {
        (self.veracity * self.phase.cos(), self.veracity * self.phase.sin())
    }
}

/// Collection of vector rails forming a network
#[derive(Debug, Clone)]
pub struct RailNetwork {
    rails: HashMap<(HexCoord, HexCoord), VectorRail>,
}

impl RailNetwork {
    pub fn new() -> Self {
        Self {
            rails: HashMap::new(),
        }
    }

    /// Add a rail between two coordinates
    pub fn add_rail(&mut self, from: HexCoord, to: HexCoord) {
        let rail = VectorRail::new(from, to);
        self.rails.insert((from, to), rail);
    }

    /// Get rail between coordinates
    pub fn get_rail(&self, from: &HexCoord, to: &HexCoord) -> Option<&VectorRail> {
        self.rails.get(&(*from, *to))
    }

    /// Get rail mutably
    pub fn get_rail_mut(&mut self, from: &HexCoord, to: &HexCoord) -> Option<&mut VectorRail> {
        self.rails.get_mut(&(*from, *to))
    }

    /// Get all rails from a coordinate
    pub fn rails_from(&self, coord: &HexCoord) -> Vec<&VectorRail> {
        self.rails
            .iter()
            .filter(|(k, _)| k.0 == *coord)
            .map(|(_, v)| v)
            .collect()
    }

    /// Number of rails
    pub fn size(&self) -> usize {
        self.rails.len()
    }

    /// Tick all rails forward
    pub fn tick_all(&mut self, dt: f64) {
        for rail in self.rails.values_mut() {
            rail.tick(dt);
        }
    }

    /// Calculate total veracity in network
    pub fn total_veracity(&self) -> f64 {
        self.rails.values().map(|r| r.veracity).sum()
    }

    /// Average veracity
    pub fn average_veracity(&self) -> f64 {
        if self.rails.is_empty() {
            return 0.0;
        }
        self.total_veracity() / self.rails.len() as f64
    }
}

impl Default for RailNetwork {
    fn default() -> Self {
        Self::new()
    }
}

/// Veracity wave equation simulator
/// Tests: wave propagation with dispersion
/// From GUTOE.md: "ω² = v² k² - λ_QG l_P² k⁴"
#[derive(Debug)]
pub struct WaveSimulator {
    pub lambda_qg: f64,   // Quantum gravity coupling
    pub velocity: f64,     // Wave velocity
    pub dt: f64,           // Time step
}

impl WaveSimulator {
    pub fn new(lambda_qg: f64, velocity: f64) -> Self {
        Self {
            lambda_qg,
            velocity,
            dt: 0.01,
        }
    }

    /// Calculate dispersion relation: ω² = v² k² - λ_QG l_P² k⁴
    pub fn dispersion(&self, k: f64, l_p: f64) -> f64 {
        let v2_k2 = self.velocity * self.velocity * k * k;
        let lambda_correction = self.lambda_qg * l_p * l_p * k * k * k * k;
        (v2_k2 - lambda_correction).sqrt()
    }

    /// Check if wave is stable (real frequency)
    pub fn is_stable(&self, k: f64, l_p: f64) -> bool {
        self.velocity * self.velocity >= self.lambda_qg * l_p * l_p * k * k
    }

    /// Find critical wave number where instability begins
    pub fn critical_k(&self, l_p: f64) -> f64 {
        if self.lambda_qg * l_p * l_p > 0.0 {
            self.velocity / (self.lambda_qg.sqrt() * l_p)
        } else {
            f64::INFINITY
        }
    }
}

/// Propagation test - simulates signal through rail network
pub fn simulate_propagation(
    network: &RailNetwork,
    start: &HexCoord,
    steps: usize,
) -> Vec<f64> {
    let mut amplitudes = Vec::with_capacity(steps);
    let mut current = *start;

    for _ in 0..steps {
        // Get outgoing rails
        let outgoing = network.rails_from(&current);
        if outgoing.is_empty() {
            amplitudes.push(0.0);
            break;
        }

        // Follow the strongest rail
        let best = outgoing.iter().max_by(|a, b| {
            a.veracity.partial_cmp(&b.veracity).unwrap()
        });

        if let Some(rail) = best {
            amplitudes.push(rail.amplitude());
            current = rail.to;
        } else {
            amplitudes.push(0.0);
            break;
        }
    }

    amplitudes
}

#[cfg(test)]
mod tests {
    use super::*;
    use gutoe_core::constants::LAMBDA_QG;

    #[test]
    fn test_rail_creation() {
        let from = HexCoord::new(0, 0);
        let to = HexCoord::new(1, 0);
        let rail = VectorRail::new(from, to);
        assert_eq!(rail.from, from);
        assert_eq!(rail.to, to);
    }

    #[test]
    fn test_wave_oscillation() {
        let mut rail = VectorRail::new(HexCoord::new(0, 0), HexCoord::new(1, 0));
        rail.frequency = 1.0;

        // Initial phase = 0, amplitude = 0
        assert!((rail.amplitude()).abs() < 0.001);

        // Tick forward
        rail.tick(0.25);
        // At t=0.25 with f=1, phase = π/2, sin(π/2) = 1
        assert!((rail.amplitude() - rail.veracity).abs() < 0.001);
    }

    #[test]
    fn test_dispersion_relation() {
        let sim = WaveSimulator::new(LAMBDA_QG, 1.0);

        // Without quantum gravity correction (λ_QG → 0), ω = vk
        let omega_classical = 1.0 * 1.0; // v * k

        // With correction, should be slightly less
        let omega_quantum = sim.dispersion(1.0, 0.1);

        assert!(omega_quantum < omega_classical);
    }

    #[test]
    fn test_critical_wave_number() {
        let sim = WaveSimulator::new(LAMBDA_QG, 1.0);
        let k_crit = sim.critical_k(0.1);

        // Critical k should be finite
        assert!(k_crit.is_finite());
        assert!(k_crit > 0.0);
    }

    #[test]
    fn test_stability_check() {
        let sim = WaveSimulator::new(LAMBDA_QG, 1.0);

        // Low k should be stable
        assert!(sim.is_stable(0.1, 0.1));

        // Very high k might be unstable
        assert!(!sim.is_stable(100.0, 0.1));
    }
}
