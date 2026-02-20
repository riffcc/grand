/*!
 * GUTOE Core - Hexagonal States (12-state system)
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

//! 12-State Hexagonal System
//!
//! From VOID-DIFFERENTIATION.md:
//! "The mathematical progression takes an unexpected turn when moving beyond binary logic.
//! Instead of progressing to four or eight states, the system evolves into a hexagonal
//! structure with twelve states (six on each face). This structure provides rotational
//! and reflective symmetry, creating a stable configuration that can support balanced,
//! multi-directional interactions."
//!
//! The 12 states form two hexagonal faces:
//! - Face A (positive): 0° → 60° → 120° → 180° → 240° → 300° → 0°
//! - Face B (negative): Same angles but negated (dual face)

use std::fmt;

/// 12-state hexagonal system
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum HexState {
    // Face A (positive/basis states) - 6 states at angles 0, 60, 120, 180, 240, 300
    A0,   // |0⟩ - base
    A60,  // |1⟩
    A120, // |2⟩
    A180, // |3⟩
    A240, // |4⟩
    A300, // |5⟩

    // Face B (negative/dual states) - 6 states negated
    B0,   // -|0⟩ - dual of A0
    B60,  // -|1⟩
    B120, // -|2⟩
    B180, // -|3⟩
    B240, // -|4⟩
    B300, // -|5⟩
}

impl HexState {
    /// Get all 12 states
    pub fn all() -> [HexState; 12] {
        [
            HexState::A0, HexState::A60, HexState::A120, HexState::A180, HexState::A240, HexState::A300,
            HexState::B0, HexState::B60, HexState::B120, HexState::B180, HexState::B240, HexState::B300,
        ]
    }

    /// Get angle in degrees (0-360)
    pub fn angle(&self) -> f64 {
        match self {
            HexState::A0 => 0.0,
            HexState::A60 => 60.0,
            HexState::A120 => 120.0,
            HexState::A180 => 180.0,
            HexState::A240 => 240.0,
            HexState::A300 => 300.0,
            HexState::B0 => 180.0,    // Negated
            HexState::B60 => 240.0,
            HexState::B120 => 300.0,
            HexState::B180 => 0.0,
            HexState::B240 => 60.0,
            HexState::B300 => 120.0,
        }
    }

    /// Get phase angle in radians
    pub fn phase(&self) -> f64 {
        self.angle() * std::f64::consts::PI / 180.0
    }

    /// Check if on positive face
    pub fn is_positive(&self) -> bool {
        matches!(self, HexState::A0 | HexState::A60 | HexState::A120 |
                       HexState::A180 | HexState::A240 | HexState::A300)
    }

    /// Check if on negative/dual face
    pub fn is_negative(&self) -> bool {
        !self.is_positive()
    }

    /// Rotate by +60 degrees (clockwise in angle space)
    pub fn rotate_cw(&self) -> HexState {
        match self {
            HexState::A0 => HexState::A300,
            HexState::A60 => HexState::A0,
            HexState::A120 => HexState::A60,
            HexState::A180 => HexState::A120,
            HexState::A240 => HexState::A180,
            HexState::A300 => HexState::A240,
            HexState::B0 => HexState::B300,
            HexState::B60 => HexState::B0,
            HexState::B120 => HexState::B60,
            HexState::B180 => HexState::B120,
            HexState::B240 => HexState::B180,
            HexState::B300 => HexState::B240,
        }
    }

    /// Rotate by -60 degrees (counter-clockwise)
    pub fn rotate_ccw(&self) -> HexState {
        match self {
            HexState::A0 => HexState::A60,
            HexState::A60 => HexState::A120,
            HexState::A120 => HexState::A180,
            HexState::A180 => HexState::A240,
            HexState::A240 => HexState::A300,
            HexState::A300 => HexState::A0,
            HexState::B0 => HexState::B60,
            HexState::B60 => HexState::B120,
            HexState::B120 => HexState::B180,
            HexState::B180 => HexState::B240,
            HexState::B240 => HexState::B300,
            HexState::B300 => HexState::B0,
        }
    }

    /// Negate (flip to dual face)
    pub fn negate(&self) -> HexState {
        match self {
            HexState::A0 => HexState::B0,
            HexState::A60 => HexState::B60,
            HexState::A120 => HexState::B120,
            HexState::A180 => HexState::B180,
            HexState::A240 => HexState::B240,
            HexState::A300 => HexState::B300,
            HexState::B0 => HexState::A0,
            HexState::B60 => HexState::A60,
            HexState::B120 => HexState::A120,
            HexState::B180 => HexState::A180,
            HexState::B240 => HexState::A240,
            HexState::B300 => HexState::A300,
        }
    }

    /// Get the complementary state (add 180°)
    pub fn complement(&self) -> HexState {
        match self {
            HexState::A0 => HexState::A180,
            HexState::A60 => HexState::A240,
            HexState::A120 => HexState::A300,
            HexState::A180 => HexState::A0,
            HexState::A240 => HexState::A60,
            HexState::A300 => HexState::A120,
            HexState::B0 => HexState::B180,
            HexState::B60 => HexState::B240,
            HexState::B120 => HexState::B300,
            HexState::B180 => HexState::B0,
            HexState::B240 => HexState::B60,
            HexState::B300 => HexState::B120,
        }
    }

    /// Convert to complex amplitude (unit circle)
    pub fn to_complex(&self) -> (f64, f64) {
        let angle = self.phase();
        (angle.cos(), angle.sin())
    }

    /// Distance to another state (in angular steps)
    pub fn distance(&self, other: &HexState) -> u32 {
        let a = self.angle();
        let b = other.angle();

        let diff = (a - b).abs();
        let dist = diff.min(360.0 - diff);

        (dist / 60.0) as u32
    }

    /// Check if orthogonal (180° apart)
    pub fn is_orthogonal(&self, other: &HexState) -> bool {
        self.distance(other) == 3
    }

    /// Check if adjacent (60° apart)
    pub fn is_adjacent(&self, other: &HexState) -> bool {
        self.distance(other) == 1
    }
}

impl fmt::Display for HexState {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            HexState::A0 => write!(f, "A0°"),
            HexState::A60 => write!(f, "A60°"),
            HexState::A120 => write!(f, "A120°"),
            HexState::A180 => write!(f, "A180°"),
            HexState::A240 => write!(f, "A240°"),
            HexState::A300 => write!(f, "A300°"),
            HexState::B0 => write!(f, "B0°"),
            HexState::B60 => write!(f, "B60°"),
            HexState::B120 => write!(f, "B120°"),
            HexState::B180 => write!(f, "B180°"),
            HexState::B240 => write!(f, "B240°"),
            HexState::B300 => write!(f, "B300°"),
        }
    }
}

/// A register of hexagonal states
#[derive(Debug, Clone)]
pub struct HexRegister {
    states: Vec<HexState>,
}

impl HexRegister {
    pub fn new(n: usize) -> Self {
        Self {
            states: vec![HexState::A0; n],
        }
    }

    pub fn len(&self) -> usize {
        self.states.len()
    }

    pub fn is_empty(&self) -> bool {
        self.states.is_empty()
    }

    pub fn get(&self, i: usize) -> Option<HexState> {
        self.states.get(i).copied()
    }

    pub fn set(&mut self, i: usize, state: HexState) {
        if i < self.states.len() {
            self.states[i] = state;
        }
    }

    /// Rotate a specific qubit clockwise
    pub fn rotate_cw(&mut self, i: usize) {
        if let Some(s) = self.states.get_mut(i) {
            *s = s.rotate_cw();
        }
    }

    /// Rotate a specific qubit counter-clockwise
    pub fn rotate_ccw(&mut self, i: usize) {
        if let Some(s) = self.states.get_mut(i) {
            *s = s.rotate_ccw();
        }
    }

    /// Negate (flip to dual face)
    pub fn negate(&mut self, i: usize) {
        if let Some(s) = self.states.get_mut(i) {
            *s = s.negate();
        }
    }
}

/// Hexagonal Hadamard - creates superposition between faces
pub struct HexHadamard;

impl HexHadamard {
    /// Create superposition between positive and negative faces
    pub fn apply(state: HexState) -> Vec<(HexState, f64)> {
        let negated = state.negate();
        vec![
            (state, 1.0 / 2.0_f64.sqrt()),
            (negated, 1.0 / 2.0_f64.sqrt()),
        ]
    }
}

/// Hexagonal Phase gate
pub struct HexPhase {
    angle: f64,
}

impl HexPhase {
    pub fn new(degrees: f64) -> Self {
        Self { angle: degrees }
    }

    pub fn S() -> Self { Self { angle: 60.0 } }
    pub fn T() -> Self { Self { angle: 30.0 } }

    pub fn apply(&self, state: HexState) -> HexState {
        let current_angle = state.angle();
        let new_angle = (current_angle + self.angle) % 360.0;

        // Find closest state
        HexState::all()
            .into_iter()
            .min_by(|a, b| {
                let dist_a = (a.angle() - new_angle).abs().min(360.0 - (a.angle() - new_angle).abs());
                let dist_b = (b.angle() - new_angle).abs().min(360.0 - (b.angle() - new_angle).abs());
                dist_a.partial_cmp(&dist_b).unwrap()
            })
            .unwrap_or(state)
    }
}

/// Time evolution as branching
/// From VOID-DIFFERENTIATION.md: "Each increment of the timer doesn't overwrite
/// previous states but adds new ones, allowing different branches to coexist independently"
#[derive(Debug, Clone)]
pub struct BranchingState {
    pub current: Vec<HexState>,
    pub history: Vec<Vec<HexState>>,
}

impl BranchingState {
    pub fn new(initial: Vec<HexState>) -> Self {
        Self {
            current: initial.clone(),
            history: vec![initial],
        }
    }

    /// Branch (fork the timeline)
    pub fn branch(&self) -> Self {
        Self {
            current: self.current.clone(),
            history: self.history.clone(),
        }
    }

    /// Advance time (add new state without erasing history)
    pub fn tick(&mut self, new_states: Vec<HexState>) {
        self.current = new_states;
        self.history.push(self.current.clone());
    }

    /// Number of branches (paths through history)
    pub fn branch_factor(&self) -> usize {
        // Each history entry could potentially branch
        // Simplified: count unique paths
        self.history.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_12_states_exist() {
        let states = HexState::all();
        assert_eq!(states.len(), 12);
    }

    #[test]
    fn test_positive_negative_split() {
        let states = HexState::all();
        let positives: Vec<_> = states.iter().filter(|s| s.is_positive()).collect();
        let negatives: Vec<_> = states.iter().filter(|s| s.is_negative()).collect();

        assert_eq!(positives.len(), 6);
        assert_eq!(negatives.len(), 6);
    }

    #[test]
    fn test_rotation_cw_6_times() {
        let mut state = HexState::A0;
        for _ in 0..6 {
            state = state.rotate_cw();
        }
        assert_eq!(state, HexState::A0);
    }

    #[test]
    fn test_rotation_ccw_6_times() {
        let mut state = HexState::A0;
        for _ in 0..6 {
            state = state.rotate_ccw();
        }
        assert_eq!(state, HexState::A0);
    }

    #[test]
    fn test_negate_twice_identity() {
        let state = HexState::A60;
        assert_eq!(state.negate().negate(), state);
    }

    #[test]
    fn test_complement_180() {
        let state = HexState::A0;
        assert_eq!(state.complement(), HexState::A180);
    }

    #[test]
    fn test_orthogonal_states() {
        assert!(HexState::A0.is_orthogonal(&HexState::A180));
        assert!(HexState::A60.is_orthogonal(&HexState::A240));
        assert!(HexState::B0.is_orthogonal(&HexState::B180));
    }

    #[test]
    fn test_adjacent_states() {
        assert!(HexState::A0.is_adjacent(&HexState::A60));
        assert!(HexState::A0.is_adjacent(&HexState::A300));
    }

    #[test]
    fn test_distance_wraps_around() {
        assert_eq!(HexState::A0.distance(&HexState::A300), 1);
    }

    #[test]
    fn test_hex_register() {
        let mut reg = HexRegister::new(3);
        assert_eq!(reg.len(), 3);
        reg.rotate_cw(0);
        assert_eq!(reg.get(0), Some(HexState::A300));
    }

    #[test]
    fn test_branching_state() {
        let state = BranchingState::new(vec![HexState::A0, HexState::A60]);
        assert_eq!(state.branch_factor(), 1);

        // Branch creates independent copy
        let branch = state.branch();
        assert_eq!(branch.current, state.current);
    }

    #[test]
    fn test_complex_amplitudes() {
        // States should lie on unit circle
        for s in HexState::all() {
            let (re, im) = s.to_complex();
            let mag = (re * re + im * im).sqrt();
            assert!((mag - 1.0).abs() < 1e-10);
        }
    }
}
