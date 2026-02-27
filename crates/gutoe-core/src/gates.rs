/*!
 * GUTOE Core - Tripartite Quantum Gates
 * Copyright (C) 2026  Wings
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

//! Tripartite quantum gates
//!
//! Implements the fundamental gate operations for the 4-state GUTOE system

use crate::states::TriState;
use num_complex::Complex64;
use std::f64::consts::PI;

/// TripartiteNOT gate - cycles through SINE → COSINE → TANGENT → SINE
/// This is analogous to Pauli-X but for the tripartite system
pub struct TripartiteNot;

impl TripartiteNot {
    pub fn apply(state: TriState) -> TriState {
        state.cycle()
    }

    pub fn matrix() -> [[Complex64; 4]; 4] {
        // Transformation matrix in order: VOID, COSINE, SINE, TANGENT
        [
            [
                Complex64::new(1.0, 0.0),
                Complex64::new(0.0, 0.0),
                Complex64::new(0.0, 0.0),
                Complex64::new(0.0, 0.0),
            ],
            [
                Complex64::new(0.0, 0.0),
                Complex64::new(0.0, 0.0),
                Complex64::new(1.0, 0.0),
                Complex64::new(0.0, 0.0),
            ],
            [
                Complex64::new(0.0, 0.0),
                Complex64::new(0.0, 0.0),
                Complex64::new(0.0, 0.0),
                Complex64::new(1.0, 0.0),
            ],
            [
                Complex64::new(0.0, 0.0),
                Complex64::new(1.0, 0.0),
                Complex64::new(0.0, 0.0),
                Complex64::new(0.0, 0.0),
            ],
        ]
    }
}

/// Hadamard-like gate - creates superposition from basis states
pub struct TripartiteHadamard;

impl TripartiteHadamard {
    /// Create superposition from COSINE or SINE
    pub fn apply(state: TriState) -> TriState {
        match state {
            TriState::COSINE | TriState::SINE => TriState::TANGENT,
            TriState::TANGENT => {
                // Hadamard on superposition - this should be identity
                TriState::TANGENT
            }
            TriState::VOID => TriState::VOID,
        }
    }
}

/// Phase gate - rotates phase of TANGENT state
pub struct PhaseGate {
    pub phase: f64,
}

impl PhaseGate {
    pub fn new(phase: f64) -> Self {
        Self { phase }
    }

    pub fn S() -> Self {
        Self { phase: PI / 2.0 }
    }

    pub fn T() -> Self {
        Self { phase: PI / 4.0 }
    }
}

/// Controlled TripartiteNOT - analogous to CNOT
pub struct CTripartiteNot;

impl CTripartiteNot {
    /// Apply controlled-NOT: if control is SINE, cycle target
    pub fn apply(control: TriState, target: TriState) -> TriState {
        if control == TriState::SINE {
            target.cycle()
        } else {
            target
        }
    }
}

/// Veracity gate - measures/collapses relationship between two states
pub struct VeracityGate {
    pub threshold: f64,
}

impl VeracityGate {
    pub fn new(threshold: f64) -> Self {
        Self { threshold }
    }

    /// Calculate veracity (coherence strength) between two states.
    ///
    /// Physical basis (hexagonal spectral geometry):
    /// - Same state = full coherence = 1.0
    /// - SINE ↔ COSINE: 60° adjacent hex modes → √3/2 ≈ 0.866 (NOT orthogonal)
    /// - TANGENT with SINE/COSINE: 0.5 (partial — ratio/slope state)
    /// - VOID: 0.0 (no connections)
    pub fn measure(a: TriState, b: TriState) -> f64 {
        match (a, b) {
            (TriState::VOID, _) | (_, TriState::VOID) => 0.0,
            (TriState::COSINE, TriState::COSINE)
            | (TriState::SINE, TriState::SINE)
            | (TriState::TANGENT, TriState::TANGENT) => 1.0,
            (TriState::COSINE, TriState::SINE) | (TriState::SINE, TriState::COSINE) => {
                3.0_f64.sqrt() / 2.0
            }
            _ => 0.5, // TANGENT-SINE, TANGENT-COSINE
        }
    }
}

/// Three-qubit Toffoli-like gate (C2NOT)
pub struct Toffoli3;

impl Toffoli3 {
    /// Apply Toffoli: cycle target only if both controls are SINE
    pub fn apply(c1: TriState, c2: TriState, target: TriState) -> TriState {
        if c1 == TriState::SINE && c2 == TriState::SINE {
            target.cycle()
        } else {
            target
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── TripartiteNot: basic cycle map ────────────────────────────────────────

    #[test]
    fn tripartite_not_cycles_correctly() {
        assert_eq!(TripartiteNot::apply(TriState::SINE), TriState::COSINE);
        assert_eq!(TripartiteNot::apply(TriState::COSINE), TriState::TANGENT);
        assert_eq!(TripartiteNot::apply(TriState::TANGENT), TriState::SINE);
        assert_eq!(TripartiteNot::apply(TriState::VOID), TriState::VOID);
    }

    // ── TripartiteNot: Z₃ structure, NOT Z₂ (not self-inverse) ──────────────

    #[test]
    fn tripartite_not_is_not_self_inverse() {
        // Classic Pauli-X: X² = I.  This gate is Z₃ so cycle² ≠ I.
        assert_ne!(
            TripartiteNot::apply(TripartiteNot::apply(TriState::COSINE)),
            TriState::COSINE,
            "TripartiteNot² = id on COSINE — that's Z₂ behavior, but cycle is Z₃"
        );
    }

    #[test]
    fn tripartite_not_has_order_3() {
        // Applying the gate 3 times must return the original state.
        for s in [TriState::SINE, TriState::COSINE, TriState::TANGENT] {
            let thrice = TripartiteNot::apply(TripartiteNot::apply(TripartiteNot::apply(s)));
            assert_eq!(
                thrice, s,
                "TripartiteNot³({s:?}) ≠ {s:?} — gate does not have order 3"
            );
        }
    }

    // ── CTripartiteNot: NOT self-inverse (Z₃, not Z₂) ────────────────────────

    #[test]
    fn c_tripartite_not_is_not_self_inverse() {
        // Classic CNOT: C² = I. CTripartiteNot uses cycle (Z₃), so C² ≠ I.
        let twice = CTripartiteNot::apply(
            TriState::SINE,
            CTripartiteNot::apply(TriState::SINE, TriState::COSINE),
        );
        assert_ne!(
            twice,
            TriState::COSINE,
            "CTripartiteNot² = identity on COSINE — requires Z₂, but cycle has order 3"
        );
    }

    #[test]
    fn c_tripartite_not_has_order_3_with_sine_control() {
        for target in [TriState::SINE, TriState::COSINE, TriState::TANGENT] {
            let t1 = CTripartiteNot::apply(TriState::SINE, target);
            let t2 = CTripartiteNot::apply(TriState::SINE, t1);
            let t3 = CTripartiteNot::apply(TriState::SINE, t2);
            assert_eq!(
                t3, target,
                "CTripartiteNot³({target:?}) ≠ {target:?} with SINE control — should have order 3"
            );
        }
    }

    #[test]
    fn c_tripartite_not_is_identity_when_control_not_sine() {
        for control in [TriState::COSINE, TriState::TANGENT, TriState::VOID] {
            for target in [TriState::SINE, TriState::COSINE, TriState::TANGENT] {
                assert_eq!(
                    CTripartiteNot::apply(control, target),
                    target,
                    "CTripartiteNot({control:?}, {target:?}) fired without SINE control"
                );
            }
        }
    }

    // ── Toffoli3: NOT self-inverse (same Z₃ reason) ───────────────────────────

    #[test]
    fn toffoli3_cycles_correctly_with_both_controls_sine() {
        assert_eq!(
            Toffoli3::apply(TriState::SINE, TriState::SINE, TriState::COSINE),
            TriState::TANGENT
        );
        assert_eq!(
            Toffoli3::apply(TriState::SINE, TriState::COSINE, TriState::COSINE),
            TriState::COSINE
        );
    }

    #[test]
    fn toffoli3_is_not_self_inverse() {
        let twice = Toffoli3::apply(
            TriState::SINE,
            TriState::SINE,
            Toffoli3::apply(TriState::SINE, TriState::SINE, TriState::COSINE),
        );
        assert_ne!(
            twice,
            TriState::COSINE,
            "Toffoli3² = identity on COSINE — requires Z₂, but cycle has order 3"
        );
    }

    #[test]
    fn toffoli3_has_order_3_when_both_controls_sine() {
        for target in [TriState::SINE, TriState::COSINE, TriState::TANGENT] {
            let t1 = Toffoli3::apply(TriState::SINE, TriState::SINE, target);
            let t2 = Toffoli3::apply(TriState::SINE, TriState::SINE, t1);
            let t3 = Toffoli3::apply(TriState::SINE, TriState::SINE, t2);
            assert_eq!(
                t3, target,
                "Toffoli3³({target:?}) ≠ {target:?} — should have order 3"
            );
        }
    }

    // ── TripartiteHadamard: NOT a valid quantum gate ───────────────────────────
    // A unitary gate must be bijective (injective).  This one collapses two
    // distinct inputs (COSINE, SINE) to the same output (TANGENT), so it
    // destroys information and cannot be unitary.

    #[test]
    fn hadamard_is_not_injective() {
        let h_cosine = TripartiteHadamard::apply(TriState::COSINE);
        let h_sine = TripartiteHadamard::apply(TriState::SINE);
        assert_eq!(
            h_cosine,
            TriState::TANGENT,
            "Hadamard(COSINE) should give TANGENT"
        );
        assert_eq!(
            h_sine,
            TriState::TANGENT,
            "Hadamard(SINE) should give TANGENT"
        );
        // Same output for two distinct inputs → not injective, not invertible, not unitary.
        assert_eq!(
            h_cosine, h_sine,
            "Hadamard maps COSINE and SINE to different states — expected both to TANGENT"
        );
    }

    #[test]
    fn hadamard_is_not_involutive() {
        // Real Hadamard: H² = I.
        // TripartiteHadamard: COSINE → TANGENT → TANGENT (stuck; never returns to COSINE).
        let twice = TripartiteHadamard::apply(TripartiteHadamard::apply(TriState::COSINE));
        assert_ne!(
            twice,
            TriState::COSINE,
            "TripartiteHadamard² = identity on COSINE — cannot be true if not injective"
        );
    }

    // ── VeracityGate vs quantum fidelity ─────────────────────────────────────
    // Standard amplitude fidelity: |⟨SINE|COSINE⟩|² = 0 (orthogonal amplitudes).
    // But veracity(SINE, COSINE) = √3/2 ≈ 0.866 (60° hexagonal adjacency).
    // Veracity is a hexagonal-geometric coherence measure, NOT standard fidelity.

    #[test]
    fn veracity_tangent_self_is_1() {
        // TANGENT with itself = full coherence, same as any non-VOID state.
        let v = VeracityGate::measure(TriState::TANGENT, TriState::TANGENT);
        assert!(
            (v - 1.0).abs() < 1e-12,
            "veracity(TANGENT, TANGENT) = {v:.4}, expected 1.0"
        );
    }

    // ── Real Python CNOT vs Rust port ─────────────────────────────────────────
    // Python TripartiteCNOT: SINE↔COSINE swap when control=SINE (self-inverse, Z₂).
    // Rust CTripartiteNot:   cycle() when control=SINE (Z₃, NOT self-inverse).

    /// Mirrors Python's tripartite_cnot.py SINE↔COSINE swap semantics.
    fn real_cnot(control: TriState, target: TriState) -> TriState {
        if control == TriState::SINE {
            match target {
                TriState::SINE => TriState::COSINE,
                TriState::COSINE => TriState::SINE,
                other => other,
            }
        } else {
            target
        }
    }

    #[test]
    fn real_python_cnot_is_self_inverse() {
        // Python CNOT is a clean Z₂ involution — two applications = identity.
        for target in [
            TriState::SINE,
            TriState::COSINE,
            TriState::TANGENT,
            TriState::VOID,
        ] {
            let twice = real_cnot(TriState::SINE, real_cnot(TriState::SINE, target));
            assert_eq!(
                twice, target,
                "real_cnot²({target:?}) = {twice:?} ≠ {target:?} — should be self-inverse"
            );
        }
    }

    #[test]
    fn rust_port_cnot_is_not_self_inverse_on_cosine() {
        // Rust CTripartiteNot uses cycle(), giving Z₃ — NOT self-inverse.
        // Python TripartiteCNOT swaps, giving Z₂ — IS self-inverse.
        let rust_twice = CTripartiteNot::apply(
            TriState::SINE,
            CTripartiteNot::apply(TriState::SINE, TriState::COSINE),
        );
        let real_twice = real_cnot(TriState::SINE, real_cnot(TriState::SINE, TriState::COSINE));

        assert_ne!(
            rust_twice,
            TriState::COSINE,
            "Rust CTripartiteNot² = identity on COSINE — that's Z₂ behavior, but cycle is Z₃"
        );
        assert_eq!(
            real_twice,
            TriState::COSINE,
            "Python real_cnot² ≠ identity on COSINE — should be self-inverse Z₂ swap"
        );
    }

    #[test]
    fn rust_port_and_python_cnot_disagree_on_cosine_with_sine_control() {
        // Rust: COSINE.cycle() = TANGENT
        // Python: COSINE swapped to SINE
        let rust_result = CTripartiteNot::apply(TriState::SINE, TriState::COSINE);
        let real_result = real_cnot(TriState::SINE, TriState::COSINE);
        assert_eq!(
            rust_result,
            TriState::TANGENT,
            "Rust CTripartiteNot(SINE, COSINE) should be TANGENT (cycle)"
        );
        assert_eq!(
            real_result,
            TriState::SINE,
            "Python real_cnot(SINE, COSINE) should be SINE (swap)"
        );
        assert_ne!(
            rust_result, real_result,
            "Rust port and Python cnot should disagree on (SINE, COSINE)"
        );
    }

    #[test]
    fn veracity_diverges_from_amplitude_fidelity_for_sine_cosine() {
        // SINE and COSINE have orthogonal amplitudes: |⟨SINE|COSINE⟩|² = 0.
        // But veracity(SINE, COSINE) = √3/2 ≈ 0.866 (60° hexagonal adjacency).
        // Veracity is NOT the standard inner-product fidelity.
        let (a0_s, a1_s) = TriState::SINE.to_amplitude();
        let (a0_c, a1_c) = TriState::COSINE.to_amplitude();
        // Standard inner-product fidelity: |⟨SINE|COSINE⟩|² = |a0_s*·a0_c + a1_s*·a1_c|²
        let inner = a0_s.conj() * a0_c + a1_s.conj() * a1_c;
        let amplitude_fidelity = inner.norm_sqr();
        let veracity = VeracityGate::measure(TriState::SINE, TriState::COSINE);

        assert!((amplitude_fidelity - 0.0).abs() < 1e-12,
            "amplitude fidelity of SINE vs COSINE should be 0.0 (orthogonal), got {amplitude_fidelity:.6}");
        assert!(
            (veracity - 3.0_f64.sqrt() / 2.0).abs() < 1e-12,
            "veracity(SINE, COSINE) should be √3/2, got {veracity:.6}"
        );
        assert!(
            veracity > 0.8,
            "veracity ({veracity:.4}) ≈ amplitude fidelity ({amplitude_fidelity:.4}): \
             they should differ significantly — veracity uses hex geometry, not inner product"
        );
    }
}
