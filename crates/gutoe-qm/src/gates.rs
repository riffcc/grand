// GUTOE QM — Unitary gates for the quantum Clifford lattice
// Copyright (C) 2026 Riff Labs, AGPL-3.0-or-later
//
// All gates preserve the norm: ‖U|ψ⟩‖ = ‖|ψ⟩‖ (unitarity).
// This is the fundamental requirement: quantum mechanics = unitary dynamics.
//
// Gates implemented:
//   hop_unitary    — beam-splitter between two spatial sites (2×2 unitary)
//   em_phase       — U(1) phase rotation from EM field (1×1 unitary)
//   z3_perm        — Z₃ permutation of Clifford basis states (16×16 unitary)
//   clifford_phase — diagonal phase from the Clifford metric (16×16 unitary)

use crate::hilbert::{SiteAmp, SpatialPsi};
use num_complex::Complex64;

// ── Spatial gates ──────────────────────────────────────────────────────────────

/// Beam-splitter between spatial sites i and j with mixing angle θ.
///
/// |ψᵢ'⟩ = cos(θ) |ψᵢ⟩ + i sin(θ) |ψⱼ⟩
/// |ψⱼ'⟩ = i sin(θ) |ψᵢ⟩ + cos(θ) |ψⱼ⟩
///
/// This is an exact 2×2 unitary. At θ = π/4: 50-50 split (Hadamard-like).
/// At θ = π/2: full swap with phase.
///
/// Physical meaning: a lepton can tunnel between two sites with amplitude
/// proportional to sin(θ), without committing to either path (superposition).
/// This is the quantum replacement for the classical "hop with probability p".
pub fn hop_unitary(psi: &mut SpatialPsi, i: usize, j: usize, theta: f64) {
    let ai = psi[i];
    let aj = psi[j];
    let c = theta.cos();
    let s = theta.sin();
    psi[i] = ai * c + aj * Complex64::new(0.0, s);
    psi[j] = ai * Complex64::new(0.0, s) + aj * c;
}

/// U(1) phase rotation from the electromagnetic field.
///
/// |ψᵢ'⟩ = e^{i q φ} |ψᵢ⟩
///
/// where q = charge and φ = Coulomb potential at site i.
///
/// Physical meaning: the lepton (charge q = −1) accumulates phase −φ at each site.
/// This is the Aharonov-Bohm effect: the PHASE (not the force) carries the EM information.
/// A path through high-φ regions accumulates MORE phase than a path through low-φ.
/// When two paths meet, their phase difference causes interference.
///
/// This replaces the classical "hop toward max-φ neighbor" with quantum phase modulation.
/// The interference pattern of the quantum lepton NATURALLY concentrates the wave
/// function in the proton's Coulomb shell — without any gradient-following.
pub fn em_phase(psi: &mut SpatialPsi, site: usize, phi: f64, charge: f64) {
    psi[site] *= Complex64::from_polar(1.0, charge * phi);
}

/// Apply EM phase at every site simultaneously.
///
/// This is the full quantum EM coupling: every lattice site gets a phase
/// proportional to its Coulomb potential × lepton charge.
/// After this gate, interference of the wave function encodes the EM field.
pub fn em_phase_all(psi: &mut SpatialPsi, phi_field: &[f64], charge: f64) {
    assert_eq!(psi.len(), phi_field.len());
    for (i, amp) in psi.iter_mut().enumerate() {
        *amp *= Complex64::from_polar(1.0, charge * phi_field[i]);
    }
}

// ── Internal (Clifford) gates ──────────────────────────────────────────────────

/// Z₃ permutation gate on the 16-dimensional Clifford Hilbert space.
///
/// This is a PERMUTATION matrix — automatically unitary.
/// Fixed points: {γ⁰ (s=2), scalar (s=1), γ¹²³ (s=15), pseudoscalar (s=16)}.
/// Quark orbits: {γ¹,γ²,γ³} and 4 other grade-2/3 triplets cycle.
///
/// In the quantum theory, this gate acts on the INTERNAL Clifford degrees
/// of freedom at a single site, rotating the color state of a quark without
/// changing its spatial position.
pub fn z3_gate(amp: &mut SiteAmp) {
    // Z₃ table: Z3_TABLE[s] gives the image of state s under Z₃
    // From gutoe_em::sim::Z3_TABLE (0-indexed here)
    const Z3: [usize; 17] = [
        0,  // VOID → VOID
        1,  // s=1 (scalar) → 1 (fixed)
        2,  // s=2 (γ⁰, lepton) → 2 (FIXED POINT)
        5,  // s=3 (γ¹) → s=5 (γ²)
        6,  // s=4 (γ⁰¹) → s=6 (γ⁰²)
        9,  // s=5 (γ²) → s=9 (γ³)
        10, // s=6 (γ⁰²) → s=10 (γ⁰³)
        13, // s=7 (γ¹²) → s=13 (γ²³)
        14, // s=8 (γ⁰¹²) → s=14 (γ⁰²³)
        3,  // s=9 (γ³) → s=3 (γ¹)
        4,  // s=10 (γ⁰³) → s=4 (γ⁰¹)
        7,  // s=11 (γ¹³) → s=7 (γ¹²)
        8,  // s=12 (γ⁰¹³) → s=8 (γ⁰¹²)
        11, // s=13 (γ²³) → s=11 (γ¹³)
        12, // s=14 (γ⁰²³) → s=12 (γ⁰¹³)
        15, // s=15 (γ¹²³) → 15 (fixed, all spatial bits = 1)
        16, // s=16 (γ⁰¹²³) → 16 (fixed, pseudoscalar)
    ];
    let old = *amp;
    for s in 0..=16 {
        amp[Z3[s]] = old[s];
    }
}

/// Diagonal phase gate from the Clifford metric.
///
/// Assigns a phase e^{iθ_s} to each basis state s based on its grade
/// and the Minkowski signature (−,+,+,+):
///
///   θ(grade-0) = 0    (scalar, no phase)
///   θ(grade-1) = −π/4 × n_timelike_bits  (timelike directions get phase −π/4)
///   θ(grade-2) = depends on temporal/spatial content
///   ...
///
/// This implements the time evolution from the kinetic term of the Dirac Hamiltonian.
/// In the continuum: e^{−iH_0 t} where H_0 = iγ^μ∂_μ
pub fn clifford_phase(amp: &mut SiteAmp, dt: f64) {
    // Phase for each basis state based on its grade and Minkowski signature
    // γ⁰ squares to -1 (timelike), γⁱ squares to +1 (spacelike)
    // The phase is: Q(state) × dt where Q is the quadratic form
    // Q(grade-0) = 0 (scalar has no kinetic energy)
    // Q(γ⁰) = -1 (timelike: imaginary mass-like term)
    // Q(γⁱ) = +1 (spacelike: real kinetic term)
    // Q(γ⁰ⁱ) = -1 (one timelike index)
    // etc.

    // For state s, mi = s - 1, count of timelike bits = (mi & 1)
    // Q(s) = (number of spacelike bits) - (number of timelike bits) * 1
    for s in 1..=16usize {
        let mi = s - 1;
        let n_timelike = mi & 1; // bit 0 = γ⁰ (timelike)
        let n_spacelike = ((mi >> 1) & 1) + ((mi >> 2) & 1) + ((mi >> 3) & 1);
        let q = n_spacelike as f64 - n_timelike as f64;
        amp[s] *= Complex64::from_polar(1.0, q * dt);
    }
}

// ── Time evolution ─────────────────────────────────────────────────────────────

/// One quantum time step for a lepton in the EM field.
///
/// Protocol:
///   1. Apply EM phase at current position (Aharonov-Bohm accumulation)
///   2. Apply beam-splitter hops to all neighbors (tunneling)
///   3. Renormalize (numerical stability)
///
/// This is the quantum replacement for the classical "hop toward max-φ" rule.
/// The lepton doesn't "decide" to move — it's everywhere simultaneously.
/// Interference between paths concentrates the wave function in the proton shell.
pub fn quantum_lepton_step(
    psi: &mut SpatialPsi,
    phi_field: &[f64],
    nbr_lists: &[Vec<usize>],
    theta: f64,
    charge: f64,
) {
    let n = psi.len();

    // 1. EM phase accumulation (Aharonov-Bohm effect)
    em_phase_all(psi, phi_field, charge);

    // 2. Hop unitaries to all neighbors
    // Apply beam-splitters site by site (Strang splitting approximation)
    for i in 0..n {
        for &j in &nbr_lists[i] {
            if i < j {
                // Apply once per pair (i < j ensures no double-counting)
                hop_unitary(psi, i, j, theta);
            }
        }
    }

    // 3. Renormalize (floating-point drift correction)
    crate::hilbert::spatial_normalize(psi);
}
