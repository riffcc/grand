// GUTOE QM — Clifford Hilbert space: 16-dimensional complex wave function per site
// Copyright (C) 2026 Riff Labs, AGPL-3.0-or-later
//
// The Clifford algebra Cl(1,3) is 16-dimensional. Each lattice site holds a
// complex amplitude for each of the 16 basis elements {1, γ⁰, γ¹, ..., γ⁰¹²³}.
//
// Classical simulation: one u8 per site (which basis element it IS)
// Quantum simulation:   [Complex; 16] per site (amplitude for EACH basis element)
//
// The Born rule: |⟨s|ψ⟩|² = amplitude[s-1].norm_sqr() = probability of measuring s

use num_complex::Complex64;

// ── Types ─────────────────────────────────────────────────────────────────────

/// Complex amplitude for one site's 16-dimensional Clifford Hilbert space.
/// amp[s] = complex amplitude for basis state s+1 (s ∈ 0..15, state = s+1).
/// amp[0] = VOID amplitude (state 0 is the vacuum).
///
/// The inner product: ⟨φ|ψ⟩ = Σ_s φ[s]* ψ[s]
/// The norm: ‖ψ‖² = Σ_s |ψ[s]|² = 1 for a normalized state
pub type SiteAmp = [Complex64; 17]; // indices 0..16 (s=0 VOID, s=1..16 Clifford)

/// Quantum wave function for a single SPATIAL degree of freedom (e.g., lepton position).
/// psi[i] = complex amplitude for finding the particle at lattice site i.
/// Σ_i |psi[i]|² = 1 after normalization.
///
/// This is a SPATIAL superposition — a particle exists at multiple sites simultaneously.
/// Different from SiteAmp which is an INTERNAL superposition over Clifford basis elements.
pub type SpatialPsi = Vec<Complex64>;

// ── Site amplitude operations ──────────────────────────────────────────────────

/// Create a site amplitude in a pure Clifford basis state s (s ∈ 0..16).
pub fn pure_state(s: usize) -> SiteAmp {
    assert!(s <= 16, "state s must be 0..16");
    let mut amp = [Complex64::new(0.0, 0.0); 17];
    amp[s] = Complex64::new(1.0, 0.0);
    amp
}

/// Inner product ⟨φ|ψ⟩ = Σ_s φ[s]* ψ[s]
pub fn inner(phi: &SiteAmp, psi: &SiteAmp) -> Complex64 {
    phi.iter().zip(psi.iter()).map(|(a, b)| a.conj() * b).sum()
}

/// Norm squared ‖ψ‖² = ⟨ψ|ψ⟩
pub fn norm_sq(amp: &SiteAmp) -> f64 {
    amp.iter().map(|a| a.norm_sqr()).sum()
}

/// Normalize in place. Returns false if the state is the zero vector.
pub fn normalize_site(amp: &mut SiteAmp) -> bool {
    let n = norm_sq(amp).sqrt();
    if n < 1e-15 {
        return false;
    }
    for a in amp.iter_mut() {
        *a /= n;
    }
    true
}

/// Born rule: probability of measuring state s (0..16).
pub fn born_prob(amp: &SiteAmp, s: usize) -> f64 {
    amp[s].norm_sqr()
}

// ── Spatial wave function operations ──────────────────────────────────────────

/// Initialize lepton at a single site (pure position eigenstate).
pub fn init_at(site: usize, n_sites: usize) -> SpatialPsi {
    let mut psi = vec![Complex64::new(0.0, 0.0); n_sites];
    psi[site] = Complex64::new(1.0, 0.0);
    psi
}

/// Initialize equal superposition over two sites (coherent 50-50 split).
pub fn init_superposition(site1: usize, site2: usize, n_sites: usize) -> SpatialPsi {
    let mut psi = vec![Complex64::new(0.0, 0.0); n_sites];
    let a = 1.0 / 2.0_f64.sqrt();
    psi[site1] = Complex64::new(a, 0.0);
    psi[site2] = Complex64::new(a, 0.0);
    psi
}

/// Norm squared of a spatial wave function. Should be 1.0 after normalize.
pub fn spatial_norm_sq(psi: &SpatialPsi) -> f64 {
    psi.iter().map(|a| a.norm_sqr()).sum()
}

/// Normalize the spatial wave function in place.
pub fn spatial_normalize(psi: &mut SpatialPsi) {
    let n = spatial_norm_sq(psi).sqrt();
    if n > 1e-15 {
        for a in psi.iter_mut() {
            *a /= n;
        }
    }
}

/// Born rule probability at site i: |ψᵢ|²
pub fn spatial_prob(psi: &SpatialPsi, site: usize) -> f64 {
    psi[site].norm_sqr()
}

/// Born rule measurement: sample a site with probability |ψᵢ|², then collapse.
pub fn measure_and_collapse(psi: &mut SpatialPsi, rng: &mut impl rand::Rng) -> usize {
    let total: f64 = spatial_norm_sq(psi);
    let r = rng.gen::<f64>() * total;
    let mut cumsum = 0.0;
    let mut site = psi.len() - 1;
    for (i, a) in psi.iter().enumerate() {
        cumsum += a.norm_sqr();
        if cumsum >= r {
            site = i;
            break;
        }
    }
    // Collapse: project onto |site⟩
    for (i, a) in psi.iter_mut().enumerate() {
        *a = if i == site {
            Complex64::new(1.0, 0.0)
        } else {
            Complex64::new(0.0, 0.0)
        };
    }
    site
}
