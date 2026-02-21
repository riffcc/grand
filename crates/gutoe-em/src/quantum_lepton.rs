// GUTOE EM — Quantum lepton: Schrödinger equation on the hex lattice
// Copyright (C) 2026 Riff Labs, AGPL-3.0-or-later
//
// Replaces the classical stochastic lepton hop with quantum mechanics.
//
// Classical Phase 2:  lepton at one site, hops to max-φ neighbour with prob em_prob
// Quantum Phase 2:    lepton as SpatialPsi = Vec<Complex>, evolves under H = -∇²_hex + V
//
// The Hamiltonian on the hex lattice (per-layer, intra-layer neighbours):
//   (Hψ)[i] = -(∇²_hex ψ)[i] + V[i] ψ[i]
//           = ψ[i] - (1/k) Σ_{j∈nbrs} ψ[j]   +   q × φ[i] × ψ[i]
//
// where k=6 (hex coordination), q=-1 (lepton charge), φ = Coulomb potential.
//
// Ground state: imaginary time evolution dψ/dτ = -H ψ → normalise → repeat.
// Converges to the lowest eigenstate of H (the bound state, if one exists).
//
// Binding energy: E = ⟨ψ|H|ψ⟩.  E < 0 → bound (hydrogen). E ≥ 0 → unbound.

use num_complex::Complex64;

use crate::config::LatticeConfig;
use crate::geometry::{mesh_neighbours, site_coords};

/// Quantum lepton wave function: complex amplitude at each lattice site.
/// Σ_i |psi[i]|² = 1 after normalization.
pub type LeptonPsi = Vec<Complex64>;

// ── Hamiltonian ────────────────────────────────────────────────────────────────

/// Apply the hex-lattice Schrödinger Hamiltonian: H = -∇²_hex + V.
///
/// (Hψ)[i] = kinetic[i] + potential[i]
///
/// kinetic[i]   = ψ[i] − (1/k) Σ_{j∈nbrs(i)} ψ[j]   [discrete Laplacian]
/// potential[i] = q × φ[i] × ψ[i]                    [Aharonov-Bohm potential]
///
/// Note: only intra-layer neighbours (mesh_neighbours never crosses layers),
/// so this is the 2D Schrödinger equation on a single hex layer.
pub fn apply_hamiltonian(
    psi: &LeptonPsi,
    phi: &[f64],
    cfg: &LatticeConfig,
    charge: f64,
) -> LeptonPsi {
    let n = psi.len();
    let mut h_psi = vec![Complex64::new(0.0, 0.0); n];

    for site in 0..n {
        let (r, c, z) = site_coords(site, cfg);
        let nbrs = mesh_neighbours(r, c, z, cfg);
        let k = nbrs.len() as f64;

        // Discrete Laplacian: ψ[i] − mean(ψ[nbrs])
        let nbr_sum: Complex64 = nbrs.iter().map(|&j| psi[j]).sum();
        let kinetic = psi[site] - nbr_sum / k;

        // Coulomb potential: V[i] = q × φ[i]
        let potential = psi[site] * Complex64::new(charge * phi[site], 0.0);

        h_psi[site] = kinetic + potential;
    }

    h_psi
}

/// Expected energy: ⟨ψ|H|ψ⟩ = Σ_i ψ[i]* (Hψ)[i].
/// Returns (total, kinetic, potential) in lattice units.
/// E < 0 → bound state (hydrogen). E ≥ 0 → unbound (scattering state).
pub fn expected_energy(
    psi: &LeptonPsi,
    phi: &[f64],
    cfg: &LatticeConfig,
    charge: f64,
) -> (f64, f64, f64) {
    let n = psi.len();
    let mut kinetic = 0.0;
    let mut potential = 0.0;

    for site in 0..n {
        let (r, c, z) = site_coords(site, cfg);
        let nbrs = mesh_neighbours(r, c, z, cfg);
        let k = nbrs.len() as f64;

        let nbr_sum: Complex64 = nbrs.iter().map(|&j| psi[j]).sum();
        let kin_term = psi[site] - nbr_sum / k;
        let pot_term = charge * phi[site] * psi[site];

        kinetic   += (psi[site].conj() * kin_term).re;
        potential += (psi[site].conj() * pot_term).re;
    }

    (kinetic + potential, kinetic, potential)
}

// ── Ground state via imaginary time evolution ──────────────────────────────────

/// One imaginary-time step: ψ → (1 − δτ H) ψ, then normalise.
///
/// Converges to the ground state: the component of ψ along each eigenstate
/// e_n decays as exp(−E_n δτ). Lowest E_n survives longest.
pub fn imaginary_time_step(
    psi: &mut LeptonPsi,
    phi: &[f64],
    cfg: &LatticeConfig,
    charge: f64,
    dtau: f64,
) {
    let h_psi = apply_hamiltonian(psi, phi, cfg, charge);
    for (i, h_i) in h_psi.iter().enumerate() {
        psi[i] -= *h_i * dtau;
    }
    normalise(psi);
}

fn normalise(psi: &mut LeptonPsi) {
    let norm: f64 = psi.iter().map(|a| a.norm_sqr()).sum::<f64>().sqrt();
    if norm > 1e-30 {
        for a in psi.iter_mut() {
            *a /= norm;
        }
    }
}

// ── Full ground state solver ───────────────────────────────────────────────────

/// Find the quantum lepton ground state in the proton's Coulomb field.
///
/// Protocol:
///   1. Start with ψ uniform over the proton shell (Gaussian centred on shell)
///   2. Evolve by imaginary time: ψ → e^{-H δτ} ψ (normalised)
///   3. Repeat until energy converges
///   4. Return (ground_state_psi, binding_energy)
///
/// If binding_energy < 0: the lepton is bound → quantum hydrogen formed.
/// If binding_energy ≥ 0: the Coulomb well is too shallow for this lattice.
pub fn quantum_hydrogen_ground_state(
    phi: &[f64],
    proton_shell_sites: &[usize],
    cfg: &LatticeConfig,
    n_iter: usize,
    dtau: f64,
) -> (LeptonPsi, f64, f64, f64) {
    let n = cfg.n_sites();

    // Initialise wave function uniform over the proton shell
    let mut psi = vec![Complex64::new(0.0, 0.0); n];
    if proton_shell_sites.is_empty() {
        // Fallback: uniform over layer 0
        let layer_sz = cfg.hex_rows * cfg.hex_cols;
        for i in 0..layer_sz {
            psi[i] = Complex64::new(1.0, 0.0);
        }
    } else {
        for &s in proton_shell_sites {
            psi[s] = Complex64::new(1.0, 0.0);
        }
    }
    normalise(&mut psi);

    let charge = -1.0_f64; // lepton

    // Imaginary time evolution
    for _ in 0..n_iter {
        imaginary_time_step(&mut psi, phi, cfg, charge, dtau);
    }

    let (e_total, e_kin, e_pot) = expected_energy(&psi, phi, cfg, charge);
    (psi, e_total, e_kin, e_pot)
}

// ── Hydrogen enrichment from Born rule ────────────────────────────────────────

/// Fraction of the lepton wave function in the proton shell vs background.
///
/// enrichment = (|ψ|² in shell / shell_size) / (|ψ|² in bg / bg_size)
///
/// enrichment > 1 → lepton concentrated near proton (quantum hydrogen).
/// enrichment = 1 → uniform (unbound, no hydrogen).
pub fn quantum_shell_enrichment(
    psi: &LeptonPsi,
    proton_sites: &std::collections::HashSet<usize>,
    cfg: &LatticeConfig,
) -> f64 {
    let n = cfg.n_sites();
    let layer_stride = cfg.layer_stride();

    // Shell = non-proton sites adjacent to any proton quark
    let mut shell = std::collections::HashSet::new();
    for &s in proton_sites {
        let (r, c, z) = site_coords(s, cfg);
        for nb in mesh_neighbours(r, c, z, cfg) {
            if !proton_sites.contains(&nb) {
                shell.insert(nb);
            }
        }
    }

    // Restrict to proton-containing layers
    let proton_layers: std::collections::HashSet<usize> =
        proton_sites.iter().map(|&s| s / layer_stride).collect();
    let bg: Vec<usize> = (0..n)
        .filter(|&s| {
            proton_layers.contains(&(s / layer_stride))
                && !proton_sites.contains(&s)
                && !shell.contains(&s)
        })
        .collect();

    let shell_sz = shell.len().max(1);
    let bg_sz = bg.len().max(1);

    let p_shell: f64 = shell.iter().map(|&s| psi[s].norm_sqr()).sum();
    let p_bg: f64 = bg.iter().map(|&s| psi[s].norm_sqr()).sum();

    let rs = p_shell / shell_sz as f64;
    let rb = p_bg / bg_sz as f64;

    if rb > 1e-30 { (rs / rb).min(20.0) } else { 20.0 }
}

// ── Tests ──────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::LatticeConfig;

    fn single_layer_cfg() -> LatticeConfig {
        LatticeConfig { layers: 1, ..Default::default() }
    }

    #[test]
    fn hamiltonian_zero_potential_zero_state_is_zero() {
        let cfg = single_layer_cfg();
        let n = cfg.n_sites();
        let psi = vec![Complex64::new(0.0, 0.0); n];
        let phi = vec![0.0; n];
        let h_psi = apply_hamiltonian(&psi, &phi, &cfg, -1.0);
        for (i, h) in h_psi.iter().enumerate() {
            assert!(h.norm() < 1e-14, "H|0⟩ ≠ 0 at site {i}: {h}");
        }
    }

    #[test]
    fn hamiltonian_uniform_state_zero_kinetic() {
        // Uniform wave function: every site has the same amplitude.
        // Kinetic term = ψ[i] − mean(ψ[nbrs]) = ψ[i] − ψ[i] = 0
        let cfg = single_layer_cfg();
        let n = cfg.n_sites();
        let a = 1.0 / (n as f64).sqrt();
        let psi = vec![Complex64::new(a, 0.0); n];
        let phi = vec![0.0; n];
        let h_psi = apply_hamiltonian(&psi, &phi, &cfg, -1.0);
        for (i, h) in h_psi.iter().enumerate() {
            assert!(
                h.norm() < 1e-12,
                "Uniform state should have zero kinetic energy at site {i}: {h}"
            );
        }
    }

    #[test]
    fn imaginary_time_preserves_norm() {
        let cfg = single_layer_cfg();
        let n = cfg.n_sites();
        let mut psi: LeptonPsi = (0..n)
            .map(|i| Complex64::new((i as f64 * 0.1).sin(), (i as f64 * 0.1).cos()))
            .collect();
        normalise(&mut psi);
        let phi = vec![0.5; n];

        for _ in 0..10 {
            imaginary_time_step(&mut psi, &phi, &cfg, -1.0, 0.05);
            let norm: f64 = psi.iter().map(|a| a.norm_sqr()).sum();
            assert!((norm - 1.0).abs() < 1e-12, "Norm drifted to {norm} during imaginary time");
        }
    }

    #[test]
    fn attractive_potential_gives_negative_energy() {
        // Coulomb well V = q × φ = (-1) × φ: if φ > 0 at the lepton site,
        // potential energy = -φ < 0. A concentrated wave function in the well
        // has negative total energy → bound state.
        let cfg = single_layer_cfg();
        let n = cfg.n_sites();

        // Point source Coulomb field: high φ at center, decaying outward
        let center = n / 2;
        let mut phi = vec![0.1; n]; // small background
        phi[center] = 2.0; // strong source at center

        // Concentrate the wave function at the center
        let mut psi = vec![Complex64::new(0.0, 0.0); n];
        psi[center] = Complex64::new(1.0, 0.0);

        let (e, ek, ep) = expected_energy(&psi, &phi, &cfg, -1.0);
        // Potential energy = (-1) × φ[center] = -2.0
        // Kinetic energy = ψ[center] - mean(ψ[nbrs]) = 1 - 0 = 1.0
        // Total energy = 1.0 - 2.0 = -1.0 < 0 → bound state
        assert!(
            ep < 0.0,
            "Attractive potential: V = {ep:.4}, expected < 0 (charge=-1, phi=+2)"
        );
        assert!(
            e < 0.0,
            "Bound state energy: E = {e:.4} (kin={ek:.4}, pot={ep:.4}), expected < 0"
        );
    }

    #[test]
    fn ground_state_localises_near_proton() {
        // With a Coulomb well (φ > 0 at center), imaginary time evolution
        // should concentrate the wave function near the well.
        let cfg = single_layer_cfg();
        let n = cfg.n_sites();
        let center = n / 2;

        // Point Coulomb field
        let mut phi = vec![0.0; n];
        phi[center] = 1.5;
        for (i, p) in phi.iter_mut().enumerate() {
            let dist = ((i as i64 - center as i64).pow(2) as f64).sqrt();
            if dist > 0.0 { *p = 0.8 / dist; }
        }

        let shell: Vec<usize> = vec![center]; // simplification: shell = center
        let (psi, e_total, _ek, _ep) = quantum_hydrogen_ground_state(
            &phi, &shell, &cfg, 200, 0.05
        );

        // The ground state should have negative energy (bound)
        assert!(
            e_total < 0.0,
            "Ground state energy = {e_total:.6}, expected < 0 (bound state)"
        );

        // Wave function should peak near the well (center)
        let p_center = psi[center].norm_sqr();
        let p_avg: f64 = psi.iter().map(|a| a.norm_sqr()).sum::<f64>() / n as f64;
        assert!(
            p_center > p_avg,
            "Ground state should be localised: P(center)={p_center:.6} < avg={p_avg:.6}"
        );

        println!("  Quantum hydrogen ground state:");
        println!("    E = {e_total:.6} (bound: E < 0 ✓)");
        println!("    P(center) = {p_center:.6}  vs  P(avg) = {p_avg:.6}");
        println!("    Localisation ratio: {:.2}×", p_center / p_avg);
    }
}
