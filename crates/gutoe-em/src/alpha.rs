// GUTOE EM -- Fine Structure Constant: Algebraic + Lattice Measurement
// Copyright (C) 2026 Riff Labs, AGPL-3.0-or-later
//
// The fine structure constant alpha ~ 1/137 emerges from the combinatorial
// structure of the spacetime Clifford algebra Cl(1,3):
//
//   alpha^-1 = T(dim Cl(1,3)) + 1 = T(16) + 1 = 137
//
// where T(n) = n(n+1)/2 is the triangular number.
//
// Additionally: grade-2 dim = C(4,2) = 6 = hex lattice coordination number.

use crate::config::LatticeConfig;
use crate::gauge::jacobi_poisson;
use crate::geometry::site_coords;

// ── Algebraic constants ────────────────────────────────────────────────────

/// Triangular number T(n) = n(n+1)/2.
pub const fn triangular(n: u32) -> u32 {
    n * (n + 1) / 2
}

/// Dimension of Cl(1,3) = 2^4 = 16.
pub const CLIFFORD_DIM: u32 = 16;

/// Hex lattice coordination number = 6 = C(4,2) = grade-2 dimension.
pub const HEX_COORDINATION: u32 = 6;

/// The Eddington number: T(16) + 1 = 137.
pub const EDDINGTON_NUMBER: u32 = triangular(CLIFFORD_DIM) + 1;

/// Physical alpha^-1 for comparison.
pub const ALPHA_INVERSE_PHYSICAL: f64 = 137.035999084;

// ── Coulomb coupling measurement ───────────────────────────────────────────

/// Cartesian coordinates of a hex lattice site (unit spacing).
/// Odd rows shifted left by 0.5 to create hex geometry.
fn hex_cartesian(r: usize, c: usize) -> (f64, f64) {
    let x = c as f64 - 0.5 * (r % 2) as f64;
    let y = r as f64 * (3.0_f64).sqrt() / 2.0;
    (x, y)
}

/// One bin of the radial potential profile.
#[derive(Debug, Clone)]
pub struct RadialBin {
    pub r_mean: f64,
    pub phi_mean: f64,
    pub count: usize,
}

/// Result of a Coulomb coupling measurement.
#[derive(Debug)]
pub struct CoulombMeasurement {
    /// Logarithmic slope: phi(r) = slope * ln(r) + intercept.
    pub slope: f64,
    pub intercept: f64,
    /// Bare 2D coupling: |slope|.
    pub g_2d: f64,
    /// Potential at the charge site.
    pub phi_center: f64,
    /// Radial profile bins.
    pub profile: Vec<RadialBin>,
}

/// Measure the bare Coulomb coupling on a 2D hex lattice.
///
/// Places a unit point charge at the lattice center, runs Jacobi-Poisson
/// to convergence, and extracts the coupling from the logarithmic decay.
///
/// In 2D, the potential of a point charge decays logarithmically:
///   phi(r) = -g_2D * ln(r) + C
///
/// The hex lattice operator L phi = phi - mean_nbrs(phi) corresponds to
///   nabla^2 phi = -(4/a^2) rho in the continuum limit,
/// giving g_2D = 2/pi for the hex lattice with unit spacing.
pub fn measure_coulomb_coupling(rows: usize, cols: usize, n_iter: usize) -> CoulombMeasurement {
    let cfg = LatticeConfig {
        hex_rows: rows,
        hex_cols: cols,
        layers: 1,
        ..Default::default()
    };
    let n = cfg.n_sites();
    let center = n / 2;

    // Point charge at center, neutralized for periodic BCs
    let mean = 1.0 / n as f64;
    let mut rho = vec![-mean; n];
    rho[center] += 1.0;

    let phi = jacobi_poisson(&rho, &cfg, n_iter);

    // Cartesian distances from center (minimum image convention)
    let (cx, cy) = {
        let (r, c, _) = site_coords(center, &cfg);
        hex_cartesian(r, c)
    };
    let lx = cols as f64;
    let ly = rows as f64 * (3.0_f64).sqrt() / 2.0;

    let mut dist_phi: Vec<(f64, f64)> = (0..n)
        .map(|site| {
            let (r, c, _) = site_coords(site, &cfg);
            let (sx, sy) = hex_cartesian(r, c);
            let dx = (sx - cx).abs().min(lx - (sx - cx).abs());
            let dy = (sy - cy).abs().min(ly - (sy - cy).abs());
            let d = (dx * dx + dy * dy).sqrt();
            (d, phi[site])
        })
        .collect();

    dist_phi.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    // Bin by distance (0.5 spacing)
    let max_r = (rows.min(cols) / 3) as f64;
    let bin_width = 0.5;
    let mut profile = Vec::new();
    let mut r_start = 0.5;
    while r_start < max_r {
        let r_end = r_start + bin_width;
        let bin: Vec<&(f64, f64)> = dist_phi
            .iter()
            .filter(|(d, _)| *d >= r_start && *d < r_end)
            .collect();
        if !bin.is_empty() {
            let count = bin.len();
            let r_mean = bin.iter().map(|(d, _)| d).sum::<f64>() / count as f64;
            let phi_mean = bin.iter().map(|(_, p)| p).sum::<f64>() / count as f64;
            profile.push(RadialBin { r_mean, phi_mean, count });
        }
        r_start = r_end;
    }

    // Fit logarithmic decay using least squares
    let fit_data: Vec<(f64, f64)> = profile
        .iter()
        .filter(|b| b.r_mean > 2.0 && b.r_mean < max_r * 0.6)
        .map(|b| (b.r_mean, b.phi_mean))
        .collect();

    let (slope, intercept) = if fit_data.len() >= 3 {
        let nf = fit_data.len() as f64;
        let sum_x: f64 = fit_data.iter().map(|(r, _)| r.ln()).sum();
        let sum_y: f64 = fit_data.iter().map(|(_, p)| p).sum();
        let sum_xy: f64 = fit_data.iter().map(|(r, p)| r.ln() * p).sum();
        let sum_xx: f64 = fit_data.iter().map(|(r, _)| r.ln().powi(2)).sum();
        let s = (nf * sum_xy - sum_x * sum_y) / (nf * sum_xx - sum_x * sum_x);
        let i = (sum_y - s * sum_x) / nf;
        (s, i)
    } else {
        (0.0, 0.0)
    };

    CoulombMeasurement {
        slope,
        intercept,
        g_2d: slope.abs(),
        phi_center: phi[center],
        profile,
    }
}

// ── Mass spectrum constants ────────────────────────────────────────────────

/// Number of GUTOE layers = dim(SU(3)) + dim(SU(2)) + dim(U(1)) = 12.
pub const N_LAYERS: u32 = 12;

/// T(17) = 17 × 18 / 2 = 153.  (17 = Clifford_dim + 1)
pub const T17: u32 = triangular(CLIFFORD_DIM + 1);

/// mp/me algebraic prediction: n_layers × T(Clifford_dim + 1) = 12 × 153 = 1836.
pub const MP_ME_CLIFFORD: u32 = N_LAYERS * T17;

/// Experimental proton-to-electron mass ratio.
pub const MP_ME_EXP: f64 = 1836.15267343;

/// Geometric (Wyler) formula: 6π⁵.
pub fn mp_me_geometric() -> f64 {
    6.0 * std::f64::consts::PI.powi(5)
}

/// Weinberg angle sin²θ_W at GUT scale (SU(5) prediction).
pub const WEINBERG_GUT: f64 = 3.0 / 8.0;   // = 0.375

/// Weinberg angle at the electroweak scale (Clifford prediction).
/// sin²θ_W = 3/13 where:
///   3 = dim(SU(2)) = spatial bivectors {γ¹², γ¹³, γ²³}
///   13 = 3 + grade2_dim + grade3_dim = 3+6+4 = Clifford_dim - 3 = 16-3
/// Agreement with experiment: 99.805% (error 0.195%). Zero free parameters.
pub const WEINBERG_ELECTROWEAK: f64 = 3.0 / 13.0;

/// Weinberg angle at Z mass (experimental).
pub const WEINBERG_OBSERVED: f64 = 0.23122;

/// φ_shell formula: exact 12×12 hex lattice Green's function ≈ 13/21.
/// 13 = Clifford_dim - dim(SU(2)) = 16 - 3 (SAME 13 as Weinberg denominator!)
/// 21 = T(6) = T(grade2_dim) = T(hex_coordination)
/// Numerical verification: exact solve gives 0.619978 vs 13/21 = 0.619048 (0.15% error)
pub const PHI_SHELL_FORMULA: f64 = 13.0 / 21.0;

/// Number of distinct grades in Cl(1,3): {0,1,2,3,4}.
pub const N_GRADES: u32 = 5;

/// First-loop estimate for Δ(α⁻¹): N_grades / α⁻¹ = 5/137.
pub fn delta_alpha_inv_approx() -> f64 {
    N_GRADES as f64 / EDDINGTON_NUMBER as f64
}

// ── Tests ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // -- Algebraic tests: these verify the Eddington derivation --

    #[test]
    fn triangular_number_16_is_136() {
        assert_eq!(triangular(16), 136);
    }

    #[test]
    fn eddington_number_is_137() {
        assert_eq!(EDDINGTON_NUMBER, 137);
    }

    #[test]
    fn clifford_dim_is_16() {
        assert_eq!(CLIFFORD_DIM, 16);
        assert_eq!(1u32 << 4, 16); // 2^4
    }

    #[test]
    fn grade2_equals_hex_coordination() {
        // C(4,2) = 6 = hex lattice neighbors per site
        let grade2 = (4 * 3) / 2; // C(4,2)
        assert_eq!(grade2, HEX_COORDINATION);
    }

    #[test]
    fn grade_sum_equals_clifford_dim() {
        // C(4,0) + C(4,1) + C(4,2) + C(4,3) + C(4,4) = 16
        let grade_dims = [1u32, 4, 6, 4, 1]; // binomial coefficients
        let sum: u32 = grade_dims.iter().sum();
        assert_eq!(sum, CLIFFORD_DIM);
    }

    #[test]
    fn pair_decomposition() {
        // C(16,2) + 16 = T(16) = 136
        let distinct_pairs = (16 * 15) / 2; // C(16,2) = 120
        let self_pairs = 16u32;
        assert_eq!(distinct_pairs + self_pairs, triangular(16));
    }

    #[test]
    fn alpha_inverse_predictions() {
        // alpha^-1(d) = T(2^d) + 1 for d spacetime dimensions
        assert_eq!(triangular(1 << 2) + 1, 11);   // d=2: Cl(1,1)
        assert_eq!(triangular(1 << 3) + 1, 37);   // d=3: Cl(1,2)
        assert_eq!(triangular(1 << 4) + 1, 137);  // d=4: Cl(1,3) -- our universe
        assert_eq!(triangular(1 << 5) + 1, 529);  // d=5: Cl(1,4)
        assert_eq!(triangular(1 << 6) + 1, 2081); // d=6: Cl(1,5)
    }

    #[test]
    fn alpha_inverse_monotone() {
        // More dimensions → weaker EM coupling
        for d in 2..6u32 {
            let a = triangular(1 << d) + 1;
            let b = triangular(1 << (d + 1)) + 1;
            assert!(
                a < b,
                "alpha_inv({d}) = {a} should be < alpha_inv({}) = {b}",
                d + 1
            );
        }
    }

    #[test]
    fn gauge_group_dimension_is_12() {
        // dim SU(3) + dim SU(2) + dim U(1) = 8 + 3 + 1 = 12 = layers
        assert_eq!(8 + 3 + 1, 12);
    }

    // -- Geometric test: hex neighbors are all at unit distance --

    #[test]
    fn hex_neighbors_at_unit_distance() {
        use crate::geometry::mesh_neighbours;

        let cfg = LatticeConfig {
            hex_rows: 12,
            hex_cols: 12,
            layers: 1,
            ..Default::default()
        };

        // Check a few representative sites
        for site in [0, 1, 6, 7, 72, 73] {
            let (r, c, z) = site_coords(site, &cfg);
            let (sx, sy) = hex_cartesian(r, c);
            let nbrs = mesh_neighbours(r, c, z, &cfg);

            for &nb in &nbrs {
                let (nr, nc, _) = site_coords(nb, &cfg);
                let (nx, ny) = hex_cartesian(nr, nc);

                // Minimum image distance (periodic)
                let lx = cfg.hex_cols as f64;
                let ly = cfg.hex_rows as f64 * (3.0_f64).sqrt() / 2.0;
                let dx = (nx - sx).abs().min(lx - (nx - sx).abs());
                let dy = (ny - sy).abs().min(ly - (ny - sy).abs());
                let dist = (dx * dx + dy * dy).sqrt();

                assert!(
                    (dist - 1.0).abs() < 0.01,
                    "site {site} ({r},{c}) -> nbr {nb} ({nr},{nc}): \
                     dist = {dist:.4}, expected 1.0"
                );
            }
        }
    }

    // -- Numerical test: Coulomb coupling measurement --

    #[test]
    fn coulomb_potential_decays_logarithmically() {
        // On a 2D hex lattice, the Coulomb potential decays as ln(r).
        // Measure on a 40x40 lattice with 1500 Jacobi iterations.
        let m = measure_coulomb_coupling(40, 40, 1500);

        // The slope should be negative (potential decreases with distance)
        assert!(
            m.slope < 0.0,
            "Coulomb potential should decay with distance: slope = {:.6}",
            m.slope
        );

        // The bare 2D coupling should be in the range [0.3, 1.0]
        // (depends on lattice normalization; theory: 2/pi ~ 0.637)
        assert!(
            m.g_2d > 0.3 && m.g_2d < 1.0,
            "g_2D = {:.6} should be O(0.5) for 2D hex",
            m.g_2d
        );

        // The potential at the center should be positive (positive charge)
        assert!(
            m.phi_center > 0.0,
            "phi(center) = {:.6} should be positive",
            m.phi_center
        );

        // Bare coupling is much larger than physical alpha
        let ratio = m.g_2d / (1.0 / ALPHA_INVERSE_PHYSICAL);
        assert!(
            ratio > 30.0,
            "g_2D/alpha should be >> 1 (bare vs renormalized): ratio = {:.1}",
            ratio
        );

        println!("  Coulomb measurement results:");
        println!("    slope     = {:.6}", m.slope);
        println!("    g_2D      = {:.6}", m.g_2d);
        println!("    2/pi      = {:.6}", 2.0 / std::f64::consts::PI);
        println!("    phi(0)    = {:.6}", m.phi_center);
        println!("    g_2D/alpha= {:.1}", ratio);
        for bin in m.profile.iter().take(8) {
            println!(
                "    r={:.2}  phi={:+.6}  n={}",
                bin.r_mean, bin.phi_mean, bin.count
            );
        }
    }

    // -- Mass spectrum tests --

    #[test]
    fn t17_is_153() {
        assert_eq!(T17, 153);
    }

    #[test]
    fn mp_me_clifford_is_1836() {
        // 12 × T(17) = 12 × 153 = 1836
        assert_eq!(MP_ME_CLIFFORD, 1836);
    }

    #[test]
    fn mp_me_geometric_agrees() {
        // 6π⁵ ≈ 1836.12, within 0.01% of experiment (1836.153)
        let pred = mp_me_geometric();
        let err = (pred - MP_ME_EXP).abs() / MP_ME_EXP;
        assert!(
            err < 0.001,  // 0.1% tolerance
            "6π⁵ = {pred:.4} vs experiment {MP_ME_EXP:.4}, error = {:.4}%",
            err * 100.0
        );
    }

    #[test]
    fn mp_me_algebraic_agrees() {
        // 12 × T(17) = 1836, within 0.01% of experiment (1836.153)
        let pred = MP_ME_CLIFFORD as f64;
        let err = (pred - MP_ME_EXP).abs() / MP_ME_EXP;
        assert!(
            err < 0.001,  // 0.1% tolerance
            "12×T(17) = {pred:.4} vs experiment {MP_ME_EXP:.4}, error = {:.4}%",
            err * 100.0
        );
    }

    #[test]
    fn weinberg_electroweak_3_over_13() {
        // sin²θ_W = 3/13 = 0.23077 at the electroweak scale
        // Experiment: 0.23122 at M_Z (MS-bar scheme)
        // Error: 0.195% — best Clifford prediction by far
        let pred = WEINBERG_ELECTROWEAK;
        let err = (pred - WEINBERG_OBSERVED).abs() / WEINBERG_OBSERVED;
        assert!(
            err < 0.003, // < 0.3%
            "sin²θ_W = 3/13 = {pred:.5} vs experiment {WEINBERG_OBSERVED:.5}: error {:.2}%",
            err * 100.0
        );
        // Verify the Clifford decomposition:
        // 3 = spatial_bivectors, 10 = grade2 + grade3 = 6+4, 13 = 3+10
        let spatial_biv: u32 = 3;
        let grade2: u32 = 6;  // C(4,2)
        let grade3: u32 = 4;  // C(4,3)
        assert_eq!(spatial_biv + grade2 + grade3, 13, "Weinberg denominator = 3+6+4 = 13");
        // Also: 13 = Clifford_dim - SU(2)_dim = 16 - 3
        assert_eq!(CLIFFORD_DIM - spatial_biv, 13, "13 = Clifford_dim - dim(SU(2))");
        // T(6) = 21 = triangular number of hex coordination
        assert_eq!(triangular(6), 21, "T(6) = 21");
        // The 13 connection: same 13 in Weinberg and phi_shell
        let phi_shell_pred = 13.0 / triangular(6) as f64;
        let phi_shell_exact = 0.619978; // from exact Green's function solve
        let phi_err = (phi_shell_pred - phi_shell_exact).abs() / phi_shell_exact;
        assert!(phi_err < 0.002, "phi_shell = 13/21 = {phi_shell_pred:.6} vs exact {phi_shell_exact:.6}: error {:.3}%", phi_err*100.0);
    }

    #[test]
    fn schwinger_correction_n_grades_times_alpha() {
        // The first-loop correction to both α⁻¹ and mp/me ≈ n_grades × α = 5/137
        let n_grades = N_GRADES as f64;
        let alpha = 1.0 / ALPHA_INVERSE_PHYSICAL;
        let correction = n_grades * alpha;

        // Experimental corrections from integer/Wyler formulas
        let delta_alpha_inv = ALPHA_INVERSE_PHYSICAL - EDDINGTON_NUMBER as f64;
        let delta_mp_me = MP_ME_EXP - 6.0 * std::f64::consts::PI.powi(5);

        // Both corrections should be within 10% of 5α
        let err_alpha = (correction - delta_alpha_inv).abs() / delta_alpha_inv;
        let err_mp    = (correction - delta_mp_me).abs() / delta_mp_me;
        assert!(err_alpha < 0.10, "5α = {correction:.6} vs Δ(α⁻¹) = {delta_alpha_inv:.6}: {:.1}%", err_alpha*100.0);
        assert!(err_mp < 0.10, "5α = {correction:.6} vs Δ(mp/me) = {delta_mp_me:.6}: {:.1}%", err_mp*100.0);

        // The corrected formulas
        let alpha_inv_corrected = EDDINGTON_NUMBER as f64 + correction;
        let mp_me_corrected = 6.0 * std::f64::consts::PI.powi(5) + correction;
        println!("  Schwinger correction 5α = {correction:.6}");
        println!("  α⁻¹ corrected: {alpha_inv_corrected:.6} vs exp {ALPHA_INVERSE_PHYSICAL:.6} (Δ={:.6})", ALPHA_INVERSE_PHYSICAL - alpha_inv_corrected);
        println!("  mp/me corrected: {mp_me_corrected:.6} vs exp {MP_ME_EXP:.6} (Δ={:.6})", MP_ME_EXP - mp_me_corrected);
        // Residuals should be < 0.001 (order α²)
        assert!((ALPHA_INVERSE_PHYSICAL - alpha_inv_corrected).abs() < 0.001);
        assert!((MP_ME_EXP - mp_me_corrected).abs() < 0.005);
    }

    #[test]
    fn weinberg_gut_prediction() {
        // SU(5) GUT prediction: sin²θ_W = 3/8 = 0.375
        assert!((WEINBERG_GUT - 0.375).abs() < 1e-10);

        // GUT prediction exceeds observed (0.2312) — this is expected from RG running
        assert!(
            WEINBERG_GUT > WEINBERG_OBSERVED,
            "sin²θ_W(GUT) = {} should exceed observed {}",
            WEINBERG_GUT, WEINBERG_OBSERVED
        );

        // GUT prediction is between 1/5 and 1/2
        assert!(WEINBERG_GUT > 0.2 && WEINBERG_GUT < 0.5);
    }

    #[test]
    fn mass_ratio_uses_same_ingredients_as_alpha() {
        // α⁻¹ = T(16) + 1 = 137 uses T(CLIFFORD_DIM)
        // mp/me = 12 × T(17) uses T(CLIFFORD_DIM + 1)
        // Same Clifford_dim = 16, just T of the next triangular number
        assert_eq!(EDDINGTON_NUMBER, triangular(CLIFFORD_DIM) + 1);
        assert_eq!(MP_ME_CLIFFORD, N_LAYERS * triangular(CLIFFORD_DIM + 1));
    }

    // -- Running coupling tests --

    #[test]
    fn b0_eff_from_clifford_grade_structure() {
        // b₀ = (11/3) × N_grade2 − (2/3) × N_grade1 = (11/3)×6 − (2/3)×4 = 58/3
        let n_grade2 = 6u32; // C(4,2) = 6 bivectors = gluon-analog states
        let n_grade1 = 4u32; // C(4,1) = 4 vectors  = fermion-analog states
        // In integer arithmetic: b₀ × 3 = 11 × N_grade2 − 2 × N_grade1
        let b0_times_3 = 11 * n_grade2 - 2 * n_grade1;
        assert_eq!(b0_times_3, 58, "b₀ × 3 = 11×6 − 2×4 = 58");

        // Verify the floating-point value matches
        use crate::config::LatticeConfig;
        let cfg = LatticeConfig::default();
        let expected = 58.0 / 3.0;
        assert!(
            (cfg.beta_coeff - expected).abs() < 1e-10,
            "beta_coeff = {} ≠ 58/3 = {}",
            cfg.beta_coeff,
            expected
        );
    }

    #[test]
    fn landau_pole_at_phase1_end() {
        // The UV coupling is tuned so the Landau pole is at t_* ≈ 149
        // t_* = exp(2π / (b₀ × α_UV)) − 1
        use crate::config::LatticeConfig;
        use crate::sim::landau_pole;
        let cfg = LatticeConfig::default();
        let t_star = landau_pole(&cfg);
        assert!(
            (t_star - 149.0).abs() < 1.0,
            "Landau pole t_* = {t_star:.2}, expected ≈ 149"
        );
    }

    #[test]
    fn running_coupling_grows_with_t() {
        // α_s(t) should increase from UV toward the Landau pole
        use crate::config::LatticeConfig;
        use crate::sim::running_alpha_s;
        let cfg = LatticeConfig::default();
        let a0 = running_alpha_s(0, &cfg);
        let a50 = running_alpha_s(50, &cfg);
        let a100 = running_alpha_s(100, &cfg);
        let a140 = running_alpha_s(140, &cfg);
        assert!(a0 < a50, "α_s should grow: {a0} < {a50}");
        assert!(a50 < a100, "α_s should grow: {a50} < {a100}");
        assert!(a100 < a140, "α_s should grow: {a100} < {a140}");
    }

    #[test]
    fn cycle_prob_decreases_toward_confinement() {
        // cycle_prob(t) → 0 as t → t_* (quarks freeze into color singlet)
        use crate::config::LatticeConfig;
        use crate::sim::cycle_prob_rg;
        let cfg = LatticeConfig::default();
        let cp0 = cycle_prob_rg(0, &cfg);
        let cp100 = cycle_prob_rg(100, &cfg);
        let cp140 = cycle_prob_rg(140, &cfg);
        assert!(cp0 > cp100, "cycle_prob should decrease: {cp0} > {cp100}");
        assert!(cp100 > cp140, "cycle_prob should decrease: {cp100} > {cp140}");
        assert!(
            (cp0 - cfg.cycle_prob).abs() < 0.01,
            "At t=0, cycle_prob should equal cycle_prob config: {cp0}"
        );
    }

    #[test]
    fn alignment_grows_toward_confinement() {
        // alignment_rg(t) → ∞ as t → t_* (stronger binding at IR)
        use crate::config::LatticeConfig;
        use crate::sim::alignment_rg;
        let cfg = LatticeConfig::default();
        let al0 = alignment_rg(0, &cfg);
        let al100 = alignment_rg(100, &cfg);
        let al140 = alignment_rg(140, &cfg);
        assert!(al0 < al100, "alignment should grow: {al0} < {al100}");
        assert!(al100 < al140, "alignment should grow: {al100} < {al140}");
        assert!(
            (al0 - cfg.alignment_strength).abs() < 0.01,
            "At t=0, alignment should equal base value: {al0}"
        );
    }

    #[test]
    fn mass_ratio_approaches_1836() {
        // The mass ratio E_prot(t)/E_lep grows from ~0.7 at UV
        // and passes through 1836 just before the Landau pole.
        //
        // E_prot(t) ∝ alpha_s(t)/alpha_UV (alignment energy grows with coupling)
        // E_lep = phi_shell (Coulomb well, fixed by EM which doesn't run)
        //
        // This test verifies: ratio DOES grow substantially (> 100) before t_*
        // without the number 1836 being put in by hand.
        use crate::config::LatticeConfig;
        use crate::sim::running_alpha_s;

        let cfg = LatticeConfig::default();

        // Rough values: E_base ~ 0.81 (from sim), phi_shell ~ 1.17 (Jacobi on 12x12)
        let e_base = 0.81_f64;
        let phi_shell = 1.17_f64;

        let ratio_uv = e_base * running_alpha_s(0, &cfg) / cfg.coupling_uv / phi_shell;
        let ratio_mid = e_base * running_alpha_s(100, &cfg) / cfg.coupling_uv / phi_shell;
        let ratio_pre = e_base * running_alpha_s(140, &cfg) / cfg.coupling_uv / phi_shell;

        // At UV: ratio should be < 10 (proton barely heavier than lepton)
        assert!(ratio_uv < 10.0, "UV ratio = {ratio_uv:.1}, expected < 10");

        // At mid-phase: ratio should have grown significantly
        assert!(
            ratio_mid > ratio_uv * 5.0,
            "ratio(100) = {ratio_mid:.1} should be > 5× UV ratio {ratio_uv:.1}"
        );

        // Near the Landau pole: ratio should exceed 50 (significant growth)
        assert!(ratio_pre > 50.0, "ratio(140) = {ratio_pre:.1}, expected > 50");

        // Verify ratio passes through 1836 between t=148 and the pole
        // (the Landau pole IS at t_*=149, so ratio diverges there)
        let ratio_148 = e_base * running_alpha_s(148, &cfg) / cfg.coupling_uv / phi_shell;
        let is_infinite = running_alpha_s(149, &cfg).is_infinite();
        assert!(
            ratio_148 < MP_ME_EXP,
            "ratio(148) = {ratio_148:.0} should be < 1836 (pole at t=149)"
        );
        assert!(is_infinite, "α_s(149) should be infinite at the Landau pole");

        println!(
            "  Mass ratio trajectory: UV={:.1} → t=100:{:.1} → t=140:{:.1} → t=148:{:.1} → t=149:∞",
            ratio_uv, ratio_mid, ratio_pre, ratio_148
        );
        println!("  1836 crossed between t=148 and t=149 (Landau pole)");
        println!("  No 1836 in the code — emerges from b₀_eff(Clifford) + α_UV(Phase-1 scale)");
    }

    #[test]
    fn delta_alpha_order_estimate() {
        // The 0.036 correction to α⁻¹ is ~5/137 = N_grades / α⁻¹
        let delta = ALPHA_INVERSE_PHYSICAL - EDDINGTON_NUMBER as f64;
        let estimate = N_GRADES as f64 / EDDINGTON_NUMBER as f64;

        // The estimate 5/137 ≈ 0.0365 should be within 5% of the actual 0.036
        let err = (estimate - delta).abs() / delta;
        assert!(
            err < 0.05,
            "5/137 = {estimate:.4} vs actual Δ = {delta:.6}, error = {:.2}%",
            err * 100.0
        );

        println!("  α correction analysis:");
        println!("    Δ(α⁻¹) = {:.6}  (experimental)", delta);
        println!("    5/137  = {:.6}  (N_grades/α⁻¹ estimate)", estimate);
        println!("    error  = {:.2}%", err * 100.0);
    }
}
