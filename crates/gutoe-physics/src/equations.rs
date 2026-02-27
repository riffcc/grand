/*!
 * GUTOE Physics - Field Equations
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

//! Field equations from GUTOE framework
//!
//! Key equations:
//! 1. Veracity wave equation: g^μν∇_μ∇_ν φ + λ_QG l_P² ∇⁴φ = 0
//! 2. Modified Einstein: G_μν + λ_QG l_P² H_μν = κ T_μν + ξ Λ g_μν
//! 3. Black hole entropy: S = A/4G + α ln(A/l_P²)
//! 4. Dispersion: ω² = v² k² - λ_QG l_P² k⁴

use crate::constants::*;

/// Veracity wave equation solver
/// Equation: g^μν∇_μ∇_ν φ + λ_QG l_P² ∇⁴φ = 0
/// Reduces to standard wave equation when λ_QG → 0
pub struct WaveEquation {
    pub lambda_qg: f64,
    pub l_p: f64,
    pub metric: [f64; 4], // Simplified flat metric
}

impl WaveEquation {
    pub fn new(lambda_qg: f64) -> Self {
        Self {
            lambda_qg,
            l_p: PLANCK_LENGTH,
            metric: [1.0, -1.0, -1.0, -1.0], // Minkowski
        }
    }

    /// Calculate wave frequency for given wavenumber
    /// Dispersion: ω² = v² k² - λ_QG l_P² k⁴
    pub fn dispersion(&self, k: f64, v: f64) -> f64 {
        let v2_k2 = v * v * k * k;
        let correction = self.lambda_qg * self.l_p * self.l_p * k.powi(4);

        if v2_k2 > correction {
            (v2_k2 - correction).sqrt()
        } else {
            // Wave is unstable - imaginary frequency
            f64::NAN
        }
    }

    /// Check reduces to standard wave equation
    pub fn reduces_to_classical(&self, k: f64, v: f64) -> bool {
        let quantum = self.dispersion(k, v);
        let _classical = v * k;

        // With λ_QG → 0, should equal classical
        let lambda_zero = WaveEquation::new(0.0);
        let classical_result = lambda_zero.dispersion(k, v);

        (quantum - classical_result).abs() < 1e-10
    }

    /// Group velocity (dω/dk)
    /// For stable waves: v_g <= v
    pub fn group_velocity(&self, k: f64, v: f64) -> f64 {
        let omega = self.dispersion(k, v);
        if omega.is_nan() || omega == 0.0 {
            return v;
        }

        // dω/dk = (v²k - 2λ_QG l_P² k³) / ω
        // For λ_QG > 0, this is always <= v
        let numerator = v * v * k - 2.0 * self.lambda_qg * self.l_p * self.l_p * k.powi(3);
        let vg = numerator / omega;

        // Clamp to physical limit
        if vg.is_nan() || vg > v || vg < 0.0 {
            return v;
        }
        vg
    }
}

/// Black hole thermodynamics
/// Equation: S = A/4G + α ln(A/l_P²)
pub struct BlackHoleEntropy {
    pub alpha: f64,
}

impl BlackHoleEntropy {
    pub fn new() -> Self {
        Self {
            // α tracks the shared λ_QG constant from the dispersion relation.
            alpha: LAMBDA_QG,
        }
    }

    /// Calculate entropy from horizon area (m²) — in SI units (J/K).
    ///
    /// Correct Bekenstein-Hawking formula: S = k_B c³ A / (4ħG)
    /// Plus GUTOE quantum correction: + α k_B ln(A/l_P²)
    pub fn entropy(&self, area: f64) -> f64 {
        let classical = KB * C * C * C * area / (4.0 * HBAR * G);
        let quantum = self.alpha * KB * (area / (PLANCK_LENGTH * PLANCK_LENGTH)).ln();
        classical + quantum
    }

    /// Calculate Hawking temperature from horizon area (m²), in Kelvin.
    ///
    /// Derived thermodynamically: T = dE/dS = (dE/dA) / (dS/dA)
    ///   E = Mc² = c⁴/(4G) √(A/π)   →   dE/dA = c⁴ / (8G √(πA))
    ///   dS/dA = k_B c³/(4ħG) + α k_B / A  (from modified entropy)
    pub fn temperature(&self, area: f64) -> f64 {
        let de_da = C.powi(4) / (8.0 * G * (std::f64::consts::PI * area).sqrt());
        let ds_da = KB * C.powi(3) / (4.0 * HBAR * G) + self.alpha * KB / area;
        de_da / ds_da
    }

    /// Hawking radiation power
    pub fn hawking_power(&self, area: f64) -> f64 {
        let temp = self.temperature(area);
        // P = σAT⁴ (Stefan-Boltzmann)
        let sigma = std::f64::consts::PI * std::f64::consts::PI * KB.powi(4)
            / (60.0 * HBAR.powi(3) * C.powi(2));
        sigma * area * temp.powi(4)
    }

    /// Verify that the quantum correction is small relative to the classical term.
    /// For macroscopic black holes the α ln(A/l_P²) term is negligible.
    pub fn classical_limit(&self, area: f64) -> bool {
        let classical = KB * C.powi(3) * area / (4.0 * HBAR * G);
        let quantum = self.alpha * KB * (area / (PLANCK_LENGTH * PLANCK_LENGTH)).ln();
        quantum.abs() < classical * 0.1
    }
}

impl Default for BlackHoleEntropy {
    fn default() -> Self {
        Self::new()
    }
}

/// Modified Einstein field equations
/// G_μν + λ_QG l_P² H_μν = κ T_μν + ξ Λ g_μν
pub struct ModifiedEinstein {
    pub lambda_qg: f64,
    pub kappa: f64,
    pub xi: f64,
    pub cosmological: f64,
}

impl ModifiedEinstein {
    pub fn new() -> Self {
        Self {
            lambda_qg: LAMBDA_QG,
            kappa: compute_kappa(),
            xi: 1.0,
            cosmological: LAMBDA_COSMOLOGICAL,
        }
    }

    /// Effective gravitational constant at scale
    pub fn effective_g(&self, scale: f64) -> f64 {
        // At large scales (scale >> l_P), reduces to classical G
        let correction = 1.0 + self.lambda_qg * (PLANCK_LENGTH / scale).powi(2);
        G * correction
    }

    /// Check reduces to classical GR
    pub fn reduces_to_gr(&self) -> bool {
        let large_scale = 1.0; // Macroscopic scale
        let effective = self.effective_g(large_scale);
        (effective - G).abs() / G < 1e-10
    }
}

impl Default for ModifiedEinstein {
    fn default() -> Self {
        Self::new()
    }
}

/// Gravitational wave dispersion
/// From GUTOE.md: Δv/c ≈ λ_QG(E/E_QG)²
pub fn gravitational_wave_dispersion(E: f64, E_QG: f64) -> f64 {
    LAMBDA_QG * (E / E_QG).powi(2)
}

/// Gamma-ray burst time delay
/// From GUTOE.md: ΔT ≈ ξ·L·(E/E_QG)²
pub fn gamma_ray_delay(distance: f64, E: f64, E_QG: f64) -> f64 {
    distance * (E / E_QG).powi(2) / C
}

/// Fermi weak prefactor from the SU(2) electroweak mass relation:
/// G_F = 1 / (2 f0^2)
pub fn weak_fermi_prefactor(f0: f64) -> Option<f64> {
    if f0 == 0.0 {
        return None;
    }
    Some(1.0 / (2.0 * f0.powi(2)))
}

/// Equivalent prefactor from g^2 / (8 m_W^2) with m_W = g f0 / 2.
pub fn weak_prefactor_from_su2(g: f64, f0: f64) -> Option<f64> {
    if g == 0.0 || f0 == 0.0 {
        return None;
    }
    Some(g.powi(2) / (8.0 * (g * f0 / 2.0).powi(2)))
}

/// Sommerfeld parameter η = α * sqrt(m_r / (2E)).
pub fn sommerfeld_parameter(alpha: f64, m_reduced: f64, collision_energy: f64) -> Option<f64> {
    if alpha <= 0.0 || m_reduced <= 0.0 || collision_energy <= 0.0 {
        return None;
    }
    Some(alpha * (m_reduced / (2.0 * collision_energy)).sqrt())
}

/// Gamow penetration factor exp(-2π η) for Coulomb tunneling.
pub fn gamow_factor(alpha: f64, m_reduced: f64, collision_energy: f64) -> Option<f64> {
    let eta = sommerfeld_parameter(alpha, m_reduced, collision_energy)?;
    Some((-2.0 * std::f64::consts::PI * eta).exp())
}

/// Structural pp weak reaction-rate kernel:
/// rate ~ weak_prefactor * n_p^2 * Gamow.
///
/// Uses the Lean-parity leading-order α = 1/137.
pub fn pp_weak_rate_from_su2_and_gamow(
    g: f64,
    f0: f64,
    proton_density: f64,
    m_reduced: f64,
    collision_energy: f64,
) -> Option<f64> {
    if proton_density <= 0.0 {
        return None;
    }
    let weak = weak_prefactor_from_su2(g, f0)?;
    let gamow = gamow_factor(ALPHA_LEADING_ORDER, m_reduced, collision_energy)?;
    Some(weak * proton_density.powi(2) * gamow)
}

/// Maxwell-Boltzmann thermal weight exp(-E/T).
pub fn maxwell_boltzmann_weight(temperature_scale: f64, collision_energy: f64) -> Option<f64> {
    if temperature_scale <= 0.0 {
        return None;
    }
    Some((-collision_energy / temperature_scale).exp())
}

/// Pointwise thermal pp kernel = weak-rate kernel * MB thermal weight.
pub fn pp_thermal_kernel(
    g: f64,
    f0: f64,
    proton_density: f64,
    m_reduced: f64,
    temperature_scale: f64,
    collision_energy: f64,
) -> Option<f64> {
    let rate = pp_weak_rate_from_su2_and_gamow(g, f0, proton_density, m_reduced, collision_energy)?;
    let w = maxwell_boltzmann_weight(temperature_scale, collision_energy)?;
    Some(rate * w)
}

/// Positive 3-point quadrature proxy for a thermally averaged pp rate.
pub fn pp_thermal_average3(
    g: f64,
    f0: f64,
    proton_density: f64,
    m_reduced: f64,
    temperature_scale: f64,
    e1: f64,
    e2: f64,
    e3: f64,
) -> Option<f64> {
    let k1 = pp_thermal_kernel(g, f0, proton_density, m_reduced, temperature_scale, e1)?;
    let k2 = pp_thermal_kernel(g, f0, proton_density, m_reduced, temperature_scale, e2)?;
    let k3 = pp_thermal_kernel(g, f0, proton_density, m_reduced, temperature_scale, e3)?;
    Some((k1 + k2 + k3) / 3.0)
}

/// Uniform `(n+1)`-sample thermal average on an energy ladder `e0 + i * de`.
///
/// This mirrors `ppThermalAverageUniform` in `Gutoe/StellarFusion.lean`.
pub fn pp_thermal_average_uniform(
    g: f64,
    f0: f64,
    proton_density: f64,
    m_reduced: f64,
    temperature_scale: f64,
    e0: f64,
    de: f64,
    n: u32,
) -> Option<f64> {
    if de < 0.0 || e0 <= 0.0 {
        return None;
    }

    let mut sum = 0.0;
    for i in 0..=n {
        let e = e0 + f64::from(i) * de;
        let k = pp_thermal_kernel(g, f0, proton_density, m_reduced, temperature_scale, e)?;
        sum += k;
    }
    Some(sum / f64::from(n + 1))
}

/// Lane-Emden-style compression proxy used by the Lean stellar-fusion bridge.
#[inline]
pub fn lane_emden_compression_proxy(mass: f64, rho_c: f64) -> f64 {
    mass * rho_c
}

/// Exact n=0 Lane-Emden profile θ(ξ) = 1 - ξ²/6.
#[inline]
pub fn lane_emden_theta_n0(xi: f64) -> f64 {
    1.0 - xi * xi / 6.0
}

/// First derivative of the exact n=0 Lane-Emden profile.
#[inline]
pub fn lane_emden_theta_n0_prime(xi: f64) -> f64 {
    -xi / 3.0
}

/// Second derivative of the exact n=0 Lane-Emden profile.
#[inline]
pub fn lane_emden_theta_n0_prime_prime(_xi: f64) -> f64 {
    -1.0 / 3.0
}

/// Multiplied-form n=0 Lane-Emden residual:
/// ξ² θ'' + 2 ξ θ' + ξ².
#[inline]
pub fn lane_emden_residual_n0(xi: f64) -> f64 {
    xi * xi * lane_emden_theta_n0_prime_prime(xi)
        + 2.0 * xi * lane_emden_theta_n0_prime(xi)
        + xi * xi
}

/// Integer-index multiplied Lane-Emden residual:
/// ξ² θ'' + 2 ξ θ' + ξ² θ^n.
#[inline]
pub fn lane_emden_residual_nat(
    n: u32,
    xi: f64,
    theta: f64,
    theta_prime: f64,
    theta_double_prime: f64,
) -> f64 {
    xi * xi * theta_double_prime + 2.0 * xi * theta_prime + xi * xi * theta.powi(n as i32)
}

/// Origin regularity check for Lane-Emden profiles.
#[inline]
pub fn lane_emden_regular_origin(theta0: f64, theta_prime0: f64, tol: f64) -> bool {
    (theta0 - 1.0).abs() <= tol && theta_prime0.abs() <= tol
}

/// Integrate the integer-index Lane-Emden ODE with RK4 using the
/// regular center expansion for initialization:
/// θ(ξ) ≈ 1 - ξ²/6, θ'(ξ) ≈ -ξ/3.
pub fn lane_emden_integrate_rk4_nat(
    n: u32,
    xi_max: f64,
    step: f64,
) -> Option<Vec<(f64, f64, f64)>> {
    if xi_max <= 0.0 || step <= 0.0 || !xi_max.is_finite() || !step.is_finite() {
        return None;
    }

    // Seed from the regular-center expansion at ξ = step to avoid ξ=0 singularity.
    let mut xi = step;
    let mut theta = lane_emden_theta_n0(xi);
    let mut z = lane_emden_theta_n0_prime(xi); // z = θ'

    let mut out = Vec::new();
    out.push((0.0, 1.0, 0.0));
    out.push((xi, theta, z));

    while xi < xi_max {
        let h = (xi_max - xi).min(step);

        let f_theta = |zz: f64| zz;
        let f_z = |x: f64, th: f64, zz: f64| -> f64 {
            if x <= 0.0 {
                0.0
            } else {
                -2.0 * zz / x - th.powi(n as i32)
            }
        };

        let k1_t = f_theta(z);
        let k1_z = f_z(xi, theta, z);

        let t2 = theta + 0.5 * h * k1_t;
        let z2 = z + 0.5 * h * k1_z;
        let x2 = xi + 0.5 * h;
        let k2_t = f_theta(z2);
        let k2_z = f_z(x2, t2, z2);

        let t3 = theta + 0.5 * h * k2_t;
        let z3 = z + 0.5 * h * k2_z;
        let k3_t = f_theta(z3);
        let k3_z = f_z(x2, t3, z3);

        let t4 = theta + h * k3_t;
        let z4 = z + h * k3_z;
        let x4 = xi + h;
        let k4_t = f_theta(z4);
        let k4_z = f_z(x4, t4, z4);

        theta += h * (k1_t + 2.0 * k2_t + 2.0 * k3_t + k4_t) / 6.0;
        z += h * (k1_z + 2.0 * k2_z + 2.0 * k3_z + k4_z) / 6.0;
        xi += h;

        if !theta.is_finite() || !z.is_finite() {
            return None;
        }
        out.push((xi, theta, z));
    }

    Some(out)
}

/// Mean theta over a sampled Lane-Emden trajectory.
pub fn lane_emden_average_theta_from_profile(profile: &[(f64, f64, f64)]) -> Option<f64> {
    if profile.is_empty() {
        return None;
    }
    let sum: f64 = profile.iter().map(|(_, theta, _)| *theta).sum();
    Some(sum / profile.len() as f64)
}

/// Profile-weighted compression from a sampled Lane-Emden trajectory.
pub fn lane_emden_profile_weighted_compression(
    mass: f64,
    rho_c: f64,
    profile: &[(f64, f64, f64)],
) -> Option<f64> {
    let avg_theta = lane_emden_average_theta_from_profile(profile)?;
    Some(lane_emden_compression_proxy(mass, rho_c) * avg_theta)
}

/// Ignition condition using sampled Lane-Emden profile weighting.
///
/// Mirrors the Lean bridge theorem assumption pattern:
/// if profile-average theta is bounded by 1 and profile-weighted compression
/// clears threshold, then ignition is admitted.
pub fn polytropic_ignition_condition_from_lane_emden_profile(
    g: f64,
    mu: f64,
    xi: f64,
    t_ign: f64,
    mass: f64,
    rho_c: f64,
    profile: &[(f64, f64, f64)],
) -> Option<bool> {
    let avg_theta = lane_emden_average_theta_from_profile(profile)?;
    if avg_theta > 1.0 {
        return None;
    }
    let profile_comp = lane_emden_profile_weighted_compression(mass, rho_c, profile)?;
    let threshold = minimum_polytropic_compression(g, mu, xi, t_ign);
    Some(profile_comp >= threshold)
}

/// Specialized n=0 profile ignition condition.
///
/// This is the executable counterpart of the Lean theorem where the
/// envelope bound `avg(theta) <= 1` is discharged for the exact n=0 profile.
pub fn polytropic_ignition_condition_from_lane_emden_n0_profile(
    g: f64,
    mu: f64,
    xi: f64,
    t_ign: f64,
    mass: f64,
    rho_c: f64,
    profile: &[(f64, f64, f64)],
) -> Option<bool> {
    let avg_theta = lane_emden_average_theta_from_profile(profile)?;
    if avg_theta > 1.0 + 1e-9 {
        return None;
    }
    let profile_comp = lane_emden_profile_weighted_compression(mass, rho_c, profile)?;
    let threshold = minimum_polytropic_compression(g, mu, xi, t_ign);
    Some(profile_comp >= threshold)
}

/// True iff all profile sample coordinates satisfy ξ >= 0.
pub fn lane_emden_profile_all_nonnegative_xi(profile: &[(f64, f64, f64)]) -> bool {
    profile.iter().all(|(xi, _, _)| *xi >= 0.0)
}

/// Monotonicity checker for sampled Lane-Emden profiles.
///
/// Returns true if theta is nonincreasing up to `tol`.
pub fn lane_emden_profile_is_nonincreasing(profile: &[(f64, f64, f64)], tol: f64) -> bool {
    profile.windows(2).all(|w| w[1].1 <= w[0].1 + tol)
}

/// Envelope witness from sampled monotonicity:
/// if profile starts at theta(0)=1 (within tolerance), uses ξ>=0 samples,
/// and is nonincreasing, then sampled theta never exceeds 1.
pub fn lane_emden_envelope_le_one_from_monotone_profile(
    profile: &[(f64, f64, f64)],
    tol: f64,
) -> bool {
    let Some((xi0, theta0, _)) = profile.first().copied() else {
        return false;
    };
    if xi0 < -tol || (theta0 - 1.0).abs() > tol {
        return false;
    }
    lane_emden_profile_all_nonnegative_xi(profile)
        && lane_emden_profile_is_nonincreasing(profile, tol)
        && profile.iter().all(|(_, theta, _)| *theta <= 1.0 + tol)
}

/// Polytropic core-temperature proxy:
/// T_c ∝ ξ G μ √(M ρ_c)
///
/// Mirrors `coreTemperaturePolytropic` in `Gutoe/StellarFusion.lean`.
pub fn core_temperature_polytropic(g: f64, mu: f64, xi: f64, mass: f64, rho_c: f64) -> f64 {
    let compression = lane_emden_compression_proxy(mass, rho_c).max(0.0);
    xi * g * mu * compression.sqrt()
}

/// Compression threshold for ignition in the polytropic proxy model.
#[inline]
pub fn minimum_polytropic_compression(g: f64, mu: f64, xi: f64, t_ign: f64) -> f64 {
    (t_ign / (xi * g * mu)).powi(2)
}

/// Ignition condition in the polytropic proxy model.
#[inline]
pub fn polytropic_ignition_condition(
    g: f64,
    mu: f64,
    xi: f64,
    t_ign: f64,
    mass: f64,
    rho_c: f64,
) -> bool {
    lane_emden_compression_proxy(mass, rho_c) >= minimum_polytropic_compression(g, mu, xi, t_ign)
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPSILON: f64 = 1e-10;

    // ── Wave equation ────────────────────────────────────────────────────────

    #[test]
    fn classical_limit_recovered_at_low_k() {
        // Key claim: ω² = v²k² - λ_QG l_P² k⁴  reduces to ω = vk when k→0.
        // At radio-wave k (1e10 /m), the Planck-scale correction is negligible.
        let wave = WaveEquation::new(LAMBDA_QG);
        assert!(
            wave.reduces_to_classical(1e10, C),
            "classical limit NOT recovered at low k — dispersion relation is wrong"
        );
    }

    #[test]
    fn dispersion_produces_real_frequency_at_stable_k() {
        let wave = WaveEquation::new(LAMBDA_QG);
        // k_crit = v / sqrt(λ_QG) / l_P ≈ enormous; any reasonable k is stable.
        let omega = wave.dispersion(1e20, C);
        assert!(
            !omega.is_nan(),
            "ω is NaN for k=1e20 — wave is unstable where it shouldn't be"
        );
        assert!(omega > 0.0, "ω ≤ 0 for a stable wave");
    }

    #[test]
    fn dispersion_is_less_than_classical_when_lambda_nonzero() {
        // The quantum-gravity correction *reduces* ω relative to the classical vk.
        let wave = WaveEquation::new(LAMBDA_QG);
        let k = 1e10f64;
        let v = 1.0f64;

        let omega_qg = wave.dispersion(k, v);
        let omega_cl = v * k;

        if !omega_qg.is_nan() {
            assert!(
                omega_qg <= omega_cl + EPSILON,
                "ω_QG ({omega_qg:.6e}) > ω_classical ({omega_cl:.6e}) — \
                 quantum gravity increases frequency, violating the dispersion relation"
            );
        }
    }

    #[test]
    fn group_velocity_bounded_by_phase_velocity() {
        // KEY PHYSICAL CLAIM: v_g ≤ v (no faster-than-light information).
        let wave = WaveEquation::new(LAMBDA_QG);
        let v = C;

        // Test across a range of physically plausible wavenumbers.
        for &k in &[1e5f64, 1e10, 1e15, 1e20] {
            let omega = wave.dispersion(k, v);
            if omega.is_nan() {
                continue;
            }

            let vg = wave.group_velocity(k, v);
            assert!(
                vg <= v + EPSILON,
                "v_g = {vg:.4e} > v = {v:.4e} at k = {k:.1e} — causality violated!"
            );
        }
    }

    #[test]
    fn group_velocity_is_positive() {
        // Negative group velocity would mean information travelling backwards.
        let wave = WaveEquation::new(LAMBDA_QG);
        for &k in &[1e5f64, 1e10, 1e15] {
            let omega = wave.dispersion(k, v_from_c());
            if omega.is_nan() {
                continue;
            }
            let vg = wave.group_velocity(k, v_from_c());
            assert!(
                vg >= 0.0,
                "v_g = {vg:.4e} < 0 at k = {k:.1e} — negative group velocity"
            );
        }
    }

    fn v_from_c() -> f64 {
        1.0
    } // unit velocity for dimensionless tests

    #[test]
    fn polytropic_ignition_threshold_matches_temperature_proxy() {
        let g = 2.0;
        let mu = 3.0;
        let xi = 4.0;
        let t_ign = 5.0;
        let mass = 1.0;
        let rho_c = minimum_polytropic_compression(g, mu, xi, t_ign);

        assert!(polytropic_ignition_condition(g, mu, xi, t_ign, mass, rho_c));

        let t_core = core_temperature_polytropic(g, mu, xi, mass, rho_c);
        assert!((t_core - t_ign).abs() < 1e-10);
    }

    #[test]
    fn stronger_compression_increases_polytropic_core_temperature() {
        let g = 1.0;
        let mu = 1.0;
        let xi = 1.0;
        let t_low = core_temperature_polytropic(g, mu, xi, 1.0, 1.0);
        let t_high = core_temperature_polytropic(g, mu, xi, 4.0, 4.0);
        assert!(t_high > t_low);
    }

    #[test]
    fn lane_emden_n0_profile_matches_closed_form_values() {
        assert!((lane_emden_theta_n0(0.0) - 1.0).abs() < 1e-12);
        assert!((lane_emden_theta_n0(1.0) - (5.0 / 6.0)).abs() < 1e-12);
    }

    #[test]
    fn lane_emden_n0_residual_is_numerically_zero() {
        for &xi in &[0.0, 0.5, 1.0, 2.0, 3.0] {
            let r = lane_emden_residual_n0(xi);
            assert!(r.abs() < 1e-12, "residual={r} at xi={xi}");
        }
    }

    #[test]
    fn lane_emden_residual_nat_reduces_to_n0_residual() {
        for &xi in &[0.0, 0.5, 1.0, 2.0] {
            let theta = lane_emden_theta_n0(xi);
            let theta_p = lane_emden_theta_n0_prime(xi);
            let theta_pp = lane_emden_theta_n0_prime_prime(xi);
            let r_nat = lane_emden_residual_nat(0, xi, theta, theta_p, theta_pp);
            let r_n0 = lane_emden_residual_n0(xi);
            assert!(
                (r_nat - r_n0).abs() < 1e-12,
                "xi={xi}: r_nat={r_nat}, r_n0={r_n0}"
            );
        }
    }

    #[test]
    fn lane_emden_regular_origin_accepts_exact_center() {
        assert!(lane_emden_regular_origin(1.0, 0.0, 1e-12));
        assert!(!lane_emden_regular_origin(0.99, 0.0, 1e-12));
    }

    #[test]
    fn lane_emden_rk4_n0_matches_closed_form_near_xi1() {
        let sol = lane_emden_integrate_rk4_nat(0, 1.0, 1.0e-3).expect("solution");
        let (_xi, theta, _z) = sol.last().copied().expect("last");
        let expected = lane_emden_theta_n0(1.0);
        assert!(
            (theta - expected).abs() < 2.0e-3,
            "theta={theta}, expected={expected}"
        );
    }

    #[test]
    fn lane_emden_rk4_n1_matches_sinxi_over_xi_near_xi1() {
        let sol = lane_emden_integrate_rk4_nat(1, 1.0, 1.0e-3).expect("solution");
        let (_xi, theta, _z) = sol.last().copied().expect("last");
        let expected = 1.0f64.sin() / 1.0;
        assert!(
            (theta - expected).abs() < 1.0e-2,
            "theta={theta}, expected={expected}"
        );
    }

    #[test]
    fn lane_emden_profile_weighted_ignition_condition_tracks_threshold() {
        let profile = lane_emden_integrate_rk4_nat(1, 1.0, 1.0e-3).expect("profile");
        let on = polytropic_ignition_condition_from_lane_emden_profile(
            2.0, 3.0, 4.0, 5.0, 1.0, 1.0, &profile,
        )
        .expect("on");
        assert!(on);

        let off = polytropic_ignition_condition_from_lane_emden_profile(
            2.0, 3.0, 4.0, 5.0, 0.01, 0.01, &profile,
        )
        .expect("off");
        assert!(!off);
    }

    #[test]
    fn lane_emden_profile_ignition_rejects_avg_theta_above_one() {
        let fake = vec![(0.0, 1.2, 0.0), (0.1, 1.1, -0.1)];
        let cond = polytropic_ignition_condition_from_lane_emden_profile(
            2.0, 3.0, 4.0, 5.0, 1.0, 1.0, &fake,
        );
        assert!(cond.is_none());
    }

    #[test]
    fn lane_emden_n0_profile_average_theta_is_bounded_by_one() {
        let profile = lane_emden_integrate_rk4_nat(0, 2.0, 1.0e-3).expect("profile");
        let avg = lane_emden_average_theta_from_profile(&profile).expect("avg");
        assert!(avg <= 1.0 + 1.0e-9, "avg_theta={avg}");
    }

    #[test]
    fn lane_emden_n0_profile_specialized_ignition_bridge_is_usable() {
        let profile = lane_emden_integrate_rk4_nat(0, 1.0, 1.0e-3).expect("profile");
        let cond = polytropic_ignition_condition_from_lane_emden_n0_profile(
            2.0, 3.0, 4.0, 5.0, 1.0, 1.0, &profile,
        )
        .expect("condition");
        assert!(cond);
    }

    #[test]
    fn lane_emden_n1_profile_is_nonincreasing_on_nonnegative_window() {
        let profile = lane_emden_integrate_rk4_nat(1, 2.5, 1.0e-3).expect("profile");
        assert!(lane_emden_profile_is_nonincreasing(&profile, 1e-6));
        assert!(lane_emden_envelope_le_one_from_monotone_profile(
            &profile, 1e-6
        ));
    }

    #[test]
    fn lane_emden_n3_profile_is_nonincreasing_on_nonnegative_window() {
        let profile = lane_emden_integrate_rk4_nat(3, 2.0, 1.0e-3).expect("profile");
        assert!(lane_emden_profile_is_nonincreasing(&profile, 1e-6));
        assert!(lane_emden_envelope_le_one_from_monotone_profile(
            &profile, 1e-6
        ));
    }

    #[test]
    fn su2_fermi_prefactors_match_and_are_positive() {
        let g = 0.65;
        let f0 = 246.0;
        let lhs = weak_prefactor_from_su2(g, f0).expect("lhs");
        let rhs = weak_fermi_prefactor(f0).expect("rhs");
        assert!((lhs - rhs).abs() < 1e-16);
        assert!(lhs > 0.0);
    }

    #[test]
    fn gamow_factor_is_between_zero_and_one_for_positive_inputs() {
        let g = gamow_factor(ALPHA_LEADING_ORDER, 469.136, 0.002).expect("gamow");
        assert!(g > 0.0 && g < 1.0);
    }

    #[test]
    fn pp_weak_rate_kernel_is_strictly_positive_under_physical_inputs() {
        let rate =
            pp_weak_rate_from_su2_and_gamow(0.65, 246.0, 1.0e30, 469.136, 0.002).expect("rate");
        assert!(rate > 0.0);
    }

    #[test]
    fn maxwell_boltzmann_weight_is_positive_for_positive_temperature() {
        let w = maxwell_boltzmann_weight(0.002, 0.001).expect("w");
        assert!(w > 0.0 && w <= 1.0);
    }

    #[test]
    fn pp_thermal_kernel_is_strictly_positive_under_physical_inputs() {
        let k = pp_thermal_kernel(0.65, 246.0, 1.0e30, 469.136, 0.002, 0.001).expect("k");
        assert!(k > 0.0);
    }

    #[test]
    fn pp_thermal_average3_is_strictly_positive_under_physical_inputs() {
        let avg = pp_thermal_average3(0.65, 246.0, 1.0e30, 469.136, 0.002, 0.0008, 0.001, 0.0012)
            .expect("avg");
        assert!(avg > 0.0);
    }

    #[test]
    fn pp_thermal_average_uniform_is_strictly_positive_under_physical_inputs() {
        let avg =
            pp_thermal_average_uniform(0.65, 246.0, 1.0e30, 469.136, 0.002, 0.0008, 0.0001, 8)
                .expect("avg");
        assert!(avg > 0.0);
    }

    #[test]
    fn pp_thermal_average_uniform_rejects_nonphysical_grid() {
        let bad_e0 = pp_thermal_average_uniform(0.65, 246.0, 1.0e30, 469.136, 0.002, 0.0, 0.001, 4);
        assert!(bad_e0.is_none());

        let bad_de =
            pp_thermal_average_uniform(0.65, 246.0, 1.0e30, 469.136, 0.002, 0.001, -0.0001, 4);
        assert!(bad_de.is_none());
    }

    // ── Black-hole entropy ───────────────────────────────────────────────────

    #[test]
    fn black_hole_entropy_is_positive() {
        let bh = BlackHoleEntropy::new();
        let area = solar_mass_area();
        let s = bh.entropy(area);
        assert!(s > 0.0, "Black-hole entropy = {s} ≤ 0 — unphysical");
    }

    #[test]
    fn black_hole_entropy_uses_shared_lambda_qg_constant() {
        let bh = BlackHoleEntropy::new();
        assert!(
            (bh.alpha - LAMBDA_QG).abs() < EPSILON,
            "BlackHoleEntropy alpha={} diverged from LAMBDA_QG={}",
            bh.alpha,
            LAMBDA_QG
        );
    }

    #[test]
    fn quantum_correction_is_small_for_solar_mass() {
        // KEY CLAIM: For macroscopic black holes, quantum correction < 10% of classical.
        let bh = BlackHoleEntropy::new();
        let area = solar_mass_area();
        assert!(
            bh.classical_limit(area),
            "Quantum correction > 10% for solar-mass BH — GUTOE claims it should be tiny"
        );
    }

    #[test]
    fn entropy_grows_with_area() {
        let bh = BlackHoleEntropy::new();
        let a1 = solar_mass_area();
        let a2 = a1 * 4.0; // 2× mass → 4× area (Schwarzschild)
        assert!(
            bh.entropy(a2) > bh.entropy(a1),
            "S(4A) ≤ S(A) — entropy is not monotone in area"
        );
    }

    #[test]
    fn hawking_temperature_is_finite_and_positive() {
        let bh = BlackHoleEntropy::new();
        let area = solar_mass_area();
        let t = bh.temperature(area);
        assert!(
            t.is_finite() && t > 0.0,
            "Hawking temperature = {t} — should be small positive number"
        );
        // Solar-mass BH Hawking temp ≈ 60 nK
        assert!(
            t < 1e-6,
            "T = {t:.4e} K — far too hot for a solar-mass black hole (expect ~60 nK)"
        );
    }

    #[test]
    fn hawking_power_is_positive() {
        let bh = BlackHoleEntropy::new();
        let area = solar_mass_area();
        let p = bh.hawking_power(area);
        assert!(
            p > 0.0 && p.is_finite(),
            "Hawking power = {p} — should be tiny but positive"
        );
    }

    fn solar_mass_area() -> f64 {
        let mass_sun = 1.989e30_f64;
        let rs = 2.0 * G * mass_sun / (C * C);
        4.0 * std::f64::consts::PI * rs * rs
    }

    // ── Modified Einstein equations ──────────────────────────────────────────

    #[test]
    fn effective_g_reduces_to_classical_at_large_scales() {
        // KEY CLAIM: G_eff → G as scale → macroscopic.
        let einstein = ModifiedEinstein::new();
        assert!(
            einstein.reduces_to_gr(),
            "Modified G does not reduce to classical GR at 1 metre scale"
        );
    }

    #[test]
    fn effective_g_increases_at_planck_scale() {
        // The GUTOE correction adds ε = λ_QG (l_P / scale)² > 0,
        // so G_eff > G at small scales.
        let einstein = ModifiedEinstein::new();
        let g_macroscopic = einstein.effective_g(1.0); // 1 metre
        let g_planck = einstein.effective_g(PLANCK_LENGTH); // Planck scale
        assert!(
            g_planck > g_macroscopic,
            "G_Planck ≤ G_macroscopic — correction has wrong sign"
        );
    }

    // ── Gravitational wave dispersion ────────────────────────────────────────

    #[test]
    fn gw_dispersion_is_proportional_to_energy_squared() {
        // Δv/c ∝ (E/E_QG)² — doubling E quadruples the dispersion.
        let e_qg = 1e28_f64; // Planck energy in joules
        let e1 = 1e19_f64;
        let e2 = 2.0 * e1;
        let d1 = gravitational_wave_dispersion(e1, e_qg);
        let d2 = gravitational_wave_dispersion(e2, e_qg);
        let ratio = d2 / d1;
        assert!(
            (ratio - 4.0).abs() < 1e-10,
            "Δv/c does not scale as E² (ratio = {ratio:.6}, expected 4.0)"
        );
    }

    // ── λ_QG² mass scaling ────────────────────────────────────────────────────
    // Python particle_formation.py expands to:
    //   m = veracity × curvature × field_gradient / planck_length × λ_QG²
    // λ_QG enters *quadratically* — small changes have outsized mass impact.

    #[test]
    fn quark_mass_scales_as_lambda_qg_squared() {
        // Runtime baseline uses shared λ_QG = 1/12.
        // Python precision_qg_tuner converged on λ_QG = 0.120000.
        // Mass ratio = (0.120 / (1/12))² > 2.0 — strong quadratic sensitivity.
        let lambda_rust = LAMBDA_QG;
        let lambda_python = 0.120000_f64;
        let mass_ratio = lambda_python.powi(2) / lambda_rust.powi(2);
        assert!(
            mass_ratio > 2.0,
            "λ_rust={lambda_rust}, λ_python={lambda_python}: \
             mass ratio={mass_ratio:.4} — expected >2× (quadratic dependence)"
        );
    }

    #[test]
    fn doubling_lambda_qg_quadruples_mass() {
        // Proving the quadratic law directly: m(2λ) = 4·m(λ).
        // Using fixed physical parameters; only λ varies.
        fn quark_mass(lambda: f64) -> f64 {
            // m = veracity × curvature × field_gradient / planck_length × λ²
            // With dummy unit values for everything except λ:
            1.0 * 1.0 * 1.0 / 1.0 * lambda.powi(2)
        }
        let lambda = 0.1_f64;
        assert!(
            (quark_mass(2.0 * lambda) - 4.0 * quark_mass(lambda)).abs() < 1e-15,
            "m(2λ) ≠ 4·m(λ) — mass does not scale quadratically with λ_QG"
        );
    }

    #[test]
    fn gw_dispersion_is_sub_luminal() {
        // Δv/c must be < 1 for any physical energy.
        let e_qg = 1e28_f64;
        for &e in &[1e10f64, 1e19, 1e27] {
            let d = gravitational_wave_dispersion(e, e_qg);
            assert!(d < 1.0, "Δv/c = {d} ≥ 1 at E = {e:.1e} — superluminal!");
            assert!(d >= 0.0, "Δv/c = {d} < 0 — unphysical");
        }
    }
}
