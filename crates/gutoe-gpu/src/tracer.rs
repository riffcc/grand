// GUTOE Geodesic Ray Tracer — null geodesics in the GUTOE Schwarzschild metric
// Copyright (C) 2026 Riff Labs, AGPL-3.0-or-later
//
// Traces photon paths through the GUTOE-corrected Schwarzschild metric, producing
// black hole images with the SC lattice corrections from Cl(1,3).
//
// ── Physics ───────────────────────────────────────────────────────────────────
//
// GUTOE correction: coordinate r → areal radius r_eff = √(r² + r_core²)
// where r_core = √C_∞ × l_P resolves the classical r = 0 singularity.
//
// Null geodesic conserved quantities (metric in equatorial plane θ = π/2):
//   E = (1 − r_s/r_eff) × ṫ            (energy)
//   L = r_eff² × φ̇                     (angular momentum)
//   b = L/E                              (impact parameter)
//
// From the null condition 0 = g_μν ẋ^μ ẋ^ν:
//   (dr/dφ)² = r²r_eff²(1/b² − V)      [H(r)]
// where V(r) = (1 − r_s/r_eff) / r_eff²
//
// Orbit ODE via d²r/dφ² = (1/2) dH/dr:
//   d²r/dφ² = r(2r²+r_c²)/b² − r + r_s·r(r²+2r_c²)/(2r_eff³)
//
// Limit r_c → 0: reduces exactly to the Schwarzschild Binet equation
//   d²u/dφ² + u = (3r_s/2)u²   (u = 1/r) ✓
//
// Photon sphere: r_eff = 3r_s/2, b_crit = 3√3r_s/2 (same as GR in areal radius)
//
// ── Integration ───────────────────────────────────────────────────────────────
//
// State: (r, p = dr/dφ), integrated with 4th-order Runge–Kutta in orbital angle φ.
// Disk crossings (z = 0) detected at φ = n×π for tilted orbits.
// For equatorial orbits (bᵧ = 0): disk hit detected on first inward crossing.
//
// ── Rendering ─────────────────────────────────────────────────────────────────
//
// Orthographic camera at infinity, edge-on by default (inclination = 90°).
// Each screen pixel (sₓ, sᵧ) maps to impact parameters (bₓ = sₓ, bᵧ = sᵧ sin θ_obs).
// Result coloured by disk temperature profile T ∝ (r_ISCO/r_eff)^(3/4).

use std::f64::consts::PI;

use crate::kerr::KerrMetric;
use crate::metric::GutoeMetric;

// ── Critical impact parameter ─────────────────────────────────────────────────

/// b_crit = (3√3/2) × r_s: impact parameter of the unstable photon sphere orbit.
///
/// Photons with b < b_crit are captured by the black hole.
/// Photons with b > b_crit escape (possibly after transiting the disk).
///
/// This value is the same as GR: the photon sphere sits at areal radius 3r_s/2
/// in both theories (the GUTOE correction shifts only the *coordinate* radius).
pub fn b_critical(r_s: f64) -> f64 {
    1.5 * 3.0_f64.sqrt() * r_s
}

// ── Orbit equation ────────────────────────────────────────────────────────────

/// d²r/dφ² for a null geodesic in the GUTOE Schwarzschild metric.
///
/// Derivation: H(r) = (dr/dφ)² = r²r_eff²(1/b² − V).
///   dH/dr = 2r(2r²+r_c²)/b² − 2r + r_s·r(r²+2r_c²)/r_eff³
///   d²r/dφ² = (1/2) dH/dr
///
/// For r_c = 0 and u = 1/r this gives d²u/dφ² + u = (3r_s/2)u² ✓
fn orbit_accel(r: f64, b: f64, r_s: f64, r_c: f64) -> f64 {
    let re2 = r * r + r_c * r_c;
    let re3 = re2 * re2.sqrt();
    r * (2.0 * r * r + r_c * r_c) / (b * b)
        - r
        + r_s * r * (r * r + 2.0 * r_c * r_c) / (2.0 * re3)
}

/// (dr/dφ)² = r²·r_eff²·(1/b² − V) where V = (1 − r_s/r_eff)/r_eff².
///
/// Positive in the classically allowed region, zero at turning points,
/// and negative (unphysical) in the forbidden zone near the photon sphere.
fn orbit_vr_sq(r: f64, b: f64, r_s: f64, r_c: f64) -> f64 {
    let re2 = r * r + r_c * r_c;
    let re = re2.sqrt();
    let f = 1.0 - r_s / re;
    r * r * re2 / (b * b) - r * r * f
}

/// 4th-order Runge–Kutta step for the system (dr/dφ = p, dp/dφ = orbit_accel(r)).
fn rk4_step(r: f64, p: f64, b: f64, r_s: f64, r_c: f64, dphi: f64) -> (f64, f64) {
    // The orbit_accel depends only on r (not p), so the system is:
    //   dr/dφ = p
    //   dp/dφ = A(r)
    let a = |ri: f64| orbit_accel(ri, b, r_s, r_c);

    let k1r = p;
    let k1p = a(r);
    let k2r = p + 0.5 * dphi * k1p;
    let k2p = a(r + 0.5 * dphi * k1r);
    let k3r = p + 0.5 * dphi * k2p;
    let k3p = a(r + 0.5 * dphi * k2r);
    let k4r = p + dphi * k3p;
    let k4p = a(r + dphi * k3r);

    (
        r + dphi * (k1r + 2.0 * k2r + 2.0 * k3r + k4r) / 6.0,
        p + dphi * (k1p + 2.0 * k2p + 2.0 * k3p + k4p) / 6.0,
    )
}

// ── Trace result ──────────────────────────────────────────────────────────────

/// Result of tracing a single null geodesic through the GUTOE metric.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TraceResult {
    /// Photon absorbed by the black hole (areal radius dropped below horizon).
    Captured,

    /// Photon intercepted the accretion disk.
    DiskHit {
        /// Areal radius r_eff at the disk crossing.
        r_eff: f64,
        /// Total orbital angle at the hit (φ = n×π for tilted orbits).
        phi_orb: f64,
        /// Which equatorial crossing this is: 1 = far side, 2 = near side again, …
        n_cross: u32,
    },

    /// Photon escaped to infinity.
    Escaped {
        /// Total orbital angle traversed (deflection ≈ phi_total − π for nearly straight).
        phi_total: f64,
    },
}

// ── Trace function ────────────────────────────────────────────────────────────

/// Trace a null geodesic through the GUTOE Schwarzschild metric.
///
/// # Geometry
///
/// The orthographic camera is at infinity, looking toward the black hole.
/// Screen pixel (sₓ, sᵧ) maps to impact parameters:
///   bₓ = sₓ             (equatorial, horizontal in screen)
///   bᵧ = sᵧ × sin_obs   (polar, vertical scaled by observer inclination)
///
/// The orbit lies in a plane inclined by angle i = arcsin(bᵧ/b) from the
/// equatorial plane. Disk crossings (z = 0) occur at orbital angle φ = n×π.
///
/// # Arguments
/// - `metric`:        GUTOE metric (r_s, l_planck)
/// - `disk_inner_re`: inner disk boundary, areal radius (≥ r_ISCO = 3 r_s)
/// - `disk_outer_re`: outer disk boundary, areal radius (e.g. 10 r_s)
/// - `bx`:            equatorial impact parameter (same units as r_s)
/// - `by`:            polar impact parameter
/// - `max_phi`:       maximum orbital angle to integrate (radians; 20π = 10 half-orbits)
/// - `dphi`:          RK4 step size (radians; 0.01 gives good accuracy)
pub fn trace_photon(
    metric: &GutoeMetric,
    disk_inner_re: f64,
    disk_outer_re: f64,
    bx: f64,
    by: f64,
    max_phi: f64,
    dphi: f64,
) -> TraceResult {
    let r_s = metric.r_s;
    let r_c = metric.r_core();
    let b = (bx * bx + by * by).sqrt();
    if b < 1e-12 {
        return TraceResult::Captured; // head-on — definitely falls in
    }

    // Deep-shadow shortcut: any photon with b ≪ b_crit is inevitably captured.
    // Avoids integrating numerically unstable orbits that plunge straight into the horizon.
    if b < 0.5 * b_critical(r_s) {
        return TraceResult::Captured;
    }

    let sin_i = by / b; // sine of orbital plane inclination from equatorial plane
    let is_equatorial = sin_i.abs() < 1e-6;

    // Start at r_start = 3b. Stability analysis: orbit_accel ≈ 2r³/b², p ≈ r²/b at large r.
    // Change per step = accel·dphi/|p| ≈ 2r·dphi/b. At r=3b: ratio = 6·dphi = 0.06 (6%).
    // This is scale-independent and numerically stable for all b > 0.5·b_crit.
    let r_start = 3.0 * b;

    // Initial radial velocity: ingoing (p < 0), magnitude from orbit equation
    let vr0_sq = orbit_vr_sq(r_start, b, r_s, r_c);
    let p_start = if vr0_sq > 0.0 { -vr0_sq.sqrt() } else { -r_start * r_start / b };

    // Capture boundary: areal radius r_eff = sqrt(r² + r_c²) drops below the horizon.
    // The event horizon is at areal radius r_s (not the coordinate radius r_horizon).
    // Using metric.r_s directly gives the correct areal comparison with re_new.
    let r_capture_re = metric.r_s * 0.99;

    let mut r = r_start;
    let mut p = p_start;
    let mut phi = 0.0_f64;
    let mut n_cross = 0_u32; // number of equatorial crossings completed (first is at φ = π)
    let mut turned = false; // photon has passed its radial turning point

    // For equatorial orbit detection: track whether we are inside the disk
    let mut in_disk_eq = false;

    let max_steps = (max_phi / dphi).ceil() as usize + 1;

    for _step in 0..max_steps {
        let (r_new, p_rk4) = rk4_step(r, p, b, r_s, r_c, dphi);
        // Enforce orbital constraint p² = orbit_vr_sq(r). Keeps integration on the
        // constraint surface, preventing centrifugal blowup at large r. Preserves the
        // sign (direction) from RK4: for b < b_crit, orbit_vr_sq > 0 everywhere so p
        // stays negative (ingoing) → correct capture. At the turning point for b > b_crit,
        // orbit_vr_sq → 0 → p → 0, RK4 sign flips to positive → turned fires correctly.
        let vr2_new = orbit_vr_sq(r_new, b, r_s, r_c).max(0.0);
        let p_new = if p_rk4 >= 0.0 { vr2_new.sqrt() } else { -vr2_new.sqrt() };
        let phi_new = phi + dphi;

        // Areal radius at the new position
        let re_new = (r_new * r_new + r_c * r_c).sqrt();

        // Capture: areal radius dropped below horizon (or BH is sub-Planckian)
        if re_new < r_capture_re || r_new < r_c * 0.01 {
            return TraceResult::Captured;
        }

        // Turning point: radial motion reverses from ingoing to outgoing
        if !turned && p < 0.0 && p_new >= 0.0 {
            turned = true;
        }

        // Escape: photon has turned and returned to the starting radius
        if turned && r_new >= r_start * 0.99 {
            return TraceResult::Escaped { phi_total: phi_new };
        }

        // ── Disk hit detection ──────────────────────────────────────────────

        if is_equatorial {
            // Equatorial orbit (sin_i ≈ 0): photon stays in z = 0.
            // Disk hit when r_eff first enters the disk zone going inward.
            let re_cur = (r * r + r_c * r_c).sqrt();
            let cur_in_disk = re_cur >= disk_inner_re && re_cur <= disk_outer_re;
            if !in_disk_eq && cur_in_disk && p < 0.0 {
                return TraceResult::DiskHit { r_eff: re_cur, phi_orb: phi, n_cross: 1 };
            }
            if !cur_in_disk {
                in_disk_eq = false;
            }
        } else {
            // Tilted orbit: disk crossings at φ = n×π (equatorial plane crossings).
            // n_cross = 0 → first crossing target is π (far side from camera).
            let target = (n_cross as f64 + 1.0) * PI;
            if phi < target && phi_new >= target {
                // Linear interpolation to find r at the exact crossing
                let t = (target - phi) / dphi;
                let r_cross = r + t * (r_new - r);
                let re_cross = (r_cross * r_cross + r_c * r_c).sqrt();
                n_cross += 1;

                if re_cross >= disk_inner_re && re_cross <= disk_outer_re {
                    return TraceResult::DiskHit {
                        r_eff: re_cross,
                        phi_orb: target,
                        n_cross,
                    };
                }
            }
        }

        r = r_new;
        p = p_new;
        phi = phi_new;
    }

    // Ran out of steps. Classify by current position.
    if r >= r_start * 0.5 {
        TraceResult::Escaped { phi_total: phi }
    } else {
        TraceResult::Captured
    }
}

/// Experimental Kerr null-geodesic tracer (exterior camera, CPU path).
///
/// Integrates Carter first-order equations in Boyer-Lindquist form using an
/// affine-like step parameter. This is intended as a physically grounded spin
/// path (not a shader warp), but still an experimental integrator.
pub fn trace_photon_kerr(
    kerr: &KerrMetric,
    disk_inner_re: f64,
    disk_outer_re: f64,
    bx: f64,
    by: f64,
    inclination_deg: f64,
    max_lambda: f64,
    dlambda: f64,
) -> TraceResult {
    let r_s = kerr.r_s;
    let (r_plus, _) = kerr.horizons();
    // Kerr image constants must use the true observer inclination.
    // The renderer passes Kerr screen coordinates as (alpha,beta) without
    // pre-applying inclination to beta, so we map with theta_obs directly.
    let theta_obs = inclination_deg.to_radians().clamp(1e-4, PI - 1e-4);
    let (xi, eta) = kerr.image_to_constants(bx, by, theta_obs);

    // Start far away so image-plane mapping approximates asymptotic observer.
    let b = (bx * bx + by * by).sqrt();
    let r_start = (40.0 * r_s).max(12.0 * b).max(20.0);

    let mut r = r_start;
    let mut theta = theta_obs;
    let mut phi = 0.0_f64;
    let mut n_cross = 0_u32;

    // Ingoing branch from observer.
    let mut sgn_r = -1.0_f64;
    // Initial polar direction from image-plane beta sign.
    // In our screen convention, +beta is toward decreasing θ (toward +z),
    // so the dθ/dλ sign is opposite beta.
    let mut sgn_th = if by >= 0.0 { -1.0 } else { 1.0 };

    let mut lambda = 0.0_f64;
    let max_steps = (max_lambda / dlambda).ceil() as usize + 1;

    for _ in 0..max_steps {
        let a = kerr.a();
        let sin_th = theta.sin();
        let sin2 = (sin_th * sin_th).max(1e-9);
        let sigma = kerr.sigma(r, theta).max(1e-12);
        let delta = kerr.delta(r);
        let p = (r * r + a * a) - a * xi;

        let rpot = kerr.radial_potential(r, xi, eta);
        if rpot <= 1e-12 {
            sgn_r = -sgn_r;
        }
        let tpot = kerr.polar_potential(theta, xi, eta);
        if tpot <= 1e-12 {
            sgn_th = -sgn_th;
        }

        let rr = rpot.max(0.0).sqrt();
        let thh = tpot.max(0.0).sqrt();

        let dr = sgn_r * rr / sigma;
        let dth = sgn_th * thh / sigma;

        // Carter azimuth equation: Σ dφ/dλ = ξ csc²θ - a + aP/Δ
        // Clamp near Δ→0 to keep integration stable at horizon edge.
        let dphi = if delta.abs() < 1e-9 {
            (xi / sin2 - a) / sigma
        } else {
            (xi / sin2 - a + a * p / delta) / sigma
        };

        let r_new = r + dlambda * dr;
        let th_new = (theta + dlambda * dth).clamp(1e-4, PI - 1e-4);
        let phi_new = phi + dlambda * dphi;

        if !r_new.is_finite() || !th_new.is_finite() || !phi_new.is_finite() {
            return TraceResult::Captured;
        }

        // Capture: crossed outer horizon.
        if r_new <= r_plus * 1.001 {
            return TraceResult::Captured;
        }

        // Escape: returned to asymptotic radius on outgoing branch.
        if sgn_r > 0.0 && r_new >= r_start * 0.995 {
            return TraceResult::Escaped { phi_total: phi_new.abs() };
        }

        // Disk crossing at θ = π/2.
        let pi2 = PI * 0.5;
        if (theta - pi2) * (th_new - pi2) <= 0.0 {
            n_cross += 1;
            let t = ((pi2 - theta) / (th_new - theta)).clamp(0.0, 1.0);
            let r_cross = r + t * (r_new - r);
            if r_cross >= disk_inner_re && r_cross <= disk_outer_re {
                return TraceResult::DiskHit {
                    r_eff: r_cross,
                    phi_orb: phi_new.abs(),
                    n_cross,
                };
            }
        }

        r = r_new;
        theta = th_new;
        phi = phi_new;
        lambda += dlambda;
        if lambda >= max_lambda {
            break;
        }
    }

    if sgn_r > 0.0 && r > 0.5 * r_start {
        TraceResult::Escaped { phi_total: phi.abs() }
    } else if r > (3.0 * r_s).max(1.2 * r_plus) {
        // Conservative timeout classification: still far outside strong-field region.
        TraceResult::Escaped { phi_total: phi.abs() }
    } else {
        TraceResult::Captured
    }
}

// ── Interior-camera trace ─────────────────────────────────────────────────────

/// Trace a null geodesic fired **outward** from a camera inside the event horizon.
///
/// # Physics
///
/// From inside the horizon (r_cam < r_horizon), the potential V(r_cam) < 0, so
/// `orbit_vr_sq(r_cam, b, …) > 0` for every finite b.  Every photon can start moving
/// outward (p_start > 0).
///
/// - `b < b_crit`:  photon clears the photon-sphere barrier → escapes to infinity
///                  (sees the outside universe, accretion disk, star field)
/// - `b > b_crit`:  photon reaches a turning point before the photon sphere,
///                  reverses, falls back → returns `Captured` (GUTOE core)
///
/// # Arguments
/// - `r_cam`: coordinate radius of the camera (must satisfy r_cam < metric.r_horizon())
pub fn trace_photon_interior(
    metric: &GutoeMetric,
    disk_inner_re: f64,
    disk_outer_re: f64,
    r_cam: f64,
    bx: f64,
    by: f64,
    max_phi: f64,
    dphi: f64,
) -> TraceResult {
    let r_s = metric.r_s;
    let r_c = metric.r_core();
    let b   = (bx * bx + by * by).sqrt();

    // Pure-radial photon: always escapes
    if b < 1e-12 {
        return TraceResult::Escaped { phi_total: max_phi };
    }

    // Inside the horizon V < 0 → orbit_vr_sq > 0 always
    let vr0_sq  = orbit_vr_sq(r_cam, b, r_s, r_c).max(0.0);
    let p_start = vr0_sq.sqrt();   // outward (positive)

    // Camera areal radius: capture returned photons once back below this level
    let re_cam  = (r_cam * r_cam + r_c * r_c).sqrt();
    let re_cap  = re_cam * 1.05;

    // Escape once past the outer disk / strong-gravity region
    let r_escape = (3.0 * b)
        .max(disk_outer_re * 1.5)
        .max(20.0 * r_s);

    let sin_i        = by / b;
    let is_equatorial = sin_i.abs() < 1e-6;

    let mut r   = r_cam;
    let mut p   = p_start;
    let mut phi = 0.0_f64;
    let mut n_cross = 0_u32;
    let mut turned  = false;

    let max_steps = (max_phi / dphi).ceil() as usize + 1;

    for _step in 0..max_steps {
        let (r_new, p_rk4) = rk4_step(r, p, b, r_s, r_c, dphi);

        let vr2_new = orbit_vr_sq(r_new, b, r_s, r_c).max(0.0);
        let p_new   = if p_rk4 >= 0.0 { vr2_new.sqrt() } else { -vr2_new.sqrt() };
        let phi_new = phi + dphi;
        let re_new  = (r_new * r_new + r_c * r_c).sqrt();

        // Turning point: radial motion reversed from outward to inward
        if !turned && p > 0.0 && p_new <= 0.0 {
            turned = true;
        }

        // After turning: capture once back below the camera's areal radius
        if turned && re_new < re_cap {
            return TraceResult::Captured;
        }

        // Escape: risen past strong-gravity region
        if !turned && r_new >= r_escape {
            return TraceResult::Escaped { phi_total: phi_new };
        }

        // Disk hit detection (same φ = nπ crossing logic as exterior trace)
        if is_equatorial {
            let re_cur = (r * r + r_c * r_c).sqrt();
            if !turned && re_cur >= disk_inner_re && re_cur <= disk_outer_re && p > 0.0 {
                return TraceResult::DiskHit { r_eff: re_cur, phi_orb: phi, n_cross: 1 };
            }
        } else {
            let target = (n_cross as f64 + 1.0) * PI;
            if phi < target && phi_new >= target {
                let t       = (target - phi) / dphi;
                let r_cross = r + t * (r_new - r);
                let re_cross = (r_cross * r_cross + r_c * r_c).sqrt();
                n_cross += 1;
                if re_cross >= disk_inner_re && re_cross <= disk_outer_re {
                    return TraceResult::DiskHit {
                        r_eff: re_cross,
                        phi_orb: target,
                        n_cross,
                    };
                }
            }
        }

        r   = r_new;
        p   = p_new;
        phi = phi_new;
    }

    // Timed out
    if r >= r_escape * 0.5 && !turned {
        TraceResult::Escaped { phi_total: phi }
    } else {
        TraceResult::Captured
    }
}

/// Trace a null geodesic fired **toward the core** from a camera inside the horizon.
///
/// This is the companion to `trace_photon_interior` (which fires outward).  It models an
/// interior observer turning around to look down at the regularized lattice core.
///
/// For future-directed null rays inside the horizon this is a plunging branch: the ray is
/// initialized with inward radial momentum and integrated until it reaches the core shell.
pub fn trace_photon_interior_core(
    metric: &GutoeMetric,
    r_cam: f64,
    bx: f64,
    by: f64,
    max_phi: f64,
    dphi: f64,
) -> TraceResult {
    let r_s = metric.r_s;
    let r_c = metric.r_core();
    let b = (bx * bx + by * by).sqrt();

    // Radial center pixel: direct plunge.
    if b < 1e-12 {
        return TraceResult::DiskHit { r_eff: r_c, phi_orb: 0.0, n_cross: 1 };
    }

    // Start at camera radius with inward radial momentum.
    let vr0_sq = orbit_vr_sq(r_cam, b, r_s, r_c).max(0.0);
    let p_start = -vr0_sq.sqrt();

    // Capture at the regularized core shell.
    let re_core_cap = (1.02 * r_c).max(r_c + 1e-9);

    let mut r = r_cam;
    let mut p = p_start;
    let mut phi = 0.0_f64;
    let max_steps = (max_phi / dphi).ceil() as usize + 1;

    for _step in 0..max_steps {
        let (r_new, p_rk4) = rk4_step(r, p, b, r_s, r_c, dphi);
        let vr2_new = orbit_vr_sq(r_new, b, r_s, r_c).max(0.0);
        let p_new = if p_rk4 >= 0.0 { vr2_new.sqrt() } else { -vr2_new.sqrt() };
        let phi_new = phi + dphi;
        let re_new = (r_new * r_new + r_c * r_c).sqrt();

        if re_new <= re_core_cap || r_new <= r_c * 0.01 {
            return TraceResult::DiskHit { r_eff: re_new, phi_orb: phi_new, n_cross: 1 };
        }

        // Defensive escape classification for numerical edge-cases.
        if p < 0.0 && p_new >= 0.0 && r_new > r_cam * 1.01 {
            return TraceResult::Escaped { phi_total: phi_new };
        }

        r = r_new;
        p = p_new;
        phi = phi_new;
    }

    TraceResult::DiskHit { r_eff: re_core_cap, phi_orb: phi, n_cross: 1 }
}

// ── Render config ─────────────────────────────────────────────────────────────

/// Configuration for rendering a black hole image.
#[derive(Debug, Clone)]
pub struct RenderConfig {
    /// Image width in pixels.
    pub width: usize,
    /// Image height in pixels.
    pub height: usize,
    /// Half-width of the image in units of r_s.
    /// e.g. fov_rs = 12.0 → the image spans ±12 r_s horizontally.
    pub fov_rs: f64,
    /// Observer inclination from the disk normal, degrees.
    /// 90° = edge-on (see the disk ring), 0° = face-on (see the disk face).
    pub inclination_deg: f64,
    /// Maximum orbital angle in radians. 20π allows ~10 half-orbits.
    pub max_phi: f64,
    /// RK4 step size in radians. 0.01 gives good accuracy and speed.
    pub dphi: f64,
}

impl Default for RenderConfig {
    fn default() -> Self {
        Self {
            width: 400,
            height: 400,
            fov_rs: 12.0,
            inclination_deg: 80.0,
            max_phi: 20.0 * PI,
            dphi: 0.01,
        }
    }
}

impl RenderConfig {
    /// Quick 200×200 preview with coarser integration.
    pub fn preview() -> Self {
        Self {
            width: 200,
            height: 200,
            fov_rs: 12.0,
            inclination_deg: 80.0,
            max_phi: 15.0 * PI,
            dphi: 0.02,
        }
    }

    /// High-quality 800×800 with fine integration.
    pub fn high_quality() -> Self {
        Self {
            width: 800,
            height: 800,
            fov_rs: 12.0,
            inclination_deg: 80.0,
            max_phi: 30.0 * PI,
            dphi: 0.005,
        }
    }
}

// ── Color functions ───────────────────────────────────────────────────────────

/// Disk temperature colour using a simplified Novikov–Thorne profile.
///
/// T_rel = (r_ISCO / r_eff)^(3/4): inner disk is hotter and brighter.
/// Higher-order images (n_cross > 1) are dimmed by 0.7^(n-1).
/// Colour gradient: deep orange (outer/cool) → yellow → white (inner/hot).
fn disk_color(r_eff: f64, r_isco: f64, n_cross: u32) -> [u8; 3] {
    let t_rel = (r_isco / r_eff).powf(0.75).clamp(0.005, 1.0);
    let fade = 0.7_f64.powi(n_cross as i32 - 1);
    let b = (t_rel * fade).clamp(0.0, 1.0); // brightness [0, 1]

    // Orange → yellow → white as brightness increases
    let r = (255.0 * b.powf(0.4)).clamp(0.0, 255.0) as u8;
    let g = (200.0 * b.powf(0.7)).clamp(0.0, 255.0) as u8;
    let bl = (120.0 * b.powf(1.8)).clamp(0.0, 255.0) as u8;
    [r, g, bl]
}

/// Dark blue-black background (night sky).
fn background_color() -> [u8; 3] {
    [5, 5, 20]
}

/// Pure black shadow.
fn shadow_color() -> [u8; 3] {
    [0, 0, 0]
}

// ── Main render function ──────────────────────────────────────────────────────

/// Render a GUTOE black hole image via backward ray-tracing.
///
/// Returns a flat `Vec<[u8; 3]>` of RGB pixels in row-major order (top to bottom).
/// Pixel (ix, iy) is at index `iy * cfg.width + ix`.
///
/// # Camera geometry
///
/// Orthographic (parallel rays) camera at inclination `cfg.inclination_deg` above
/// the disk plane. Screen x-axis aligns with the equatorial direction; screen y-axis
/// is tilted toward the disk normal.
///
/// Impact parameters:
///   bx = sx = (ix − (W−1)/2) × scale × r_s
///   by = sy × sin(inclination)
///
/// For edge-on (90°): by = sy, gives the full disk ring structure.
/// For face-on (0°): by = 0, all orbits equatorial, shows circular shadow only.
pub fn render(
    metric: &GutoeMetric,
    disk_inner_rs: f64,
    disk_outer_rs: f64,
    cfg: &RenderConfig,
) -> Vec<[u8; 3]> {
    let r_s = metric.r_s;
    let disk_inner_re = disk_inner_rs * r_s;
    let disk_outer_re = disk_outer_rs * r_s;
    let r_isco = 3.0 * r_s; // areal ISCO radius

    let sin_inc = cfg.inclination_deg.to_radians().sin();

    // Pixel size in units of r_s
    let scale = 2.0 * cfg.fov_rs * r_s / cfg.width as f64;

    let mut pixels = vec![[0u8; 3]; cfg.width * cfg.height];

    for iy in 0..cfg.height {
        for ix in 0..cfg.width {
            // Screen coordinates centred at (0, 0), y-axis upward
            let sx = (ix as f64 - 0.5 * (cfg.width as f64 - 1.0)) * scale;
            let sy = (0.5 * (cfg.height as f64 - 1.0) - iy as f64) * scale;

            // Impact parameters: bx horizontal (equatorial), by vertical (polar)
            let bx = sx;
            let by = sy * sin_inc;

            let result = trace_photon(
                metric,
                disk_inner_re,
                disk_outer_re,
                bx,
                by,
                cfg.max_phi,
                cfg.dphi,
            );

            pixels[iy * cfg.width + ix] = match result {
                TraceResult::Captured => shadow_color(),
                TraceResult::DiskHit { r_eff, n_cross, .. } => {
                    disk_color(r_eff, r_isco, n_cross)
                }
                TraceResult::Escaped { .. } => background_color(),
            };
        }
    }

    pixels
}

// ── PPM output ────────────────────────────────────────────────────────────────

/// Encode pixels as a binary PPM (P6) byte vector.
///
/// Write the result to a file with `std::fs::write(path, write_ppm(...))`.
pub fn write_ppm(pixels: &[[u8; 3]], width: usize, height: usize) -> Vec<u8> {
    let header = format!("P6\n{width} {height}\n255\n");
    let mut out = Vec::with_capacity(header.len() + pixels.len() * 3);
    out.extend_from_slice(header.as_bytes());
    for px in pixels {
        out.extend_from_slice(px);
    }
    out
}

/// Encode pixels as ASCII PPM (P3) string (useful for debugging small images).
pub fn write_ppm_ascii(pixels: &[[u8; 3]], width: usize, height: usize) -> String {
    let mut s = format!("P3\n{width} {height}\n255\n");
    for px in pixels {
        s.push_str(&format!("{} {} {}\n", px[0], px[1], px[2]));
    }
    s
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kerr::KerrMetric;
    use crate::metric::GutoeMetric;

    const R_S: f64 = 1.0; // unit Schwarzschild radius

    // ── Orbit equation consistency ──────────────────────────────────────────

    #[test]
    fn b_critical_is_three_sqrt_three_half_rs() {
        let b = b_critical(R_S);
        let expected = 1.5 * 3.0_f64.sqrt() * R_S;
        assert!(
            (b - expected).abs() < 1e-12,
            "b_crit = {b:.10}, expected 3√3/2 × r_s = {expected:.10}"
        );
        // Numerical check: 3√3/2 ≈ 2.598076
        assert!((b - 2.598_076_211).abs() < 1e-6);
    }

    #[test]
    fn orbit_accel_schwarzschild_limit() {
        // For r_c = 0: d²r/dφ² = 2r³/b² − r + r_s/2
        let r = 5.0 * R_S;
        let b = 4.0 * R_S;
        let a = orbit_accel(r, b, R_S, 0.0);
        let expected = 2.0 * r * r * r / (b * b) - r + R_S / 2.0;
        assert!(
            (a - expected).abs() < 1e-10,
            "Schwarzschild limit: got {a:.10}, expected {expected:.10}"
        );
    }

    #[test]
    fn orbit_vr_sq_zero_at_photon_sphere() {
        // For r_c = 0, b = b_crit, r = 3r_s/2: (dr/dφ)² = 0 (photon sphere)
        let b = b_critical(R_S);
        let r_ph = 1.5 * R_S;
        let vr2 = orbit_vr_sq(r_ph, b, R_S, 0.0);
        assert!(
            vr2.abs() < 1e-10,
            "(dr/dφ)² at photon sphere = {vr2:.2e}, expected ≈ 0"
        );
    }

    #[test]
    fn orbit_accel_zero_at_photon_sphere() {
        // At the photon sphere with b = b_crit: d²r/dφ² = 0 (circular orbit)
        let b = b_critical(R_S);
        let r_ph = 1.5 * R_S;
        let a = orbit_accel(r_ph, b, R_S, 0.0);
        assert!(
            a.abs() < 1e-10,
            "d²r/dφ² at photon sphere = {a:.2e}, expected 0"
        );
    }

    #[test]
    fn gutoe_orbit_accel_differs_from_gr() {
        // r_core = √C_∞ ≈ 0.739 l_P modifies the orbit dynamics
        let r_c = 0.739; // r_core for r_s = l_P = 1
        let b = b_critical(R_S);
        let r = 1.5 * R_S; // GR photon sphere radius
        let a_gr = orbit_accel(r, b, R_S, 0.0);
        let a_gut = orbit_accel(r, b, R_S, r_c);
        assert!(
            (a_gut - a_gr).abs() > 1e-6,
            "GUTOE accel {a_gut:.6} should differ from GR {a_gr:.6} at r_c = {r_c}"
        );
    }

    // ── Capture and escape ──────────────────────────────────────────────────

    #[test]
    fn head_on_photon_is_captured() {
        let m = GutoeMetric::planck_units(R_S);
        let r = trace_photon(&m, 3.0, 10.0, 0.0, 0.0, 20.0 * PI, 0.01);
        assert_eq!(r, TraceResult::Captured, "b=0 must be captured");
    }

    #[test]
    fn subcritical_photon_is_captured() {
        // No disk (inner > outer) so the disk does not intercept the photon before it is captured.
        let m = GutoeMetric::planck_units(R_S);
        let b_crit = b_critical(R_S);
        let r = trace_photon(&m, 1.0, 0.0, 0.9 * b_crit, 0.0, 30.0 * PI, 0.01);
        assert_eq!(r, TraceResult::Captured, "b = 0.9 b_crit must be captured");
    }

    #[test]
    fn supercritical_photon_escapes() {
        // b = 10.5 r_s >> b_crit ≈ 2.6 r_s: turning point ≈ 10 r_s, escapes to infinity.
        // No disk (inner > outer) so we test escape physics cleanly.
        let m = GutoeMetric::planck_units(R_S);
        let r = trace_photon(&m, 1.0, 0.0, 10.5 * R_S, 0.0, 20.0 * PI, 0.01);
        match r {
            TraceResult::Escaped { .. } => {}
            _ => panic!("b = 10.5 r_s must escape, got {r:?}"),
        }
    }

    #[test]
    fn photon_sphere_boundary() {
        // b = 0.995 b_crit → captured
        // b = 1.1  b_crit → escapes or hits disk (definitely not captured forever)
        let m = GutoeMetric::planck_units(R_S);
        let b_crit = b_critical(R_S);

        // No disk (inner > outer): test pure capture/escape, no disk interception.
        let below = trace_photon(&m, 1.0, 0.0, 0.995 * b_crit, 0.0, 40.0 * PI, 0.005);
        assert_eq!(below, TraceResult::Captured, "b = 0.995 b_crit must be captured");

        let above = trace_photon(&m, 1.0, 0.0, 1.1 * b_crit, 0.0, 40.0 * PI, 0.005);
        match above {
            TraceResult::Captured => panic!("b = 1.1 b_crit must NOT be captured, got Captured"),
            _ => {} // Escaped or DiskHit both acceptable
        }
    }

    // ── Disk hits ───────────────────────────────────────────────────────────

    #[test]
    fn equatorial_photon_hits_disk() {
        // bx = 5 r_s is between ISCO (3 r_s) and disk outer (10 r_s)
        // The photon enters the disk outer edge on the way inward
        let m = GutoeMetric::planck_units(R_S);
        let r = trace_photon(&m, 3.0, 10.0, 5.0 * R_S, 0.0, 30.0 * PI, 0.01);
        match r {
            TraceResult::DiskHit { r_eff, .. } => {
                assert!(
                    r_eff >= 3.0 * R_S && r_eff <= 10.01 * R_S,
                    "disk hit at r_eff = {r_eff:.4}, expected in [3, 10] r_s"
                );
            }
            _ => panic!("equatorial b=5 r_s must hit disk, got {r:?}"),
        }
    }

    #[test]
    fn tilted_photon_hits_far_side_disk() {
        // bx = 3.5 r_s, by = 1 r_s: tilted orbit, b ≈ 3.6 r_s ≈ 1.4 b_crit.
        // The orbit bends ~90° total, so at φ = π (far side) r_eff ≈ 5.3 r_s — inside
        // the disk [3, 10].  (bx = 5 doesn't work: nearly-straight orbit exits the disk
        // at φ ≈ π/2 and the equatorial-plane crossing at φ = π is at r ≈ 40 r_s.)
        let m = GutoeMetric::planck_units(R_S);
        let r = trace_photon(&m, 3.0, 10.0, 3.5 * R_S, 1.0 * R_S, 30.0 * PI, 0.01);
        match r {
            TraceResult::DiskHit { r_eff, n_cross, .. } => {
                assert!(
                    r_eff >= 3.0 * R_S && r_eff <= 10.01 * R_S,
                    "far-side disk hit at r_eff = {r_eff:.4}, expected in [3, 10] r_s"
                );
                // The first crossing after start is n_cross = 1 (far side, φ = π)
                assert!(n_cross >= 1, "must be at least the first crossing");
            }
            _ => panic!("tilted bx=5 by=1 must hit disk, got {r:?}"),
        }
    }

    // ── Render ──────────────────────────────────────────────────────────────

    #[test]
    fn render_preview_pixel_count() {
        let m = GutoeMetric::planck_units(R_S);
        let cfg = RenderConfig::preview();
        let pixels = render(&m, 3.0, 10.0, &cfg);
        assert_eq!(
            pixels.len(),
            cfg.width * cfg.height,
            "render must return exactly width×height pixels"
        );
    }

    #[test]
    fn render_shadow_center_is_black() {
        // The center pixel (bx ≈ 0, by ≈ 0) is deep inside the shadow → [0, 0, 0]
        let m = GutoeMetric::planck_units(R_S);
        let cfg = RenderConfig {
            width: 100,
            height: 100,
            fov_rs: 12.0,
            inclination_deg: 90.0,
            max_phi: 20.0 * PI,
            dphi: 0.02,
        };
        let pixels = render(&m, 3.0, 10.0, &cfg);
        // Closest-to-centre pixel: (50, 50). b ≈ 0.24 r_s << b_crit ≈ 2.6 r_s → shadow.
        let px = pixels[50 * 100 + 50];
        assert_eq!(px, [0, 0, 0], "center pixel must be black (shadow), got {px:?}");
    }

    #[test]
    fn render_has_non_black_pixels() {
        // A correctly rendered image must have bright pixels (disk + escape photons)
        let m = GutoeMetric::planck_units(R_S);
        let cfg = RenderConfig {
            width: 80,
            height: 80,
            fov_rs: 12.0,
            inclination_deg: 80.0,
            max_phi: 20.0 * PI,
            dphi: 0.02,
        };
        let pixels = render(&m, 3.0, 10.0, &cfg);
        let non_black = pixels.iter().filter(|&&p| p != [0u8, 0, 0]).count();
        assert!(
            non_black > 0,
            "rendered image must have at least some non-black pixels (disk + background)"
        );
    }

    // ── Colour ──────────────────────────────────────────────────────────────

    #[test]
    fn disk_color_inner_brighter_than_outer() {
        let r_isco = 3.0 * R_S;
        let inner = disk_color(r_isco * 1.01, r_isco, 1);
        let outer = disk_color(r_isco * 3.0, r_isco, 1);
        let inner_lum: u32 = inner.iter().map(|&x| x as u32).sum();
        let outer_lum: u32 = outer.iter().map(|&x| x as u32).sum();
        assert!(
            inner_lum > outer_lum,
            "inner disk (lum={inner_lum}) must be brighter than outer (lum={outer_lum})"
        );
    }

    #[test]
    fn secondary_image_dimmer_than_primary() {
        let r_isco = 3.0 * R_S;
        let r_eff = 5.0 * R_S;
        let primary = disk_color(r_eff, r_isco, 1);
        let secondary = disk_color(r_eff, r_isco, 2);
        let p_lum: u32 = primary.iter().map(|&x| x as u32).sum();
        let s_lum: u32 = secondary.iter().map(|&x| x as u32).sum();
        assert!(
            p_lum > s_lum,
            "primary image (lum={p_lum}) must be brighter than secondary (lum={s_lum})"
        );
    }

    // ── PPM output ──────────────────────────────────────────────────────────

    #[test]
    fn write_ppm_valid_header_and_size() {
        let pixels = vec![[255u8, 0, 0]; 4]; // 2×2 all-red image
        let ppm = write_ppm(&pixels, 2, 2);
        // Check header prefix
        assert!(
            ppm.starts_with(b"P6\n"),
            "PPM must start with 'P6\\n'"
        );
        // Total size: "P6\n2 2\n255\n" = 11 bytes + 4 pixels × 3 = 23 bytes
        let expected_len = b"P6\n2 2\n255\n".len() + 4 * 3;
        assert_eq!(ppm.len(), expected_len, "PPM byte count mismatch");
        // Last 12 bytes are pixel data: RRGGBB × 4 = [255,0,0] × 4
        let pixel_bytes = &ppm[ppm.len() - 12..];
        assert_eq!(pixel_bytes, &[255, 0, 0, 255, 0, 0, 255, 0, 0, 255, 0, 0]);
    }

    #[test]
    fn write_ppm_ascii_valid() {
        let pixels = vec![[10u8, 20, 30]];
        let ppm = write_ppm_ascii(&pixels, 1, 1);
        assert!(ppm.starts_with("P3\n1 1\n255\n"));
        assert!(ppm.contains("10 20 30"));
    }

    #[test]
    fn kerr_tracer_large_b_escapes() {
        let k = KerrMetric::new(1.0, 0.6).expect("valid");
        // No disk interception: inner > outer.
        let r = trace_photon_kerr(&k, 2.0, 1.0, 12.0, 0.0, 70.0, 80.0, 0.01);
        match r {
            TraceResult::Escaped { .. } => {}
            _ => panic!("large-b Kerr ray should escape, got {r:?}"),
        }
    }

    #[test]
    fn kerr_tracer_captures_small_b() {
        let k = KerrMetric::new(1.0, 0.9).expect("valid");
        // No disk interception: inner > outer.
        let r = trace_photon_kerr(&k, 2.0, 1.0, 0.2, 0.0, 75.0, 80.0, 0.01);
        assert_eq!(r, TraceResult::Captured);
    }
}
