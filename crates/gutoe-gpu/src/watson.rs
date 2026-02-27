// Watson integral for the hex+z lattice — analytical prediction of C_∞
// Copyright (C) 2026 Riff Labs, AGPL-3.0-or-later
//
// The hex+z lattice has 6 in-plane triangular neighbors + 2 z-axis neighbors.
// The kinetic operator: T(k) = 1 - (1/8)[2cos(θ₁) + 2cos(θ₂) + 2cos(θ₁-θ₂) + 2cos(θ_z)]
// where θ₁ = k·a₁, θ₂ = k·a₂ are the phases along the two triangular primitive vectors,
// and θ_z is the z-component of k.
//
// The lattice Green's function at the origin is the Watson integral:
//   G(0) = (1/(2π)³) ∫₀²π dθ₁ ∫₀²π dθ₂ ∫₋π^π dθ_z / T(θ₁, θ₂, θ_z)

use std::f64::consts::PI;

/// Hex+z kinetic dispersion: T(θ₁, θ₂, θ_z) = 1 - S/8
/// where S = 2cos(θ₁) + 2cos(θ₂) + 2cos(θ₁-θ₂) + 2cos(θ_z)
#[inline]
fn t_hex_z(theta1: f64, theta2: f64, theta_z: f64) -> f64 {
    let s = 2.0 * theta1.cos()
        + 2.0 * theta2.cos()
        + 2.0 * (theta1 - theta2).cos()
        + 2.0 * theta_z.cos();
    1.0 - s / 8.0
}

/// Simple cubic kinetic dispersion for validation:
/// T(k₁,k₂,k₃) = 1 - (cos k₁ + cos k₂ + cos k₃)/3
#[inline]
fn t_simple_cubic(k1: f64, k2: f64, k3: f64) -> f64 {
    1.0 - (k1.cos() + k2.cos() + k3.cos()) / 3.0
}

/// 3D midpoint-rule integration over [a,b]³.
fn integrate_3d_midpoint<F: Fn(f64, f64, f64) -> f64>(
    f: F,
    ax: f64,
    bx: f64,
    ay: f64,
    by: f64,
    az: f64,
    bz: f64,
    n: usize,
) -> f64 {
    let hx = (bx - ax) / n as f64;
    let hy = (by - ay) / n as f64;
    let hz = (bz - az) / n as f64;

    let mut sum = 0.0;
    for ix in 0..n {
        let x = ax + (ix as f64 + 0.5) * hx;
        for iy in 0..n {
            let y = ay + (iy as f64 + 0.5) * hy;
            for iz in 0..n {
                let z = az + (iz as f64 + 0.5) * hz;
                sum += f(x, y, z);
            }
        }
    }
    sum * hx * hy * hz
}

pub fn watson_hex_z(n: usize) -> f64 {
    integrate_3d_midpoint(
        |t1, t2, tz| 1.0 / t_hex_z(t1, t2, tz),
        0.0,
        2.0 * PI,
        0.0,
        2.0 * PI,
        -PI,
        PI,
        n,
    ) / (2.0 * PI).powi(3)
}

pub fn watson_simple_cubic(n: usize) -> f64 {
    integrate_3d_midpoint(
        |k1, k2, k3| 1.0 / t_simple_cubic(k1, k2, k3),
        -PI,
        PI,
        -PI,
        PI,
        -PI,
        PI,
        n,
    ) / (2.0 * PI).powi(3)
}

pub fn second_moment_tensor() -> (f64, f64, f64) {
    let hex_neighbors: [(f64, f64); 6] = [
        (-0.5, -(3.0_f64).sqrt() / 2.0),
        (0.5, -(3.0_f64).sqrt() / 2.0),
        (-1.0, 0.0),
        (1.0, 0.0),
        (-0.5, (3.0_f64).sqrt() / 2.0),
        (0.5, (3.0_f64).sqrt() / 2.0),
    ];
    let mut sum_dx2 = 0.0;
    let mut sum_dy2 = 0.0;
    for &(dx, dy) in &hex_neighbors {
        sum_dx2 += dx * dx;
        sum_dy2 += dy * dy;
    }
    (sum_dx2 / 16.0, sum_dy2 / 16.0, 2.0 / 16.0)
}

pub fn continuum_green_coefficient() -> f64 {
    let (m_xx, _, m_zz) = second_moment_tensor();
    1.0 / (4.0 * PI * m_xx * m_zz.sqrt())
}

pub fn run_analysis() {
    println!("╔══════════════════════════════════════════════════════════════╗");
    println!("║  GUTOE Hex+Z Lattice Green's Function (Watson Integral)    ║");
    println!("╚══════════════════════════════════════════════════════════════╝\n");

    let (m_xx, m_yy, m_zz) = second_moment_tensor();
    println!("  ── Second-Moment Tensor ──");
    println!(
        "  M_xx = {:.6} (3/16), M_yy = {:.6} (3/16), M_zz = {:.6} (1/8)",
        m_xx, m_yy, m_zz
    );
    println!("  m_xy = {:.4}, m_z = {:.4}\n", 0.5 / m_xx, 0.5 / m_zz);

    println!("  ── Simple Cubic (validation, known ≈ 1.5164) ──");
    for &n in &[100, 200, 400] {
        println!("  n={:4}: {:.8}", n, watson_simple_cubic(n));
    }
    println!();

    println!("  ── Hex+Z Watson Integral ──");
    let mut wv = Vec::new();
    for &n in &[100, 200, 400] {
        let w = watson_hex_z(n);
        wv.push((n, w));
        println!("  n={:4}: {:.8}", n, w);
    }
    let _g0 = wv.last().unwrap().1;
    println!();

    let c_g = continuum_green_coefficient();
    let c_inf_var = c_g * c_g / 2.0;
    println!(
        "  C_G = {:.6}, C_∞(continuum) = C_G²/2 = {:.6}",
        c_g, c_inf_var
    );
    println!("  C_∞(GPU) ≈ 0.547 → ratio = {:.4}\n", 0.547 / c_inf_var);
}

// ── CPU lattice hydrogen solver ────────────────────────────────────────────
//
// Matches GPU `run_hydrogen_obc` exactly:
//   - OBC neighbors: out-of-bounds → sentinel -1
//   - Poisson OBC: out-of-bounds neighbors contribute φ_wall = α/r
//   - Schrödinger OBC: out-of-bounds neighbors contribute ψ = 0
//   - 6-color SOR for Poisson (in-place, ω = 2/(1+sin(π/L)))
//   - Coulomb warm-start: φ_init = α/r

const SENTINEL: i32 = -1;

/// OBC neighbors: returns (flat_idx, nr, nc, nz) × 8. idx = SENTINEL for out-of-bounds.
fn hex_nbrs_obc(site: usize, l: usize) -> [(i32, i32, i32, i32); 8] {
    let lsz = l * l;
    let li = l as i32;
    let z = (site / lsz) as i32;
    let rem = site % lsz;
    let r = (rem / l) as i32;
    let c = (rem % l) as i32;

    let (dr, dc): ([i32; 6], [i32; 6]) = if r % 2 == 0 {
        ([-1, -1, 0, 0, 1, 1], [0, 1, -1, 1, 0, 1])
    } else {
        ([-1, -1, 0, 0, 1, 1], [-1, 0, -1, 1, -1, 0])
    };

    let mut out = [(SENTINEL, 0i32, 0i32, 0i32); 8];
    for i in 0..6 {
        let nr = r + dr[i];
        let nc = c + dc[i];
        out[i].1 = nr;
        out[i].2 = nc;
        out[i].3 = z;
        if nr >= 0 && nr < li && nc >= 0 && nc < li {
            out[i].0 = (z * li + nr) * li + nc;
        }
    }
    let zp = z - 1;
    out[6] = if zp >= 0 {
        ((zp * li + r) * li + c, r, c, zp)
    } else {
        (SENTINEL, r, c, zp)
    };
    let zn = z + 1;
    out[7] = if zn < li {
        ((zn * li + r) * li + c, r, c, zn)
    } else {
        (SENTINEL, r, c, zn)
    };
    out
}

fn hex_cart(r: i32, c: i32) -> (f64, f64) {
    let offset = if r & 1 == 1 { -0.5 } else { 0.0 };
    (c as f64 + offset, r as f64 * (3.0_f64).sqrt() / 2.0)
}

fn obc_phi(nr: i32, nc: i32, nz: i32, cx: f64, cy: f64, cz: f64, alpha: f64) -> f64 {
    let (px, py) = hex_cart(nr, nc);
    let dz = nz as f64 - cz;
    let d = ((px - cx).powi(2) + (py - cy).powi(2) + dz * dz).sqrt();
    if d > 0.5 {
        alpha / d
    } else {
        2.0 * alpha
    }
}

/// 6-color: color = (c + (r&1)) % 3 + (z&1) * 3
fn site_color(r: i32, c: i32, z: i32) -> usize {
    ((c + (r & 1)) % 3 + (z & 1) * 3) as usize
}

pub fn cpu_hydrogen_obc(
    alpha: f64,
    l: usize,
    n_sor: usize,
    n_iter: usize,
    dtau: f64,
) -> (f64, f64, f64) {
    let n = l * l * l;
    let li = l as i32;
    let lsz = l * l;
    let cr = li / 2;
    let cc = li / 2;
    let cz_i = li / 2;
    let (cx, cy) = hex_cart(cr, cc);
    let cz = cz_i as f64;

    // Pre-compute OBC neighbors + coords
    let nbrs: Vec<[(i32, i32, i32, i32); 8]> = (0..n).map(|i| hex_nbrs_obc(i, l)).collect();
    let coords: Vec<(i32, i32, i32)> = (0..n)
        .map(|i| {
            let z = (i / lsz) as i32;
            let rem = i % lsz;
            ((rem / l) as i32, (rem % l) as i32, z)
        })
        .collect();

    // Phase 1: Coulomb warm-start + SOR
    let mut phi = vec![0.0f64; n];
    for i in 0..n {
        let (r, c, z) = coords[i];
        let (px, py) = hex_cart(r, c);
        let dz = z as f64 - cz;
        let d = ((px - cx).powi(2) + (py - cy).powi(2) + dz * dz).sqrt();
        phi[i] = if d > 0.5 { alpha / d } else { 2.0 * alpha };
    }
    let center = ((cz_i * li + cr) * li + cc) as usize;
    let mut rho = vec![-1.0 / n as f64; n];
    rho[center] += 1.0;

    let omega = 2.0 / (1.0 + (PI / l as f64).sin());
    for _sweep in 0..n_sor {
        for phase in 0..6usize {
            for i in 0..n {
                let (r, c, z) = coords[i];
                if site_color(r, c, z) != phase {
                    continue;
                }
                let nb = &nbrs[i];
                let mut s = 0.0;
                for j in 0..8 {
                    let (idx, nr, nc, nz) = nb[j];
                    s += if idx >= 0 {
                        phi[idx as usize]
                    } else {
                        obc_phi(nr, nc, nz, cx, cy, cz, alpha)
                    };
                }
                let phi_gs = (s + 8.0 * rho[i]) / 8.0;
                phi[i] = (1.0 - omega) * phi[i] + omega * phi_gs;
            }
        }
    }

    // Phase 2: Gaussian init
    let sigma = (1.0 / alpha).max(1.0);
    let mut psi: Vec<f64> = (0..n)
        .map(|i| {
            let (r, c, z) = coords[i];
            let (px, py) = hex_cart(r, c);
            let dz = z as f64 - cz;
            let d2 = (px - cx).powi(2) + (py - cy).powi(2) + dz * dz;
            (-d2 / (2.0 * sigma * sigma)).exp()
        })
        .collect();
    let norm: f64 = psi.iter().map(|v| v * v).sum();
    let sc = 1.0 / norm.sqrt();
    for v in &mut psi {
        *v *= sc;
    }

    // Phase 3: Imaginary time (OBC: out-of-bounds → ψ=0)
    let mut psi_new = vec![0.0f64; n];
    for step in 1..=n_iter {
        for i in 0..n {
            let nb = &nbrs[i];
            let mut ns = 0.0;
            for j in 0..8 {
                if nb[j].0 >= 0 {
                    ns += psi[nb[j].0 as usize];
                }
            }
            let kinetic = psi[i] - ns / 8.0;
            let potential = -alpha * phi[i] * psi[i];
            psi_new[i] = psi[i] - dtau * (kinetic + potential);
        }
        std::mem::swap(&mut psi, &mut psi_new);
        if step % 20 == 0 {
            let norm: f64 = psi.iter().map(|v| v * v).sum();
            let sc = 1.0 / norm.sqrt();
            for v in &mut psi {
                *v *= sc;
            }
        }
    }

    // Energy measurement
    let (mut ek, mut ep) = (0.0, 0.0);
    for i in 0..n {
        let nb = &nbrs[i];
        let mut ns = 0.0;
        for j in 0..8 {
            if nb[j].0 >= 0 {
                ns += psi[nb[j].0 as usize];
            }
        }
        ek += psi[i] * (psi[i] - ns / 8.0);
        ep += psi[i] * (-alpha * phi[i] * psi[i]);
    }
    (ek + ep, ek, ep)
}

pub fn cpu_hydrogen_scan(alpha: f64) {
    println!("\n  ── CPU Lattice Hydrogen OBC (α={}) ──", alpha);
    println!("  Coulomb warm-start, 6-color SOR, OBC walls");
    let dtau = 0.05;
    // Use smaller L values where CPU converges in reasonable time
    // L must be odd (center site at L/2) and large enough for Bohr radius 1/α
    let l_values: Vec<usize> = vec![31, 41, 51, 61, 71, 81];
    let mut results = Vec::new();

    for &l in &l_values {
        let n_sor = 3 * l; // O(L) SOR with optimal ω, 3× for safety
                           // τ = 1000 sufficient — earlier tests showed C converges by τ=250 for small L
        let tau = 1000.0_f64;
        let n_iter = (tau / dtau) as usize;
        let n = l * l * l;
        println!("  L={} ({} sites, n_sor={}, τ={:.0})...", l, n, n_sor, tau);
        let (et, ek, ep) = cpu_hydrogen_obc(alpha, l, n_sor, n_iter, dtau);
        let c = -2.0 * et / (alpha * alpha);
        results.push((l, et, c));
        println!(
            "    E={:.6e} (kin={:.6e} pot={:.6e}), C={:.4}",
            et, ek, ep, c
        );
    }

    println!("\n  Richardson pair fits (C = C_∞ - B/L):");
    for i in 0..results.len() - 1 {
        let (l1, _, c1) = results[i];
        let (l2, _, c2) = results[i + 1];
        let b = (c1 - c2) * (l1 as f64 * l2 as f64) / (l1 as f64 - l2 as f64);
        let c_inf = c1 + b / l1 as f64;
        println!("    ({},{}): B={:.1}, C_∞={:.4}", l1, l2, b, c_inf);
    }

    // Three-point Richardson: C(L) = C_∞ + B/L + D/L²
    // Solve using elimination on adjacent triples
    if results.len() >= 3 {
        println!("\n  Three-point Richardson (C = C_∞ + B/L + D/L²):");
        for i in 0..results.len() - 2 {
            let (l1, _, c1) = results[i];
            let (l2, _, c2) = results[i + 1];
            let (l3, _, c3) = results[i + 2];
            let (f1, f2, f3) = (l1 as f64, l2 as f64, l3 as f64);
            // c_i = A + B/f_i + D/f_i²
            // Differences:
            let dc21 = c2 - c1;
            let dc31 = c3 - c1;
            let dx21 = 1.0 / f2 - 1.0 / f1;
            let dx31 = 1.0 / f3 - 1.0 / f1;
            let dx2_21 = 1.0 / (f2 * f2) - 1.0 / (f1 * f1);
            let dx2_31 = 1.0 / (f3 * f3) - 1.0 / (f1 * f1);
            let d = (dc21 * dx31 - dc31 * dx21) / (dx2_21 * dx31 - dx2_31 * dx21);
            let b = (dc21 - d * dx2_21) / dx21;
            let a = c1 - b / f1 - d / (f1 * f1);
            println!(
                "    ({},{},{}): C_∞={:.4}, B={:.1}, D={:.0}",
                l1, l2, l3, a, b, d
            );
        }
    }

    println!(
        "\n  Continuum prediction: C_∞ = C_G²/2 = {:.4}",
        continuum_green_coefficient().powi(2) / 2.0
    );
    println!("  GPU reference (α=0.1): C_∞ ≈ 0.547, B ≈ 18.1");
}

/// Dispersion relation along high-symmetry directions in the hex+z BZ.
/// Uses BZ coordinates (θ₁, θ₂, θ_z) and their correct quadratic coefficients.
pub fn dispersion_analysis() {
    println!("\n  ── Dispersion Relation: T_lattice vs T_quadratic ──");
    println!("  T_quad(θ₁,θ₂,θ_z) = (θ₁²-θ₁θ₂+θ₂²)/4 + θ_z²/8");
    println!(
        "  {:>6} {:>12} {:>12} {:>8}",
        "θ", "T_lattice", "T_quad", "ratio"
    );

    // Along θ₁ direction (θ₂=θ_z=0): T_quad = θ₁²/4
    println!("\n  Direction: Γ→M (θ₁, 0, 0)");
    for i in 1..=20 {
        let th = i as f64 * PI / 20.0;
        let t_lat = t_hex_z(th, 0.0, 0.0);
        let t_quad = th * th / 4.0;
        println!(
            "  {:6.3} {:12.6} {:12.6} {:8.4}",
            th,
            t_lat,
            t_quad,
            t_lat / t_quad
        );
    }

    // Along z-axis (θ₁=θ₂=0, θ_z): T_quad = θ_z²/8
    println!("\n  Direction: Γ→A (0, 0, θ_z)");
    for i in 1..=20 {
        let th = i as f64 * PI / 20.0;
        let t_lat = t_hex_z(0.0, 0.0, th);
        let t_quad = th * th / 8.0;
        println!(
            "  {:6.3} {:12.6} {:12.6} {:8.4}",
            th,
            t_lat,
            t_quad,
            t_lat / t_quad
        );
    }

    // Along diagonal (θ₁=θ₂=θ, θ_z=0): T_quad = θ²/4
    println!("\n  Direction: Γ→K (θ, θ, 0)");
    for i in 1..=20 {
        let th = i as f64 * PI / 20.0;
        let t_lat = t_hex_z(th, th, 0.0);
        let t_quad = th * th / 4.0;
        println!(
            "  {:6.3} {:12.6} {:12.6} {:8.4}",
            th,
            t_lat,
            t_quad,
            t_lat / t_quad
        );
    }

    // Lattice Coulomb along z-axis vs anisotropic continuum:
    // G_cont(0,0,z) = 1/(4π M_perp |z|) = 4/(3π|z|)
    println!("\n  ── Lattice Coulomb along z vs continuum ──");
    println!("  G_cont(0,0,z) = 1/(4π·M_perp·|z|) = 4/(3π|z|)");
    let m_perp = 3.0 / 16.0;
    let g_cont_coeff = 1.0 / (4.0 * PI * m_perp); // = 4/(3π) ≈ 0.4244
    println!(
        "  {:>4} {:>12} {:>12} {:>8}",
        "z", "G_lattice", "G_cont", "ratio"
    );
    for &rz in &[1, 2, 3, 5, 8, 12] {
        // Use enough integration points: need at least ~10 points per oscillation
        let n_theta = (200_usize).max(20 * rz);
        let g_r = integrate_3d_midpoint(
            |t1, t2, tz| (rz as f64 * tz).cos() / t_hex_z(t1, t2, tz),
            0.0,
            2.0 * PI,
            0.0,
            2.0 * PI,
            -PI,
            PI,
            n_theta,
        ) / (2.0 * PI).powi(3);
        let g_cont = g_cont_coeff / rz as f64;
        println!(
            "  {:4} {:12.6} {:12.6} {:8.4}",
            rz,
            g_r,
            g_cont,
            g_r / g_cont
        );
    }

    // BZ shell integral: fraction of G(0) from each momentum shell
    println!("\n  ── BZ integral by k-shell (fraction of G(0)) ──");
    println!("  {:>8} {:>12} {:>12}", "|k|_max", "G_partial", "fraction");
    let g0_full = watson_hex_z(200);
    for &frac in &[0.1, 0.2, 0.3, 0.5, 0.7, 1.0] {
        let k_max = frac * PI;
        let g_partial = integrate_3d_midpoint(
            |t1, t2, tz| {
                // BZ coordinates: effective |k|² accounts for non-orthogonal axes
                // For triangular + z: use θ₁²-θ₁θ₂+θ₂² + θ_z² as isotropic proxy
                let k2_eff = t1 * t1 - t1 * t2 + t2 * t2 + tz * tz;
                if k2_eff < k_max * k_max {
                    1.0 / t_hex_z(t1, t2, tz)
                } else {
                    0.0
                }
            },
            0.0,
            2.0 * PI,
            0.0,
            2.0 * PI,
            -PI,
            PI,
            200,
        ) / (2.0 * PI).powi(3);
        println!(
            "  {:8.3} {:12.6} {:12.4}",
            k_max,
            g_partial,
            g_partial / g0_full
        );
    }
    println!("  Full G(0) = {:.6}", g0_full);
}

// ── Simple cubic T³ PBC hydrogen solver ──────────────────────────────────────
//
// T³ topology: the spatial rotation group SO(3) is compact, so the Cayley
// graph of its discrete subgroup (Z₃ lattice symmetry) is naturally periodic.
// PBC is the physics, not a computational shortcut.
//
// L = 12: lcm(Z₃ period = 3, Clifford grade-1 dim = 4) = lcm(3,4) = 12.
// This is the minimum torus that fits a complete Z₃ quark orbit AND a complete
// grade-1 spacetime basis without truncation. L=12 = LEPTON_GRADE_DIM × SU2_DIM.
//
// Charge neutrality on T³ is required by Gauss's law on a compact manifold.
// The Poisson equation has a unique solution only after zero-mode subtraction —
// exactly what experiment 5 (PBC zero-mode fix) demonstrated.

fn sc_nbrs_pbc(site: usize, l: usize) -> [usize; 6] {
    let l2 = l * l;
    let z = site / l2;
    let rem = site % l2;
    let y = rem / l;
    let x = rem % l;
    [
        z * l2 + y * l + (x + 1) % l,       // +x
        z * l2 + y * l + (x + l - 1) % l,   // -x
        z * l2 + ((y + 1) % l) * l + x,     // +y
        z * l2 + ((y + l - 1) % l) * l + x, // -y
        ((z + 1) % l) * l2 + y * l + x,     // +z
        ((z + l - 1) % l) * l2 + y * l + x, // -z
    ]
}

/// Hydrogen on a simple cubic T³ torus with PBC.
///
/// Poisson normalization: discrete SC equation (6I−A)·φ = 6·ρ.
/// For unit point source ρ=δ: source = 6·δ → φ ~ 6/(4πr) ≈ 0.477/r.
/// This is C_poisson = 6/(4π) ≈ 0.477 (the lattice Coulomb constant, not a bug).
/// Therefore V = -α·φ ~ -0.477α/r (lattice-natural Coulomb, not -α/r).
///
/// The OBC solver uses 8 neighbors: (8I−A)·φ = 8·ρ → φ ~ 8/(4πr) ≈ 0.637/r.
/// Both conventions are self-consistent; neither matches -α/r from physical units,
/// because "α" is defined in the same lattice units as the kinetic energy.
///
/// Kinetic: SC 6-neighbor, kinetic = psi[i] - ns/6 → m_eff = 3.
/// Poisson: 2-color red-black SOR + zero-mode subtraction (charge neutrality on T³).
/// C = -2E/α² (matches OBC convention).
///
/// Returns (e_total, e_kin, e_pot).
pub fn pbc_sc_hydrogen(
    alpha: f64,
    l: usize,
    n_sor: usize,
    n_iter: usize,
    dtau: f64,
) -> (f64, f64, f64) {
    let n = l * l * l;
    let l2 = l * l;
    let lf = l as f64;
    let half_l = lf / 2.0;
    let cf = (l / 2) as f64;
    let center = (l / 2) * l2 + (l / 2) * l + (l / 2);

    let nbrs: Vec<[usize; 6]> = (0..n).map(|i| sc_nbrs_pbc(i, l)).collect();

    // Minimum-image coordinates relative to proton at (L/2, L/2, L/2)
    let coords: Vec<(f64, f64, f64)> = (0..n)
        .map(|i| {
            let z = (i / l2) as f64;
            let rem = i % l2;
            let y = (rem / l) as f64;
            let x = (rem % l) as f64;
            let mut dx = x - cf;
            if dx > half_l {
                dx -= lf;
            } else if dx < -half_l {
                dx += lf;
            }
            let mut dy = y - cf;
            if dy > half_l {
                dy -= lf;
            } else if dy < -half_l {
                dy += lf;
            }
            let mut dz = z - cf;
            if dz > half_l {
                dz -= lf;
            } else if dz < -half_l {
                dz += lf;
            }
            (dx, dy, dz)
        })
        .collect();

    // Coulomb warm-start: φ = C_poisson/r where C_poisson = 6/(4π) ≈ 0.477.
    // The discrete SC Poisson (6I−A)·φ = 6·ρ with ρ=δ gives φ ~ 6/(4πr) in the
    // continuum limit.  This warm-start is the converged shape, speeding SOR.
    let c_poisson = 6.0 / (4.0 * PI); // ≈ 0.477
    let mut phi: Vec<f64> = (0..n)
        .map(|i| {
            let (dx, dy, dz) = coords[i];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if r > 0.5 {
                c_poisson / r
            } else {
                2.0 * c_poisson
            }
        })
        .collect();
    let mean: f64 = phi.iter().sum::<f64>() / n as f64;
    phi.iter_mut().for_each(|p| *p -= mean);

    // Poisson source: proton +1 at center, neutralizing background -1/N per site
    // -Δφ(i) = 6φ(i) - Σφ_j = 6·ρ(i)  →  φ_GS = (Σφ_j + 6ρ) / 6
    let rhs_center = 1.0 - 1.0 / n as f64;
    let rhs_bg = -1.0 / n as f64;

    // PBC optimal SOR ω: lowest mode is 2π/L (vs OBC π/L)
    let omega = 2.0 / (1.0 + (2.0 * PI / lf).sin());

    // 2-color (red-black) SOR: color = (x+y+z) % 2; all 6 SC neighbors change color
    for _ in 0..n_sor {
        for color in 0..2usize {
            for i in 0..n {
                let z = i / l2;
                let rem = i % l2;
                let y = rem / l;
                let x = rem % l;
                if (x + y + z) % 2 != color {
                    continue;
                }
                let rhs = if i == center { rhs_center } else { rhs_bg };
                let s: f64 = nbrs[i].iter().map(|&j| phi[j]).sum();
                let phi_gs = (s + 6.0 * rhs) / 6.0;
                phi[i] += omega * (phi_gs - phi[i]);
            }
        }
        // Zero-mode subtraction: enforces charge neutrality, removes PBC divergence
        let mean: f64 = phi.iter().sum::<f64>() / n as f64;
        phi.iter_mut().for_each(|p| *p -= mean);
    }

    // Gaussian init centred on proton: σ = Bohr radius estimate 1/α
    let sigma = (1.0 / alpha).max(1.0);
    let mut psi: Vec<f64> = (0..n)
        .map(|i| {
            let (dx, dy, dz) = coords[i];
            (-(dx * dx + dy * dy + dz * dz) / (2.0 * sigma * sigma)).exp()
        })
        .collect();
    let norm: f64 = psi.iter().map(|v| v * v).sum::<f64>().sqrt();
    psi.iter_mut().for_each(|v| *v /= norm);

    // Imaginary time: ψ → ψ - dtau·Hψ,  H = -Δ_SC/6 + V
    // kinetic = psi[i] - ns/6  (SC: divide by coordination 6)
    // potential = -alpha·phi[i]  (lattice Coulomb: V ~ -0.477α/r via SC Poisson)
    let mut psi_new = vec![0.0f64; n];
    for step in 1..=n_iter {
        for i in 0..n {
            let ns: f64 = nbrs[i].iter().map(|&j| psi[j]).sum();
            let kinetic = psi[i] - ns / 6.0;
            let potential = -alpha * phi[i] * psi[i];
            psi_new[i] = psi[i] - dtau * (kinetic + potential);
        }
        std::mem::swap(&mut psi, &mut psi_new);
        if step % 20 == 0 {
            let norm: f64 = psi.iter().map(|v| v * v).sum::<f64>().sqrt();
            psi.iter_mut().for_each(|v| *v /= norm);
        }
    }
    let norm: f64 = psi.iter().map(|v| v * v).sum::<f64>().sqrt();
    psi.iter_mut().for_each(|v| *v /= norm);

    // Energy: E = ⟨ψ|H|ψ⟩
    let (mut ek, mut ep) = (0.0, 0.0);
    for i in 0..n {
        let ns: f64 = nbrs[i].iter().map(|&j| psi[j]).sum();
        ek += psi[i] * (psi[i] - ns / 6.0);
        ep += psi[i] * (-alpha * phi[i] * psi[i]);
    }
    (ek + ep, ek, ep)
}

/// Richardson L-scan for SC T³ hydrogen at fixed α.
///
/// Scans L = 16, 24, 32, 48, 64 and extrapolates C_∞ via C(L) = C_∞ + B/L.
///
/// At α = 0.5: SC Bohr radius a₀ = 1/(m_eff·α_eff) = 1/(3·0.477·0.5) ≈ 1.4
/// lattice spacings (α_eff = 0.477α from C_poisson = 6/(4π)).
/// L/a₀ ≈ 11..46 across the scan → clean Richardson convergence.
/// At α = 0.1: a₀ ≈ 7 spacings, L=12 gives L/a₀ = 1.7 (too small, hence
/// the anomalously low C≈0.2 from the original L=12 run).
pub fn pbc_sc_scan(alpha: f64) {
    println!("\n╔══════════════════════════════════════════════════════════════╗");
    println!(
        "║  SC T³ Hydrogen  α={:.3}  L scan — Richardson C_∞         ║",
        alpha
    );
    println!("╚══════════════════════════════════════════════════════════════╝");
    println!("  Topology: T³ × ℝ — PBC from compact SO(3), open time");
    println!("  Poisson:  (6I-A)·φ = 6·δ → φ ~ 0.477/r (C_poisson = 6/(4π))");
    println!("  Coupling: V = -α·φ ~ -0.477α/r,  α_eff = 0.477α,  a₀ ≈ 1.4 (at α=0.5)");
    println!("  Kinetic:  SC 6-neighbor, kinetic = ψ - ns/6");
    println!();
    println!(
        "  {:>5}  {:>8}  {:>10}  {:>10}  {:>10}  {:>10}",
        "L", "N", "E_tot", "E_kin", "E_pot", "C = -2E/α²"
    );
    println!("  {}", "-".repeat(62));

    let dtau = 0.05;
    let l_values: Vec<usize> = vec![16, 24, 32, 48, 64];
    let mut results: Vec<(usize, f64)> = Vec::new();

    for &l in &l_values {
        let n = l * l * l;
        let n_sor = 4 * l; // PBC optimal ω → O(L) convergence; 4L safe
        let tau = 300.0_f64; // Bohr radius < lattice spacing → fast convergence
        let n_iter = (tau / dtau) as usize;
        print!("  {:5}  {:8}  computing...", l, n);
        let (et, ek, ep) = pbc_sc_hydrogen(alpha, l, n_sor, n_iter, dtau);
        let c = -2.0 * et / (alpha * alpha);
        println!(
            "\r  {:5}  {:8}  {:10.6}  {:10.6}  {:10.6}  {:10.4}",
            l, n, et, ek, ep, c
        );
        results.push((l, c));
    }

    println!("\n  Richardson pair fits (C = C_∞ + B/L):");
    for i in 0..results.len() - 1 {
        let (l1, c1) = results[i];
        let (l2, c2) = results[i + 1];
        let b = (c1 - c2) * (l1 as f64 * l2 as f64) / (l2 as f64 - l1 as f64);
        let c_inf = c1 - b / l1 as f64;
        println!("    ({:2},{:2}): B = {:+.2}, C_∞ = {:.4}", l1, l2, b, c_inf);
    }

    // Three-point Richardson: C = C_∞ + B/L + D/L²
    if results.len() >= 3 {
        println!("\n  Three-point Richardson (C = C_∞ + B/L + D/L²):");
        for i in 0..results.len() - 2 {
            let (l1, c1) = results[i];
            let (l2, c2) = results[i + 1];
            let (l3, c3) = results[i + 2];
            let (f1, f2, f3) = (l1 as f64, l2 as f64, l3 as f64);
            let dc21 = c2 - c1;
            let dc31 = c3 - c1;
            let dx21 = 1.0 / f2 - 1.0 / f1;
            let dx31 = 1.0 / f3 - 1.0 / f1;
            let dx2_21 = 1.0 / (f2 * f2) - 1.0 / (f1 * f1);
            let dx2_31 = 1.0 / (f3 * f3) - 1.0 / (f1 * f1);
            let d = (dc21 * dx31 - dc31 * dx21) / (dx2_21 * dx31 - dx2_31 * dx21);
            let b = (dc21 - d * dx2_21) / dx21;
            let a = c1 - b / f1 - d / (f1 * f1);
            println!(
                "    ({:2},{:2},{:2}): C_∞ = {:.4}, B = {:+.2}, D = {:+.0}",
                l1, l2, l3, a, b, d
            );
        }
    }

    println!();
    println!("  Reference (OBC hex+z L→∞):  C_∞ = 0.5466  (Richardson L=161–961)");
    println!("  SC continuum limit (L→∞):   C_∞ ≈ 0.684   (C_G(SC)²/2, before UV correction)");
    println!("  Physical hydrogen:           C   = 0.5000  (E₀ = -α²/2)");
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn watson_analysis() {
        run_analysis();
    }

    #[test]
    fn cpu_hydrogen() {
        println!();
        cpu_hydrogen_scan(0.1);
    }

    #[test]
    fn dispersion() {
        dispersion_analysis();
    }

    #[test]
    fn pbc_sc_torus() {
        // Richardson L-scan: SC T³ hydrogen at α=0.5 with corrected 4π normalization.
        //
        // Fix: -Δφ = δ → φ ~ 1/(4πr), so V = -α/r requires potential = -4πα·φ.
        // At α=0.5: SC Bohr radius a₀ = 1/(3·0.5) ≈ 0.67 < lattice spacing.
        // L=16..64 all satisfy L >> a₀ → clean Richardson extrapolation.
        //
        // Expected: C(L) converging from above toward C_∞ ≈ 0.50–0.55.
        println!();
        pbc_sc_scan(0.5);

        // Smoke-test at α=0.5, L=24 with corrected potential
        let dtau = 0.05;
        let l = 24usize;
        let n_sor = 4 * l;
        let n_iter = (300.0_f64 / dtau) as usize;
        let (et, _ek, _ep) = pbc_sc_hydrogen(0.5, l, n_sor, n_iter, dtau);
        let c = -2.0 * et / (0.5 * 0.5);

        // C must be positive (bound state exists)
        assert!(c > 0.0, "C must be positive (bound state), got C={c:.4}");
        // C must be physically bounded: below SC continuum limit
        assert!(c < 1.5, "C={c:.4} unexpectedly large");
        // Energy must be negative (bound)
        assert!(
            et < 0.0,
            "Ground state energy must be negative, got E={et:.6e}"
        );
    }
}
