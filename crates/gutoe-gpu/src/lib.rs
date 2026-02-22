// GUTOE GPU — Accelerated 3D Schrödinger solver
// Copyright (C) 2026 Riff Labs, AGPL-3.0-or-later
//
// Backends:
//   CPU:  always available, wraps gutoe-em::quantum_lepton
//   CUDA: --features cuda  (NVIDIA sm_86 = RTX 3070 Ti default)
//   ROCm: --features rocm  (AMD via HIP — same kernel source)
//
// To compile:
//   NVIDIA: CUDA_ARCH=sm_86 cargo build -p gutoe-gpu --features cuda
//   AMD:    cargo build -p gutoe-gpu --features rocm
//
// The 3070 Ti benchmark at α=0.1, L=144:
//   2.985M sites × 5000 Jacobi + 20000 imaginary-time steps
//   Memory: ~200 MB (fits in 12 GB)
//   Expected runtime: ~15s (FP64 memory-bandwidth limited)

pub mod watson;
pub mod speculative;

use gutoe_em::quantum_lepton::{bohr_test_3d, BohrResult};

// ── Backend enum ──────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Backend { Cpu, #[cfg(feature = "cuda")] Cuda, #[cfg(feature = "rocm")] Rocm }

pub fn detect_backend() -> Backend {
    #[cfg(feature = "rocm")] { return Backend::Rocm; }
    #[cfg(feature = "cuda")] { return Backend::Cuda; }
    #[allow(unreachable_code)]
    Backend::Cpu
}

// ── Solver config ─────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct SolverConfig {
    pub alpha:     f64,
    pub l:         usize,   // L×L×L cube
    pub n_jacobi:  usize,
    pub n_iter:    usize,
    pub dtau:      f64,
    pub backend:   Backend,
}

impl SolverConfig {
    /// RTX 3070 Ti: α=0.1, L=144 — Bohr radius (10) fits, ~200 MB VRAM
    pub fn bohr_3070ti() -> Self {
        Self { alpha: 0.1, l: 144, n_jacobi: 5000, n_iter: 20_000, dtau: 0.001,
               backend: detect_backend() }
    }
    /// Fast CPU test: α=0.5, L=12
    pub fn fast_test(backend: Backend) -> Self {
        Self { alpha: 0.5, l: 12, n_jacobi: 200, n_iter: 5_000, dtau: 0.01, backend }
    }
    /// Multi-point scan for Bohr convergence
    pub fn scan_point(alpha: f64, backend: Backend) -> Self {
        let l = ((6.0 / alpha) as usize).max(12).min(300);
        let n_iter = ((3.0 / (alpha * alpha)) as usize).min(50_000).max(2_000);
        let dtau = 0.003 * alpha.min(1.0);
        let n_jacobi = 1000;
        Self { alpha, l, n_jacobi, n_iter, dtau, backend }
    }
}

// ── Main solver ───────────────────────────────────────────────────────────────

pub fn solve_hydrogen_3d(cfg: &SolverConfig) -> BohrResult {
    match cfg.backend {
        Backend::Cpu                     => solve_cpu(cfg),
        #[cfg(feature = "cuda")] Backend::Cuda => solve_gpu(cfg),
        #[cfg(feature = "rocm")] Backend::Rocm => solve_gpu(cfg),
    }
}

fn solve_cpu(cfg: &SolverConfig) -> BohrResult {
    bohr_test_3d(cfg.alpha, cfg.l, cfg.l, cfg.n_jacobi, cfg.n_iter, cfg.dtau)
}

// ── GPU FFI (shared between CUDA and ROCm — same symbol, same ABI) ───────────

#[cfg(any(feature = "cuda", feature = "rocm"))]
extern "C" {
    fn run_hydrogen_cuda(
        alpha_em:    f64,
        hex_rows:    i32, hex_cols: i32, layers: i32,
        n_jacobi:    i32, n_iter:   i32, renorm_every: i32,
        dtau:        f64,
        out_e_total: *mut f64,
        out_e_kin:   *mut f64,
        out_e_pot:   *mut f64,
    );
    /// Open-boundary variant: small box, φ=α/r and ψ=0 at walls.
    /// Equivalent to L=∞ periodic for localised states. Memory ∝ (2R+1)³
    /// where R ≈ 8/α — independent of the physical box size.
    fn run_hydrogen_obc(
        alpha_em:    f64,
        hex_rows:    i32, hex_cols: i32, layers: i32,
        n_jacobi:    i32, n_iter:   i32, renorm_every: i32,
        dtau:        f64,
        out_e_total: *mut f64,
        out_e_kin:   *mut f64,
        out_e_pot:   *mut f64,
    );
    /// Periodic-boundary variant: SOR Poisson with PBC (cold start),
    /// periodic Hamiltonian. Includes Madelung-like image corrections.
    fn run_hydrogen_pbc(
        alpha_em:    f64,
        hex_rows:    i32, hex_cols: i32, layers: i32,
        n_jacobi:    i32, n_iter:   i32, renorm_every: i32,
        dtau:        f64,
        out_e_total: *mut f64,
        out_e_kin:   *mut f64,
        out_e_pot:   *mut f64,
    );
}

#[cfg(any(feature = "cuda", feature = "rocm"))]
fn solve_gpu(cfg: &SolverConfig) -> BohrResult {
    let (mut e_total, mut e_kin, mut e_pot) = (0.0_f64, 0.0_f64, 0.0_f64);
    unsafe {
        run_hydrogen_cuda(
            cfg.alpha,
            cfg.l as i32, cfg.l as i32, cfg.l as i32,
            cfg.n_jacobi as i32, cfg.n_iter as i32, 10_i32,
            cfg.dtau,
            &mut e_total, &mut e_kin, &mut e_pot,
        );
    }
    let bohr_3d = -cfg.alpha * cfg.alpha / 2.0;
    BohrResult { alpha: cfg.alpha, l: cfg.l, e_total, e_kin, e_pot,
                 bohr_3d, ratio: e_total / bohr_3d }
}

/// Open-boundary Schrödinger solver — CoW memory: allocates only (2R+1)³ sites.
/// For α=0.1: R≈82, box=165³≈4.5M sites regardless of physical L.
/// For α=0.01: R≈802, box=1605³ — scale R with 1/α.
#[cfg(any(feature = "cuda", feature = "rocm"))]
pub fn solve_hydrogen_obc(alpha: f64, n_jacobi: usize, n_iter: usize, dtau: f64,
                           _backend: Backend) -> BohrResult {
    // Active radius: 8 Bohr radii (exp(-8)≈3e-4, energy error negligible)
    let r_active = ((8.0 / alpha).ceil() as usize).max(8);
    let l = 2 * r_active + 1;   // small cube, open BC
    let (mut e_total, mut e_kin, mut e_pot) = (0.0_f64, 0.0_f64, 0.0_f64);
    unsafe {
        run_hydrogen_obc(
            alpha,
            l as i32, l as i32, l as i32,
            n_jacobi as i32, n_iter as i32, 10_i32,
            dtau,
            &mut e_total, &mut e_kin, &mut e_pot,
        );
    }
    let bohr_3d = -alpha * alpha / 2.0;
    BohrResult { alpha, l, e_total, e_kin, e_pot, bohr_3d, ratio: e_total / bohr_3d }
}

// ── Bohr convergence scan ─────────────────────────────────────────────────────

/// Scan α values to show exp_3d converging to 2.0 (Bohr formula).
/// On GPU: each point takes ~15s (3070 Ti) or ~0.05s (395+), so a
/// 10-point scan is feasible as a single benchmark run.
pub fn bohr_convergence_scan(alphas: &[f64], backend: Backend) {
    println!("GUTOE Bohr convergence: E₀ ∝ α^n, n → 2.0 as L → ∞");
    println!("{:>8}  {:>6}  {:>10}  {:>10}  {:>8}", "α", "L", "E_total", "−α²/2", "ratio");
    println!("{:>8}  {:>6}  {:>10}  {:>10}  {:>8}", "─", "─", "─", "─", "─");

    let mut prev_e: Option<(f64, f64)> = None;
    for &alpha in alphas {
        let cfg = SolverConfig::scan_point(alpha, backend);
        let r = solve_hydrogen_3d(&cfg);
        println!("{:>8.4}  {:>6}  {:>10.6}  {:>10.6}  {:>8.3}",
            alpha, cfg.l, r.e_total, r.bohr_3d, r.ratio);

        if let Some((a0, e0)) = prev_e {
            let exp = (e0.abs() / r.e_total.abs()).ln() / (a0 / alpha).ln();
            println!("         exp[{a0:.3}→{alpha:.3}] = {exp:.3}  (Bohr target: 2.0)");
        }
        prev_e = Some((alpha, r.e_total));
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cpu_backend_gives_bound_state() {
        let cfg = SolverConfig::fast_test(Backend::Cpu);
        let r = solve_hydrogen_3d(&cfg);
        assert!(r.e_total < 0.0,
            "CPU: α={:.2} L={} E={:.6} must be bound", cfg.alpha, cfg.l, r.e_total);
        println!("  CPU: α={:.2} L={} E={:.6} C={:.3}", cfg.alpha, cfg.l, r.e_total, r.ratio);
    }

    #[test]
    fn backend_detection() {
        let b = detect_backend();
        println!("  Active backend: {b:?}");
        let _ = b;
    }

    /// Run the GPU Bohr test if CUDA/ROCm feature is enabled.
    /// cargo test -p gutoe-gpu --features cuda -- --nocapture
    #[test]
    #[cfg(any(feature = "cuda", feature = "rocm"))]
    fn gpu_bohr_test() {
        let cfg = SolverConfig::bohr_3070ti();
        println!("  Running on {:?}: α={} L={}³ N={}", cfg.backend, cfg.alpha, cfg.l, cfg.l.pow(3));
        let r = solve_hydrogen_3d(&cfg);
        println!("  E_kin={:.6} E_pot={:.6} E_total={:.6} C={:.3}",
            r.e_kin, r.e_pot, r.e_total, r.ratio);
        println!("  Bohr −α²/2 = {:.6}", r.bohr_3d);
        assert!(r.e_total < 0.0, "GPU: must be bound at α=0.1 L=144");
    }

    /// Convergence scan: C = E₀/(−α²/2) vs lattice size.
    /// Run with: cargo test -p gutoe-gpu --features rocm --release -- bohr_convergence_scan --nocapture
    #[test]
    #[cfg(any(feature = "cuda", feature = "rocm"))]
    fn bohr_convergence_scan() {
        let backend = detect_backend();
        // Scan increasing L at fixed α=0.1 (a₀=10 lattice units).
        // Keep iterations fixed to bound runtime; convergence is geometric
        // so we can extrapolate C_∞ from the falling series.
        // L=144 → L/a₀=14.4, L=192 → 19.2, L=240 → 24.0, L=288 → 28.8
        // Scale iterations down for large L to stay within ROCm iGPU limits.
        // C(L) trend is what matters; we don't need full ground-state convergence.
        let configs = [
            (0.1_f64, 336_usize, 5_000, 20_000, 0.001),  // sanity: should work
            (0.1_f64, 432_usize, 2_000,  8_000, 0.001),
            (0.1_f64, 576_usize, 2_000,  8_000, 0.001),
            (0.1_f64, 720_usize, 1_000,  4_000, 0.001),
            (0.1_f64, 960_usize, 1_000,  4_000, 0.001),
        ];
        println!("\n  GUTOE Bohr convergence scan — backend: {backend:?}");
        println!("  {:>6}  {:>8}  {:>10}  {:>10}  {:>8}  {:>8}",
                 "L", "L/a₀", "E_total", "−α²/2", "C", "sites");
        println!("  {:>6}  {:>8}  {:>10}  {:>10}  {:>8}  {:>8}",
                 "─", "─", "─", "─", "─", "─");
        for (alpha, l, n_jacobi, n_iter, dtau) in configs {
            let cfg = SolverConfig { alpha, l, n_jacobi, n_iter, dtau, backend };
            let n = l * l * l;
            println!("  L={l} ({n} sites)…");
            let r = solve_hydrogen_3d(&cfg);
            println!("  {:>6}  {:>8.2}  {:>10.6}  {:>10.6}  {:>8.4}  {:>8}",
                     l, (l as f64) * alpha, r.e_total, r.bohr_3d, r.ratio, n);
            assert!(r.e_total < 0.0, "must be bound at L={l}");
        }
    }

    /// OBC α-scan: scan α at fixed L/a₀ ≈ 16 (L = 2⌈8/α⌉+1) with proper Jacobi convergence.
    /// n_jacobi ≈ 0.93×L² ensures <1% Jacobi residual (convergence rate ≈ π²/(2L²) per step).
    /// Box-size correction: C_∞(α) ≈ C_box(α) + A/L, A ≈ 19 (from L-scan at α=0.1).
    /// At fixed L/a₀=16: correction = 19α/16 ≈ 1.19α — significant at large α.
    ///
    /// cargo test -p gutoe-gpu --features rocm --release -- bohr_obc_scan --nocapture
    #[test]
    #[cfg(any(feature = "cuda", feature = "rocm"))]
    fn bohr_obc_scan() {
        let backend = detect_backend();
        // Scan α: smaller α = finer discretisation, C_∞(α→0) = lattice Bohr constant.
        // Jacobi convergence rate ≈ π²/(2L²) per step; need n_jacobi ≈ L² for 1% φ accuracy.
        // Imaginary-time: need τ = n_iter×dtau >> 1/ΔE where ΔE ≈ 3α²×C/8 ≈ 0.12α².
        // dtau=0.05: stable (dtau × max_eigenvalue ≈ 0.05×2 = 0.1 << 1).
        // n_jacobi ∝ L² = (16/α)²; n_iter ∝ 1/α² so τ∝1/α² ≈ 2/ΔE.
        let alphas: &[(f64, usize, usize, f64)] = &[
            //  α    n_sor    n_iter   dtau    L=16/α  ω≈2-2π/L  theory n_sor×2
            (0.30,    200,   4_000, 0.05),  // L=55,  ω=1.889,  n≈ 80
            (0.20,    300,   9_000, 0.05),  // L=81,  ω=1.925,  n≈120
            (0.10,    500,  34_000, 0.05),  // L=161, ω=1.962,  n≈200
            (0.07,    700,  70_000, 0.05),  // L=231, ω=1.973,  n≈270
            (0.05,   1000, 134_000, 0.05),  // L=321, ω=1.980,  n≈345
        ];
        println!("\n  GUTOE OBC α-scan (open BC, L/a₀≈16, n_jac∝L²) — backend: {backend:?}");
        println!("  {:>6}  {:>6}  {:>8}  {:>10}  {:>10}  {:>8}",
                 "α", "L_box", "sites", "E_total", "−α²/2", "C");
        println!("  {:>6}  {:>6}  {:>8}  {:>10}  {:>10}  {:>8}",
                 "─", "─", "─", "─", "─", "─");
        for &(alpha, n_jacobi, n_iter, dtau) in alphas {
            let r_active = ((8.0 / alpha).ceil() as usize).max(8);
            let l = 2 * r_active + 1;
            let n = l * l * l;
            println!("  α={alpha:.2} → L_box={l} ({n} sites)…");
            let res = solve_hydrogen_obc(alpha, n_jacobi, n_iter, dtau, backend);
            println!("  {:>6.2}  {:>6}  {:>8}  {:>10.6}  {:>10.6}  {:>8.4}",
                     alpha, l, n, res.e_total, res.bohr_3d, res.ratio);
            assert!(res.e_total < 0.0, "OBC: must be bound at α={alpha}");
        }
        // C_∞(α) ≈ C_box(α) + 19/L where L=16/α; from L-scan: C_∞(0.1) ≈ 0.550.
    }

    /// Single-point L=961: confirm Richardson extrapolation C_∞(α=0.1)=0.547.
    /// Predicted C(961) = 0.547 − 18.2/961 = 0.528.
    /// 961³=888M sites, ~36 GB peak (d_rho freed after Poisson) — fits tealc's 96 GiB GPU.
    /// Optimised: n_sor=2000 (theory ~700, 3× safety), n_iter=10000 (τ=500, 3.4× gap).
    ///
    /// cargo test -p gutoe-gpu --features rocm --release -- bohr_obc_l960 --nocapture
    #[test]
    #[cfg(any(feature = "cuda", feature = "rocm"))]
    fn bohr_obc_l960() {
        let alpha = 0.10_f64;
        let a0 = 1.0 / alpha;
        let backend = detect_backend();
        let bohr = -alpha * alpha / 2.0;

        let l = 961_usize;   // odd = centred proton; 961³=888M, ~36 GB peak
        let n_sor = 1_000;   // convergence diagnostic: 5× more SOR than previous run
        let n_iter = 30_000; // τ=1500, matching alpha-scan convergence depth
        let dtau = 0.05;
        let n = l * l * l;

        println!("\n  GUTOE OBC L=961 (α={alpha}, n_sor={n_sor}, n_iter={n_iter}, dtau={dtau}) — {backend:?}");
        println!("  L={l} ({n} sites, {:.1} GB φ+ψ)", n as f64 * 40.0 / 1e9);
        let (mut et, mut ek, mut ep) = (0.0_f64, 0.0_f64, 0.0_f64);
        unsafe {
            run_hydrogen_obc(
                alpha,
                l as i32, l as i32, l as i32,
                n_sor as i32, n_iter as i32, 10_i32,
                dtau,
                &mut et, &mut ek, &mut ep,
            );
        }
        let c = et / bohr;
        let c_pred = 0.547 - 18.2 / l as f64;
        println!("  L={l}  L/a₀={:.1}  E_total={et:.6}  −α²/2={bohr:.6}  C={c:.4}",
                 l as f64 / a0);
        println!("  Richardson prediction: C={c_pred:.4}  (Δ={:.4})", c - c_pred);
        assert!(et < 0.0, "OBC: must be bound at L={l}");
    }

    /// OBC L-scan: fix α=0.1, vary L to isolate finite-size vs α-dependence.
    /// n_jacobi ∝ L² ensures 1% Jacobi convergence regardless of box size.
    /// n_iter fixed so all points have the same imaginary-time τ=1700.
    ///
    /// cargo test -p gutoe-gpu --features rocm --release -- bohr_obc_lscan --nocapture
    #[test]
    #[cfg(any(feature = "cuda", feature = "rocm"))]
    fn bohr_obc_lscan() {
        let alpha = 0.10_f64;
        let a0 = 1.0 / alpha;   // Bohr radius = 10 lattice sites
        let backend = detect_backend();
        let bohr = -alpha * alpha / 2.0;

        // n_jacobi ≈ 0.93×L² for 1% Jacobi residual (rate ≈ π²/(2L²) per step)
        // n_iter = 34000, dtau = 0.05 → τ = 1700 for all (same imaginary time)
        let configs: &[(usize, usize, usize, f64)] = &[
            //   L    n_sor    n_iter   dtau    ω=2/(1+sin(π/L))  theory n_sor×2
            ( 161,    500,  34_000, 0.05),  // L/a₀=16, ω=1.962,  n≈200
            ( 241,    700,  34_000, 0.05),  // L/a₀=24, ω=1.974,  n≈270
            ( 321,   1000,  34_000, 0.05),  // L/a₀=32, ω=1.980,  n≈345
            ( 481,   1500,  34_000, 0.05),  // L/a₀=48, ω=1.987,  n≈950
            ( 960,   2000,  17_000, 0.05),  // L/a₀=96, ω=1.993,  n≈1050 (τ=850)
        ];
        println!("\n  GUTOE OBC L-scan (α={alpha}, open BC, n_iter=34K, dtau=0.05) — {backend:?}");
        println!("  {:>6}  {:>8}  {:>10}  {:>10}  {:>8}",
                 "L", "L/a₀", "E_total", "−α²/2", "C");
        println!("  {:>6}  {:>8}  {:>10}  {:>10}  {:>8}",
                 "─", "─", "─", "─", "─");
        for &(l, n_jacobi, n_iter, dtau) in configs {
            let n = l * l * l;
            println!("  L={l} ({n} sites, n_jacobi={n_jacobi})…");
            let (mut et, mut ek, mut ep) = (0.0_f64, 0.0_f64, 0.0_f64);
            unsafe {
                run_hydrogen_obc(
                    alpha,
                    l as i32, l as i32, l as i32,
                    n_jacobi as i32, n_iter as i32, 10_i32,
                    dtau,
                    &mut et, &mut ek, &mut ep,
                );
            }
            let c = et / bohr;
            println!("  {:>6}  {:>8.1}  {:>10.6}  {:>10.6}  {:>8.4}",
                     l, (l as f64) / a0, et, bohr, c);
            assert!(et < 0.0, "OBC L-scan: must be bound at L={l}");
        }
    }

    /// OBC L-scan at α=0.2: extract C_∞(0.2) via Richardson extrapolation.
    /// a₀=5 sites, everything converges fast: τ∝1/α²=25, need ~2000 steps.
    /// If C_∞(0.2) ≈ 0.547 → universal lattice constant (same as α=0.1).
    /// If C_∞(0.2) ≠ 0.547 → α-dependent, not universal.
    ///
    /// cargo test -p gutoe-gpu --features cuda --release -- bohr_obc_lscan_a02 --nocapture
    #[test]
    #[cfg(any(feature = "cuda", feature = "rocm"))]
    fn bohr_obc_lscan_a02() {
        let alpha = 0.20_f64;
        let a0 = 1.0 / alpha;   // Bohr radius = 5 lattice sites
        let backend = detect_backend();
        let bohr = -alpha * alpha / 2.0;

        // τ = n_iter × dtau.  ΔE ≈ 3α²C/8 ≈ 0.008.  Need τ >> 1/ΔE = 125.
        // τ = 500 gives e^{-ΔEτ} = e^{-4} ≈ 2% excited-state contamination.
        let configs: &[(usize, usize, usize, f64)] = &[
            //   L    n_sor    n_iter   dtau
            (  81,    200,  10_000, 0.05),  // L/a₀=16.2
            ( 121,    300,  10_000, 0.05),  // L/a₀=24.2
            ( 161,    400,  10_000, 0.05),  // L/a₀=32.2
            ( 241,    600,  10_000, 0.05),  // L/a₀=48.2
        ];
        println!("\n  GUTOE OBC L-scan α=0.2 (a₀=5) — {backend:?}");
        println!("  {:>6}  {:>8}  {:>10}  {:>10}  {:>8}",
                 "L", "L/a₀", "E_total", "−α²/2", "C");
        println!("  {:>6}  {:>8}  {:>10}  {:>10}  {:>8}",
                 "─", "─", "─", "─", "─");
        let mut results: Vec<(usize, f64)> = Vec::new();
        for &(l, n_sor, n_iter, dtau) in configs {
            let n = l * l * l;
            println!("  L={l} ({n} sites, n_sor={n_sor})…");
            let (mut et, mut ek, mut ep) = (0.0_f64, 0.0_f64, 0.0_f64);
            unsafe {
                run_hydrogen_obc(
                    alpha,
                    l as i32, l as i32, l as i32,
                    n_sor as i32, n_iter as i32, 10_i32,
                    dtau,
                    &mut et, &mut ek, &mut ep,
                );
            }
            let c = et / bohr;
            println!("  {:>6}  {:>8.1}  {:>10.6}  {:>10.6}  {:>8.4}",
                     l, (l as f64) / a0, et, bohr, c);
            assert!(et < 0.0, "OBC L-scan α=0.2: must be bound at L={l}");
            results.push((l, c));
        }
        // Richardson extrapolation: C(L) = C_∞ - B/L
        // From last two points: B = (C2-C1)/( 1/L1 - 1/L2 ), C_∞ = C2 + B/L2
        println!("\n  Richardson pair fits:");
        for i in 1..results.len() {
            let (l1, c1) = results[i-1];
            let (l2, c2) = results[i];
            let b = (c2 - c1) / (1.0/l1 as f64 - 1.0/l2 as f64);
            let c_inf = c2 + b / l2 as f64;
            println!("    ({l1},{l2}): B={b:.1}, C_∞={c_inf:.4}");
        }
    }

    /// OBC L-scan at α=0.07: a₀=14.3 sites. Finer discretization probe.
    ///
    /// cargo test -p gutoe-gpu --features rocm --release -- bohr_obc_lscan_a007 --nocapture
    #[test]
    #[cfg(any(feature = "cuda", feature = "rocm"))]
    fn bohr_obc_lscan_a007() {
        let alpha = 0.07_f64;
        let a0 = 1.0 / alpha;
        let backend = detect_backend();
        let bohr = -alpha * alpha / 2.0;

        // ΔE ≈ 3α²C/8 ≈ 0.001.  τ=500 → e^{-0.5} ≈ 60% — need more τ.
        // τ=2500 → e^{-2.5} ≈ 8%.  τ=5000 → e^{-5} ≈ 0.7%.
        let configs: &[(usize, usize, usize, f64)] = &[
            //   L    n_sor    n_iter   dtau
            ( 231,    600,  70_000, 0.05),  // L/a₀=16.2, τ=3500
            ( 341,    900,  70_000, 0.05),  // L/a₀=23.9
            ( 461,   1200,  70_000, 0.05),  // L/a₀=32.3
        ];
        println!("\n  GUTOE OBC L-scan α=0.07 (a₀=14.3) — {backend:?}");
        println!("  {:>6}  {:>8}  {:>10}  {:>10}  {:>8}",
                 "L", "L/a₀", "E_total", "−α²/2", "C");
        println!("  {:>6}  {:>8}  {:>10}  {:>10}  {:>8}",
                 "─", "─", "─", "─", "─");
        let mut results: Vec<(usize, f64)> = Vec::new();
        for &(l, n_sor, n_iter, dtau) in configs {
            let n = l * l * l;
            println!("  L={l} ({n} sites, n_sor={n_sor})…");
            let (mut et, mut ek, mut ep) = (0.0_f64, 0.0_f64, 0.0_f64);
            unsafe {
                run_hydrogen_obc(
                    alpha,
                    l as i32, l as i32, l as i32,
                    n_sor as i32, n_iter as i32, 10_i32,
                    dtau,
                    &mut et, &mut ek, &mut ep,
                );
            }
            let c = et / bohr;
            println!("  {:>6}  {:>8.1}  {:>10.6}  {:>10.6}  {:>8.4}",
                     l, (l as f64) / a0, et, bohr, c);
            assert!(et < 0.0, "OBC L-scan α=0.07: must be bound at L={l}");
            results.push((l, c));
        }
        println!("\n  Richardson pair fits:");
        for i in 1..results.len() {
            let (l1, c1) = results[i-1];
            let (l2, c2) = results[i];
            let b = (c2 - c1) / (1.0/l1 as f64 - 1.0/l2 as f64);
            let c_inf = c2 + b / l2 as f64;
            println!("    ({l1},{l2}): B={b:.1}, C_∞={c_inf:.4}");
        }
    }

    /// OBC L-scan at α=0.05: a₀=20 sites. Finest discretization we can do fast.
    ///
    /// cargo test -p gutoe-gpu --features rocm --release -- bohr_obc_lscan_a005 --nocapture
    #[test]
    #[cfg(any(feature = "cuda", feature = "rocm"))]
    fn bohr_obc_lscan_a005() {
        let alpha = 0.05_f64;
        let a0 = 1.0 / alpha;
        let backend = detect_backend();
        let bohr = -alpha * alpha / 2.0;

        // ΔE ≈ 3α²C/8 ≈ 0.0005.  Need τ >> 2000.
        // τ=6700 → e^{-3.35} ≈ 3.5%.
        let configs: &[(usize, usize, usize, f64)] = &[
            //   L    n_sor     n_iter   dtau
            ( 321,   1000,  134_000, 0.05),  // L/a₀=16.1, τ=6700
            ( 481,   1500,  134_000, 0.05),  // L/a₀=24.1
            ( 641,   2000,  134_000, 0.05),  // L/a₀=32.1
        ];
        println!("\n  GUTOE OBC L-scan α=0.05 (a₀=20) — {backend:?}");
        println!("  {:>6}  {:>8}  {:>10}  {:>10}  {:>8}",
                 "L", "L/a₀", "E_total", "−α²/2", "C");
        println!("  {:>6}  {:>8}  {:>10}  {:>10}  {:>8}",
                 "─", "─", "─", "─", "─");
        let mut results: Vec<(usize, f64)> = Vec::new();
        for &(l, n_sor, n_iter, dtau) in configs {
            let n = l * l * l;
            println!("  L={l} ({n} sites, n_sor={n_sor})…");
            let (mut et, mut ek, mut ep) = (0.0_f64, 0.0_f64, 0.0_f64);
            unsafe {
                run_hydrogen_obc(
                    alpha,
                    l as i32, l as i32, l as i32,
                    n_sor as i32, n_iter as i32, 10_i32,
                    dtau,
                    &mut et, &mut ek, &mut ep,
                );
            }
            let c = et / bohr;
            println!("  {:>6}  {:>8.1}  {:>10.6}  {:>10.6}  {:>8.4}",
                     l, (l as f64) / a0, et, bohr, c);
            assert!(et < 0.0, "OBC L-scan α=0.05: must be bound at L={l}");
            results.push((l, c));
        }
        println!("\n  Richardson pair fits:");
        for i in 1..results.len() {
            let (l1, c1) = results[i-1];
            let (l2, c2) = results[i];
            let b = (c2 - c1) / (1.0/l1 as f64 - 1.0/l2 as f64);
            let c_inf = c2 + b / l2 as f64;
            println!("    ({l1},{l2}): B={b:.1}, C_∞={c_inf:.4}");
        }
    }

    /// OBC L-scan at α=0.3: fastest scan, a₀=3.3 sites.
    /// Third independent C_∞ measurement.
    ///
    /// cargo test -p gutoe-gpu --features cuda --release -- bohr_obc_lscan_a03 --nocapture
    #[test]
    #[cfg(any(feature = "cuda", feature = "rocm"))]
    fn bohr_obc_lscan_a03() {
        let alpha = 0.30_f64;
        let a0 = 1.0 / alpha;
        let backend = detect_backend();
        let bohr = -alpha * alpha / 2.0;

        // ΔE ≈ 3α²C/8 ≈ 0.018.  τ=500 → e^{-9} ≈ 0.01%.  Very converged.
        let configs: &[(usize, usize, usize, f64)] = &[
            //   L    n_sor    n_iter   dtau
            (  55,    150,  10_000, 0.05),  // L/a₀=16.5
            (  81,    200,  10_000, 0.05),  // L/a₀=24.3
            ( 121,    300,  10_000, 0.05),  // L/a₀=36.3
            ( 161,    400,  10_000, 0.05),  // L/a₀=48.3
        ];
        println!("\n  GUTOE OBC L-scan α=0.3 (a₀=3.3) — {backend:?}");
        println!("  {:>6}  {:>8}  {:>10}  {:>10}  {:>8}",
                 "L", "L/a₀", "E_total", "−α²/2", "C");
        println!("  {:>6}  {:>8}  {:>10}  {:>10}  {:>8}",
                 "─", "─", "─", "─", "─");
        let mut results: Vec<(usize, f64)> = Vec::new();
        for &(l, n_sor, n_iter, dtau) in configs {
            let n = l * l * l;
            println!("  L={l} ({n} sites, n_sor={n_sor})…");
            let (mut et, mut ek, mut ep) = (0.0_f64, 0.0_f64, 0.0_f64);
            unsafe {
                run_hydrogen_obc(
                    alpha,
                    l as i32, l as i32, l as i32,
                    n_sor as i32, n_iter as i32, 10_i32,
                    dtau,
                    &mut et, &mut ek, &mut ep,
                );
            }
            let c = et / bohr;
            println!("  {:>6}  {:>8.1}  {:>10.6}  {:>10.6}  {:>8.4}",
                     l, (l as f64) / a0, et, bohr, c);
            assert!(et < 0.0, "OBC L-scan α=0.3: must be bound at L={l}");
            results.push((l, c));
        }
        println!("\n  Richardson pair fits:");
        for i in 1..results.len() {
            let (l1, c1) = results[i-1];
            let (l2, c2) = results[i];
            let b = (c2 - c1) / (1.0/l1 as f64 - 1.0/l2 as f64);
            let c_inf = c2 + b / l2 as f64;
            println!("    ({l1},{l2}): B={b:.1}, C_∞={c_inf:.4}");
        }
    }

    /// Boundary diagnostic: OBC vs PBC vs analytical correction at α=0.1.
    ///
    /// Tests all three hypotheses in one run:
    /// 1. Boundary |ψ|² — is the wavefunction hitting the walls?
    /// 2. PBC comparison — does removing grounded walls change C?
    /// 3. Analytical correction C + 2/(αL) — does image-charge formula work?
    ///
    /// cargo test -p gutoe-gpu --features cuda --release -- bohr_boundary_diagnostics --nocapture
    #[test]
    #[cfg(any(feature = "cuda", feature = "rocm"))]
    fn bohr_boundary_diagnostics() {
        let alpha = 0.10_f64;
        let a0 = 1.0 / alpha;
        let backend = detect_backend();
        let bohr = -alpha * alpha / 2.0;

        let configs: &[(usize, usize, usize, f64)] = &[
            //   L    n_sor    n_iter   dtau
            ( 161,    500,  34_000, 0.05),
            ( 241,    700,  34_000, 0.05),
            ( 321,   1000,  34_000, 0.05),
        ];

        println!("\n  ═══ GUTOE Boundary Diagnostics α={alpha} (a₀={a0}) — {backend:?} ═══");

        // ── Part 1: OBC with boundary |ψ|² ──
        println!("\n  ── OBC (grounded walls: ψ=0, φ=α/r) ──");
        println!("  {:>6}  {:>8}  {:>8}  {:>12}  {:>8}",
                 "L", "C_raw", "C_corr", "bnd_|ψ|²", "2/(αL)");
        let mut obc_results: Vec<(usize, f64)> = Vec::new();
        for &(l, n_sor, n_iter, dtau) in configs {
            let (mut et, mut ek, mut ep) = (0.0_f64, 0.0_f64, 0.0_f64);
            println!("  L={l} ({} sites, n_sor={n_sor})…", (l as u64).pow(3));
            unsafe {
                run_hydrogen_obc(
                    alpha,
                    l as i32, l as i32, l as i32,
                    n_sor as i32, n_iter as i32, 10_i32,
                    dtau,
                    &mut et, &mut ek, &mut ep,
                );
            }
            let c_raw = et / bohr;
            let correction = 2.0 / (alpha * l as f64);
            let c_corr = c_raw + correction;
            // boundary |ψ|² is printed by C code; we just print the rest
            println!("  {:>6}  {:>8.4}  {:>8.4}  {:>12}  {:>8.4}",
                     l, c_raw, c_corr, "(see above)", correction);
            obc_results.push((l, c_raw));
        }
        println!("\n  OBC Richardson pair fits:");
        for i in 1..obc_results.len() {
            let (l1, c1) = obc_results[i-1];
            let (l2, c2) = obc_results[i];
            let b = (c2 - c1) / (1.0/l1 as f64 - 1.0/l2 as f64);
            let c_inf = c2 + b / l2 as f64;
            println!("    ({l1},{l2}): B={b:.1}, C_∞={c_inf:.4}");
        }

        // ── Part 2: PBC (periodic boundaries, no walls) ──
        println!("\n  ── PBC (periodic: no walls, Madelung images) ──");
        println!("  {:>6}  {:>8}", "L", "C_pbc");
        let mut pbc_results: Vec<(usize, f64)> = Vec::new();
        for &(l, n_sor, n_iter, dtau) in configs {
            let (mut et, mut ek, mut ep) = (0.0_f64, 0.0_f64, 0.0_f64);
            // PBC needs more SOR sweeps (cold start, no Coulomb warm-start)
            let n_sor_pbc = n_sor * 3;
            println!("  L={l} ({} sites, n_sor={n_sor_pbc})…", (l as u64).pow(3));
            unsafe {
                run_hydrogen_pbc(
                    alpha,
                    l as i32, l as i32, l as i32,
                    n_sor_pbc as i32, n_iter as i32, 10_i32,
                    dtau,
                    &mut et, &mut ek, &mut ep,
                );
            }
            let c_pbc = et / bohr;
            println!("  {:>6}  {:>8.4}", l, c_pbc);
            pbc_results.push((l, c_pbc));
        }
        println!("\n  PBC Richardson pair fits:");
        for i in 1..pbc_results.len() {
            let (l1, c1) = pbc_results[i-1];
            let (l2, c2) = pbc_results[i];
            let b = (c2 - c1) / (1.0/l1 as f64 - 1.0/l2 as f64);
            let c_inf = c2 + b / l2 as f64;
            println!("    ({l1},{l2}): B={b:.1}, C_∞={c_inf:.4}");
        }

        // ── Part 3: Summary ──
        println!("\n  ═══ SUMMARY ═══");
        println!("  If OBC boundary |ψ|² is ~0: wavefunction doesn't reach walls, B is NOT image-charge.");
        println!("  If PBC C_∞ ≈ OBC C_∞: the 1/L correction is NOT from grounded-wall images.");
        println!("  If C_corr is constant across L: B = 2/(αL) is the exact correction.\n");
    }
}
