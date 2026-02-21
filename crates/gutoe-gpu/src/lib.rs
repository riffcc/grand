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

use gutoe_em::quantum_lepton::{bohr_test_3d, BohrResult};

// ── Backend enum ──────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Backend { Cpu, #[cfg(feature = "cuda")] Cuda, #[cfg(feature = "rocm")] Rocm }

pub fn detect_backend() -> Backend {
    #[cfg(feature = "rocm")] return Backend::Rocm;
    #[cfg(feature = "cuda")] return Backend::Cuda;
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
}
