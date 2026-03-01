//! Topology constraint gate.
//!
//! Fails (exit code 1) if overdetermined ratio-lock identities drift beyond
//! tolerance.

use std::process::ExitCode;

fn env_f64(name: &str, default: f64) -> f64 {
    std::env::var(name)
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(default)
}

fn main() -> ExitCode {
    let tol = env_f64("GUTOE_TOPOLOGY_GATE_TOL", 1e-12).abs();

    let branching = 3.0_f64;
    let void_frac = 3.0_f64 / 16.0_f64;
    let eta = 4.0_f64 / 6.0_f64;
    let infra = 16.0_f64 / 6.0_f64;

    let g = branching * void_frac * eta * infra;
    let residual = (g - 1.0).abs();

    let inferred_infra = 1.0 / (branching * void_frac * eta);
    let inferred_eta = 1.0 / (branching * void_frac * infra);
    let inferred_void = 1.0 / (branching * eta * infra);
    let inferred_branching = 1.0 / (void_frac * eta * infra);

    let errs = [
        residual,
        (inferred_infra - infra).abs(),
        (inferred_eta - eta).abs(),
        (inferred_void - void_frac).abs(),
        (inferred_branching - branching).abs(),
    ];
    let max_err = errs.iter().copied().fold(0.0_f64, f64::max);

    println!("[topology_constraint_gate]");
    println!("tol={:.12e}", tol);
    println!("gain={:.12e} residual={:.12e}", g, residual);
    println!(
        "errors: infra={:.12e} eta={:.12e} void={:.12e} branching={:.12e}",
        (inferred_infra - infra).abs(),
        (inferred_eta - eta).abs(),
        (inferred_void - void_frac).abs(),
        (inferred_branching - branching).abs()
    );
    println!("max_err={:.12e}", max_err);

    if max_err <= tol {
        println!("status=PASS");
        ExitCode::SUCCESS
    } else {
        println!("status=FAIL");
        ExitCode::from(1)
    }
}

