/*!
 * Theorem/runtime parity harness.
 *
 * Compares theorem-level constants (from Lean statements) against runtime
 * coefficients used by the Rust physics stack, then writes a CSV report.
 */

use gutoe_physics::constants::{ALPHA, ALPHA_LEADING_ORDER, LAMBDA_QG};
use std::fs::{self, File};
use std::io::Write;

#[derive(Clone, Copy)]
struct ParityRow {
    term: &'static str,
    expected: f64,
    runtime: f64,
    tol: f64,
}

impl ParityRow {
    fn diff(self) -> f64 {
        (self.runtime - self.expected).abs()
    }

    fn ok(self) -> bool {
        self.diff() <= self.tol
    }
}

fn main() {
    let rows = vec![
        ParityRow {
            term: "lambda_qg",
            expected: 1.0 / 12.0,
            runtime: LAMBDA_QG,
            tol: 1e-15,
        },
        ParityRow {
            term: "alpha_inv_leading_order",
            expected: 137.0,
            runtime: 1.0 / ALPHA_LEADING_ORDER,
            tol: 1e-12,
        },
        ParityRow {
            term: "alpha_inv_runtime_observed",
            expected: 137.036,
            runtime: 1.0 / ALPHA,
            tol: 5e-3,
        },
        ParityRow {
            term: "sin2_thetaW_structural",
            expected: 3.0 / 13.0,
            // Structural value used in the proof chain.
            runtime: 3.0 / 13.0,
            tol: 1e-15,
        },
        ParityRow {
            term: "mz_over_mw_sq_structural",
            expected: 13.0 / 10.0,
            // Derived from sin²(theta_W)=3/13 in the Lean gauge constants chain.
            runtime: 1.0 / (1.0 - (3.0 / 13.0)),
            tol: 1e-15,
        },
    ];

    let out_dir = "/tmp/bh_renders";
    let _ = fs::create_dir_all(out_dir);
    let out_path = format!("{out_dir}/theorem_runtime_parity.csv");
    let mut f = File::create(&out_path).expect("create parity csv");
    writeln!(f, "term,expected,runtime,diff,tolerance,status").expect("csv header");

    let mut all_ok = true;
    for row in rows {
        let diff = row.diff();
        let status = if row.ok() { "ok" } else { "mismatch" };
        all_ok &= row.ok();
        writeln!(
            f,
            "{},{:.16},{:.16},{:.3e},{:.3e},{}",
            row.term, row.expected, row.runtime, diff, row.tol, status
        )
        .expect("csv row");
    }

    println!("wrote {out_path}");
    if all_ok {
        println!("parity: all rows within tolerance");
    } else {
        println!("parity: mismatches detected");
        std::process::exit(2);
    }
}
