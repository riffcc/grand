//! Overdetermined topology ratio probe.
//!
//! Checks the closure:
//!   G = branching * void * eta * infra = 1
//! and solves each factor from the other three to verify overdetermination.

use serde_json::json;
use std::fs;
use std::path::PathBuf;

fn main() {
    let out_dir = std::env::var("GUTOE_TOPOLOGY_OVERDETERMINATION_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/topology_overdetermination_probe".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    // Structural counts/invariants
    let branching = 3.0_f64; // |Z3|
    let void_frac = 3.0_f64 / 16.0_f64; // 3/16
    let grade1 = 4.0_f64;
    let grade2 = 6.0_f64;
    let total = 16.0_f64;
    let eta = grade1 / grade2; // 2/3
    let infra = total / grade2; // 8/3
    let even_supp = 8.0_f64 / 16.0_f64; // 1/2

    let g = branching * void_frac * eta * infra;

    // Infer each factor from unit closure and the other three.
    let inferred_infra = 1.0 / (branching * void_frac * eta);
    let inferred_eta = 1.0 / (branching * void_frac * infra);
    let inferred_void = 1.0 / (branching * eta * infra);
    let inferred_branching = 1.0 / (void_frac * eta * infra);

    let r16_over_9 = eta * infra;
    let r3_over_8 = branching * void_frac * eta;
    let r3_over_2 = branching * void_frac * infra;

    let residual = (g - 1.0).abs();

    let payload = json!({
      "structural_counts": {
        "branching_z3": branching,
        "void_fraction": void_frac,
        "grade1_count": grade1,
        "grade2_count": grade2,
        "total_basis": total,
        "even_suppression": even_supp
      },
      "ratios": {
        "eta_grade1_over_grade2": eta,
        "infra_total_over_grade2": infra,
        "eta_times_infra": r16_over_9,
        "branching_void_eta": r3_over_8,
        "branching_void_infra": r3_over_2
      },
      "closure": {
        "gain": g,
        "target": 1.0,
        "abs_residual": residual,
        "is_closed": residual < 1e-12
      },
      "overdetermination": {
        "inferred_infra": inferred_infra,
        "inferred_eta": inferred_eta,
        "inferred_void": inferred_void,
        "inferred_branching": inferred_branching,
        "infra_abs_error": (inferred_infra - infra).abs(),
        "eta_abs_error": (inferred_eta - eta).abs(),
        "void_abs_error": (inferred_void - void_frac).abs(),
        "branching_abs_error": (inferred_branching - branching).abs()
      },
      "interpretation": "any one ratio is forced by the other three under unit closure"
    });

    let txt_path = out.join("topology_overdetermination_probe.txt");
    let json_path = out.join("topology_overdetermination_probe.json");

    let mut txt = String::new();
    txt.push_str("[topology_overdetermination_probe]\n");
    txt.push_str(&format!(
        "G = b*void*eta*infra = {:.12e} (residual {:.12e})\n",
        g, residual
    ));
    txt.push_str(&format!(
        "b={:.12e}, void={:.12e}, eta={:.12e}, infra={:.12e}\n",
        branching, void_frac, eta, infra
    ));
    txt.push_str(&format!(
        "inferred: infra={:.12e}, eta={:.12e}, void={:.12e}, branching={:.12e}\n",
        inferred_infra, inferred_eta, inferred_void, inferred_branching
    ));

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(&json_path, serde_json::to_string_pretty(&payload).expect("json")).expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
}

