//! GRAND-216: Holonomy diagnostics from lattice parallel transport.

use gutoe_em::{
    sample_holonomy_diagnostics, u1_geometric_phase, u1_phase_composition_residual, LatticeConfig,
    RestrictedHolonomySignature, Su2Links,
};
use gutoe_physics::StandardModelDynamicsMap;
use rand::rngs::StdRng;
use rand::SeedableRng;
use serde_json::json;
use std::fs;
use std::path::PathBuf;

const WILSON_RESIDUAL_EPS: f64 = 1.0e-12;
const PHASE_RESIDUAL_EPS: f64 = 1.0e-12;

fn env_parse<T: std::str::FromStr>(key: &str, default: T) -> T {
    std::env::var(key)
        .ok()
        .and_then(|s| s.parse::<T>().ok())
        .unwrap_or(default)
}

fn main() {
    let out_dir = std::env::var("GUTOE_HOLONOMY_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/holonomy".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let cfg = LatticeConfig {
        hex_rows: env_parse("HOLO_ROWS", 18usize),
        hex_cols: env_parse("HOLO_COLS", 18usize),
        layers: env_parse("HOLO_LAYERS", 1usize),
        ..Default::default()
    };
    let seed: u64 = env_parse("HOLO_SEED", 137_u64);
    let beta: f64 = env_parse("HOLO_BETA", 1.0_f64);
    let eps: f64 = env_parse("HOLO_EPS", 0.45_f64);
    let n_therm: usize = env_parse("HOLO_THERM", 250usize);
    let max_samples: usize = env_parse("HOLO_SAMPLES", 64usize);

    let mut rng = StdRng::seed_from_u64(seed);
    let mut links = Su2Links::hot_start(&mut rng, &cfg);
    for _ in 0..n_therm {
        links.metropolis_sweep(&mut rng, beta, eps, &cfg);
    }
    let diagnostics = sample_holonomy_diagnostics(&links, &cfg, max_samples);

    let signature = RestrictedHolonomySignature::from_clifford_z3();
    let sm = StandardModelDynamicsMap::from_clifford_z3();
    let restricted_holonomy_ok = signature.recovers_sm()
        && signature.su3_generators == sm.su3_generators
        && signature.su2_generators == sm.su2_generators
        && signature.u1_generators == sm.u1_generators
        && signature.total_generators == sm.total_gauge_generators;

    let wilson_bridge_ok = diagnostics.max_wilson_residual_abs < WILSON_RESIDUAL_EPS;

    // Berry/geometric phase witness: exact U(1) composition and unitarity.
    let theta = std::f64::consts::TAU / signature.total_generators as f64;
    let phase = u1_geometric_phase(theta);
    let phase_unitarity_residual = (phase.norm() - 1.0).abs();
    let phase_composition_residual = u1_phase_composition_residual(theta, -theta);
    let geometric_phase_ok = phase_unitarity_residual < PHASE_RESIDUAL_EPS
        && phase_composition_residual < PHASE_RESIDUAL_EPS;

    let passes_all = restricted_holonomy_ok && wilson_bridge_ok && geometric_phase_ok;

    let txt_path = out.join("holonomy_report.txt");
    let json_path = out.join("holonomy_report.json");

    let mut txt = String::new();
    txt.push_str("[meta]\n");
    txt.push_str("lane = GRAND-216_holonomy\n");
    txt.push_str(&format!(
        "lattice = {}x{}x{}\nseed = {}\nbeta = {:.6}\neps = {:.6}\ntherm_sweeps = {}\nmax_samples = {}\n\n",
        cfg.hex_rows, cfg.hex_cols, cfg.layers, seed, beta, eps, n_therm, max_samples
    ));

    txt.push_str("[restricted_holonomy]\n");
    txt.push_str(&format!("z3_order = {}\n", signature.z3_order));
    txt.push_str(&format!("u1_generators = {}\n", signature.u1_generators));
    txt.push_str(&format!("su2_generators = {}\n", signature.su2_generators));
    txt.push_str(&format!("su3_generators = {}\n", signature.su3_generators));
    txt.push_str(&format!("total_generators = {}\n", signature.total_generators));
    txt.push_str(&format!(
        "runtime_map_total_generators = {}\n\n",
        sm.total_gauge_generators
    ));

    txt.push_str("[holonomy_samples]\n");
    txt.push_str(&format!("sample_count = {}\n", diagnostics.samples.len()));
    txt.push_str(&format!(
        "mean_trace_over_2 = {:.12e}\n",
        diagnostics.mean_trace_over_2
    ));
    txt.push_str(&format!(
        "mean_class_angle_rad = {:.12e}\n",
        diagnostics.mean_class_angle_rad
    ));
    txt.push_str(&format!(
        "max_wilson_residual_abs = {:.12e}\n",
        diagnostics.max_wilson_residual_abs
    ));
    for (idx, s) in diagnostics.samples.iter().take(12).enumerate() {
        txt.push_str(&format!(
            "sample_{idx:02}: tri=({},{},{}), tr/2={:.12e}, theta={:.12e}, wilson_resid={:.12e}\n",
            s.i, s.j, s.k, s.trace_over_2, s.class_angle_rad, s.wilson_residual_abs
        ));
    }
    txt.push('\n');

    txt.push_str("[geometric_phase]\n");
    txt.push_str(&format!("theta = {:.12e}\n", theta));
    txt.push_str(&format!("phase_re = {:.12e}\n", phase.re));
    txt.push_str(&format!("phase_im = {:.12e}\n", phase.im));
    txt.push_str(&format!(
        "phase_unitarity_residual = {:.12e}\n",
        phase_unitarity_residual
    ));
    txt.push_str(&format!(
        "phase_composition_residual = {:.12e}\n\n",
        phase_composition_residual
    ));

    txt.push_str("[gate]\n");
    txt.push_str(&format!("restricted_holonomy_ok = {}\n", restricted_holonomy_ok));
    txt.push_str(&format!("wilson_bridge_ok = {}\n", wilson_bridge_ok));
    txt.push_str(&format!("geometric_phase_ok = {}\n", geometric_phase_ok));
    txt.push_str(&format!("passes_all = {}\n", passes_all));

    let payload = json!({
        "meta": {
            "lane": "GRAND-216_holonomy",
            "lattice": {
                "hex_rows": cfg.hex_rows,
                "hex_cols": cfg.hex_cols,
                "layers": cfg.layers,
            },
            "sampling": {
                "seed": seed,
                "beta": beta,
                "eps": eps,
                "therm_sweeps": n_therm,
                "max_samples": max_samples,
            },
            "thresholds": {
                "wilson_residual_eps": WILSON_RESIDUAL_EPS,
                "phase_residual_eps": PHASE_RESIDUAL_EPS,
            }
        },
        "restricted_holonomy": {
            "signature": signature,
            "runtime_map": {
                "u1_generators": sm.u1_generators,
                "su2_generators": sm.su2_generators,
                "su3_generators": sm.su3_generators,
                "total_gauge_generators": sm.total_gauge_generators,
            }
        },
        "holonomy_samples": diagnostics,
        "geometric_phase": {
            "theta": theta,
            "phase_re": phase.re,
            "phase_im": phase.im,
            "phase_unitarity_residual": phase_unitarity_residual,
            "phase_composition_residual": phase_composition_residual,
        },
        "gate": {
            "restricted_holonomy_ok": restricted_holonomy_ok,
            "wilson_bridge_ok": wilson_bridge_ok,
            "geometric_phase_ok": geometric_phase_ok,
            "passes_all": passes_all,
        }
    });

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(
        &json_path,
        serde_json::to_string_pretty(&payload).expect("serialize json"),
    )
    .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "holonomy: pass={} max_wilson_resid={:.3e} phase_comp_resid={:.3e}",
        passes_all, diagnostics.max_wilson_residual_abs, phase_composition_residual
    );
}
