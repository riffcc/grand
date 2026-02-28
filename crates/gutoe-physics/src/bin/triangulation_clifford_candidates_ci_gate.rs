//! CI gate for candidate Cl(1,3) reconstructions of triangulated constants.

use gutoe_em::triangulate_neutrino_from_splittings;
use std::fs::{self, File};
use std::io::Write;
use std::process;

const SOLAR_DM21_TARGET_EV2: f64 = 7.53e-5;
const ATMOSPHERIC_DM32_TARGET_EV2: f64 = 2.453e-3;
const SIN2_THETA_W_MZ_TARGET: f64 = 0.23122;

const D: f64 = 16.0;
const SU2: f64 = 3.0;
const GRADE2: f64 = 6.0;
const LATTICE_SHIFT: f64 = GRADE2 + 1.0; // 7
const COMPLEMENT: f64 = D - SU2; // 13
const TOTAL_GAUGE: f64 = 12.0;
const ALPHA_INV: f64 = 137.0;
const T16: f64 = 136.0;

const P_REL_MAX: f64 = 2.0e-6;
const KAPPA_REL_MAX: f64 = 1.0e-6;
const EW_REL_MAX: f64 = 1.0e-7;

fn rel_err(obs: f64, tgt: f64) -> f64 {
    if tgt.abs() < 1.0e-30 {
        0.0
    } else {
        (obs - tgt) / tgt
    }
}

fn main() {
    let out_dir = std::env::var("GUTOE_TRIANG_CLIFFORD_CAND_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let json_path = format!("{out_dir}/triangulation_clifford_candidates_ci_gate.json");

    let tri = triangulate_neutrino_from_splittings(SOLAR_DM21_TARGET_EV2, ATMOSPHERIC_DM32_TARGET_EV2);
    let p_target = tri.p_triangulated;
    let kappa_target = tri.kappa_geo;
    let alpha = 1.0 / ALPHA_INV;
    let ew_target = (SIN2_THETA_W_MZ_TARGET - (3.0 / 13.0)) / (alpha * alpha);

    let p_candidate = ALPHA_INV / 10.0 - 1.0 / (LATTICE_SHIFT * TOTAL_GAUGE);
    let kappa_candidate =
        (60.0 / 11.0) * ((D + SU2) / SU2 + 1.0 / (GRADE2 * GRADE2) + 1.0 / (LATTICE_SHIFT * COMPLEMENT * T16));
    let ew_candidate = D / 2.0 + GRADE2 / COMPLEMENT - 1.0 / (LATTICE_SHIFT * T16);

    let p_rel = rel_err(p_candidate, p_target);
    let kappa_rel = rel_err(kappa_candidate, kappa_target);
    let ew_rel = rel_err(ew_candidate, ew_target);

    let p_ok = p_rel.abs() <= P_REL_MAX;
    let kappa_ok = kappa_rel.abs() <= KAPPA_REL_MAX;
    let ew_ok = ew_rel.abs() <= EW_REL_MAX;
    let overall_pass = p_ok && kappa_ok && ew_ok;

    let mut json = File::create(&json_path).expect("create candidate gate json");
    writeln!(
        json,
        "{{\n  \"overall_pass\": {},\n  \"windows\": {{\"p_rel_max\": {:.3e}, \"kappa_rel_max\": {:.3e}, \"ew_rel_max\": {:.3e}}},\n  \"targets\": {{\"p\": {:.12}, \"kappa\": {:.12}, \"ew_coeff\": {:.12}}},\n  \"candidates\": {{\"p\": {:.12}, \"kappa\": {:.12}, \"ew_coeff\": {:.12}}},\n  \"residuals\": {{\"p_rel\": {:.12e}, \"kappa_rel\": {:.12e}, \"ew_rel\": {:.12e}}},\n  \"checks\": {{\"p_ok\": {}, \"kappa_ok\": {}, \"ew_ok\": {}}}\n}}",
        if overall_pass { "true" } else { "false" },
        P_REL_MAX,
        KAPPA_REL_MAX,
        EW_REL_MAX,
        p_target,
        kappa_target,
        ew_target,
        p_candidate,
        kappa_candidate,
        ew_candidate,
        p_rel,
        kappa_rel,
        ew_rel,
        if p_ok { "true" } else { "false" },
        if kappa_ok { "true" } else { "false" },
        if ew_ok { "true" } else { "false" }
    )
    .expect("write candidate gate json");

    println!(
        "triangulation_clifford_candidates_ci_gate: pass={} p_rel={:.3e} kappa_rel={:.3e} ew_rel={:.3e}",
        overall_pass, p_rel, kappa_rel, ew_rel
    );
    println!("wrote {json_path}");

    if !overall_pass {
        eprintln!(
            "FAIL: p_ok={} kappa_ok={} ew_ok={} p_rel={:.3e} kappa_rel={:.3e} ew_rel={:.3e}",
            p_ok, kappa_ok, ew_ok, p_rel, kappa_rel, ew_rel
        );
        process::exit(2);
    }
}
