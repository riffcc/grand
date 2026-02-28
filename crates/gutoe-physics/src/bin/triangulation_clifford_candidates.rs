//! Candidate Cl(1,3) reconstructions for triangulated constants.
//!
//! This lane evaluates low-complexity structural formulas against the frozen
//! triangulated targets `(p_ratio, kappa_geo, ew_coeff_required)`.

use gutoe_em::triangulate_neutrino_from_splittings;

const SOLAR_DM21_TARGET_EV2: f64 = 7.53e-5;
const ATMOSPHERIC_DM32_TARGET_EV2: f64 = 2.453e-3;
const SIN2_THETA_W_MZ_TARGET: f64 = 0.23122;

// Shared Cl(1,3) counts.
const D: f64 = 16.0;
const SU2: f64 = 3.0;
const GRADE1: f64 = 4.0;
const GRADE2: f64 = 6.0;
const LATTICE_SHIFT: f64 = GRADE2 + 1.0; // 7
const COMPLEMENT: f64 = D - SU2; // 13
const TOTAL_GAUGE: f64 = 12.0;
const ALPHA_INV: f64 = 137.0;
const T16: f64 = 136.0; // triangular(16)

fn rel_err(obs: f64, tgt: f64) -> f64 {
    if tgt.abs() < 1e-30 {
        0.0
    } else {
        (obs - tgt) / tgt
    }
}

fn main() {
    let tri = triangulate_neutrino_from_splittings(SOLAR_DM21_TARGET_EV2, ATMOSPHERIC_DM32_TARGET_EV2);
    let p_target = tri.p_triangulated;
    let kappa_target = tri.kappa_geo;

    // Reference from frozen EW target:
    // coeff_required = (sin2_target - 3/13) / alpha^2
    let alpha = 1.0 / ALPHA_INV;
    let ew_coeff_target = (SIN2_THETA_W_MZ_TARGET - (3.0 / 13.0)) / (alpha * alpha);

    // Candidate A: p = α⁻¹/(|grade1|+|grade2|) - 1/(lattice_shift * total_gauge)
    let p_candidate = ALPHA_INV / (GRADE1 + GRADE2) - 1.0 / (LATTICE_SHIFT * TOTAL_GAUGE);

    // Candidate B:
    // κ = (60/11) * [ (d+|SU2|)/|SU2| + 1/|grade2|² + 1/(lattice_shift * complement * T16) ]
    let kappa_candidate =
        (60.0 / 11.0) * ((D + SU2) / SU2 + 1.0 / (GRADE2 * GRADE2) + 1.0 / (LATTICE_SHIFT * COMPLEMENT * T16));

    // Candidate C: ew_coeff = d/2 + grade2/complement - 1/(lattice_shift * T16)
    let ew_coeff_candidate = D / 2.0 + GRADE2 / COMPLEMENT - 1.0 / (LATTICE_SHIFT * T16);

    let p_rel = rel_err(p_candidate, p_target);
    let kappa_rel = rel_err(kappa_candidate, kappa_target);
    let ew_rel = rel_err(ew_coeff_candidate, ew_coeff_target);

    println!("triangulation_clifford_candidates");
    println!("targets:");
    println!("  p_ratio_target = {:.12}", p_target);
    println!("  kappa_target = {:.12}", kappa_target);
    println!("  ew_coeff_target = {:.12}", ew_coeff_target);
    println!();
    println!("candidates:");
    println!("  p_candidate = {:.12}   rel_err = {:+.3e}", p_candidate, p_rel);
    println!(
        "  kappa_candidate = {:.12}   rel_err = {:+.3e}",
        kappa_candidate, kappa_rel
    );
    println!(
        "  ew_coeff_candidate = {:.12}   rel_err = {:+.3e}",
        ew_coeff_candidate, ew_rel
    );
    println!();
    println!("formulas:");
    println!("  p = 137/10 - 1/(7*12)");
    println!("  kappa = (60/11) * (19/3 + 1/36 + 1/(7*13*136))");
    println!("  ew_coeff = 8 + 6/13 - 1/(7*136)");
}
