//! Tiny nonzero neutrino-mass closure report.
//!
//! Uses texture eigenvalue ordering from Cl(1,3) and a structural suppression
//! scale from the electron lane to produce absolute eV-scale masses.

use gutoe_em::alpha::ALPHA_INVERSE_PHYSICAL;
use gutoe_em::{electron_mass_from_proton_anchor, neutrino_hierarchy_prediction, neutrino_texture_eigenvalues};
use std::fs::{self, File};
use std::io::Write;

fn main() {
    let out_dir = std::env::var("GUTOE_NEUTRINO_TINY_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/neutrino_tiny_mass_report".to_string());
    let _ = fs::create_dir_all(&out_dir);

    // Texture lane (dimensionless shape/eigen-ordering).
    let mut raw = neutrino_texture_eigenvalues().map(|x| x.abs());
    raw.sort_by(|a, b| a.total_cmp(b));
    let hierarchy = neutrino_hierarchy_prediction();

    // Structural absolute scale:
    // m_scale = m_e * α^4 * (60/11), using shared dark/visible geometric factor.
    let alpha = 1.0 / ALPHA_INVERSE_PHYSICAL;
    let me_ev = electron_mass_from_proton_anchor() * 1.0e6;
    let m_scale_ev = me_ev * alpha.powi(4) * (60.0 / 11.0);

    // Normalize to the largest texture eigenvalue, then apply structural scale.
    let raw_max = raw[2];
    let m1_ev = m_scale_ev * (raw[0] / raw_max);
    let m2_ev = m_scale_ev * (raw[1] / raw_max);
    let m3_ev = m_scale_ev;

    let sum_ev = m1_ev + m2_ev + m3_ev;
    let dm21_ev2 = m2_ev * m2_ev - m1_ev * m1_ev;
    let dm31_ev2 = m3_ev * m3_ev - m1_ev * m1_ev;

    let katrin_cap_ev = 0.8;
    let cosmology_sum_cap_ev = 0.12;

    let nonzero_ok = m1_ev > 0.0 && m2_ev > 0.0 && m3_ev > 0.0;
    let tiny_ok = m3_ev < katrin_cap_ev && sum_ev < cosmology_sum_cap_ev;
    let normal_ok = hierarchy == "normal";
    let overall_pass = nonzero_ok && tiny_ok && normal_ok;

    let txt_path = format!("{out_dir}/neutrino_tiny_mass_report.txt");
    let json_path = format!("{out_dir}/neutrino_tiny_mass_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "alpha = {:.12}", alpha).expect("write");
    writeln!(txt, "electron_mass_anchor_ev = {:.12e}", me_ev).expect("write");
    writeln!(txt, "mass_scale_ev = {:.12e}", m_scale_ev).expect("write");
    writeln!(txt, "hierarchy_prediction = {}", hierarchy).expect("write");
    writeln!(txt, "").expect("write");
    writeln!(txt, "[masses_ev]").expect("write");
    writeln!(txt, "m1 = {:.12e}", m1_ev).expect("write");
    writeln!(txt, "m2 = {:.12e}", m2_ev).expect("write");
    writeln!(txt, "m3 = {:.12e}", m3_ev).expect("write");
    writeln!(txt, "sum = {:.12e}", sum_ev).expect("write");
    writeln!(txt, "delta_m21_ev2 = {:.12e}", dm21_ev2).expect("write");
    writeln!(txt, "delta_m31_ev2 = {:.12e}", dm31_ev2).expect("write");
    writeln!(txt, "").expect("write");
    writeln!(txt, "[bounds]").expect("write");
    writeln!(txt, "katrin_cap_ev = {:.6}", katrin_cap_ev).expect("write");
    writeln!(txt, "cosmology_sum_cap_ev = {:.6}", cosmology_sum_cap_ev).expect("write");
    writeln!(txt, "nonzero_ok = {}", nonzero_ok).expect("write");
    writeln!(txt, "tiny_ok = {}", tiny_ok).expect("write");
    writeln!(txt, "normal_ok = {}", normal_ok).expect("write");
    writeln!(txt, "overall_pass = {}", overall_pass).expect("write");

    let mut json = File::create(&json_path).expect("create json");
    writeln!(
        json,
        "{{\n  \"alpha\": {:.12},\n  \"electron_mass_anchor_ev\": {:.12e},\n  \"mass_scale_ev\": {:.12e},\n  \"hierarchy_prediction\": \"{}\",\n  \"masses_ev\": {{\"m1\": {:.12e}, \"m2\": {:.12e}, \"m3\": {:.12e}, \"sum\": {:.12e}}},\n  \"delta_m2_ev2\": {{\"dm21\": {:.12e}, \"dm31\": {:.12e}}},\n  \"bounds\": {{\"katrin_cap_ev\": {:.6}, \"cosmology_sum_cap_ev\": {:.6}, \"nonzero_ok\": {}, \"tiny_ok\": {}, \"normal_ok\": {}, \"overall_pass\": {}}}\n}}",
        alpha,
        me_ev,
        m_scale_ev,
        hierarchy,
        m1_ev,
        m2_ev,
        m3_ev,
        sum_ev,
        dm21_ev2,
        dm31_ev2,
        katrin_cap_ev,
        cosmology_sum_cap_ev,
        if nonzero_ok { "true" } else { "false" },
        if tiny_ok { "true" } else { "false" },
        if normal_ok { "true" } else { "false" },
        if overall_pass { "true" } else { "false" }
    )
    .expect("write json");

    println!("wrote {txt_path}");
    println!("wrote {json_path}");
    println!(
        "neutrino tiny masses: overall_pass={} (m3={:.3e} eV, sum={:.3e} eV)",
        overall_pass, m3_ev, sum_ev
    );
}

