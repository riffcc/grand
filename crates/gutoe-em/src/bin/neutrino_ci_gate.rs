//! GRAND-342 neutrino hierarchy + mass-character CI gate.
//!
//! Hard-gates:
//! - hierarchy prediction is normal ordering
//! - texture mass-character lane is structurally Dirac
//! - absolute tiny-mass transduction remains nonzero and below observational caps

use gutoe_em::alpha::ALPHA_INVERSE_PHYSICAL;
use gutoe_em::{
    electron_mass_from_proton_anchor, neutrino_dirac_majorana_prediction, neutrino_hierarchy_prediction,
    neutrino_majorana_symmetry_residual, neutrino_texture_eigenvalues,
};
use std::fs::{self, File};
use std::io::Write;
use std::process;

fn main() {
    let out_dir =
        std::env::var("GUTOE_NEUTRINO_GATE_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let json_path = format!("{out_dir}/neutrino_ci_gate.json");
    let mut json = File::create(&json_path).expect("create gate json");

    let m_tex = neutrino_texture_eigenvalues();
    let hierarchy = neutrino_hierarchy_prediction();
    let mass_character = neutrino_dirac_majorana_prediction();
    let majorana_symmetry_residual = neutrino_majorana_symmetry_residual();

    // Absolute tiny-mass transduction lane (same as neutrino_tiny_mass_report).
    let mut raw = m_tex.map(|x| x.abs());
    raw.sort_by(|a, b| a.total_cmp(b));
    let alpha = 1.0 / ALPHA_INVERSE_PHYSICAL;
    let me_ev = electron_mass_from_proton_anchor() * 1.0e6;
    let m_scale_ev = me_ev * alpha.powi(4) * (60.0 / 11.0);
    let raw_max = raw[2].max(1.0e-18);
    let m1_ev = m_scale_ev * (raw[0] / raw_max);
    let m2_ev = m_scale_ev * (raw[1] / raw_max);
    let m3_ev = m_scale_ev;
    let sum_ev = m1_ev + m2_ev + m3_ev;

    let dm21 = m_tex[1] - m_tex[0];
    let dm31 = m_tex[2] - m_tex[0];

    let katrin_cap_ev = 0.8;
    let cosmology_sum_cap_ev = 0.12;
    let majorana_symmetry_max = 1.0e-12;

    let hierarchy_ok = hierarchy == "normal";
    let mass_character_ok = mass_character == "dirac";
    let nonzero_ok = m1_ev > 0.0 && m2_ev > 0.0 && m3_ev > 0.0;
    let tiny_ok = m3_ev < katrin_cap_ev && sum_ev < cosmology_sum_cap_ev;
    let majorana_excluded_ok = majorana_symmetry_residual > majorana_symmetry_max;
    let overall_pass =
        hierarchy_ok && mass_character_ok && nonzero_ok && tiny_ok && majorana_excluded_ok;

    writeln!(
        json,
        "{{\n  \"overall_pass\": {},\n  \"windows\": {{\"katrin_cap_ev\": {:.12}, \"cosmology_sum_cap_ev\": {:.12}, \"majorana_symmetry_max\": {:.12e}}},\n  \"texture_lane\": {{\"m1\": {:.12e}, \"m2\": {:.12e}, \"m3\": {:.12e}, \"delta_m21\": {:.12e}, \"delta_m31\": {:.12e}, \"hierarchy_prediction\": \"{}\", \"mass_character_prediction\": \"{}\", \"majorana_symmetry_residual\": {:.12e}}},\n  \"absolute_lane\": {{\"alpha\": {:.12}, \"electron_mass_anchor_ev\": {:.12e}, \"mass_scale_ev\": {:.12e}, \"m1_ev\": {:.12e}, \"m2_ev\": {:.12e}, \"m3_ev\": {:.12e}, \"sum_ev\": {:.12e}}},\n  \"checks\": {{\"hierarchy_ok\": {}, \"mass_character_ok\": {}, \"majorana_excluded_ok\": {}, \"nonzero_ok\": {}, \"tiny_ok\": {}}}\n}}",
        if overall_pass { "true" } else { "false" },
        katrin_cap_ev,
        cosmology_sum_cap_ev,
        majorana_symmetry_max,
        m_tex[0],
        m_tex[1],
        m_tex[2],
        dm21,
        dm31,
        hierarchy,
        mass_character,
        majorana_symmetry_residual,
        alpha,
        me_ev,
        m_scale_ev,
        m1_ev,
        m2_ev,
        m3_ev,
        sum_ev,
        if hierarchy_ok { "true" } else { "false" },
        if mass_character_ok { "true" } else { "false" },
        if majorana_excluded_ok { "true" } else { "false" },
        if nonzero_ok { "true" } else { "false" },
        if tiny_ok { "true" } else { "false" }
    )
    .expect("write gate json");

    println!(
        "neutrino gate: pass={} (hierarchy={}, type={}, majorana_resid={:.3e}, m3={:.3e} eV, sum={:.3e} eV)",
        overall_pass, hierarchy, mass_character, majorana_symmetry_residual, m3_ev, sum_ev
    );
    println!("wrote {json_path}");

    if !overall_pass {
        eprintln!(
            "FAIL: hierarchy_ok={} mass_character_ok={} majorana_excluded_ok={} nonzero_ok={} tiny_ok={}",
            hierarchy_ok, mass_character_ok, majorana_excluded_ok, nonzero_ok, tiny_ok
        );
        process::exit(2);
    }
}
