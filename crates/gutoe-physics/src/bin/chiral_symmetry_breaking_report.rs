//! GRAND-126: chiral symmetry breaking structural report.

use gutoe_physics::{
    chiral_explicit_breaking_alpha, confinement_chiral_link_strength, pion_mass_proxy,
    pion_mass_sq_from_explicit_breaking, pion_mass_sq_proxy, pseudo_goldstone_ratio,
    quark_condensate_proxy,
};
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

fn main() {
    let out_dir = std::env::var("GUTOE_CHIRAL_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/chiral_symmetry_breaking".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let condensate = quark_condensate_proxy();
    let alpha = chiral_explicit_breaking_alpha();
    let pion_sq = pion_mass_sq_proxy();
    let pion = pion_mass_proxy();
    let ratio = pseudo_goldstone_ratio();
    let chiral_limit_pion_sq = pion_mass_sq_from_explicit_breaking(0.0);
    let link = confinement_chiral_link_strength();

    let condensate_nonzero = condensate < 0.0;
    let pseudo_goldstone_ok = (ratio - alpha).abs() < 1.0e-15;
    let chiral_limit_ok = chiral_limit_pion_sq.abs() < 1.0e-18;
    let confinement_link_ok = link > 0.0;
    let passes_all =
        condensate_nonzero && pseudo_goldstone_ok && chiral_limit_ok && confinement_link_ok;

    let txt_path = out.join("chiral_symmetry_breaking_report.txt");
    let json_path = out.join("chiral_symmetry_breaking_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[chiral_symmetry_breaking_structural]").expect("write");
    writeln!(txt, "quark_condensate_proxy = {:.12e}", condensate).expect("write");
    writeln!(txt, "explicit_breaking_alpha = {:.12e}", alpha).expect("write");
    writeln!(txt, "pion_mass_sq_proxy = {:.12e}", pion_sq).expect("write");
    writeln!(txt, "pion_mass_proxy = {:.12e}", pion).expect("write");
    writeln!(txt, "pseudo_goldstone_ratio = {:.12e}", ratio).expect("write");
    writeln!(txt, "chiral_limit_pion_mass_sq = {:.12e}", chiral_limit_pion_sq).expect("write");
    writeln!(txt, "confinement_chiral_link_strength = {:.12e}", link).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[gate]").expect("write");
    writeln!(txt, "condensate_nonzero = {}", condensate_nonzero).expect("write");
    writeln!(txt, "pseudo_goldstone_ok = {}", pseudo_goldstone_ok).expect("write");
    writeln!(txt, "chiral_limit_ok = {}", chiral_limit_ok).expect("write");
    writeln!(txt, "confinement_link_ok = {}", confinement_link_ok).expect("write");
    writeln!(txt, "passes_all = {}", passes_all).expect("write");

    let mut json = File::create(&json_path).expect("create json");
    writeln!(
        json,
        "{{\n  \"structural\": {{\"quark_condensate_proxy\": {:.12e}, \"explicit_breaking_alpha\": {:.12e}, \"pion_mass_sq_proxy\": {:.12e}, \"pion_mass_proxy\": {:.12e}, \"pseudo_goldstone_ratio\": {:.12e}, \"chiral_limit_pion_mass_sq\": {:.12e}, \"confinement_chiral_link_strength\": {:.12e}}},\n  \"gate\": {{\"condensate_nonzero\": {}, \"pseudo_goldstone_ok\": {}, \"chiral_limit_ok\": {}, \"confinement_link_ok\": {}, \"passes_all\": {}}}\n}}",
        condensate,
        alpha,
        pion_sq,
        pion,
        ratio,
        chiral_limit_pion_sq,
        link,
        condensate_nonzero,
        pseudo_goldstone_ok,
        chiral_limit_ok,
        confinement_link_ok,
        passes_all
    )
    .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "chiral_symmetry_breaking: condensate={:.6e}, pion_sq={:.6e}, ratio={:.6e}, passes_all={}",
        condensate, pion_sq, ratio, passes_all
    );
}
