//! GRAND-368: homochirality parity-violation energy split report.

use gutoe_physics::{
    amino_acid_enantiomer_split_ev, amino_backbone_nuclear_factor, amino_backbone_parity_factor,
    amino_backbone_parity_proxy, handedness_energy_shift_ev, preferred_amino_handedness,
    rydberg_energy_structural_ev, weak_electron_scale_dimensionless, weak_nuclear_charge,
    Handedness, CHIRAL_PROJECTION_FACTOR, NITROGEN_N, NITROGEN_Z, OXYGEN_N, OXYGEN_Z,
    WEAK_GAUGE_FRACTION,
};
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

fn main() {
    let out_dir = std::env::var("GUTOE_HOMO_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/homochirality_report".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let qn = weak_nuclear_charge(NITROGEN_Z, NITROGEN_N);
    let qo = weak_nuclear_charge(OXYGEN_Z, OXYGEN_N);
    let nuclear = amino_backbone_nuclear_factor();
    let parity = amino_backbone_parity_factor();
    let proxy = amino_backbone_parity_proxy();
    let weak_scale = weak_electron_scale_dimensionless();
    let e_ryd = rydberg_energy_structural_ev();
    let split = amino_acid_enantiomer_split_ev();
    let e_left = handedness_energy_shift_ev(Handedness::Left);
    let e_right = handedness_energy_shift_ev(Handedness::Right);

    let txt_path = out.join("homochirality_report.txt");
    let json_path = out.join("homochirality_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[homochirality_structural]").expect("write");
    writeln!(txt, "weak_gauge_fraction = {:.12}", WEAK_GAUGE_FRACTION).expect("write");
    writeln!(txt, "chiral_projection_factor = {:.12}", CHIRAL_PROJECTION_FACTOR).expect("write");
    writeln!(txt, "QW_N14 = {:.12}", qn).expect("write");
    writeln!(txt, "QW_O16 = {:.12}", qo).expect("write");
    writeln!(txt, "backbone_nuclear_factor = {:.12}", nuclear).expect("write");
    writeln!(txt, "backbone_parity_factor = {:.12}", parity).expect("write");
    writeln!(txt, "alpha_suppressed_proxy = {:.12e}", proxy).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[scale]").expect("write");
    writeln!(txt, "weak_electron_scale_GF_me2 = {:.12e}", weak_scale).expect("write");
    writeln!(txt, "rydberg_structural_ev = {:.12}", e_ryd).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[prediction]").expect("write");
    writeln!(txt, "deltaE_LR_ev = {:.12e}", split).expect("write");
    writeln!(txt, "left_shift_ev = {:.12e}", e_left).expect("write");
    writeln!(txt, "right_shift_ev = {:.12e}", e_right).expect("write");
    writeln!(txt, "preferred_handedness = {:?}", preferred_amino_handedness()).expect("write");

    let mut json = File::create(&json_path).expect("create json");
    writeln!(
        json,
        "{{\n  \"structural\": {{\"weak_gauge_fraction\": {:.12}, \"chiral_projection_factor\": {:.12}, \"QW_N14\": {:.12}, \"QW_O16\": {:.12}, \"backbone_nuclear_factor\": {:.12}, \"backbone_parity_factor\": {:.12}, \"alpha_suppressed_proxy\": {:.12e}}},\n  \"scale\": {{\"weak_electron_scale_GF_me2\": {:.12e}, \"rydberg_structural_ev\": {:.12}}},\n  \"prediction\": {{\"deltaE_LR_ev\": {:.12e}, \"left_shift_ev\": {:.12e}, \"right_shift_ev\": {:.12e}, \"preferred_handedness\": \"{:?}\"}}\n}}",
        WEAK_GAUGE_FRACTION,
        CHIRAL_PROJECTION_FACTOR,
        qn,
        qo,
        nuclear,
        parity,
        proxy,
        weak_scale,
        e_ryd,
        split,
        e_left,
        e_right,
        preferred_amino_handedness(),
    )
    .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "Homochirality: ΔE_LR={:.3e} eV, left_shift={:.3e}, right_shift={:.3e}, preferred={:?}",
        split,
        e_left,
        e_right,
        preferred_amino_handedness()
    );
}

