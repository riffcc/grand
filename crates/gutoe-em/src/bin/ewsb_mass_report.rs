//! EWSB/Higgs mass-sector closure report from structural Cl(1,3) inputs.

use gutoe_em::weak::{
    electroweak_vev_from_fermi, electroweak_vev_from_lattice_order_parameter,
    electron_mass_from_proton_anchor, higgs_mass_from_vev, normalized_higgs_order_parameter,
    w_mass_from_vev_and_alpha, w_z_mass_ratio, z_mass_from_vev_and_alpha, ALPHA_EW_MZ,
    EWSB_SCALE_FACTOR, HIGGS_CRITICAL_VOID_FRACTION, HIGGS_QUARTIC_LAMBDA, PROTON_MASS_ANCHOR_MEV,
    VEV_OVER_PROTON,
};
use std::fs::{self, File};
use std::io::Write;

fn main() {
    let g_f = 1.166_378_7e-5;
    let f0_vac = 1.0;
    let order = normalized_higgs_order_parameter(f0_vac);

    let v_fermi = electroweak_vev_from_fermi(g_f);
    let m_w_fermi = w_mass_from_vev_and_alpha(v_fermi, ALPHA_EW_MZ);
    let m_z_fermi = z_mass_from_vev_and_alpha(v_fermi, ALPHA_EW_MZ);
    let m_h_fermi = higgs_mass_from_vev(v_fermi);

    let v_lattice = electroweak_vev_from_lattice_order_parameter(f0_vac);
    let m_w_lattice = w_mass_from_vev_and_alpha(v_lattice, ALPHA_EW_MZ);
    let m_z_lattice = z_mass_from_vev_and_alpha(v_lattice, ALPHA_EW_MZ);
    let m_h_lattice = higgs_mass_from_vev(v_lattice);

    let w_ref = 80.377_f64;
    let z_ref = 91.1876_f64;
    let h_ref = 125.25_f64;

    let dw_fermi = m_w_fermi - w_ref;
    let dz_fermi = m_z_fermi - z_ref;
    let dh_fermi = m_h_fermi - h_ref;

    let dw_lattice = m_w_lattice - w_ref;
    let dz_lattice = m_z_lattice - z_ref;
    let dh_lattice = m_h_lattice - h_ref;

    let out_dir = "/tmp/bh_renders";
    let _ = fs::create_dir_all(out_dir);
    let txt_path = format!("{out_dir}/ewsb_mass_report.txt");
    let json_path = format!("{out_dir}/ewsb_mass_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[ewsb_structural]").expect("write");
    writeln!(txt, "higgs_lambda = {:.12}", HIGGS_QUARTIC_LAMBDA).expect("write");
    writeln!(txt, "critical_void_fraction = {:.12}", HIGGS_CRITICAL_VOID_FRACTION).expect("write");
    writeln!(txt, "sin2_theta_w = {:.12}", 3.0 / 13.0).expect("write");
    writeln!(txt, "w_over_z_ratio = {:.12}", w_z_mass_ratio()).expect("write");
    writeln!(txt, "ewsb_scale_factor = {:.12}", EWSB_SCALE_FACTOR).expect("write");
    writeln!(txt, "proton_mass_anchor_mev = {:.12}", PROTON_MASS_ANCHOR_MEV).expect("write");
    writeln!(txt, "electron_mass_anchor_mev = {:.12}", electron_mass_from_proton_anchor()).expect("write");
    writeln!(txt, "vev_over_proton = {:.12}", VEV_OVER_PROTON).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[mass_sector_fermi_branch]").expect("write");
    writeln!(txt, "g_f = {:.12e}", g_f).expect("write");
    writeln!(txt, "alpha_mz = {:.12}", ALPHA_EW_MZ).expect("write");
    writeln!(txt, "vev = {:.12}", v_fermi).expect("write");
    writeln!(txt, "m_w = {:.12}", m_w_fermi).expect("write");
    writeln!(txt, "m_z = {:.12}", m_z_fermi).expect("write");
    writeln!(txt, "m_h = {:.12}", m_h_fermi).expect("write");
    writeln!(txt, "delta_m_w = {:.12}", dw_fermi).expect("write");
    writeln!(txt, "delta_m_z = {:.12}", dz_fermi).expect("write");
    writeln!(txt, "delta_m_h = {:.12}", dh_fermi).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[mass_sector_lattice_branch]").expect("write");
    writeln!(txt, "f0 = {:.12}", f0_vac).expect("write");
    writeln!(txt, "order = {:.12}", order).expect("write");
    writeln!(txt, "alpha_mz = {:.12}", ALPHA_EW_MZ).expect("write");
    writeln!(txt, "vev = {:.12}", v_lattice).expect("write");
    writeln!(txt, "m_w = {:.12}", m_w_lattice).expect("write");
    writeln!(txt, "m_z = {:.12}", m_z_lattice).expect("write");
    writeln!(txt, "m_h = {:.12}", m_h_lattice).expect("write");
    writeln!(txt, "delta_m_w = {:.12}", dw_lattice).expect("write");
    writeln!(txt, "delta_m_z = {:.12}", dz_lattice).expect("write");
    writeln!(txt, "delta_m_h = {:.12}", dh_lattice).expect("write");

    let mut json = File::create(&json_path).expect("create json");
    writeln!(
        json,
        "{{\n  \"higgs_lambda\": {:.12},\n  \"critical_void_fraction\": {:.12},\n  \"sin2_theta_w\": {:.12},\n  \"w_over_z_ratio\": {:.12},\n  \"ewsb_scale_factor\": {:.12},\n  \"proton_mass_anchor_mev\": {:.12},\n  \"electron_mass_anchor_mev\": {:.12},\n  \"vev_over_proton\": {:.12},\n  \"fermi_branch\": {{\n    \"inputs\": {{ \"g_f\": {:.12e}, \"alpha_mz\": {:.12} }},\n    \"masses\": {{\n      \"vev\": {:.12},\n      \"m_w\": {:.12},\n      \"m_z\": {:.12},\n      \"m_h\": {:.12},\n      \"delta_m_w\": {:.12},\n      \"delta_m_z\": {:.12},\n      \"delta_m_h\": {:.12}\n    }}\n  }},\n  \"lattice_branch\": {{\n    \"inputs\": {{ \"f0\": {:.12}, \"order\": {:.12}, \"alpha_mz\": {:.12} }},\n    \"masses\": {{\n      \"vev\": {:.12},\n      \"m_w\": {:.12},\n      \"m_z\": {:.12},\n      \"m_h\": {:.12},\n      \"delta_m_w\": {:.12},\n      \"delta_m_z\": {:.12},\n      \"delta_m_h\": {:.12}\n    }}\n  }}\n}}",
        HIGGS_QUARTIC_LAMBDA,
        HIGGS_CRITICAL_VOID_FRACTION,
        3.0 / 13.0,
        w_z_mass_ratio(),
        EWSB_SCALE_FACTOR,
        PROTON_MASS_ANCHOR_MEV,
        electron_mass_from_proton_anchor(),
        VEV_OVER_PROTON,
        g_f,
        ALPHA_EW_MZ,
        v_fermi,
        m_w_fermi,
        m_z_fermi,
        m_h_fermi,
        dw_fermi,
        dz_fermi,
        dh_fermi,
        f0_vac,
        order,
        ALPHA_EW_MZ,
        v_lattice,
        m_w_lattice,
        m_z_lattice,
        m_h_lattice,
        dw_lattice,
        dz_lattice,
        dh_lattice,
    )
    .expect("write json");

    println!("wrote {txt_path}");
    println!("wrote {json_path}");
    println!(
        "Fermi branch:   v={:.3} GeV, m_W={:.3}, m_Z={:.3}, m_H={:.3}",
        v_fermi,
        m_w_fermi,
        m_z_fermi,
        m_h_fermi
    );
    println!(
        "Lattice branch: v={:.3} GeV, m_W={:.3}, m_Z={:.3}, m_H={:.3}",
        v_lattice,
        m_w_lattice,
        m_z_lattice,
        m_h_lattice
    );
    println!(
        "λ_H={:.6}, f_c={:.6}, v/mp={:.6}",
        HIGGS_QUARTIC_LAMBDA,
        HIGGS_CRITICAL_VOID_FRACTION,
        VEV_OVER_PROTON,
    );
}
