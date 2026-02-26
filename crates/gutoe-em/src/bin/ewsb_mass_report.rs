//! EWSB/Higgs mass-sector closure report from structural Cl(1,3) inputs.

use gutoe_em::weak::{
    electroweak_vev_from_fermi, higgs_mass_from_vev, w_mass_from_vev_and_alpha,
    w_z_mass_ratio, z_mass_from_vev_and_alpha, ALPHA_EW_MZ, HIGGS_CRITICAL_VOID_FRACTION,
    HIGGS_QUARTIC_LAMBDA,
};
use std::fs::{self, File};
use std::io::Write;

fn main() {
    let g_f = 1.166_378_7e-5;
    let v = electroweak_vev_from_fermi(g_f);
    let m_w = w_mass_from_vev_and_alpha(v, ALPHA_EW_MZ);
    let m_z = z_mass_from_vev_and_alpha(v, ALPHA_EW_MZ);
    let m_h = higgs_mass_from_vev(v);

    let w_ref = 80.377_f64;
    let z_ref = 91.1876_f64;
    let h_ref = 125.25_f64;

    let dw = m_w - w_ref;
    let dz = m_z - z_ref;
    let dh = m_h - h_ref;

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
    writeln!(txt).expect("write");
    writeln!(txt, "[mass_sector]").expect("write");
    writeln!(txt, "g_f = {:.12e}", g_f).expect("write");
    writeln!(txt, "alpha_mz = {:.12}", ALPHA_EW_MZ).expect("write");
    writeln!(txt, "vev = {:.12}", v).expect("write");
    writeln!(txt, "m_w = {:.12}", m_w).expect("write");
    writeln!(txt, "m_z = {:.12}", m_z).expect("write");
    writeln!(txt, "m_h = {:.12}", m_h).expect("write");
    writeln!(txt, "delta_m_w = {:.12}", dw).expect("write");
    writeln!(txt, "delta_m_z = {:.12}", dz).expect("write");
    writeln!(txt, "delta_m_h = {:.12}", dh).expect("write");

    let mut json = File::create(&json_path).expect("create json");
    writeln!(
        json,
        "{{\n  \"higgs_lambda\": {:.12},\n  \"critical_void_fraction\": {:.12},\n  \"sin2_theta_w\": {:.12},\n  \"w_over_z_ratio\": {:.12},\n  \"inputs\": {{ \"g_f\": {:.12e}, \"alpha_mz\": {:.12} }},\n  \"masses\": {{\n    \"vev\": {:.12},\n    \"m_w\": {:.12},\n    \"m_z\": {:.12},\n    \"m_h\": {:.12},\n    \"delta_m_w\": {:.12},\n    \"delta_m_z\": {:.12},\n    \"delta_m_h\": {:.12}\n  }}\n}}",
        HIGGS_QUARTIC_LAMBDA,
        HIGGS_CRITICAL_VOID_FRACTION,
        3.0 / 13.0,
        w_z_mass_ratio(),
        g_f,
        ALPHA_EW_MZ,
        v,
        m_w,
        m_z,
        m_h,
        dw,
        dz,
        dh,
    )
    .expect("write json");

    println!("wrote {txt_path}");
    println!("wrote {json_path}");
    println!(
        "λ_H={:.6}, f_c={:.6}, v={:.3} GeV, m_W={:.3}, m_Z={:.3}, m_H={:.3}",
        HIGGS_QUARTIC_LAMBDA,
        HIGGS_CRITICAL_VOID_FRACTION,
        v,
        m_w,
        m_z,
        m_h
    );
}
