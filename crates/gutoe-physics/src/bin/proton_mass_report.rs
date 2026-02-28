use gutoe_physics::StandardModelDynamicsMap;
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

const ELECTRON_MASS_MEV_OBS: f64 = 0.510_998_950;
const PROTON_MASS_MEV_OBS: f64 = 938.272_088_16;

fn triangular(n: u32) -> u32 {
    n * (n + 1) / 2
}

fn main() {
    let out_dir = std::env::var("GUTOE_PROTON_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/proton_mass_report".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let sm = StandardModelDynamicsMap::from_clifford_z3();
    let mp_me_struct = sm.total_gauge_generators * triangular(sm.clifford_dim + 1);

    // Anchor route A (current implemented route): electron mass anchor -> proton prediction.
    let proton_pred_from_e = ELECTRON_MASS_MEV_OBS * mp_me_struct as f64;
    let proton_rel_err = (proton_pred_from_e - PROTON_MASS_MEV_OBS) / PROTON_MASS_MEV_OBS;

    // Reverse route B: proton anchor -> electron prediction.
    let electron_pred_from_p = PROTON_MASS_MEV_OBS / mp_me_struct as f64;
    let electron_rel_err = (electron_pred_from_p - ELECTRON_MASS_MEV_OBS) / ELECTRON_MASS_MEV_OBS;

    let txt_path = out.join("proton_mass_report.txt");
    let json_path = out.join("proton_mass_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[structural]").expect("write");
    writeln!(txt, "clifford_dim = {}", sm.clifford_dim).expect("write");
    writeln!(txt, "total_gauge_generators = {}", sm.total_gauge_generators).expect("write");
    writeln!(txt, "triangular(clifford_dim+1) = {}", triangular(sm.clifford_dim + 1)).expect("write");
    writeln!(txt, "mp_me_struct = {}", mp_me_struct).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[route_a_electron_anchor]").expect("write");
    writeln!(txt, "electron_mass_mev_obs = {:.12}", ELECTRON_MASS_MEV_OBS).expect("write");
    writeln!(txt, "proton_mass_mev_pred = {:.12}", proton_pred_from_e).expect("write");
    writeln!(txt, "proton_mass_mev_obs = {:.12}", PROTON_MASS_MEV_OBS).expect("write");
    writeln!(txt, "proton_rel_error = {:.12e}", proton_rel_err).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[route_b_proton_anchor]").expect("write");
    writeln!(txt, "proton_mass_mev_obs = {:.12}", PROTON_MASS_MEV_OBS).expect("write");
    writeln!(txt, "electron_mass_mev_pred = {:.12}", electron_pred_from_p).expect("write");
    writeln!(txt, "electron_mass_mev_obs = {:.12}", ELECTRON_MASS_MEV_OBS).expect("write");
    writeln!(txt, "electron_rel_error = {:.12e}", electron_rel_err).expect("write");

    let mut json = File::create(&json_path).expect("create json");
    writeln!(
        json,
        "{{\n  \"structural\": {{\"clifford_dim\": {}, \"total_gauge_generators\": {}, \"triangular_clifford_dim_plus_1\": {}, \"mp_me_struct\": {}}},\n  \"route_a_electron_anchor\": {{\"electron_mass_mev_obs\": {:.12}, \"proton_mass_mev_pred\": {:.12}, \"proton_mass_mev_obs\": {:.12}, \"proton_rel_error\": {:.12e}}},\n  \"route_b_proton_anchor\": {{\"proton_mass_mev_obs\": {:.12}, \"electron_mass_mev_pred\": {:.12}, \"electron_mass_mev_obs\": {:.12}, \"electron_rel_error\": {:.12e}}}\n}}",
        sm.clifford_dim,
        sm.total_gauge_generators,
        triangular(sm.clifford_dim + 1),
        mp_me_struct,
        ELECTRON_MASS_MEV_OBS,
        proton_pred_from_e,
        PROTON_MASS_MEV_OBS,
        proton_rel_err,
        PROTON_MASS_MEV_OBS,
        electron_pred_from_p,
        ELECTRON_MASS_MEV_OBS,
        electron_rel_err
    )
    .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "mp/me={} | proton_pred={:.6} MeV vs {:.6} (rel err {:.3e})",
        mp_me_struct, proton_pred_from_e, PROTON_MASS_MEV_OBS, proton_rel_err
    );
}
