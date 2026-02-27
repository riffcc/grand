//! GRAND-353 uncertainty propagation report across upstream cosmology lanes.

use gutoe_physics::{
    evaluate_uncertainty, DistributionSummary, UncertaintyAssumptions, UniverseAssumptions,
    UniverseWindows,
};
use std::fs::{self, File};
use std::io::Write;

fn write_summary_line(txt: &mut File, name: &str, s: DistributionSummary) {
    writeln!(
        txt,
        "{}: mean={:.9e}, std={:.9e}, p05={:.9e}, p50={:.9e}, p95={:.9e}, min={:.9e}, max={:.9e}",
        name, s.mean, s.std, s.p05, s.p50, s.p95, s.min, s.max
    )
    .expect("write summary");
}

fn main() {
    let ua = UncertaintyAssumptions::default();
    let s = evaluate_uncertainty(
        UniverseAssumptions::default(),
        UniverseWindows::default(),
        ua,
    );

    let out_dir =
        std::env::var("GUTOE_UNCERTAINTY_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let txt_path = format!("{out_dir}/uncertainty_report.txt");
    let json_path = format!("{out_dir}/uncertainty_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[sampling]").expect("write");
    writeln!(txt, "requested_samples = {}", s.requested_samples).expect("write");
    writeln!(txt, "valid_samples = {}", s.valid_samples).expect("write");
    writeln!(txt, "pass_fraction = {:.9}", s.pass_fraction).expect("write");
    writeln!(
        txt,
        "inflation_pass_fraction = {:.9}",
        s.inflation_pass_fraction
    )
    .expect("write");
    writeln!(
        txt,
        "baryogenesis_pass_fraction = {:.9}",
        s.baryogenesis_pass_fraction
    )
    .expect("write");
    writeln!(txt, "bbn_pass_fraction = {:.9}", s.bbn_pass_fraction).expect("write");
    writeln!(txt, "dark_pass_fraction = {:.9}", s.dark_pass_fraction).expect("write");
    writeln!(
        txt,
        "transfer_pass_fraction = {:.9}",
        s.transfer_pass_fraction
    )
    .expect("write");
    writeln!(
        txt,
        "microphysics_pass_fraction = {:.9}",
        s.microphysics_pass_fraction
    )
    .expect("write");
    writeln!(
        txt,
        "background_pass_fraction = {:.9}",
        s.background_pass_fraction
    )
    .expect("write");

    writeln!(txt).expect("write");
    writeln!(txt, "[distributions]").expect("write");
    write_summary_line(&mut txt, "n_s", s.n_s);
    write_summary_line(&mut txt, "A_s", s.a_s);
    write_summary_line(&mut txt, "eta10", s.eta10);
    write_summary_line(&mut txt, "dm_fraction", s.dm_fraction);
    write_summary_line(&mut txt, "H0_km_s_mpc", s.h0_km_s_mpc);
    write_summary_line(&mut txt, "age_gyr", s.age_gyr);
    write_summary_line(&mut txt, "r_s_drag_mpc", s.rs_drag_mpc);
    write_summary_line(&mut txt, "theta_star_rad", s.theta_star_rad);
    write_summary_line(&mut txt, "l_peak1", s.l_peak1);
    write_summary_line(&mut txt, "l_peak2", s.l_peak2);
    write_summary_line(&mut txt, "Yp_network", s.yp_network);
    write_summary_line(&mut txt, "D_H_network", s.dh_network);
    write_summary_line(&mut txt, "z_visibility_peak", s.z_visibility_peak);

    let mut json = File::create(&json_path).expect("create json");
    writeln!(
        json,
        "{{\n  \"sampling\": {{\"requested_samples\": {}, \"valid_samples\": {}, \"pass_fraction\": {:.9}, \"inflation_pass_fraction\": {:.9}, \"baryogenesis_pass_fraction\": {:.9}, \"bbn_pass_fraction\": {:.9}, \"dark_pass_fraction\": {:.9}, \"transfer_pass_fraction\": {:.9}, \"microphysics_pass_fraction\": {:.9}, \"background_pass_fraction\": {:.9}}},",
        s.requested_samples,
        s.valid_samples,
        s.pass_fraction,
        s.inflation_pass_fraction,
        s.baryogenesis_pass_fraction,
        s.bbn_pass_fraction,
        s.dark_pass_fraction,
        s.transfer_pass_fraction,
        s.microphysics_pass_fraction,
        s.background_pass_fraction,
    )
    .expect("write json");

    macro_rules! write_dist {
        ($name:literal, $d:expr, $comma:expr) => {
            writeln!(
                json,
                "  \"{}\": {{\"mean\": {:.9e}, \"std\": {:.9e}, \"p05\": {:.9e}, \"p50\": {:.9e}, \"p95\": {:.9e}, \"min\": {:.9e}, \"max\": {:.9e}}}{}",
                $name,
                $d.mean,
                $d.std,
                $d.p05,
                $d.p50,
                $d.p95,
                $d.min,
                $d.max,
                $comma,
            )
            .expect("write dist");
        };
    }

    writeln!(json, "  \"distributions\": {{").expect("write");
    write_dist!("n_s", s.n_s, ",");
    write_dist!("a_s", s.a_s, ",");
    write_dist!("eta10", s.eta10, ",");
    write_dist!("dm_fraction", s.dm_fraction, ",");
    write_dist!("h0_km_s_mpc", s.h0_km_s_mpc, ",");
    write_dist!("age_gyr", s.age_gyr, ",");
    write_dist!("rs_drag_mpc", s.rs_drag_mpc, ",");
    write_dist!("theta_star_rad", s.theta_star_rad, ",");
    write_dist!("l_peak1", s.l_peak1, ",");
    write_dist!("l_peak2", s.l_peak2, ",");
    write_dist!("yp_network", s.yp_network, ",");
    write_dist!("dh_network", s.dh_network, ",");
    write_dist!("z_visibility_peak", s.z_visibility_peak, "");
    writeln!(json, "  }}").expect("write");
    writeln!(json, "}}").expect("write");

    println!("wrote {txt_path}");
    println!("wrote {json_path}");
    println!(
        "Uncertainty: valid={} pass={:.3} H0(p50)={:.3} theta*(p50)={:.6e} Yp(p50)={:.5}",
        s.valid_samples, s.pass_fraction, s.h0_km_s_mpc.p50, s.theta_star_rad.p50, s.yp_network.p50,
    );
}
