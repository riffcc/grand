//! GRAND-350: end-to-end universe simulation report from assembled GUTOE lanes.

use gutoe_physics::{evaluate_universe_gate, UniverseAssumptions, UniverseWindows, SEC_PER_YEAR};
use std::fs::{self, File};
use std::io::Write;

fn main() {
    let assumptions = UniverseAssumptions::default();
    let windows = UniverseWindows::default();
    let score = evaluate_universe_gate(assumptions, windows);

    let out_dir =
        std::env::var("GUTOE_UNIVERSE_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let txt_path = format!("{out_dir}/universe_sim_report.txt");
    let json_path = format!("{out_dir}/universe_sim_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[universe_assembly]").expect("write");
    writeln!(txt, "lambda_full = {:.12e}", score.lambda_full).expect("write");
    writeln!(txt, "H0_km_s_mpc = {:.9}", score.h0_km_s_mpc).expect("write");
    writeln!(txt, "H0_rel_error = {:.12}", score.h0_rel_error).expect("write");
    writeln!(txt, "omega_b0 = {:.12}", score.omega_b0).expect("write");
    writeln!(txt, "omega_dm0 = {:.12}", score.omega_dm0).expect("write");
    writeln!(txt, "omega_m0 = {:.12}", score.omega_m0).expect("write");
    writeln!(txt, "omega_r0 = {:.12}", score.omega_r0).expect("write");
    writeln!(txt, "omega_lambda0 = {:.12}", score.omega_lambda0).expect("write");
    writeln!(txt, "age_gyr = {:.9}", score.age_gyr).expect("write");
    writeln!(
        txt,
        "recombination_age_kyr = {:.6}",
        score.recombination_age_kyr
    )
    .expect("write");
    writeln!(txt, "bbn_age_seconds = {:.6}", score.bbn_age_seconds).expect("write");
    writeln!(txt, "r_s_drag_mpc = {:.9}", score.transfer.rs_drag_mpc).expect("write");
    writeln!(
        txt,
        "theta_star_rad = {:.9e}",
        score.transfer.theta_star_rad
    )
    .expect("write");
    writeln!(txt, "l_peak1 = {:.6}", score.transfer.l_peak1).expect("write");
    writeln!(txt, "l_peak2 = {:.6}", score.transfer.l_peak2).expect("write");
    writeln!(txt, "micro_yp = {:.9}", score.microphysics.yp_network).expect("write");
    writeln!(txt, "micro_dh = {:.12e}", score.microphysics.dh_network).expect("write");
    writeln!(
        txt,
        "micro_z_visibility_peak = {:.6}",
        score.microphysics.z_visibility_peak
    )
    .expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[pipeline_gates]").expect("write");
    writeln!(txt, "inflation_ok = {}", score.inflation_ok).expect("write");
    writeln!(txt, "baryogenesis_ok = {}", score.baryogenesis_ok).expect("write");
    writeln!(txt, "bbn_ok = {}", score.bbn_ok).expect("write");
    writeln!(txt, "microphysics_ok = {}", score.microphysics_ok).expect("write");
    writeln!(
        txt,
        "dark_matter_unified_ok = {}",
        score.dark_matter_unified_ok
    )
    .expect("write");
    writeln!(txt, "transfer_ok = {}", score.transfer_ok).expect("write");
    writeln!(txt, "h0_ok = {}", score.h0_ok).expect("write");
    writeln!(txt, "age_ok = {}", score.age_ok).expect("write");
    writeln!(txt, "recombination_ok = {}", score.recombination_ok).expect("write");
    writeln!(txt, "bbn_timing_ok = {}", score.bbn_timing_ok).expect("write");
    writeln!(
        txt,
        "passes_early_universe = {}",
        score.passes_early_universe()
    )
    .expect("write");
    writeln!(
        txt,
        "passes_late_universe = {}",
        score.passes_late_universe()
    )
    .expect("write");
    writeln!(txt, "passes_all = {}", score.passes_all()).expect("write");
    writeln!(txt).expect("write");
    writeln!(txt, "[epochs]").expect("write");
    for e in &score.epochs {
        writeln!(
            txt,
            "{}: z={:.6e}, age={:.6} yr, T={:.6e} K, H={:.6} km/s/Mpc, Ω_r={:.6}, Ω_m={:.6}, Ω_Λ={:.6}",
            e.name,
            e.z,
            e.age_seconds / SEC_PER_YEAR,
            e.temperature_k,
            e.h_km_s_mpc,
            e.omega_r,
            e.omega_m,
            e.omega_lambda,
        )
        .expect("write");
    }
    writeln!(txt).expect("write");
    writeln!(txt, "[history_samples]").expect("write");
    for (i, r) in score.history.iter().enumerate() {
        if i % 16 == 0 || i + 1 == score.history.len() {
            writeln!(
                txt,
                "idx={:03}, z={:.6e}, age={:.6} yr, T={:.6e} K, H={:.6} km/s/Mpc, Ω_r={:.6}, Ω_m={:.6}, Ω_Λ={:.6}",
                i,
                r.z,
                r.age_seconds / SEC_PER_YEAR,
                r.temperature_k,
                r.h_km_s_mpc,
                r.omega_r,
                r.omega_m,
                r.omega_lambda,
            )
            .expect("write");
        }
    }

    let mut json = File::create(&json_path).expect("create json");
    writeln!(
        json,
        "{{\n  \"assumptions\": {{\"omega_r0\": {:.12}, \"omega_k0\": {:.12}, \"h0_ref_km_s_mpc\": {:.9}}},\n  \"windows\": {{\"h0_rel_error_max\": {:.12}, \"age_gyr_min\": {:.9}, \"age_gyr_max\": {:.9}, \"recombination_age_kyr_min\": {:.9}, \"recombination_age_kyr_max\": {:.9}, \"bbn_age_sec_min\": {:.9}, \"bbn_age_sec_max\": {:.9}}},\n  \"score\": {{\"lambda_full\": {:.12e}, \"h0_km_s_mpc\": {:.9}, \"h0_rel_error\": {:.12}, \"omega_b0\": {:.12}, \"omega_dm0\": {:.12}, \"omega_m0\": {:.12}, \"omega_r0\": {:.12}, \"omega_k0\": {:.12}, \"omega_lambda0\": {:.12}, \"age_gyr\": {:.9}, \"recombination_age_kyr\": {:.9}, \"bbn_age_seconds\": {:.9}, \"r_s_drag_mpc\": {:.9}, \"theta_star_rad\": {:.9e}, \"l_peak1\": {:.9}, \"l_peak2\": {:.9}, \"micro_yp\": {:.9}, \"micro_dh\": {:.12e}, \"micro_z_visibility_peak\": {:.9}, \"inflation_ok\": {}, \"baryogenesis_ok\": {}, \"bbn_ok\": {}, \"microphysics_ok\": {}, \"dark_matter_unified_ok\": {}, \"transfer_ok\": {}, \"h0_ok\": {}, \"age_ok\": {}, \"recombination_ok\": {}, \"bbn_timing_ok\": {}, \"passes_early_universe\": {}, \"passes_late_universe\": {}, \"passes_all\": {}}},",
        assumptions.omega_r0,
        assumptions.omega_k0,
        assumptions.h0_ref_km_s_mpc,
        windows.h0_rel_error_max,
        windows.age_gyr_min,
        windows.age_gyr_max,
        windows.recombination_age_kyr_min,
        windows.recombination_age_kyr_max,
        windows.bbn_age_sec_min,
        windows.bbn_age_sec_max,
        score.lambda_full,
        score.h0_km_s_mpc,
        score.h0_rel_error,
        score.omega_b0,
        score.omega_dm0,
        score.omega_m0,
        score.omega_r0,
        score.omega_k0,
        score.omega_lambda0,
        score.age_gyr,
        score.recombination_age_kyr,
        score.bbn_age_seconds,
        score.transfer.rs_drag_mpc,
        score.transfer.theta_star_rad,
        score.transfer.l_peak1,
        score.transfer.l_peak2,
        score.microphysics.yp_network,
        score.microphysics.dh_network,
        score.microphysics.z_visibility_peak,
        score.inflation_ok,
        score.baryogenesis_ok,
        score.bbn_ok,
        score.microphysics_ok,
        score.dark_matter_unified_ok,
        score.transfer_ok,
        score.h0_ok,
        score.age_ok,
        score.recombination_ok,
        score.bbn_timing_ok,
        score.passes_early_universe(),
        score.passes_late_universe(),
        score.passes_all(),
    )
    .expect("write");

    writeln!(
        json,
        "  \"transfer\": {{\"h\": {:.9}, \"omega_b_h2\": {:.9}, \"omega_m_h2\": {:.9}, \"z_drag\": {:.6}, \"z_recomb\": {:.6}, \"rs_drag_mpc\": {:.9}, \"dm_recomb_mpc\": {:.9}, \"theta_star_rad\": {:.9e}, \"acoustic_scale_la\": {:.9}, \"l_peak1\": {:.9}, \"l_peak2\": {:.9}, \"growth_z0\": {:.9}, \"growth_z1\": {:.9}, \"pk_pivot_z0\": {:.12e}, \"pk_pivot_z1\": {:.12e}, \"rs_rel_error\": {:.9}, \"theta_star_rel_error\": {:.9}, \"l1_rel_error\": {:.9}, \"l2_rel_error\": {:.9}, \"rs_ok\": {}, \"theta_star_ok\": {}, \"l1_ok\": {}, \"l2_ok\": {}, \"transfer_positive_ok\": {}, \"passes_all\": {}}},",
        score.transfer.h,
        score.transfer.omega_b_h2,
        score.transfer.omega_m_h2,
        score.transfer.z_drag,
        score.transfer.z_recomb,
        score.transfer.rs_drag_mpc,
        score.transfer.dm_recomb_mpc,
        score.transfer.theta_star_rad,
        score.transfer.acoustic_scale_la,
        score.transfer.l_peak1,
        score.transfer.l_peak2,
        score.transfer.growth_z0,
        score.transfer.growth_z1,
        score.transfer.pk_pivot_z0,
        score.transfer.pk_pivot_z1,
        score.transfer.rs_rel_error,
        score.transfer.theta_star_rel_error,
        score.transfer.l1_rel_error,
        score.transfer.l2_rel_error,
        score.transfer.rs_ok,
        score.transfer.theta_star_ok,
        score.transfer.l1_ok,
        score.transfer.l2_ok,
        score.transfer.transfer_positive_ok,
        score.transfer.passes_all(),
    )
    .expect("write");
    writeln!(
        json,
        "  \"microphysics\": {{\"yp_network\": {:.9}, \"dh_network\": {:.12e}, \"he3h_network\": {:.12e}, \"bbn_freezeout_seconds\": {:.9}, \"z_visibility_peak\": {:.9}, \"tau_recomb\": {:.9e}, \"x_e_final\": {:.9e}, \"yp_ok\": {}, \"dh_ok\": {}, \"recombination_ok\": {}, \"opacity_positive_ok\": {}, \"passes_all\": {}}},",
        score.microphysics.yp_network,
        score.microphysics.dh_network,
        score.microphysics.he3h_network,
        score.microphysics.bbn_freezeout_seconds,
        score.microphysics.z_visibility_peak,
        score.microphysics.tau_recomb,
        score.microphysics.x_e_final,
        score.microphysics.yp_ok,
        score.microphysics.dh_ok,
        score.microphysics.recombination_ok,
        score.microphysics.opacity_positive_ok,
        score.microphysics.passes_all(),
    )
    .expect("write");
    writeln!(json, "  \"epochs\": [").expect("write");
    for (i, e) in score.epochs.iter().enumerate() {
        writeln!(
            json,
            "    {{\"name\":\"{}\",\"z\":{:.9e},\"age_seconds\":{:.9e},\"temperature_k\":{:.9e},\"h_km_s_mpc\":{:.9},\"omega_r\":{:.9},\"omega_m\":{:.9},\"omega_lambda\":{:.9}}}{}",
            e.name,
            e.z,
            e.age_seconds,
            e.temperature_k,
            e.h_km_s_mpc,
            e.omega_r,
            e.omega_m,
            e.omega_lambda,
            if i + 1 == score.epochs.len() { "" } else { "," }
        )
        .expect("write");
    }
    writeln!(json, "  ],").expect("write");

    writeln!(json, "  \"history\": [").expect("write");
    for (i, r) in score.history.iter().enumerate() {
        writeln!(
            json,
            "    {{\"z\":{:.9e},\"age_seconds\":{:.9e},\"temperature_k\":{:.9e},\"h_km_s_mpc\":{:.9},\"omega_r\":{:.9},\"omega_m\":{:.9},\"omega_lambda\":{:.9}}}{}",
            r.z,
            r.age_seconds,
            r.temperature_k,
            r.h_km_s_mpc,
            r.omega_r,
            r.omega_m,
            r.omega_lambda,
            if i + 1 == score.history.len() { "" } else { "," }
        )
        .expect("write");
    }
    writeln!(json, "  ]\n}}").expect("write");

    println!("wrote {txt_path}");
    println!("wrote {json_path}");
    println!(
        "Universe sim: H0={:.4} km/s/Mpc, age={:.4} Gyr, BBN={:.2} s, recombination={:.2} kyr, r_s={:.2} Mpc, l1={:.2}, pass={}",
        score.h0_km_s_mpc,
        score.age_gyr,
        score.bbn_age_seconds,
        score.recombination_age_kyr,
        score.transfer.rs_drag_mpc,
        score.transfer.l_peak1,
        score.passes_all()
    );
}
