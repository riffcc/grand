//! MS molecular-mimicry report.
//!
//! Outputs a reduced mechanistic scorecard for:
//! - TCR self-epitope vs mimic epitope binding energies
//! - misrecognition risk index
//! - therapy-effect proxy (ocrelizumab-like + natalizumab-like)
//! - targeted blocker feasibility

use gutoe_physics::{
    default_ms_mimicry_input, default_natalizumab_proxy, default_ocrelizumab_proxy,
    default_targeted_blocker_input, evaluate_interface_energy, evaluate_molecular_mimicry,
    evaluate_targeted_blocker, evaluate_therapy_effect,
};
use serde_json::json;
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

fn main() {
    let mimicry_input = default_ms_mimicry_input();
    let self_iface =
        evaluate_interface_energy(mimicry_input.self_epitope_electro, mimicry_input.self_epitope_residual);
    let mimic_iface = evaluate_interface_energy(
        mimicry_input.mimic_epitope_electro,
        mimicry_input.mimic_epitope_residual,
    );

    let mimicry = evaluate_molecular_mimicry(mimicry_input);
    let therapy = evaluate_therapy_effect(
        mimicry.misrecognition_risk_index,
        default_ocrelizumab_proxy(),
        default_natalizumab_proxy(),
    );
    let blocker =
        evaluate_targeted_blocker(mimicry.activation_excess_kj_mol, default_targeted_blocker_input());

    let out_dir = std::env::var("GUTOE_MS_AUTOIMMUNE_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ms_autoimmunity".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);
    let txt_path = out.join("ms_autoimmunity_report.txt");
    let json_path = out.join("ms_autoimmunity_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[ms_autoimmunity_mimicry]").expect("write");
    writeln!(txt, "self_binding_kj_mol = {:.9}", mimicry.self_binding_kj_mol).expect("write");
    writeln!(txt, "mimic_binding_kj_mol = {:.9}", mimicry.mimic_binding_kj_mol).expect("write");
    writeln!(txt, "mimicry_gap_kj_mol = {:.9}", mimicry.mimicry_gap_kj_mol).expect("write");
    writeln!(
        txt,
        "activation_excess_kj_mol = {:.9}",
        mimicry.activation_excess_kj_mol
    )
    .expect("write");
    writeln!(
        txt,
        "misrecognition_risk_index = {:.9}",
        mimicry.misrecognition_risk_index
    )
    .expect("write");
    writeln!(
        txt,
        "combined_therapy_reduction_fraction = {:.9}",
        therapy.relative_drive_reduction_fraction
    )
    .expect("write");
    writeln!(txt, "residual_drive_index = {:.9}", therapy.residual_drive_index).expect("write");
    writeln!(txt, "blocker_required_occupancy = {:.9}", blocker.required_occupancy_fraction)
        .expect("write");
    writeln!(txt, "blocker_feasible = {}", blocker.feasible_at_given_concentration).expect("write");

    let payload = json!({
        "meta": {
            "lane": "ms_molecular_mimicry_reduced",
            "note": "mechanistic simulation scaffold, not clinical recommendation"
        },
        "self_interface": {
            "qed_hbond_kj_mol": self_iface.qed_hbond_kj_mol,
            "qed_polar_kj_mol": self_iface.qed_polar_kj_mol,
            "qed_ionic_kj_mol": self_iface.qed_ionic_kj_mol,
            "qed_total_kj_mol": self_iface.qed_total_kj_mol,
            "residual_total_kj_mol": self_iface.residual_total_kj_mol,
            "total_delta_g_kj_mol": self_iface.total_delta_g_kj_mol
        },
        "mimic_interface": {
            "qed_hbond_kj_mol": mimic_iface.qed_hbond_kj_mol,
            "qed_polar_kj_mol": mimic_iface.qed_polar_kj_mol,
            "qed_ionic_kj_mol": mimic_iface.qed_ionic_kj_mol,
            "qed_total_kj_mol": mimic_iface.qed_total_kj_mol,
            "residual_total_kj_mol": mimic_iface.residual_total_kj_mol,
            "total_delta_g_kj_mol": mimic_iface.total_delta_g_kj_mol
        },
        "mimicry_score": {
            "tolerance_threshold_kj_mol": mimicry_input.tolerance_threshold_kj_mol,
            "self_binding_kj_mol": mimicry.self_binding_kj_mol,
            "mimic_binding_kj_mol": mimicry.mimic_binding_kj_mol,
            "mimicry_gap_kj_mol": mimicry.mimicry_gap_kj_mol,
            "activation_excess_kj_mol": mimicry.activation_excess_kj_mol,
            "overlap_score": mimicry.overlap_score,
            "activation_score": mimicry.activation_score,
            "misrecognition_risk_index": mimicry.misrecognition_risk_index
        },
        "therapy": {
            "baseline_drive_index": therapy.baseline_drive_index,
            "ocrelizumab_like": {
                "occupancy_fraction": therapy.ocrelizumab.occupancy_fraction,
                "effective_drive_reduction_fraction": therapy.ocrelizumab.effective_drive_reduction_fraction
            },
            "natalizumab_like": {
                "occupancy_fraction": therapy.natalizumab.occupancy_fraction,
                "effective_drive_reduction_fraction": therapy.natalizumab.effective_drive_reduction_fraction
            },
            "residual_drive_index": therapy.residual_drive_index,
            "relative_drive_reduction_fraction": therapy.relative_drive_reduction_fraction
        },
        "targeted_blocker": {
            "occupancy_fraction": blocker.occupancy_fraction,
            "achieved_energy_shift_kj_mol": blocker.achieved_energy_shift_kj_mol,
            "required_energy_shift_kj_mol": blocker.required_energy_shift_kj_mol,
            "required_occupancy_fraction": blocker.required_occupancy_fraction,
            "feasible_at_given_concentration": blocker.feasible_at_given_concentration
        }
    });

    fs::write(
        &json_path,
        serde_json::to_string_pretty(&payload).expect("serialize"),
    )
    .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "ms_autoimmunity: risk={:.3}, mimicry_gap={:.3} kJ/mol, activation_excess={:.3} kJ/mol, therapy_reduction={:.3}, blocker_feasible={}",
        mimicry.misrecognition_risk_index,
        mimicry.mimicry_gap_kj_mol,
        mimicry.activation_excess_kj_mol,
        therapy.relative_drive_reduction_fraction,
        blocker.feasible_at_given_concentration
    );
}
