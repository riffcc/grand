//! Cardiovascular drug-target binding thermodynamics report.
//!
//! Focus pair:
//! - Drug: atorvastatin
//! - Target: HMG-CoA reductase (HMGCR)

use gutoe_physics::{
    decompose_non_electrostatic_residual, evaluate_atorvastatin_hmgcr_binding,
    BindingBenchmarkInput, ElectrostaticProxyInput, ResidualProxyInput,
};
use serde_json::json;
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

fn env_f64(key: &str, default: f64) -> f64 {
    std::env::var(key)
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(default)
}

fn main() {
    let benchmark = BindingBenchmarkInput {
        ki_nanomolar: env_f64("GUTOE_CARDIO_KI_NM", BindingBenchmarkInput::default().ki_nanomolar),
        temperature_k: env_f64(
            "GUTOE_CARDIO_TEMP_K",
            BindingBenchmarkInput::default().temperature_k,
        ),
    };
    let proxy = ElectrostaticProxyInput {
        ionic_contact_count: env_f64("GUTOE_CARDIO_IONIC_CONTACTS", 1.0),
        ionic_distance_nm: env_f64("GUTOE_CARDIO_IONIC_DISTANCE_NM", 0.30),
        ionic_dielectric: env_f64("GUTOE_CARDIO_IONIC_DIELECTRIC", 28.0),
        hbond_contact_count: env_f64("GUTOE_CARDIO_HBOND_CONTACTS", 5.0),
        hbond_charge_product: env_f64("GUTOE_CARDIO_HBOND_CHARGE_PRODUCT", 0.20),
        hbond_distance_nm: env_f64("GUTOE_CARDIO_HBOND_DISTANCE_NM", 0.29),
        hbond_dielectric: env_f64("GUTOE_CARDIO_HBOND_DIELECTRIC", 24.0),
    };
    let score = evaluate_atorvastatin_hmgcr_binding(benchmark, proxy);
    let residual_proxy = ResidualProxyInput {
        effective_hydrophobic_area_a2: env_f64("GUTOE_CARDIO_RESID_HYDROPHOBIC_AREA_A2", 225.0),
        hydrophobic_coeff_kj_per_a2: env_f64("GUTOE_CARDIO_RESID_HYDROPHOBIC_COEFF", 0.046),
        aromatic_contact_count: env_f64("GUTOE_CARDIO_RESID_AROMATIC_CONTACTS", 2.5),
        aromatic_contact_stabilization_kj: env_f64("GUTOE_CARDIO_RESID_AROMATIC_KJ", 1.10),
        released_water_count: env_f64("GUTOE_CARDIO_RESID_RELEASED_WATERS", 3.0),
        water_release_stabilization_kj: env_f64("GUTOE_CARDIO_RESID_WATER_KJ", 1.05),
        constrained_rotatable_bonds: env_f64("GUTOE_CARDIO_RESID_CONSTRAINED_ROTORS", 6.0),
        conformational_entropy_penalty_per_rotor_kj: env_f64(
            "GUTOE_CARDIO_RESID_ROTOR_PENALTY_KJ",
            0.78,
        ),
        polar_desolvated_contact_count: env_f64("GUTOE_CARDIO_RESID_POLAR_DESOLV_CONTACTS", 3.0),
        polar_desolvation_penalty_kj: env_f64("GUTOE_CARDIO_RESID_POLAR_DESOLV_KJ", 0.45),
        ligand_strain_penalty_kj: env_f64("GUTOE_CARDIO_RESID_STRAIN_KJ", 0.60),
    };
    let residual = decompose_non_electrostatic_residual(score.residual_required_kj_mol, residual_proxy);

    let out_dir = std::env::var("GUTOE_CARDIO_BINDING_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/cardiovascular_binding".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let txt_path = out.join("cardiovascular_binding_report.txt");
    let json_path = out.join("cardiovascular_binding_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[cardiovascular_binding]").expect("write");
    writeln!(txt, "pair = atorvastatin_hmgcr").expect("write");
    writeln!(txt, "ki_nanomolar = {:.9}", benchmark.ki_nanomolar).expect("write");
    writeln!(txt, "temperature_k = {:.6}", benchmark.temperature_k).expect("write");
    writeln!(txt, "experimental_delta_g_kj_mol = {:.9}", score.experimental_delta_g_kj_mol)
        .expect("write");
    writeln!(txt, "qed_ionic_floor_kj_mol = {:.9}", score.qed_ionic_floor_kj_mol).expect("write");
    writeln!(txt, "qed_hbond_floor_kj_mol = {:.9}", score.qed_hbond_floor_kj_mol).expect("write");
    writeln!(txt, "qed_floor_total_kj_mol = {:.9}", score.qed_floor_total_kj_mol).expect("write");
    writeln!(
        txt,
        "residual_required_kj_mol = {:.9}",
        score.residual_required_kj_mol
    )
    .expect("write");
    writeln!(
        txt,
        "explained_fraction_of_abs_delta_g = {:.9}",
        score.explained_fraction_of_abs_delta_g
    )
    .expect("write");
    writeln!(
        txt,
        "residual_hydrophobic_kj_mol = {:.9}",
        residual.hydrophobic_stabilization_kj_mol
    )
    .expect("write");
    writeln!(
        txt,
        "residual_aromatic_packing_kj_mol = {:.9}",
        residual.aromatic_packing_stabilization_kj_mol
    )
    .expect("write");
    writeln!(
        txt,
        "residual_water_release_kj_mol = {:.9}",
        residual.water_release_stabilization_kj_mol
    )
    .expect("write");
    writeln!(
        txt,
        "residual_conformational_entropy_penalty_kj_mol = {:.9}",
        residual.conformational_entropy_penalty_kj_mol
    )
    .expect("write");
    writeln!(
        txt,
        "residual_polar_desolvation_penalty_kj_mol = {:.9}",
        residual.polar_desolvation_penalty_kj_mol
    )
    .expect("write");
    writeln!(
        txt,
        "residual_ligand_strain_penalty_kj_mol = {:.9}",
        residual.ligand_strain_penalty_kj_mol
    )
    .expect("write");
    writeln!(
        txt,
        "residual_modeled_total_kj_mol = {:.9}",
        residual.modeled_residual_total_kj_mol
    )
    .expect("write");
    writeln!(
        txt,
        "residual_target_kj_mol = {:.9}",
        residual.target_residual_kj_mol
    )
    .expect("write");
    writeln!(txt, "residual_closure_error_kj_mol = {:.9}", residual.closure_error_kj_mol)
        .expect("write");

    let payload = json!({
        "meta": {
            "lane": "cardiovascular_binding_transduction",
            "pair": "atorvastatin_hmgcr",
            "assumptions": [
                "experimental benchmark is Ki->ΔG standard-state conversion",
                "QED floor uses α ħ c /(ε r) contact energetics",
                "residual term explicitly represents non-electrostatic and unresolved structural contributions"
            ]
        },
        "benchmark": {
            "ki_nanomolar": benchmark.ki_nanomolar,
            "temperature_k": benchmark.temperature_k
        },
        "electrostatic_proxy": {
            "ionic_contact_count": proxy.ionic_contact_count,
            "ionic_distance_nm": proxy.ionic_distance_nm,
            "ionic_dielectric": proxy.ionic_dielectric,
            "hbond_contact_count": proxy.hbond_contact_count,
            "hbond_charge_product": proxy.hbond_charge_product,
            "hbond_distance_nm": proxy.hbond_distance_nm,
            "hbond_dielectric": proxy.hbond_dielectric
        },
        "residual_proxy": {
            "effective_hydrophobic_area_a2": residual_proxy.effective_hydrophobic_area_a2,
            "hydrophobic_coeff_kj_per_a2": residual_proxy.hydrophobic_coeff_kj_per_a2,
            "aromatic_contact_count": residual_proxy.aromatic_contact_count,
            "aromatic_contact_stabilization_kj": residual_proxy.aromatic_contact_stabilization_kj,
            "released_water_count": residual_proxy.released_water_count,
            "water_release_stabilization_kj": residual_proxy.water_release_stabilization_kj,
            "constrained_rotatable_bonds": residual_proxy.constrained_rotatable_bonds,
            "conformational_entropy_penalty_per_rotor_kj": residual_proxy.conformational_entropy_penalty_per_rotor_kj,
            "polar_desolvated_contact_count": residual_proxy.polar_desolvated_contact_count,
            "polar_desolvation_penalty_kj": residual_proxy.polar_desolvation_penalty_kj,
            "ligand_strain_penalty_kj": residual_proxy.ligand_strain_penalty_kj
        },
        "decomposition_kj_mol": {
            "experimental_delta_g": score.experimental_delta_g_kj_mol,
            "qed_ionic_floor": score.qed_ionic_floor_kj_mol,
            "qed_hbond_floor": score.qed_hbond_floor_kj_mol,
            "qed_floor_total": score.qed_floor_total_kj_mol,
            "residual_required": score.residual_required_kj_mol
        },
        "residual_breakdown_kj_mol": {
            "hydrophobic_stabilization": residual.hydrophobic_stabilization_kj_mol,
            "aromatic_packing_stabilization": residual.aromatic_packing_stabilization_kj_mol,
            "water_release_stabilization": residual.water_release_stabilization_kj_mol,
            "conformational_entropy_penalty": residual.conformational_entropy_penalty_kj_mol,
            "polar_desolvation_penalty": residual.polar_desolvation_penalty_kj_mol,
            "ligand_strain_penalty": residual.ligand_strain_penalty_kj_mol,
            "modeled_residual_total": residual.modeled_residual_total_kj_mol,
            "target_residual": residual.target_residual_kj_mol,
            "closure_error": residual.closure_error_kj_mol
        },
        "diagnostics": {
            "explained_fraction_of_abs_delta_g": score.explained_fraction_of_abs_delta_g,
            "residual_absolute_error_kj_mol": residual.closure_error_kj_mol.abs()
        }
    });

    fs::write(
        &json_path,
        serde_json::to_string_pretty(&payload).expect("serialize report"),
    )
    .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "cardiovascular_binding: ΔG_exp={:.3} kJ/mol, QED_floor={:.3}, residual_target={:.3}, residual_modeled={:.3}, resid_err={:.3}, explained={:.3}",
        score.experimental_delta_g_kj_mol,
        score.qed_floor_total_kj_mol,
        score.residual_required_kj_mol,
        residual.modeled_residual_total_kj_mol,
        residual.closure_error_kj_mol,
        score.explained_fraction_of_abs_delta_g
    );
}
