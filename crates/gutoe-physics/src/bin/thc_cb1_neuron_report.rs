//! THC-CB1 neuronal simulation report.
//!
//! Emits:
//! - binding thermodynamics and decomposition
//! - non-electrostatic residual breakdown
//! - concentration sweep for receptor occupancy and neuronal suppression

use gutoe_physics::{
    decompose_thc_cb1_non_electrostatic_residual, evaluate_thc_cb1_binding,
    simulate_thc_cb1_neuron_response, NeuronCouplingInput, ThcCb1BindingInput,
    ThcElectrostaticProxyInput, ThcResidualProxyInput,
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

fn env_sweep() -> Vec<f64> {
    if let Ok(raw) = std::env::var("GUTOE_THC_SWEEP_NM") {
        let xs: Vec<f64> = raw
            .split(',')
            .filter_map(|s| s.trim().parse::<f64>().ok())
            .filter(|x| x.is_finite() && *x >= 0.0)
            .collect();
        if !xs.is_empty() {
            return xs;
        }
    }
    vec![0.0, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0]
}

fn nearest_point(points: &[gutoe_physics::NeuronResponsePoint], target_nm: f64) -> Option<gutoe_physics::NeuronResponsePoint> {
    points.iter().copied().min_by(|a, b| {
        let da = (a.concentration_nanomolar - target_nm).abs();
        let db = (b.concentration_nanomolar - target_nm).abs();
        da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
    })
}

fn main() {
    let binding = ThcCb1BindingInput {
        ki_nanomolar: env_f64("GUTOE_THC_KI_NM", ThcCb1BindingInput::default().ki_nanomolar),
        temperature_k: env_f64("GUTOE_THC_TEMP_K", ThcCb1BindingInput::default().temperature_k),
    };
    let electro = ThcElectrostaticProxyInput {
        hbond_contact_count: env_f64("GUTOE_THC_HBOND_CONTACTS", 1.4),
        hbond_charge_product: env_f64("GUTOE_THC_HBOND_CHARGE_PRODUCT", 0.16),
        hbond_distance_nm: env_f64("GUTOE_THC_HBOND_DISTANCE_NM", 0.30),
        hbond_dielectric: env_f64("GUTOE_THC_HBOND_DIELECTRIC", 28.0),
        polar_dipole_contact_count: env_f64("GUTOE_THC_POLAR_CONTACTS", 1.2),
        polar_dipole_charge_product: env_f64("GUTOE_THC_POLAR_CHARGE_PRODUCT", 0.10),
        polar_dipole_distance_nm: env_f64("GUTOE_THC_POLAR_DISTANCE_NM", 0.34),
        polar_dipole_dielectric: env_f64("GUTOE_THC_POLAR_DIELECTRIC", 32.0),
    };
    let residual_proxy = ThcResidualProxyInput {
        effective_hydrophobic_area_a2: env_f64("GUTOE_THC_RESID_HYDROPHOBIC_AREA_A2", 700.0),
        hydrophobic_coeff_kj_per_a2: env_f64("GUTOE_THC_RESID_HYDROPHOBIC_COEFF", 0.052),
        aromatic_contact_count: env_f64("GUTOE_THC_RESID_AROMATIC_CONTACTS", 3.0),
        aromatic_contact_stabilization_kj: env_f64("GUTOE_THC_RESID_AROMATIC_KJ", 1.55),
        released_water_count: env_f64("GUTOE_THC_RESID_RELEASED_WATERS", 2.8),
        water_release_stabilization_kj: env_f64("GUTOE_THC_RESID_WATER_KJ", 1.20),
        constrained_rotatable_bonds: env_f64("GUTOE_THC_RESID_CONSTRAINED_ROTORS", 5.5),
        conformational_entropy_penalty_per_rotor_kj: env_f64(
            "GUTOE_THC_RESID_ROTOR_PENALTY_KJ",
            0.60,
        ),
        polar_desolvated_contact_count: env_f64("GUTOE_THC_RESID_POLAR_DESOLV_CONTACTS", 1.4),
        polar_desolvation_penalty_kj: env_f64("GUTOE_THC_RESID_POLAR_DESOLV_KJ", 0.55),
        ligand_strain_penalty_kj: env_f64("GUTOE_THC_RESID_STRAIN_KJ", 0.75),
    };
    let coupling = NeuronCouplingInput {
        intrinsic_efficacy: env_f64("GUTOE_THC_EFFICACY", 0.55),
        max_release_inhibition_fraction: env_f64("GUTOE_THC_MAX_RELEASE_INHIBITION", 0.75),
        max_firing_suppression_fraction: env_f64("GUTOE_THC_MAX_FIRING_SUPPRESSION", 0.45),
        hill_coefficient: env_f64("GUTOE_THC_HILL", 1.0),
        baseline_release_probability: env_f64("GUTOE_THC_BASELINE_RELEASE_P", 0.35),
        baseline_firing_rate_hz: env_f64("GUTOE_THC_BASELINE_FIRING_HZ", 8.0),
    };

    let score = evaluate_thc_cb1_binding(binding, electro);
    let residual =
        decompose_thc_cb1_non_electrostatic_residual(score.residual_required_kj_mol, residual_proxy);
    let sweep = env_sweep();
    let points = simulate_thc_cb1_neuron_response(binding, coupling, &sweep);

    let p10 = nearest_point(&points, 10.0);
    let p30 = nearest_point(&points, 30.0);
    let p100 = nearest_point(&points, 100.0);

    let out_dir = std::env::var("GUTOE_THC_CB1_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/thc_cb1_neuron".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);
    let txt_path = out.join("thc_cb1_neuron_report.txt");
    let json_path = out.join("thc_cb1_neuron_report.json");
    let csv_path = out.join("thc_cb1_neuron_sweep.csv");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[thc_cb1_neuron]").expect("write");
    writeln!(txt, "ki_nanomolar = {:.9}", binding.ki_nanomolar).expect("write");
    writeln!(txt, "temperature_k = {:.6}", binding.temperature_k).expect("write");
    writeln!(txt, "experimental_delta_g_kj_mol = {:.9}", score.experimental_delta_g_kj_mol)
        .expect("write");
    writeln!(txt, "qed_floor_total_kj_mol = {:.9}", score.qed_floor_total_kj_mol).expect("write");
    writeln!(txt, "residual_required_kj_mol = {:.9}", score.residual_required_kj_mol)
        .expect("write");
    writeln!(
        txt,
        "explained_fraction_of_abs_delta_g = {:.9}",
        score.explained_fraction_of_abs_delta_g
    )
    .expect("write");
    writeln!(
        txt,
        "residual_modeled_total_kj_mol = {:.9}",
        residual.modeled_residual_total_kj_mol
    )
    .expect("write");
    writeln!(txt, "residual_closure_error_kj_mol = {:.9}", residual.closure_error_kj_mol)
        .expect("write");
    if let Some(p) = p10 {
        writeln!(
            txt,
            "at_10nM: occupancy={:.6}, release_p={:.6}, firing_hz={:.6}",
            p.occupancy_fraction, p.release_probability, p.firing_rate_hz
        )
        .expect("write");
    }
    if let Some(p) = p30 {
        writeln!(
            txt,
            "at_30nM: occupancy={:.6}, release_p={:.6}, firing_hz={:.6}",
            p.occupancy_fraction, p.release_probability, p.firing_rate_hz
        )
        .expect("write");
    }
    if let Some(p) = p100 {
        writeln!(
            txt,
            "at_100nM: occupancy={:.6}, release_p={:.6}, firing_hz={:.6}",
            p.occupancy_fraction, p.release_probability, p.firing_rate_hz
        )
        .expect("write");
    }

    let mut csv = String::from(
        "concentration_nM,occupancy_fraction,effective_activation_fraction,release_probability,firing_rate_hz,relative_firing_scale\n",
    );
    for p in &points {
        csv.push_str(&format!(
            "{:.9},{:.9},{:.9},{:.9},{:.9},{:.9}\n",
            p.concentration_nanomolar,
            p.occupancy_fraction,
            p.effective_activation_fraction,
            p.release_probability,
            p.firing_rate_hz,
            p.relative_firing_scale
        ));
    }
    fs::write(&csv_path, csv).expect("write csv");

    let payload = json!({
        "meta": {
            "lane": "thc_cb1_neuron_reduced_biophysics",
            "scope": [
                "thc_cb1_binding_thermodynamics",
                "qed_floor_vs_residual",
                "cb1_occupancy_to_neuronal_suppression_curve"
            ]
        },
        "binding_input": {
            "ki_nanomolar": binding.ki_nanomolar,
            "temperature_k": binding.temperature_k
        },
        "binding_decomposition_kj_mol": {
            "experimental_delta_g": score.experimental_delta_g_kj_mol,
            "qed_hbond_floor": score.qed_hbond_floor_kj_mol,
            "qed_polar_floor": score.qed_polar_floor_kj_mol,
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
        "neuron_coupling_input": {
            "intrinsic_efficacy": coupling.intrinsic_efficacy,
            "max_release_inhibition_fraction": coupling.max_release_inhibition_fraction,
            "max_firing_suppression_fraction": coupling.max_firing_suppression_fraction,
            "hill_coefficient": coupling.hill_coefficient,
            "baseline_release_probability": coupling.baseline_release_probability,
            "baseline_firing_rate_hz": coupling.baseline_firing_rate_hz
        },
        "sweep_points": points.iter().map(|p| json!({
            "concentration_nM": p.concentration_nanomolar,
            "occupancy_fraction": p.occupancy_fraction,
            "effective_activation_fraction": p.effective_activation_fraction,
            "release_probability": p.release_probability,
            "firing_rate_hz": p.firing_rate_hz,
            "relative_firing_scale": p.relative_firing_scale
        })).collect::<Vec<_>>(),
        "diagnostics": {
            "explained_fraction_of_abs_delta_g": score.explained_fraction_of_abs_delta_g,
            "residual_absolute_error_kj_mol": residual.closure_error_kj_mol.abs()
        }
    });

    fs::write(
        &json_path,
        serde_json::to_string_pretty(&payload).expect("serialize json"),
    )
    .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", csv_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "thc_cb1_neuron: ΔG_exp={:.3} kJ/mol, QED_floor={:.3}, residual={:.3}, explained={:.3}, residual_err={:.3}",
        score.experimental_delta_g_kj_mol,
        score.qed_floor_total_kj_mol,
        score.residual_required_kj_mol,
        score.explained_fraction_of_abs_delta_g,
        residual.closure_error_kj_mol.abs()
    );
}
