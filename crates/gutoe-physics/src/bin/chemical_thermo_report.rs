//! Periodic-table thermodynamics report from the chemical proxy lane.

use gutoe_physics::{
    phase_from_gibbs, predict_element_thermo, predict_element_thermo_coupled_with_diagnostics,
    scan_nuclear_chart, AVOGADRO, ChemicalFamily, MatterState, NucleusRecord, R_GAS_J_MOL_K,
    ScanConfig, P_REF_PA,
};
use serde_json::json;
use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

fn family_name(f: ChemicalFamily) -> &'static str {
    match f {
        ChemicalFamily::Alkali => "alkali",
        ChemicalFamily::AlkalineEarth => "alkaline_earth",
        ChemicalFamily::Transition => "transition",
        ChemicalFamily::PostTransition => "post_transition",
        ChemicalFamily::Metalloid => "metalloid",
        ChemicalFamily::Nonmetal => "nonmetal",
        ChemicalFamily::Halogen => "halogen",
        ChemicalFamily::NobleGas => "noble_gas",
        ChemicalFamily::Lanthanide => "lanthanide",
        ChemicalFamily::Actinide => "actinide",
    }
}

fn state_name(s: MatterState) -> &'static str {
    match s {
        MatterState::Solid => "solid",
        MatterState::Liquid => "liquid",
        MatterState::Gas => "gas",
    }
}

fn representative_isotopes(records: &[NucleusRecord], z_max: u16) -> BTreeMap<u16, u16> {
    let mut best: BTreeMap<u16, (&NucleusRecord, f64)> = BTreeMap::new();
    for r in records.iter().filter(|r| r.z >= 1 && r.z <= z_max) {
        let score = r.stability_score;
        let cur = best.get(&r.z).map(|(_, s)| *s).unwrap_or(f64::NEG_INFINITY);
        if score > cur {
            best.insert(r.z, (r, score));
        }
    }
    (1..=z_max)
        .map(|z| {
            let a = best
                .get(&z)
                .map(|(r, _)| r.a)
                .unwrap_or((2.5 * z as f64).round() as u16);
            (z, a)
        })
        .collect()
}

fn env_bool(name: &str, default: bool) -> bool {
    match std::env::var(name) {
        Ok(v) => match v.to_ascii_lowercase().as_str() {
            "1" | "true" | "yes" | "on" => true,
            "0" | "false" | "no" | "off" => false,
            _ => default,
        },
        Err(_) => default,
    }
}

fn main() {
    let out_dir = std::env::var("GUTOE_CHEM_THERMO_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/chemical_thermodynamics".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let z_max = std::env::var("GUTOE_CHEM_Z_MAX")
        .ok()
        .and_then(|v| v.parse::<u16>().ok())
        .unwrap_or(140);
    let use_coupled = env_bool("GUTOE_CHEM_THERMO_COUPLED", true);
    let scan_cfg = ScanConfig {
        z_min: 1,
        z_max,
        n_min: 1,
        n_max: 260,
        ..ScanConfig::default()
    };
    let records = scan_nuclear_chart(scan_cfg);
    let rep = representative_isotopes(&records, z_max);

    let mut csv = String::from(
        "Z,A,family,period,state_298k_1bar,state_273k_1bar_stp,state_1000k_1bar,state_298k_100bar,molar_mass_g_mol,atomic_mass_u,mass_per_atom_kg,atomic_radius_pm,molar_volume_cm3_mol,volume_per_atom_ang3,molar_volume_stp_l_mol_if_gas,density_g_cm3,cohesive_energy_ev,debye_temperature_k,latent_fusion_kj_mol,latent_vaporization_kj_mol,melting_temperature_k,boiling_temperature_k,vapor_pressure_298k_pa,bulk_modulus_gpa,thermal_expansion_1_per_k,cp_solid_j_mol_k,cp_liquid_j_mol_k,cp_gas_j_mol_k,coupled_mode,scf_iterations,scf_residual,valence_electrons_scf,scf_atomic_radius_pm,scf_ie_ev,scf_ea_ev,scf_chi_ev,scf_hardness_ev,coupled_radius_pm,cohesive_frontier_proxy_ev,delta_density_vs_proxy_g_cm3,delta_melting_vs_proxy_k,delta_boiling_vs_proxy_k\n",
    );

    let mut solids = 0_usize;
    let mut liquids = 0_usize;
    let mut gases = 0_usize;
    let mut solids_stp = 0_usize;
    let mut liquids_stp = 0_usize;
    let mut gases_stp = 0_usize;
    let mut solids_1000k = 0_usize;
    let mut liquids_1000k = 0_usize;
    let mut gases_1000k = 0_usize;
    let mut solids_100bar = 0_usize;
    let mut liquids_100bar = 0_usize;
    let mut gases_100bar = 0_usize;
    let mut rows_json = Vec::new();
    let mut abs_density_delta_sum = 0.0;
    let mut abs_melting_delta_sum = 0.0;
    let mut abs_boiling_delta_sum = 0.0;

    for z in 1_u16..=z_max {
        let a = rep[&z];
        let base = predict_element_thermo(z, a);
        let (p, scf_iterations, scf_residual, scf_valence_e, scf_atomic_radius_pm, scf_ie, scf_ea, scf_chi, scf_hardness, coupled_radius_pm, cohesive_frontier_proxy_ev) =
            if use_coupled {
                let (pred, diag) = predict_element_thermo_coupled_with_diagnostics(z, a);
                (
                    pred,
                    diag.scf_iterations,
                    diag.scf_residual,
                    diag.valence_electrons,
                    diag.scf_atomic_radius_pm,
                    diag.scf_ionization_energy_ev,
                    diag.scf_electron_affinity_ev,
                    diag.scf_electronegativity_mulliken_ev,
                    diag.scf_chemical_hardness_ev,
                    diag.coupled_radius_pm,
                    diag.cohesive_frontier_proxy_ev,
                )
            } else {
                (
                    base,
                    0,
                    f64::NAN,
                    0,
                    f64::NAN,
                    f64::NAN,
                    f64::NAN,
                    f64::NAN,
                    f64::NAN,
                    f64::NAN,
                    f64::NAN,
                )
            };
        let delta_density = p.density_g_cm3 - base.density_g_cm3;
        let delta_melting = p.melting_temperature_k - base.melting_temperature_k;
        let delta_boiling = p.boiling_temperature_k - base.boiling_temperature_k;
        abs_density_delta_sum += delta_density.abs();
        abs_melting_delta_sum += delta_melting.abs();
        abs_boiling_delta_sum += delta_boiling.abs();
        let state_stp = phase_from_gibbs(
            p.latent_fusion_kj_mol,
            p.latent_vaporization_kj_mol,
            p.melting_temperature_k,
            p.boiling_temperature_k,
            273.15,
            P_REF_PA,
        );
        let state_1000k = phase_from_gibbs(
            p.latent_fusion_kj_mol,
            p.latent_vaporization_kj_mol,
            p.melting_temperature_k,
            p.boiling_temperature_k,
            1000.0,
            P_REF_PA,
        );
        let state_100bar = phase_from_gibbs(
            p.latent_fusion_kj_mol,
            p.latent_vaporization_kj_mol,
            p.melting_temperature_k,
            p.boiling_temperature_k,
            298.15,
            1.0e7,
        );
        match p.ambient_state_298k {
            MatterState::Solid => solids += 1,
            MatterState::Liquid => liquids += 1,
            MatterState::Gas => gases += 1,
        }
        match state_stp {
            MatterState::Solid => solids_stp += 1,
            MatterState::Liquid => liquids_stp += 1,
            MatterState::Gas => gases_stp += 1,
        }
        match state_1000k {
            MatterState::Solid => solids_1000k += 1,
            MatterState::Liquid => liquids_1000k += 1,
            MatterState::Gas => gases_1000k += 1,
        }
        match state_100bar {
            MatterState::Solid => solids_100bar += 1,
            MatterState::Liquid => liquids_100bar += 1,
            MatterState::Gas => gases_100bar += 1,
        }
        let mass_per_atom_kg = (p.molar_mass_g_mol / 1000.0) / AVOGADRO;
        let volume_per_atom_ang3 = (p.molar_volume_cm3_mol / AVOGADRO) * 1.0e24;
        let molar_volume_stp_l_mol_if_gas = (R_GAS_J_MOL_K * 273.15 / P_REF_PA) * 1000.0;
        csv.push_str(&format!(
            "{},{},{},{},{},{},{},{},{:.6},{:.6},{:.9},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.9},{:.6},{:.6},{:.6},{},{},{:.9},{},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6}\n",
            z,
            a,
            family_name(p.family),
            p.period,
            state_name(p.ambient_state_298k),
            state_name(state_stp),
            state_name(state_1000k),
            state_name(state_100bar),
            p.molar_mass_g_mol,
            p.molar_mass_g_mol,
            mass_per_atom_kg,
            p.atomic_radius_pm,
            p.molar_volume_cm3_mol,
            volume_per_atom_ang3,
            molar_volume_stp_l_mol_if_gas,
            p.density_g_cm3,
            p.cohesive_energy_ev_per_atom,
            p.debye_temperature_k,
            p.latent_fusion_kj_mol,
            p.latent_vaporization_kj_mol,
            p.melting_temperature_k,
            p.boiling_temperature_k,
            p.vapor_pressure_pa_298k,
            p.bulk_modulus_gpa,
            p.thermal_expansion_1_per_k,
            p.cp_solid_j_mol_k,
            p.cp_liquid_j_mol_k,
            p.cp_gas_j_mol_k,
            use_coupled,
            scf_iterations,
            scf_residual,
            scf_valence_e,
            scf_atomic_radius_pm,
            scf_ie,
            scf_ea,
            scf_chi,
            scf_hardness,
            coupled_radius_pm,
            cohesive_frontier_proxy_ev,
            delta_density,
            delta_melting,
            delta_boiling
        ));
        rows_json.push(json!({
            "z": z,
            "a": a,
            "family": family_name(p.family),
            "period": p.period,
            "state_298k_1bar": state_name(p.ambient_state_298k),
            "state_273k_1bar_stp": state_name(state_stp),
            "state_1000k_1bar": state_name(state_1000k),
            "state_298k_100bar": state_name(state_100bar),
            "molar_mass_g_mol": p.molar_mass_g_mol,
            "atomic_mass_u": p.molar_mass_g_mol,
            "mass_per_atom_kg": mass_per_atom_kg,
            "atomic_radius_pm": p.atomic_radius_pm,
            "molar_volume_cm3_mol": p.molar_volume_cm3_mol,
            "volume_per_atom_ang3": volume_per_atom_ang3,
            "molar_volume_stp_l_mol_if_gas": molar_volume_stp_l_mol_if_gas,
            "density_g_cm3": p.density_g_cm3,
            "cohesive_energy_ev_per_atom": p.cohesive_energy_ev_per_atom,
            "debye_temperature_k": p.debye_temperature_k,
            "latent_fusion_kj_mol": p.latent_fusion_kj_mol,
            "latent_vaporization_kj_mol": p.latent_vaporization_kj_mol,
            "melting_temperature_k": p.melting_temperature_k,
            "boiling_temperature_k": p.boiling_temperature_k,
            "vapor_pressure_298k_pa": p.vapor_pressure_pa_298k,
            "bulk_modulus_gpa": p.bulk_modulus_gpa,
            "thermal_expansion_1_per_k": p.thermal_expansion_1_per_k,
            "cp_solid_j_mol_k": p.cp_solid_j_mol_k,
            "cp_liquid_j_mol_k": p.cp_liquid_j_mol_k,
            "cp_gas_j_mol_k": p.cp_gas_j_mol_k,
            "coupled_mode": use_coupled,
            "coupling_deltas_vs_proxy": {
                "density_g_cm3": delta_density,
                "melting_temperature_k": delta_melting,
                "boiling_temperature_k": delta_boiling
            },
            "scf_diagnostics": {
                "iterations": scf_iterations,
                "residual": scf_residual,
                "valence_electrons": scf_valence_e,
                "atomic_radius_pm": scf_atomic_radius_pm,
                "ionization_energy_ev": scf_ie,
                "electron_affinity_ev": scf_ea,
                "electronegativity_mulliken_ev": scf_chi,
                "chemical_hardness_ev": scf_hardness,
                "coupled_radius_pm": coupled_radius_pm,
                "cohesive_frontier_proxy_ev": cohesive_frontier_proxy_ev
            }
        }));
    }

    let superheavy_focus: Vec<_> = [112_u16, 114, 118, 120, 126, z_max]
        .iter()
        .copied()
        .filter(|&z| z >= 1 && z <= z_max)
        .collect::<std::collections::BTreeSet<_>>()
        .into_iter()
        .map(|z| {
            let a = rep.get(&z).copied().unwrap_or((2.5 * z as f64).round() as u16);
            let p = if use_coupled {
                predict_element_thermo_coupled_with_diagnostics(z, a).0
            } else {
                predict_element_thermo(z, a)
            };
            let state_stp = phase_from_gibbs(
                p.latent_fusion_kj_mol,
                p.latent_vaporization_kj_mol,
                p.melting_temperature_k,
                p.boiling_temperature_k,
                273.15,
                P_REF_PA,
            );
            json!({
                "z": z,
                "a": a,
                "state_298k": state_name(p.ambient_state_298k),
                "state_273k_stp": state_name(state_stp),
                "melting_temperature_k": p.melting_temperature_k,
                "boiling_temperature_k": p.boiling_temperature_k,
                "vapor_pressure_298k_pa": p.vapor_pressure_pa_298k,
                "density_g_cm3": p.density_g_cm3,
                "cohesive_energy_ev_per_atom": p.cohesive_energy_ev_per_atom,
                "bulk_modulus_gpa": p.bulk_modulus_gpa
            })
        })
        .collect();

    let report = json!({
        "meta": {
            "scan_bounds": {"z_min": 1, "z_max": z_max, "n_min": 1, "n_max": 260},
            "lane": if use_coupled { "chemical_thermodynamics_coupled_abinitio" } else { "chemical_thermodynamics_proxy" },
            "coupled_mode": use_coupled,
            "assumptions": [
                "hydrogenic radius and cohesive-energy proxies with family-dependent factors",
                "latent heats from cohesive energy via fixed family fractions",
                "Trouton-like and Richards-like entropy rules for boiling/melting transduction",
                "optional coupling to atomic SCF frontier descriptors (IE/EA/chi/hardness/radius)"
            ]
        },
        "summary": {
            "elements_modeled": z_max,
            "mean_abs_delta_vs_proxy": {
                "density_g_cm3": abs_density_delta_sum / z_max as f64,
                "melting_temperature_k": abs_melting_delta_sum / z_max as f64,
                "boiling_temperature_k": abs_boiling_delta_sum / z_max as f64
            },
            "ambient_state_counts_298k": {
                "solid": solids,
                "liquid": liquids,
                "gas": gases
            },
            "state_counts_273k_1bar_stp": {
                "solid": solids_stp,
                "liquid": liquids_stp,
                "gas": gases_stp
            },
            "state_counts_1000k_1bar": {
                "solid": solids_1000k,
                "liquid": liquids_1000k,
                "gas": gases_1000k
            },
            "state_counts_298k_100bar": {
                "solid": solids_100bar,
                "liquid": liquids_100bar,
                "gas": gases_100bar
            }
        },
        "superheavy_extrapolation": superheavy_focus,
        "elements": rows_json
    });

    let txt_path = out.join("chemical_thermo_report.txt");
    let csv_path = out.join("chemical_thermo_report.csv");
    let json_path = out.join("chemical_thermo_report.json");

    fs::write(&csv_path, csv).expect("write csv");
    fs::write(
        &json_path,
        serde_json::to_string_pretty(&report).expect("serialize report"),
    )
    .expect("write json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[chemical_thermodynamics]").expect("write");
    writeln!(txt, "elements_modeled = {}", z_max).expect("write");
    writeln!(txt, "ambient_solid_count = {}", solids).expect("write");
    writeln!(txt, "ambient_liquid_count = {}", liquids).expect("write");
    writeln!(txt, "ambient_gas_count = {}", gases).expect("write");
    writeln!(txt, "stp_solid_count = {}", solids_stp).expect("write");
    writeln!(txt, "stp_liquid_count = {}", liquids_stp).expect("write");
    writeln!(txt, "stp_gas_count = {}", gases_stp).expect("write");
    writeln!(txt, "state_1000k_1bar_counts = ({},{},{})", solids_1000k, liquids_1000k, gases_1000k)
        .expect("write");
    writeln!(txt, "state_298k_100bar_counts = ({},{},{})", solids_100bar, liquids_100bar, gases_100bar)
        .expect("write");
    writeln!(txt, "coupled_mode = {}", use_coupled).expect("write");
    writeln!(
        txt,
        "mean_abs_delta_vs_proxy(density,melting,boiling)=({:.9},{:.9},{:.9})",
        abs_density_delta_sum / z_max as f64,
        abs_melting_delta_sum / z_max as f64,
        abs_boiling_delta_sum / z_max as f64
    )
    .expect("write");
    writeln!(txt, "superheavy_focus = dynamic").expect("write");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", csv_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "chemical_thermo: elements={}, coupled={}, state_298K@1bar(s/l/g)=({}/{}/{}), state_273K_STP@1bar=({}/{}/{}), state_1000K@1bar=({}/{}/{}), state_298K@100bar=({}/{}/{})",
        z_max, use_coupled, solids, liquids, gases, solids_stp, liquids_stp, gases_stp, solids_1000k, liquids_1000k, gases_1000k, solids_100bar, liquids_100bar, gases_100bar
    );
}
