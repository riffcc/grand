//! Molecular ab-initio report (RHF/UHF + MP2 where applicable).

use gutoe_physics::{benchmark_molecules, optimize_molecule_geometry, run_molecular_ab_initio};
use serde_json::json;
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

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

fn env_usize(name: &str, default: usize) -> usize {
    std::env::var(name)
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(default)
}

fn env_f64(name: &str, default: f64) -> f64 {
    std::env::var(name)
        .ok()
        .and_then(|v| v.parse::<f64>().ok())
        .unwrap_or(default)
}

fn main() {
    let out_dir = std::env::var("GUTOE_MOL_ABINITIO_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/molecular_ab_initio".to_string());
    let run_opt = env_bool("GUTOE_MOL_ABINITIO_OPTIMIZE", false);
    let opt_iter = env_usize("GUTOE_MOL_ABINITIO_OPT_MAX_ITER", 10);
    let opt_step = env_f64("GUTOE_MOL_ABINITIO_OPT_STEP", 0.08);

    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let mut csv = String::from(
        "molecule,multiplicity,method,basis_functions,electron_count,alpha_electrons,beta_electrons,electron_pairs,scf_iterations,scf_residual,s2_expectation,nuclear_repulsion_hartree,electronic_energy_hartree,total_energy_hartree,mp2_correlation_hartree,total_energy_mp2_hartree,homo_energy_ev,lumo_energy_ev,homo_lumo_gap_ev,dipole_debye,opt_performed,opt_converged,opt_iterations,opt_grad_norm_hartree_per_angstrom,opt_total_energy_hartree,mulliken_charges,orbital_energies_ev\n",
    );

    let mut rows = Vec::new();
    let mut total_ok = 0usize;
    let mut total_opt = 0usize;
    let mut total_opt_conv = 0usize;

    for m in benchmark_molecules() {
        match run_molecular_ab_initio(m.clone()) {
            Ok(r) => {
                total_ok += 1;

                let mut opt_performed = false;
                let mut opt_converged = false;
                let mut opt_iterations = 0usize;
                let mut opt_grad_norm = f64::NAN;
                let mut opt_total_energy = f64::NAN;

                if run_opt {
                    opt_performed = true;
                    total_opt += 1;
                    if let Ok(opt) = optimize_molecule_geometry(m.clone(), opt_iter, opt_step) {
                        if opt.converged {
                            total_opt_conv += 1;
                        }
                        opt_converged = opt.converged;
                        opt_iterations = opt.iterations;
                        opt_grad_norm = opt.final_gradient_norm_hartree_per_angstrom;
                        opt_total_energy = opt.final_result.total_energy_hartree;
                    }
                }

                let charges = r
                    .mulliken_charges
                    .iter()
                    .map(|v| format!("{:.6}", v))
                    .collect::<Vec<_>>()
                    .join(";");
                let orbitals = r
                    .orbital_energies_ev
                    .iter()
                    .map(|v| format!("{:.6}", v))
                    .collect::<Vec<_>>()
                    .join(";");

                csv.push_str(&format!(
                    "{},{},{},{},{},{},{},{},{},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{},{},{},{:.9},{:.9},\"{}\",\"{}\"\n",
                    r.name,
                    r.spin_multiplicity,
                    r.method,
                    r.basis_functions,
                    r.electron_count,
                    r.alpha_electrons,
                    r.beta_electrons,
                    r.electron_pairs,
                    r.scf_iterations,
                    r.scf_residual,
                    r.s2_expectation,
                    r.nuclear_repulsion_hartree,
                    r.electronic_energy_hartree,
                    r.total_energy_hartree,
                    r.mp2_correlation_hartree,
                    r.total_energy_mp2_hartree,
                    r.homo_energy_ev,
                    r.lumo_energy_ev,
                    r.homo_lumo_gap_ev,
                    r.dipole_debye,
                    opt_performed,
                    opt_converged,
                    opt_iterations,
                    opt_grad_norm,
                    opt_total_energy,
                    charges,
                    orbitals
                ));

                rows.push(json!({
                    "molecule": r.name,
                    "multiplicity": r.spin_multiplicity,
                    "method": r.method,
                    "basis_functions": r.basis_functions,
                    "electron_count": r.electron_count,
                    "alpha_electrons": r.alpha_electrons,
                    "beta_electrons": r.beta_electrons,
                    "electron_pairs": r.electron_pairs,
                    "scf_iterations": r.scf_iterations,
                    "scf_residual": r.scf_residual,
                    "s2_expectation": r.s2_expectation,
                    "nuclear_repulsion_hartree": r.nuclear_repulsion_hartree,
                    "electronic_energy_hartree": r.electronic_energy_hartree,
                    "total_energy_hartree": r.total_energy_hartree,
                    "mp2_correlation_hartree": r.mp2_correlation_hartree,
                    "total_energy_mp2_hartree": r.total_energy_mp2_hartree,
                    "homo_energy_ev": r.homo_energy_ev,
                    "lumo_energy_ev": r.lumo_energy_ev,
                    "homo_lumo_gap_ev": r.homo_lumo_gap_ev,
                    "dipole_debye": r.dipole_debye,
                    "optimization": {
                        "performed": opt_performed,
                        "converged": opt_converged,
                        "iterations": opt_iterations,
                        "gradient_norm_hartree_per_angstrom": opt_grad_norm,
                        "total_energy_hartree": opt_total_energy
                    },
                    "mulliken_charges": r.mulliken_charges,
                    "orbital_energies_ev": r.orbital_energies_ev,
                }));
            }
            Err(e) => {
                rows.push(json!({
                    "molecule": m.name,
                    "error": e.to_string(),
                }));
            }
        }
    }

    let report = json!({
        "meta": {
            "lane": "molecular_ab_initio_rhf_uhf_mp2",
            "note": "compact Gaussian AO basis (s-type primitives); RHF for singlet closed-shell, UHF for open-shell, MP2 applied on RHF branch",
            "optimization": {
                "enabled": run_opt,
                "max_iter": opt_iter,
                "step": opt_step,
            }
        },
        "summary": {
            "benchmarks_attempted": benchmark_molecules().len(),
            "benchmarks_solved": total_ok,
            "optimizations_attempted": total_opt,
            "optimizations_converged": total_opt_conv,
        },
        "molecules": rows
    });

    let txt_path = out.join("molecular_ab_initio_report.txt");
    let csv_path = out.join("molecular_ab_initio_report.csv");
    let json_path = out.join("molecular_ab_initio_report.json");

    fs::write(&csv_path, csv).expect("write csv");
    fs::write(&json_path, serde_json::to_string_pretty(&report).expect("serialize report"))
        .expect("write json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[molecular_ab_initio]").expect("write");
    writeln!(txt, "benchmarks_attempted = {}", benchmark_molecules().len()).expect("write");
    writeln!(txt, "benchmarks_solved = {}", total_ok).expect("write");
    writeln!(txt, "optimization_enabled = {}", run_opt).expect("write");
    writeln!(txt, "optimizations_attempted = {}", total_opt).expect("write");
    writeln!(txt, "optimizations_converged = {}", total_opt_conv).expect("write");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", csv_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "molecular_ab_initio: solved {}/{} benchmark molecules (opt enabled: {})",
        total_ok,
        benchmark_molecules().len(),
        run_opt
    );
}
