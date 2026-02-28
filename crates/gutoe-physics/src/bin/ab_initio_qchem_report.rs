//! Ab-initio atomic quantum-chemistry report (SCF lane).

use gutoe_physics::{predict_atomic_scf, scan_nuclear_chart, NucleusRecord, ScanConfig};
use serde_json::json;
use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

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

fn main() {
    let out_dir = std::env::var("GUTOE_ABINITIO_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ab_initio_qchem".to_string());
    let z_max = std::env::var("GUTOE_ABINITIO_Z_MAX")
        .ok()
        .and_then(|v| v.parse::<u16>().ok())
        .unwrap_or(140);

    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

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
        "Z,A,electron_count,valence_electrons,scf_iterations,scf_residual,total_electronic_energy_ev,homo_energy_ev,lumo_energy_ev,ionization_energy_ev,electron_affinity_ev,electronegativity_mulliken_ev,chemical_hardness_ev,chemical_softness_inv_ev,atomic_radius_pm,covalent_radius_pm,polarizability_a0_cubed,electron_configuration,frontier_orbitals\n",
    );
    let mut orbital_csv =
        String::from("Z,A,orbital_n,orbital_l,label,occupation,zeff,energy_ev,mean_radius_pm\n");

    let mut rows_json = Vec::new();
    let mut ie_sum = 0.0;
    let mut chi_sum = 0.0;
    let mut radius_sum = 0.0;

    for z in 1_u16..=z_max {
        let a = rep[&z];
        let p = predict_atomic_scf(z, a);
        ie_sum += p.ionization_energy_ev;
        chi_sum += p.electronegativity_mulliken_ev;
        radius_sum += p.atomic_radius_pm;

        csv.push_str(&format!(
            "{},{},{},{},{},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},\"{}\",\"{}\"\n",
            p.z,
            p.a,
            p.electron_count,
            p.valence_electrons,
            p.scf_iterations,
            p.scf_residual,
            p.total_electronic_energy_ev,
            p.homo_energy_ev,
            p.lumo_energy_ev,
            p.ionization_energy_ev,
            p.electron_affinity_ev,
            p.electronegativity_mulliken_ev,
            p.chemical_hardness_ev,
            p.chemical_softness_inv_ev,
            p.atomic_radius_pm,
            p.covalent_radius_pm,
            p.polarizability_a0_cubed,
            p.electron_configuration,
            p.frontier_orbitals,
        ));

        for o in &p.orbitals {
            if o.occupation == 0 {
                continue;
            }
            orbital_csv.push_str(&format!(
                "{},{},{},{},{},{},{:.9},{:.9},{:.9}\n",
                p.z,
                p.a,
                o.n,
                o.l,
                format!("{}{}", o.n, match o.l {
                    0 => 's',
                    1 => 'p',
                    2 => 'd',
                    3 => 'f',
                    4 => 'g',
                    5 => 'h',
                    6 => 'i',
                    7 => 'k',
                    _ => '?',
                }),
                o.occupation,
                o.zeff,
                o.energy_ev,
                o.mean_radius_pm
            ));
        }

        rows_json.push(json!({
            "z": p.z,
            "a": p.a,
            "electron_count": p.electron_count,
            "valence_electrons": p.valence_electrons,
            "scf_iterations": p.scf_iterations,
            "scf_residual": p.scf_residual,
            "total_electronic_energy_ev": p.total_electronic_energy_ev,
            "homo_energy_ev": p.homo_energy_ev,
            "lumo_energy_ev": p.lumo_energy_ev,
            "ionization_energy_ev": p.ionization_energy_ev,
            "electron_affinity_ev": p.electron_affinity_ev,
            "electronegativity_mulliken_ev": p.electronegativity_mulliken_ev,
            "chemical_hardness_ev": p.chemical_hardness_ev,
            "chemical_softness_inv_ev": p.chemical_softness_inv_ev,
            "atomic_radius_pm": p.atomic_radius_pm,
            "covalent_radius_pm": p.covalent_radius_pm,
            "polarizability_a0_cubed": p.polarizability_a0_cubed,
            "electron_configuration": p.electron_configuration,
            "frontier_orbitals": p.frontier_orbitals
        }));
    }

    let n = z_max as f64;
    let summary = json!({
        "elements_modeled": z_max,
        "mean_ionization_energy_ev": ie_sum / n,
        "mean_electronegativity_ev": chi_sum / n,
        "mean_atomic_radius_pm": radius_sum / n,
    });

    let report = json!({
        "meta": {
            "lane": "ab_initio_atomic_scf",
            "scan_bounds": {"z_min": 1, "z_max": z_max, "n_min": 1, "n_max": 260},
            "model": {
                "type": "spherical_atomic_scf",
                "constants": ["Rydberg", "Bohr radius"],
                "filling": "Madelung order",
                "frontier_mapping": "Koopmans-like HOMO/LUMO descriptors"
            }
        },
        "summary": summary,
        "elements": rows_json
    });

    let txt_path = out.join("ab_initio_qchem_report.txt");
    let csv_path = out.join("ab_initio_qchem_report.csv");
    let json_path = out.join("ab_initio_qchem_report.json");
    let orbital_path = out.join("ab_initio_qchem_orbitals.csv");

    fs::write(&csv_path, csv).expect("write csv");
    fs::write(&orbital_path, orbital_csv).expect("write orbital csv");
    fs::write(
        &json_path,
        serde_json::to_string_pretty(&report).expect("serialize report"),
    )
    .expect("write json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[ab_initio_qchem]").expect("write");
    writeln!(txt, "elements_modeled = {}", z_max).expect("write");
    writeln!(txt, "mean_ionization_energy_ev = {:.9}", ie_sum / n).expect("write");
    writeln!(txt, "mean_electronegativity_ev = {:.9}", chi_sum / n).expect("write");
    writeln!(txt, "mean_atomic_radius_pm = {:.9}", radius_sum / n).expect("write");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", csv_path.display());
    println!("wrote {}", orbital_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "ab_initio_qchem: elements={}, mean_IE={:.6} eV, mean_chi={:.6} eV, mean_radius={:.6} pm",
        z_max,
        ie_sum / n,
        chi_sum / n,
        radius_sum / n
    );
}
