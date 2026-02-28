//! External benchmark pass for unified element table.
//!
//! Compares unified table predictions against external periodic references
//! (default: PubChem periodic table JSON export prepared as CSV) and emits
//! per-element error rows + aggregate scorecard.

use anyhow::{anyhow, Context, Result};
use std::collections::BTreeMap;
use std::env;
use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};

#[derive(Clone, Debug)]
struct RefRow {
    z: u16,
    symbol: String,
    name: String,
    state_298k: Option<String>,
    density_g_cm3: Option<f64>,
    melting_k: Option<f64>,
    boiling_k: Option<f64>,
    ionization_energy_ev: Option<f64>,
}

#[derive(Clone, Debug)]
struct BenchCounters {
    n_phase: usize,
    n_phase_match: usize,
    n_density: usize,
    n_density_stateaware: usize,
    n_melting: usize,
    n_boiling: usize,
    n_ionization: usize,
    n_density_condensed: usize,
    sum_abs_density: f64,
    sum_pct_density: f64,
    sum_abs_density_stateaware: f64,
    sum_pct_density_stateaware: f64,
    sum_abs_density_condensed: f64,
    sum_pct_density_condensed: f64,
    sum_abs_melting: f64,
    sum_pct_melting: f64,
    sum_abs_boiling: f64,
    sum_pct_boiling: f64,
    sum_abs_ionization: f64,
    sum_pct_ionization: f64,
    red_phase: usize,
    red_density: usize,
    red_density_stateaware: usize,
    red_melting: usize,
    red_boiling: usize,
    red_ionization: usize,
    red_any: usize,
}

fn gas_molecularity_guess(z: u16) -> f64 {
    match z {
        1 | 7 | 8 | 9 | 17 => 2.0, // H2, N2, O2, F2, Cl2
        _ => 1.0,
    }
}

fn parse_csv_line(line: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut cur = String::new();
    let mut in_quotes = false;
    let mut chars = line.chars().peekable();
    while let Some(ch) = chars.next() {
        match ch {
            '"' => {
                if in_quotes && matches!(chars.peek(), Some('"')) {
                    cur.push('"');
                    let _ = chars.next();
                } else {
                    in_quotes = !in_quotes;
                }
            }
            ',' if !in_quotes => {
                out.push(cur);
                cur = String::new();
            }
            _ => cur.push(ch),
        }
    }
    out.push(cur);
    out
}

fn read_csv_rows(path: &Path) -> Result<Vec<BTreeMap<String, String>>> {
    let content = fs::read_to_string(path)
        .with_context(|| format!("failed to read CSV {}", path.display()))?;
    let mut lines = content.lines();
    let header_line = lines
        .next()
        .ok_or_else(|| anyhow!("CSV has no header: {}", path.display()))?;
    let headers = parse_csv_line(header_line);
    if headers.is_empty() {
        return Err(anyhow!("CSV has empty header: {}", path.display()));
    }

    let mut rows = Vec::new();
    for line in lines {
        if line.trim().is_empty() {
            continue;
        }
        let fields = parse_csv_line(line);
        let mut row = BTreeMap::new();
        for (i, h) in headers.iter().enumerate() {
            let v = fields.get(i).cloned().unwrap_or_default();
            row.insert(h.clone(), v);
        }
        rows.push(row);
    }
    Ok(rows)
}

fn parse_u16(v: Option<&String>) -> Option<u16> {
    v.and_then(|s| s.trim().parse::<u16>().ok())
}

fn parse_f64(v: Option<&String>) -> Option<f64> {
    v.and_then(|s| {
        let t = s.trim();
        if t.is_empty() || t.eq_ignore_ascii_case("nan") {
            None
        } else {
            t.parse::<f64>().ok()
        }
    })
}

fn parse_state(v: Option<&String>) -> Option<String> {
    let s = v?.trim().to_ascii_lowercase();
    if s.is_empty() {
        return None;
    }
    if s.contains("solid") {
        Some("solid".to_string())
    } else if s.contains("liquid") {
        Some("liquid".to_string())
    } else if s.contains("gas") {
        Some("gas".to_string())
    } else {
        None
    }
}

fn env_u16(name: &str, default: u16) -> u16 {
    env::var(name)
        .ok()
        .and_then(|v| v.parse::<u16>().ok())
        .unwrap_or(default)
}

fn main() -> Result<()> {
    let unified_path = env::var("GUTOE_UNIFIED_TABLE")
        .unwrap_or_else(|_| "/tmp/nuclear_chart/element_unified_algebra_table.csv".to_string());
    let reference_path = env::var("GUTOE_REFERENCE_TABLE").unwrap_or_else(|_| {
        "crates/gutoe-physics/data/periodic_pubchem_reference.csv".to_string()
    });
    let out_dir =
        env::var("GUTOE_BENCH_OUT").unwrap_or_else(|_| "/tmp/nuclear_chart".to_string());
    let z_min = env_u16("GUTOE_BENCH_Z_MIN", 1);
    let z_max = env_u16("GUTOE_BENCH_Z_MAX", 94);

    fs::create_dir_all(&out_dir)?;
    let out = PathBuf::from(out_dir);

    let unified_rows = read_csv_rows(Path::new(&unified_path))?;
    let reference_rows = read_csv_rows(Path::new(&reference_path))?;

    let mut unified_by_z: BTreeMap<u16, BTreeMap<String, String>> = BTreeMap::new();
    for r in unified_rows {
        if let Some(z) = parse_u16(r.get("Z")) {
            unified_by_z.insert(z, r);
        }
    }

    let mut refs_by_z: BTreeMap<u16, RefRow> = BTreeMap::new();
    for r in reference_rows {
        let Some(z) = parse_u16(r.get("z")) else {
            continue;
        };
        refs_by_z.insert(
            z,
            RefRow {
                z,
                symbol: r.get("symbol").cloned().unwrap_or_default(),
                name: r.get("name").cloned().unwrap_or_default(),
                state_298k: parse_state(r.get("state_298k")),
                density_g_cm3: parse_f64(r.get("density_g_cm3")),
                melting_k: parse_f64(r.get("melting_k")),
                boiling_k: parse_f64(r.get("boiling_k")),
                ionization_energy_ev: parse_f64(r.get("ionization_energy_ev")),
            },
        );
    }

    // Red-light thresholds (can be tuned).
    const RED_DENSITY_PCT: f64 = 20.0;
    const RED_MELTING_PCT: f64 = 20.0;
    const RED_BOILING_PCT: f64 = 20.0;
    const RED_IONIZATION_PCT: f64 = 10.0;

    let mut c = BenchCounters {
        n_phase: 0,
        n_phase_match: 0,
        n_density: 0,
        n_density_stateaware: 0,
        n_melting: 0,
        n_boiling: 0,
        n_ionization: 0,
        n_density_condensed: 0,
        sum_abs_density: 0.0,
        sum_pct_density: 0.0,
        sum_abs_density_stateaware: 0.0,
        sum_pct_density_stateaware: 0.0,
        sum_abs_density_condensed: 0.0,
        sum_pct_density_condensed: 0.0,
        sum_abs_melting: 0.0,
        sum_pct_melting: 0.0,
        sum_abs_boiling: 0.0,
        sum_pct_boiling: 0.0,
        sum_abs_ionization: 0.0,
        sum_pct_ionization: 0.0,
        red_phase: 0,
        red_density: 0,
        red_density_stateaware: 0,
        red_melting: 0,
        red_boiling: 0,
        red_ionization: 0,
        red_any: 0,
    };

    let mut out_csv = String::from(
        "z,symbol,name,state_pred,state_ref,state_match,density_pred,density_pred_stateaware,density_ref,density_abs_err,density_pct_err,density_red,density_stateaware_abs_err,density_stateaware_pct_err,density_stateaware_red,melting_pred_k,melting_ref_k,melting_abs_err,melting_pct_err,melting_red,boiling_pred_k,boiling_ref_k,boiling_abs_err,boiling_pct_err,boiling_red,ionization_pred_ev,ionization_ref_ev,ionization_abs_err,ionization_pct_err,ionization_red,red_count,red_any\n",
    );

    for z in z_min..=z_max {
        let Some(pred) = unified_by_z.get(&z) else {
            continue;
        };
        let Some(reference) = refs_by_z.get(&z) else {
            continue;
        };

        let state_pred = parse_state(pred.get("chem_state_298k_1bar"));
        let state_ref = reference.state_298k.clone();
        let is_condensed_ref = !matches!(state_ref.as_deref(), Some("gas"));
        let state_match = match (&state_pred, &state_ref) {
            (Some(a), Some(b)) => a == b,
            _ => false,
        };
        let mut red_count = 0usize;
        if let (Some(_), Some(_)) = (&state_pred, &state_ref) {
            c.n_phase += 1;
            if state_match {
                c.n_phase_match += 1;
            } else {
                c.red_phase += 1;
                red_count += 1;
            }
        }

        let density_pred = parse_f64(pred.get("chem_density_g_cm3"));
        let density_pred_stateaware = {
            let molar_mass = parse_f64(pred.get("chem_molar_mass_g_mol"));
            match (molar_mass, is_condensed_ref) {
                (Some(m), false) => Some((m * gas_molecularity_guess(z)) / (22.413_97 * 1000.0)),
                _ => density_pred,
            }
        };
        let density_ref = reference.density_g_cm3;
        let (density_abs_err, density_pct_err, density_red) = match (density_pred, density_ref) {
            (Some(p), Some(r)) if r.abs() > 0.0 => {
                c.n_density += 1;
                let ae = (p - r).abs();
                let pe = 100.0 * ae / r.abs();
                c.sum_abs_density += ae;
                c.sum_pct_density += pe;
                if is_condensed_ref {
                    c.n_density_condensed += 1;
                    c.sum_abs_density_condensed += ae;
                    c.sum_pct_density_condensed += pe;
                }
                let red = pe > RED_DENSITY_PCT;
                if red {
                    c.red_density += 1;
                    red_count += 1;
                }
                (Some(ae), Some(pe), red)
            }
            _ => (None, None, false),
        };
        let (density_stateaware_abs_err, density_stateaware_pct_err, density_stateaware_red) =
            match (density_pred_stateaware, density_ref) {
                (Some(p), Some(r)) if r.abs() > 0.0 => {
                    c.n_density_stateaware += 1;
                    let ae = (p - r).abs();
                    let pe = 100.0 * ae / r.abs();
                    c.sum_abs_density_stateaware += ae;
                    c.sum_pct_density_stateaware += pe;
                    let red = pe > RED_DENSITY_PCT;
                    if red {
                        c.red_density_stateaware += 1;
                    }
                    (Some(ae), Some(pe), red)
                }
                _ => (None, None, false),
            };

        let melting_pred = parse_f64(pred.get("chem_melting_temperature_k"));
        let melting_ref = reference.melting_k;
        let (melting_abs_err, melting_pct_err, melting_red) = match (melting_pred, melting_ref) {
            (Some(p), Some(r)) if r.abs() > 0.0 => {
                c.n_melting += 1;
                let ae = (p - r).abs();
                let pe = 100.0 * ae / r.abs();
                c.sum_abs_melting += ae;
                c.sum_pct_melting += pe;
                let red = pe > RED_MELTING_PCT;
                if red {
                    c.red_melting += 1;
                    red_count += 1;
                }
                (Some(ae), Some(pe), red)
            }
            _ => (None, None, false),
        };

        let boiling_pred = parse_f64(pred.get("chem_boiling_temperature_k"));
        let boiling_ref = reference.boiling_k;
        let (boiling_abs_err, boiling_pct_err, boiling_red) = match (boiling_pred, boiling_ref) {
            (Some(p), Some(r)) if r.abs() > 0.0 => {
                c.n_boiling += 1;
                let ae = (p - r).abs();
                let pe = 100.0 * ae / r.abs();
                c.sum_abs_boiling += ae;
                c.sum_pct_boiling += pe;
                let red = pe > RED_BOILING_PCT;
                if red {
                    c.red_boiling += 1;
                    red_count += 1;
                }
                (Some(ae), Some(pe), red)
            }
            _ => (None, None, false),
        };

        let ion_pred = parse_f64(pred.get("scf_ionization_energy_ev"));
        let ion_ref = reference.ionization_energy_ev;
        let (ion_abs_err, ion_pct_err, ion_red) = match (ion_pred, ion_ref) {
            (Some(p), Some(r)) if r.abs() > 0.0 => {
                c.n_ionization += 1;
                let ae = (p - r).abs();
                let pe = 100.0 * ae / r.abs();
                c.sum_abs_ionization += ae;
                c.sum_pct_ionization += pe;
                let red = pe > RED_IONIZATION_PCT;
                if red {
                    c.red_ionization += 1;
                    red_count += 1;
                }
                (Some(ae), Some(pe), red)
            }
            _ => (None, None, false),
        };

        if red_count > 0 {
            c.red_any += 1;
        }

        out_csv.push_str(&format!(
            "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n",
            reference.z,
            reference.symbol,
            reference.name,
            state_pred.clone().unwrap_or_default(),
            state_ref.clone().unwrap_or_default(),
            state_match,
            density_pred.map(|v| format!("{v:.6}")).unwrap_or_default(),
            density_pred_stateaware
                .map(|v| format!("{v:.6}"))
                .unwrap_or_default(),
            density_ref.map(|v| format!("{v:.6}")).unwrap_or_default(),
            density_abs_err.map(|v| format!("{v:.6}")).unwrap_or_default(),
            density_pct_err.map(|v| format!("{v:.6}")).unwrap_or_default(),
            density_red,
            density_stateaware_abs_err
                .map(|v| format!("{v:.6}"))
                .unwrap_or_default(),
            density_stateaware_pct_err
                .map(|v| format!("{v:.6}"))
                .unwrap_or_default(),
            density_stateaware_red,
            melting_pred.map(|v| format!("{v:.6}")).unwrap_or_default(),
            melting_ref.map(|v| format!("{v:.6}")).unwrap_or_default(),
            melting_abs_err.map(|v| format!("{v:.6}")).unwrap_or_default(),
            melting_pct_err.map(|v| format!("{v:.6}")).unwrap_or_default(),
            melting_red,
            boiling_pred.map(|v| format!("{v:.6}")).unwrap_or_default(),
            boiling_ref.map(|v| format!("{v:.6}")).unwrap_or_default(),
            boiling_abs_err.map(|v| format!("{v:.6}")).unwrap_or_default(),
            boiling_pct_err.map(|v| format!("{v:.6}")).unwrap_or_default(),
            boiling_red,
            ion_pred.map(|v| format!("{v:.6}")).unwrap_or_default(),
            ion_ref.map(|v| format!("{v:.6}")).unwrap_or_default(),
            ion_abs_err.map(|v| format!("{v:.6}")).unwrap_or_default(),
            ion_pct_err.map(|v| format!("{v:.6}")).unwrap_or_default(),
            ion_red,
            red_count,
            red_count > 0
        ));
    }

    let phase_accuracy = if c.n_phase > 0 {
        c.n_phase_match as f64 / c.n_phase as f64
    } else {
        0.0
    };
    let density_mae = if c.n_density > 0 {
        c.sum_abs_density / c.n_density as f64
    } else {
        0.0
    };
    let density_mape = if c.n_density > 0 {
        c.sum_pct_density / c.n_density as f64
    } else {
        0.0
    };
    let density_condensed_mae = if c.n_density_condensed > 0 {
        c.sum_abs_density_condensed / c.n_density_condensed as f64
    } else {
        0.0
    };
    let density_condensed_mape = if c.n_density_condensed > 0 {
        c.sum_pct_density_condensed / c.n_density_condensed as f64
    } else {
        0.0
    };
    let density_stateaware_mae = if c.n_density_stateaware > 0 {
        c.sum_abs_density_stateaware / c.n_density_stateaware as f64
    } else {
        0.0
    };
    let density_stateaware_mape = if c.n_density_stateaware > 0 {
        c.sum_pct_density_stateaware / c.n_density_stateaware as f64
    } else {
        0.0
    };
    let melting_mae = if c.n_melting > 0 {
        c.sum_abs_melting / c.n_melting as f64
    } else {
        0.0
    };
    let melting_mape = if c.n_melting > 0 {
        c.sum_pct_melting / c.n_melting as f64
    } else {
        0.0
    };
    let boiling_mae = if c.n_boiling > 0 {
        c.sum_abs_boiling / c.n_boiling as f64
    } else {
        0.0
    };
    let boiling_mape = if c.n_boiling > 0 {
        c.sum_pct_boiling / c.n_boiling as f64
    } else {
        0.0
    };
    let ion_mae = if c.n_ionization > 0 {
        c.sum_abs_ionization / c.n_ionization as f64
    } else {
        0.0
    };
    let ion_mape = if c.n_ionization > 0 {
        c.sum_pct_ionization / c.n_ionization as f64
    } else {
        0.0
    };

    let summary_json = format!(
        concat!(
            "{{\n",
            "  \"meta\": {{\n",
            "    \"unified_table\": \"{}\",\n",
            "    \"reference_table\": \"{}\",\n",
            "    \"z_range\": {{\"min\": {}, \"max\": {}}},\n",
            "    \"reference_source\": \"PubChem periodic table JSON CSV export\"\n",
            "  }},\n",
            "  \"phase\": {{\"n\": {}, \"accuracy\": {:.9}, \"red\": {}}},\n",
            "  \"density_g_cm3\": {{\"n\": {}, \"mae\": {:.9}, \"mape_pct\": {:.9}, \"red\": {}}},\n",
            "  \"density_g_cm3_stateaware\": {{\"n\": {}, \"mae\": {:.9}, \"mape_pct\": {:.9}, \"red\": {}}},\n",
            "  \"density_g_cm3_condensed_only\": {{\"n\": {}, \"mae\": {:.9}, \"mape_pct\": {:.9}}},\n",
            "  \"melting_k\": {{\"n\": {}, \"mae\": {:.9}, \"mape_pct\": {:.9}, \"red\": {}}},\n",
            "  \"boiling_k\": {{\"n\": {}, \"mae\": {:.9}, \"mape_pct\": {:.9}, \"red\": {}}},\n",
            "  \"ionization_energy_ev\": {{\"n\": {}, \"mae\": {:.9}, \"mape_pct\": {:.9}, \"red\": {}}},\n",
            "  \"elements_with_any_red\": {}\n",
            "}}\n"
        ),
        unified_path,
        reference_path,
        z_min,
        z_max,
        c.n_phase,
        phase_accuracy,
        c.red_phase,
        c.n_density,
        density_mae,
        density_mape,
        c.red_density,
        c.n_density_stateaware,
        density_stateaware_mae,
        density_stateaware_mape,
        c.red_density_stateaware,
        c.n_density_condensed,
        density_condensed_mae,
        density_condensed_mape,
        c.n_melting,
        melting_mae,
        melting_mape,
        c.red_melting,
        c.n_boiling,
        boiling_mae,
        boiling_mape,
        c.red_boiling,
        c.n_ionization,
        ion_mae,
        ion_mape,
        c.red_ionization,
        c.red_any
    );

    let txt_summary = format!(
        concat!(
            "[element_unified_external_benchmark]\n",
            "z_range = {}..{}\n",
            "reference = {}\n",
            "phase_accuracy = {:.6} (n={}, red={})\n",
            "density_mae_g_cm3 = {:.6} (mape_pct={:.6}, n={}, red={})\n",
            "density_mae_g_cm3_stateaware = {:.6} (mape_pct={:.6}, n={}, red={})\n",
            "density_mae_g_cm3_condensed_only = {:.6} (mape_pct={:.6}, n={})\n",
            "melting_mae_k = {:.6} (mape_pct={:.6}, n={}, red={})\n",
            "boiling_mae_k = {:.6} (mape_pct={:.6}, n={}, red={})\n",
            "ionization_mae_ev = {:.6} (mape_pct={:.6}, n={}, red={})\n",
            "elements_with_any_red = {}\n"
        ),
        z_min,
        z_max,
        reference_path,
        phase_accuracy,
        c.n_phase,
        c.red_phase,
        density_mae,
        density_mape,
        c.n_density,
        c.red_density,
        density_stateaware_mae,
        density_stateaware_mape,
        c.n_density_stateaware,
        c.red_density_stateaware,
        density_condensed_mae,
        density_condensed_mape,
        c.n_density_condensed,
        melting_mae,
        melting_mape,
        c.n_melting,
        c.red_melting,
        boiling_mae,
        boiling_mape,
        c.n_boiling,
        c.red_boiling,
        ion_mae,
        ion_mape,
        c.n_ionization,
        c.red_ionization,
        c.red_any
    );

    let csv_path = out.join("element_unified_external_benchmark.csv");
    let json_path = out.join("element_unified_external_benchmark.json");
    let txt_path = out.join("element_unified_external_benchmark.txt");

    fs::write(&csv_path, out_csv)?;
    fs::write(&json_path, summary_json)?;
    let mut txt = File::create(&txt_path)?;
    txt.write_all(txt_summary.as_bytes())?;

    println!("wrote {}", csv_path.display());
    println!("wrote {}", json_path.display());
    println!("wrote {}", txt_path.display());

    Ok(())
}
