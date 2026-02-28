//! Empirical calibration lane for coupled chemical thermodynamics.
//!
//! This calibrates a small, constrained coefficient set against external density
//! references and reports train-vs-holdout metrics for structural generalization:
//! - Period holdout: train periods 1..4, hold 5..7
//! - Block holdout A: train s/p, hold d/f
//! - Block holdout B: train d/f, hold s/p

use anyhow::{anyhow, Context, Result};
use gutoe_physics::{
    family_of_z, period_of_z, prefetch_element_thermo_coupled,
    predict_element_thermo_coupled_from_prefetch_calibrated, ChemicalFamily,
    ChemicalThermoCalibration, CoupledThermoPrefetch,
};
use std::collections::BTreeMap;
use std::env;
use std::fs;
use std::path::{Path, PathBuf};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Block {
    S,
    P,
    D,
    F,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum RefState {
    Solid,
    Liquid,
    Gas,
}

#[derive(Clone, Debug)]
struct Sample {
    z: u16,
    period: u8,
    block: Block,
    state_ref: RefState,
    density_ref: f64,
    prefetch: CoupledThermoPrefetch,
}

#[derive(Clone, Copy, Debug)]
struct Metrics {
    n: usize,
    density_mae: f64,
    density_mape_pct: f64,
    phase_accuracy: f64,
}

#[derive(Clone, Copy, Debug)]
struct Scenario {
    name: &'static str,
    is_train: fn(&Sample) -> bool,
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

fn parse_state(v: Option<&String>) -> Option<RefState> {
    let s = v?.trim().to_ascii_lowercase();
    if s.contains("solid") {
        Some(RefState::Solid)
    } else if s.contains("liquid") {
        Some(RefState::Liquid)
    } else if s.contains("gas") {
        Some(RefState::Gas)
    } else {
        None
    }
}

fn block_of_z(z: u16) -> Block {
    match family_of_z(z) {
        ChemicalFamily::Alkali | ChemicalFamily::AlkalineEarth => Block::S,
        ChemicalFamily::Transition => Block::D,
        ChemicalFamily::Lanthanide | ChemicalFamily::Actinide => Block::F,
        ChemicalFamily::PostTransition
        | ChemicalFamily::Metalloid
        | ChemicalFamily::Nonmetal
        | ChemicalFamily::Halogen
        | ChemicalFamily::NobleGas => Block::P,
    }
}

fn gas_molecularity_guess(z: u16) -> f64 {
    match z {
        1 | 7 | 8 | 9 | 17 => 2.0,
        _ => 1.0,
    }
}

fn representative_mass_number(z: u16) -> u16 {
    (2.5 * z as f64).round() as u16
}

fn density_pred_stateaware(z: u16, molar_mass_g_mol: f64, density_pred: f64, state_ref: RefState) -> f64 {
    match state_ref {
        RefState::Gas => (molar_mass_g_mol * gas_molecularity_guess(z)) / (22.413_97 * 1000.0),
        _ => density_pred,
    }
}

fn eval_subset<F>(samples: &[Sample], cal: ChemicalThermoCalibration, pred: F) -> Metrics
where
    F: Fn(&Sample) -> bool,
{
    let mut n = 0usize;
    let mut sum_abs = 0.0;
    let mut sum_pct = 0.0;
    let mut phase_n = 0usize;
    let mut phase_match = 0usize;

    for s in samples.iter().filter(|s| pred(s)) {
        let p = predict_element_thermo_coupled_from_prefetch_calibrated(&s.prefetch, cal).0;
        let dp = density_pred_stateaware(s.z, p.molar_mass_g_mol, p.density_g_cm3, s.state_ref);
        let ae = (dp - s.density_ref).abs();
        let pe = if s.density_ref.abs() > 0.0 {
            100.0 * ae / s.density_ref.abs()
        } else {
            0.0
        };
        n += 1;
        sum_abs += ae;
        sum_pct += pe;

        let state_pred = match p.ambient_state_298k {
            gutoe_physics::MatterState::Solid => RefState::Solid,
            gutoe_physics::MatterState::Liquid => RefState::Liquid,
            gutoe_physics::MatterState::Gas => RefState::Gas,
        };
        phase_n += 1;
        if state_pred == s.state_ref {
            phase_match += 1;
        }
    }

    Metrics {
        n,
        density_mae: if n > 0 { sum_abs / n as f64 } else { 0.0 },
        density_mape_pct: if n > 0 { sum_pct / n as f64 } else { 0.0 },
        phase_accuracy: if phase_n > 0 {
            phase_match as f64 / phase_n as f64
        } else {
            0.0
        },
    }
}

fn cal_to_vec(cal: ChemicalThermoCalibration) -> [f64; 13] {
    [
        cal.pack_p_void_coef,
        cal.pack_d_gain_coef,
        cal.pack_f_gain_coef,
        cal.pack_open_d_mult,
        cal.pack_closed_d_mult,
        cal.pack_f_core_mult,
        cal.radius_p_gain,
        cal.radius_closed_d_mult,
        cal.radius_open_d_mult,
        cal.radius_f_core_mult,
        cal.radius_actinide_mult,
        cal.radius_lower_actinide,
        cal.radius_lower_transition_fcore,
    ]
}

fn vec_to_cal(v: [f64; 13]) -> ChemicalThermoCalibration {
    ChemicalThermoCalibration {
        pack_p_void_coef: v[0],
        pack_d_gain_coef: v[1],
        pack_f_gain_coef: v[2],
        pack_open_d_mult: v[3],
        pack_closed_d_mult: v[4],
        pack_f_core_mult: v[5],
        radius_p_gain: v[6],
        radius_closed_d_mult: v[7],
        radius_open_d_mult: v[8],
        radius_f_core_mult: v[9],
        radius_actinide_mult: v[10],
        radius_lower_actinide: v[11],
        radius_lower_transition_fcore: v[12],
    }
}

fn bounds() -> [(f64, f64); 13] {
    [
        (0.0, 1.2),
        (0.0, 0.8),
        (0.0, 0.8),
        (0.8, 1.6),
        (0.4, 1.0),
        (0.8, 1.5),
        (0.0, 0.8),
        (0.8, 1.3),
        (0.5, 1.05),
        (0.5, 1.0),
        (0.4, 1.0),
        (0.2, 0.9),
        (0.3, 0.9),
    ]
}

fn clamp_vec(mut v: [f64; 13]) -> [f64; 13] {
    for (i, (lo, hi)) in bounds().iter().enumerate() {
        v[i] = v[i].clamp(*lo, *hi);
    }
    v
}

fn objective_train(samples: &[Sample], cal: ChemicalThermoCalibration, is_train: fn(&Sample) -> bool) -> f64 {
    eval_subset(samples, cal, is_train).density_mae
}

fn optimize_for_scenario(samples: &[Sample], scenario: Scenario, init: ChemicalThermoCalibration) -> ChemicalThermoCalibration {
    let mut best_v = cal_to_vec(init);
    let mut best_score = objective_train(samples, vec_to_cal(best_v), scenario.is_train);
    let b = bounds();

    for pass in 0..8 {
        let base_frac = 0.24 * 0.5_f64.powi(pass);
        let mut improved = false;

        for i in 0..best_v.len() {
            let span = b[i].1 - b[i].0;
            let step = (span * base_frac).max(1.0e-4);
            for dir in [-1.0, 1.0] {
                let mut cand = best_v;
                cand[i] += dir * step;
                cand = clamp_vec(cand);
                let score = objective_train(samples, vec_to_cal(cand), scenario.is_train);
                if score + 1.0e-9 < best_score {
                    best_score = score;
                    best_v = cand;
                    improved = true;
                }
            }
        }

        if !improved && pass >= 4 {
            break;
        }
    }

    vec_to_cal(best_v)
}

fn train_period_le4(s: &Sample) -> bool {
    s.period <= 4
}

fn train_sp_block(s: &Sample) -> bool {
    matches!(s.block, Block::S | Block::P)
}

fn train_df_block(s: &Sample) -> bool {
    matches!(s.block, Block::D | Block::F)
}

fn main() -> Result<()> {
    let reference_path = env::var("GUTOE_REFERENCE_TABLE")
        .unwrap_or_else(|_| "crates/gutoe-physics/data/periodic_pubchem_reference.csv".to_string());
    let out_dir = env::var("GUTOE_CHEM_CAL_OUT")
        .unwrap_or_else(|_| "/tmp/nuclear_chart".to_string());
    let z_min = env::var("GUTOE_BENCH_Z_MIN")
        .ok()
        .and_then(|v| v.parse::<u16>().ok())
        .unwrap_or(1);
    let z_max = env::var("GUTOE_BENCH_Z_MAX")
        .ok()
        .and_then(|v| v.parse::<u16>().ok())
        .unwrap_or(94);

    fs::create_dir_all(&out_dir)?;
    let out = PathBuf::from(out_dir);

    let rows = read_csv_rows(Path::new(&reference_path))?;
    let mut samples = Vec::new();
    for r in rows {
        let Some(z) = parse_u16(r.get("z")) else {
            continue;
        };
        if z < z_min || z > z_max {
            continue;
        }
        let Some(density_ref) = parse_f64(r.get("density_g_cm3")) else {
            continue;
        };
        let Some(state_ref) = parse_state(r.get("state_298k")) else {
            continue;
        };
        samples.push(Sample {
            z,
            period: period_of_z(z),
            block: block_of_z(z),
            state_ref,
            density_ref,
            prefetch: prefetch_element_thermo_coupled(z, representative_mass_number(z)),
        });
    }

    if samples.is_empty() {
        return Err(anyhow!("no calibration samples loaded"));
    }

    let scenarios = [
        Scenario {
            name: "period_holdout_train_p1_p4_hold_p5_p7",
            is_train: train_period_le4,
        },
        Scenario {
            name: "block_holdout_train_sp_hold_df",
            is_train: train_sp_block,
        },
        Scenario {
            name: "block_holdout_train_df_hold_sp",
            is_train: train_df_block,
        },
    ];

    let init = ChemicalThermoCalibration::default();

    let mut txt = String::new();
    txt.push_str("[chemical_thermo_calibrate]\n");
    txt.push_str(&format!("z_range = {}..{}\n", z_min, z_max));
    txt.push_str(&format!("reference = {}\n", reference_path));
    txt.push_str(&format!("samples = {}\n\n", samples.len()));

    let mut json_rows = Vec::new();

    for sc in scenarios {
        let trained = optimize_for_scenario(&samples, sc, init);

        let train_base = eval_subset(&samples, init, sc.is_train);
        let hold_base = eval_subset(&samples, init, |s| !((sc.is_train)(s)));
        let train_fit = eval_subset(&samples, trained, sc.is_train);
        let hold_fit = eval_subset(&samples, trained, |s| !((sc.is_train)(s)));

        txt.push_str(&format!("scenario = {}\n", sc.name));
        txt.push_str(&format!(
            "  baseline_train: n={} density_mae={:.6} mape_pct={:.6} phase_acc={:.6}\n",
            train_base.n, train_base.density_mae, train_base.density_mape_pct, train_base.phase_accuracy
        ));
        txt.push_str(&format!(
            "  baseline_holdout: n={} density_mae={:.6} mape_pct={:.6} phase_acc={:.6}\n",
            hold_base.n, hold_base.density_mae, hold_base.density_mape_pct, hold_base.phase_accuracy
        ));
        txt.push_str(&format!(
            "  fitted_train: n={} density_mae={:.6} mape_pct={:.6} phase_acc={:.6}\n",
            train_fit.n, train_fit.density_mae, train_fit.density_mape_pct, train_fit.phase_accuracy
        ));
        txt.push_str(&format!(
            "  fitted_holdout: n={} density_mae={:.6} mape_pct={:.6} phase_acc={:.6}\n",
            hold_fit.n, hold_fit.density_mae, hold_fit.density_mape_pct, hold_fit.phase_accuracy
        ));
        txt.push_str(&format!(
            "  generalization_gap_holdout_minus_train = {:.6}\n\n",
            hold_fit.density_mae - train_fit.density_mae
        ));

        json_rows.push(format!(
            concat!(
                "{{",
                "\"name\":\"{}\",",
                "\"train_n\":{},\"holdout_n\":{},",
                "\"baseline_train_mae\":{:.9},\"baseline_holdout_mae\":{:.9},",
                "\"fitted_train_mae\":{:.9},\"fitted_holdout_mae\":{:.9},",
                "\"fitted_train_phase_acc\":{:.9},\"fitted_holdout_phase_acc\":{:.9},",
                "\"generalization_gap\":{:.9},",
                "\"coefficients\":{{",
                "\"pack_p_void_coef\":{:.9},",
                "\"pack_d_gain_coef\":{:.9},",
                "\"pack_f_gain_coef\":{:.9},",
                "\"pack_open_d_mult\":{:.9},",
                "\"pack_closed_d_mult\":{:.9},",
                "\"pack_f_core_mult\":{:.9},",
                "\"radius_p_gain\":{:.9},",
                "\"radius_closed_d_mult\":{:.9},",
                "\"radius_open_d_mult\":{:.9},",
                "\"radius_f_core_mult\":{:.9},",
                "\"radius_actinide_mult\":{:.9},",
                "\"radius_lower_actinide\":{:.9},",
                "\"radius_lower_transition_fcore\":{:.9}",
                "}}",
                "}}"
            ),
            sc.name,
            train_fit.n,
            hold_fit.n,
            train_base.density_mae,
            hold_base.density_mae,
            train_fit.density_mae,
            hold_fit.density_mae,
            train_fit.phase_accuracy,
            hold_fit.phase_accuracy,
            hold_fit.density_mae - train_fit.density_mae,
            trained.pack_p_void_coef,
            trained.pack_d_gain_coef,
            trained.pack_f_gain_coef,
            trained.pack_open_d_mult,
            trained.pack_closed_d_mult,
            trained.pack_f_core_mult,
            trained.radius_p_gain,
            trained.radius_closed_d_mult,
            trained.radius_open_d_mult,
            trained.radius_f_core_mult,
            trained.radius_actinide_mult,
            trained.radius_lower_actinide,
            trained.radius_lower_transition_fcore,
        ));
    }

    let txt_path = out.join("chemical_thermo_calibration_report.txt");
    let json_path = out.join("chemical_thermo_calibration_report.json");

    fs::write(&txt_path, txt)?;
    let json = format!(
        "{{\n  \"meta\": {{\"reference\": \"{}\", \"z_min\": {}, \"z_max\": {}, \"samples\": {}}},\n  \"scenarios\": [\n    {}\n  ]\n}}\n",
        reference_path,
        z_min,
        z_max,
        samples.len(),
        json_rows.join(",\n    ")
    );
    fs::write(&json_path, json)?;

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());

    Ok(())
}
