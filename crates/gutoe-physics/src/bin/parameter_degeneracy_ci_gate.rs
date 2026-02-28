//! GRAND-354 CI gate: cosmology parameter-degeneracy audit.
//!
//! Produces a machine-readable degree-of-freedom inventory and Jacobian-based
//! degeneracy summary for the assembled cosmology lane.

use gutoe_em::{ckm_from_textures, PMNS_THETA23_ALPHA2_COEFF_STRUCTURAL};
use gutoe_physics::{
    eta_baryon_from_jarlskog, evaluate_inflation_gate, evaluate_transfer_gate,
    lambda_cosmological_full_candidate, leptogenesis_multiplier, primordial_deuterium_ratio,
    primordial_helium4_mass_fraction, InflationWindows, TransferAssumptions, TransferWindows, C,
    DARK_TO_VISIBLE_GEOMETRIC_RATIO, OMEGA_BARYON_OBS,
};
use serde_json::{json, Value};
use std::fs::{self, File};
use std::io::Write;

const METER_PER_MPC: f64 = 3.085_677_581_491_367e22;

#[derive(Debug, Clone, Copy)]
struct ParameterInventory {
    name: &'static str,
    category: &'static str,
    provenance: &'static str,
    value: f64,
    free: bool,
}

#[derive(Debug, Clone, Copy)]
struct Knob {
    name: &'static str,
    category: &'static str,
    provenance: &'static str,
    value: f64,
    rel_step: f64,
    abs_step: f64,
    free: bool,
}

#[derive(Debug, Clone, Copy)]
struct KnobState {
    pmns_theta23_alpha2_c: f64,
    leptogenesis_pmns_gain: f64,
    omega_r0: f64,
    omega_k0: f64,
}

#[derive(Debug, Clone, Copy)]
struct Outputs {
    h0_km_s_mpc: f64,
    omega_m0: f64,
    omega_lambda0: f64,
    rs_drag_mpc: f64,
    theta_star_rad: f64,
    l_peak1: f64,
    l_peak2: f64,
    yp: f64,
    dh: f64,
    eta_b: f64,
}

impl Outputs {
    fn rows(self) -> [(&'static str, f64); 10] {
        [
            ("h0_km_s_mpc", self.h0_km_s_mpc),
            ("omega_m0", self.omega_m0),
            ("omega_lambda0", self.omega_lambda0),
            ("rs_drag_mpc", self.rs_drag_mpc),
            ("theta_star_rad", self.theta_star_rad),
            ("l_peak1", self.l_peak1),
            ("l_peak2", self.l_peak2),
            ("yp", self.yp),
            ("dh", self.dh),
            ("eta_b", self.eta_b),
        ]
    }
}

fn baseline_state() -> KnobState {
    KnobState {
        pmns_theta23_alpha2_c: PMNS_THETA23_ALPHA2_COEFF_STRUCTURAL,
        leptogenesis_pmns_gain: 1.0,
        omega_r0: 9.0e-5,
        omega_k0: 0.0,
    }
}

fn analyze_knobs(s: KnobState) -> Vec<Knob> {
    vec![
        Knob {
            name: "pmns_theta23_alpha2_c",
            category: "tunable_runtime",
            provenance: "baryogenesis::leptogenesis_multiplier input",
            value: s.pmns_theta23_alpha2_c,
            rel_step: 1.0e-3,
            abs_step: 1.0e-6,
            free: true,
        },
        Knob {
            name: "leptogenesis_pmns_gain",
            category: "tunable_runtime",
            provenance: "baryogenesis PMNS coupling gain",
            value: s.leptogenesis_pmns_gain,
            rel_step: 1.0e-3,
            abs_step: 1.0e-6,
            free: true,
        },
        Knob {
            name: "omega_r0",
            category: "fixed_assumption",
            provenance: "UniverseAssumptions baseline radiation density",
            value: s.omega_r0,
            rel_step: 1.0e-3,
            abs_step: 1.0e-8,
            free: false,
        },
        Knob {
            name: "omega_k0",
            category: "fixed_assumption",
            provenance: "UniverseAssumptions spatial curvature anchor",
            value: s.omega_k0,
            rel_step: 1.0e-3,
            abs_step: 1.0e-5,
            free: false,
        },
    ]
}

fn parameter_inventory(s: KnobState, infl_n_s: f64, infl_a_s: f64) -> Vec<ParameterInventory> {
    vec![
        ParameterInventory {
            name: "alpha_inverse",
            category: "derived_constant",
            provenance: "Cl(1,3) triangular closure",
            value: 137.0,
            free: false,
        },
        ParameterInventory {
            name: "dark_to_visible_ratio",
            category: "derived_constant",
            provenance: "DARK_TO_VISIBLE_GEOMETRIC_RATIO",
            value: DARK_TO_VISIBLE_GEOMETRIC_RATIO,
            free: false,
        },
        ParameterInventory {
            name: "omega_b0",
            category: "observational_anchor",
            provenance: "OMEGA_BARYON_OBS",
            value: OMEGA_BARYON_OBS,
            free: false,
        },
        ParameterInventory {
            name: "n_s",
            category: "derived_constant",
            provenance: "inflation structural lane",
            value: infl_n_s,
            free: false,
        },
        ParameterInventory {
            name: "a_s",
            category: "derived_constant",
            provenance: "inflation structural lane",
            value: infl_a_s,
            free: false,
        },
        ParameterInventory {
            name: "pmns_theta23_alpha2_c",
            category: "tunable_runtime",
            provenance: "PMNS theta23 correction coefficient",
            value: s.pmns_theta23_alpha2_c,
            free: true,
        },
        ParameterInventory {
            name: "leptogenesis_pmns_gain",
            category: "tunable_runtime",
            provenance: "PMNS coupling gain",
            value: s.leptogenesis_pmns_gain,
            free: true,
        },
        ParameterInventory {
            name: "omega_r0",
            category: "fixed_assumption",
            provenance: "UniverseAssumptions::default",
            value: s.omega_r0,
            free: false,
        },
        ParameterInventory {
            name: "omega_k0",
            category: "fixed_assumption",
            provenance: "UniverseAssumptions::default",
            value: s.omega_k0,
            free: false,
        },
        ParameterInventory {
            name: "h0_ref_km_s_mpc",
            category: "observational_anchor",
            provenance: "UniverseAssumptions comparison target",
            value: 67.4,
            free: false,
        },
    ]
}

fn evaluate_outputs(s: KnobState) -> Option<Outputs> {
    let inflation = evaluate_inflation_gate(InflationWindows::default());
    let ckm = ckm_from_textures();

    let eta_b = eta_baryon_from_jarlskog(ckm.jarlskog)
        * leptogenesis_multiplier(s.pmns_theta23_alpha2_c, s.leptogenesis_pmns_gain);
    let eta10 = eta_b * 1.0e10;

    let yp = primordial_helium4_mass_fraction(eta10);
    let dh = primordial_deuterium_ratio(eta10);

    let omega_b0 = OMEGA_BARYON_OBS;
    let omega_dm0 = omega_b0 * DARK_TO_VISIBLE_GEOMETRIC_RATIO;
    let omega_m0 = omega_b0 + omega_dm0;
    let omega_lambda0 = 1.0 - omega_m0 - s.omega_r0 - s.omega_k0;
    if omega_lambda0 <= 0.0 {
        return None;
    }

    let lambda_full = lambda_cosmological_full_candidate();
    if lambda_full <= 0.0 {
        return None;
    }
    let h0_s_inv = C * (lambda_full / (3.0 * omega_lambda0)).sqrt();
    let h0_km_s_mpc = h0_s_inv * METER_PER_MPC / 1_000.0;

    let transfer = evaluate_transfer_gate(
        TransferAssumptions {
            h0_km_s_mpc,
            omega_b0,
            omega_m0,
            omega_r0: s.omega_r0,
            omega_k0: s.omega_k0,
            omega_lambda0,
            n_s: inflation.n_s,
            a_s: inflation.a_s,
        },
        TransferWindows::default(),
    );

    Some(Outputs {
        h0_km_s_mpc,
        omega_m0,
        omega_lambda0,
        rs_drag_mpc: transfer.rs_drag_mpc,
        theta_star_rad: transfer.theta_star_rad,
        l_peak1: transfer.l_peak1,
        l_peak2: transfer.l_peak2,
        yp,
        dh,
        eta_b,
    })
}

fn perturb_state(mut s: KnobState, name: &str, delta: f64) -> KnobState {
    match name {
        "pmns_theta23_alpha2_c" => s.pmns_theta23_alpha2_c += delta,
        "leptogenesis_pmns_gain" => s.leptogenesis_pmns_gain += delta,
        "omega_r0" => s.omega_r0 += delta,
        "omega_k0" => s.omega_k0 += delta,
        _ => {}
    }
    s
}

fn sensitivity_matrix(base_state: KnobState, knobs: &[Knob], base_outputs: Outputs) -> Option<Vec<Vec<f64>>> {
    let out_rows = base_outputs.rows();
    let mut matrix = vec![vec![0.0; knobs.len()]; out_rows.len()];

    for (j, k) in knobs.iter().enumerate() {
        let scale = if k.value.abs() > 0.0 {
            k.value.abs()
        } else {
            k.abs_step
        };
        let delta = if k.value.abs() > 0.0 {
            k.rel_step * k.value.abs()
        } else {
            k.abs_step
        };

        let plus = evaluate_outputs(perturb_state(base_state, k.name, delta))?;
        let minus = evaluate_outputs(perturb_state(base_state, k.name, -delta))?;
        let plus_rows = plus.rows();
        let minus_rows = minus.rows();

        for i in 0..out_rows.len() {
            let y0 = out_rows[i].1;
            let yp = plus_rows[i].1;
            let ym = minus_rows[i].1;
            if !(y0.is_finite() && yp.is_finite() && ym.is_finite()) || y0.abs() <= f64::EPSILON {
                matrix[i][j] = f64::NAN;
                continue;
            }
            let dy_dk = (yp - ym) / (2.0 * delta);
            matrix[i][j] = dy_dk * (scale / y0);
        }
    }

    Some(matrix)
}

fn gram_matrix(a: &[Vec<f64>]) -> Vec<Vec<f64>> {
    if a.is_empty() {
        return Vec::new();
    }
    let n = a[0].len();
    let mut g = vec![vec![0.0; n]; n];
    for row in a {
        for i in 0..n {
            let ai = if row[i].is_finite() { row[i] } else { 0.0 };
            for j in 0..n {
                let aj = if row[j].is_finite() { row[j] } else { 0.0 };
                g[i][j] += ai * aj;
            }
        }
    }
    g
}

fn jacobi_eigenvalues(mut a: Vec<Vec<f64>>) -> Vec<f64> {
    let n = a.len();
    if n == 0 {
        return Vec::new();
    }
    let max_iter = 128 * n * n;
    let tol = 1.0e-12;

    for _ in 0..max_iter {
        let mut p = 0usize;
        let mut q = 0usize;
        let mut max_off = 0.0f64;

        for i in 0..n {
            for j in (i + 1)..n {
                let v = a[i][j].abs();
                if v > max_off {
                    max_off = v;
                    p = i;
                    q = j;
                }
            }
        }

        if max_off < tol {
            break;
        }

        let app = a[p][p];
        let aqq = a[q][q];
        let apq = a[p][q];
        let phi = 0.5 * (2.0 * apq).atan2(aqq - app);
        let c = phi.cos();
        let s = phi.sin();

        for k in 0..n {
            if k == p || k == q {
                continue;
            }
            let akp = a[k][p];
            let akq = a[k][q];
            let new_kp = c * akp - s * akq;
            let new_kq = s * akp + c * akq;
            a[k][p] = new_kp;
            a[p][k] = new_kp;
            a[k][q] = new_kq;
            a[q][k] = new_kq;
        }

        let app_new = c * c * app - 2.0 * s * c * apq + s * s * aqq;
        let aqq_new = s * s * app + 2.0 * s * c * apq + c * c * aqq;
        a[p][p] = app_new;
        a[q][q] = aqq_new;
        a[p][q] = 0.0;
        a[q][p] = 0.0;
    }

    (0..n).map(|i| a[i][i]).collect()
}

fn sv_summary(a: &[Vec<f64>]) -> (Vec<f64>, usize, f64) {
    let g = gram_matrix(a);
    let mut evals = jacobi_eigenvalues(g)
        .into_iter()
        .map(|e| if e > 0.0 { e.sqrt() } else { 0.0 })
        .collect::<Vec<_>>();
    evals.sort_by(|x, y| y.total_cmp(x));

    if evals.is_empty() {
        return (evals, 0, f64::INFINITY);
    }

    let max_sv = evals[0].max(1.0e-16);
    let rank_tol = max_sv * 1.0e-6;
    let rank = evals.iter().filter(|&&v| v > rank_tol).count();

    let min_nonzero = evals
        .iter()
        .rev()
        .copied()
        .find(|&v| v > rank_tol)
        .unwrap_or(0.0);
    let cond = if min_nonzero > 0.0 {
        max_sv / min_nonzero
    } else {
        f64::INFINITY
    };

    (evals, rank, cond)
}

fn submatrix_columns(a: &[Vec<f64>], cols: &[usize]) -> Vec<Vec<f64>> {
    let mut out = Vec::with_capacity(a.len());
    for row in a {
        let mut r = Vec::with_capacity(cols.len());
        for &c in cols {
            r.push(*row.get(c).unwrap_or(&f64::NAN));
        }
        out.push(r);
    }
    out
}

fn max_abs_for_outputs_and_cols(
    a: &[Vec<f64>],
    output_rows: &[usize],
    cols: &[usize],
) -> f64 {
    let mut max_v = 0.0f64;
    for &ri in output_rows {
        for &cj in cols {
            if let Some(row) = a.get(ri) {
                if let Some(v) = row.get(cj) {
                    if v.is_finite() {
                        max_v = max_v.max(v.abs());
                    }
                }
            }
        }
    }
    max_v
}

fn main() {
    let out_dir = std::env::var("GUTOE_DEGENERACY_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let txt_path = format!("{out_dir}/parameter_degeneracy_ci_gate.txt");
    let json_path = format!("{out_dir}/parameter_degeneracy_ci_gate.json");

    let base_state = baseline_state();
    let inflation = evaluate_inflation_gate(InflationWindows::default());
    let inventory = parameter_inventory(base_state, inflation.n_s, inflation.a_s);
    let knobs = analyze_knobs(base_state);

    let base_outputs = match evaluate_outputs(base_state) {
        Some(v) => v,
        None => {
            eprintln!("degeneracy gate: failed to evaluate baseline outputs");
            std::process::exit(2);
        }
    };

    let matrix = match sensitivity_matrix(base_state, &knobs, base_outputs) {
        Some(m) => m,
        None => {
            eprintln!("degeneracy gate: failed to evaluate sensitivity matrix");
            std::process::exit(2);
        }
    };

    let (sv_all, rank_all, cond_all) = sv_summary(&matrix);

    let tunable_cols = knobs
        .iter()
        .enumerate()
        .filter_map(|(i, k)| if k.free { Some(i) } else { None })
        .collect::<Vec<_>>();
    let tunable_matrix = submatrix_columns(&matrix, &tunable_cols);
    let (sv_tunable, rank_tunable, cond_tunable) = sv_summary(&tunable_matrix);

    // Transfer-like outputs (H0, Ω terms, r_s, θ*, l1, l2) should not be
    // controllable by leptogenesis-only tunables.
    let transfer_rows = vec![0usize, 1, 2, 3, 4, 5, 6];
    let tunable_to_transfer_max = max_abs_for_outputs_and_cols(&matrix, &transfer_rows, &tunable_cols);

    // Tunables should still move baryogenesis observables.
    let baryo_rows = vec![7usize, 8, 9];
    let tunable_to_baryo_max = max_abs_for_outputs_and_cols(&matrix, &baryo_rows, &tunable_cols);

    let free_count = inventory.iter().filter(|p| p.free).count();
    let hidden_reencoding_risk = tunable_to_transfer_max > 1.0e-4 && rank_tunable >= 2;

    let mut warnings = Vec::new();
    if hidden_reencoding_risk {
        warnings.push("tunable parameters significantly couple into transfer/CMB outputs".to_string());
    }
    if tunable_to_baryo_max <= 1.0e-8 {
        warnings.push("tunable parameters do not move baryogenesis outputs (dead knobs)".to_string());
    }
    if rank_all < 2 {
        warnings.push("global sensitivity rank too low for a meaningful audit".to_string());
    }

    let overall_pass = free_count <= 2
        && tunable_to_transfer_max <= 1.0e-4
        && tunable_to_baryo_max > 1.0e-8
        && !hidden_reencoding_risk;

    let output_names = base_outputs
        .rows()
        .iter()
        .map(|(n, _)| *n)
        .collect::<Vec<_>>();

    let matrix_json = output_names
        .iter()
        .enumerate()
        .map(|(i, name)| {
            let mut row_obj = serde_json::Map::new();
            row_obj.insert("output".to_string(), Value::String((*name).to_string()));
            for (j, k) in knobs.iter().enumerate() {
                row_obj.insert(k.name.to_string(), json!(matrix[i][j]));
            }
            Value::Object(row_obj)
        })
        .collect::<Vec<_>>();

    let inventory_json = inventory
        .iter()
        .map(|p| {
            json!({
                "name": p.name,
                "category": p.category,
                "provenance": p.provenance,
                "value": p.value,
                "free": p.free
            })
        })
        .collect::<Vec<_>>();

    let knobs_json = knobs
        .iter()
        .map(|k| {
            json!({
                "name": k.name,
                "category": k.category,
                "provenance": k.provenance,
                "value": k.value,
                "rel_step": k.rel_step,
                "abs_step": k.abs_step,
                "free": k.free
            })
        })
        .collect::<Vec<_>>();

    let verdict = if hidden_reencoding_risk {
        "risk_detected"
    } else {
        "no_hidden_lcdm_reencoding_in_core_lane"
    };

    let report = json!({
        "overall_pass": overall_pass,
        "verdict": verdict,
        "counts": {
            "inventory_total": inventory.len(),
            "free_parameters": free_count,
            "analyzed_knobs": knobs.len(),
            "tunable_knobs": tunable_cols.len(),
            "outputs": output_names.len()
        },
        "inventory": inventory_json,
        "analyzed_knobs": knobs_json,
        "baseline_outputs": {
            "h0_km_s_mpc": base_outputs.h0_km_s_mpc,
            "omega_m0": base_outputs.omega_m0,
            "omega_lambda0": base_outputs.omega_lambda0,
            "rs_drag_mpc": base_outputs.rs_drag_mpc,
            "theta_star_rad": base_outputs.theta_star_rad,
            "l_peak1": base_outputs.l_peak1,
            "l_peak2": base_outputs.l_peak2,
            "yp": base_outputs.yp,
            "dh": base_outputs.dh,
            "eta_b": base_outputs.eta_b
        },
        "sensitivity": {
            "definition": "dimensionless relative sensitivity: (d output / d knob) * (scale_knob / output)",
            "rows": matrix_json
        },
        "linear_algebra": {
            "all_knobs": {
                "singular_values": sv_all,
                "rank": rank_all,
                "condition_number": cond_all
            },
            "tunable_only": {
                "singular_values": sv_tunable,
                "rank": rank_tunable,
                "condition_number": cond_tunable
            }
        },
        "hidden_reencoding_checks": {
            "tunable_to_transfer_max_abs_sensitivity": tunable_to_transfer_max,
            "tunable_to_baryo_max_abs_sensitivity": tunable_to_baryo_max,
            "hidden_lcdm_reencoding_risk": hidden_reencoding_risk
        },
        "analysis_fit_knob_presence": {
            "notes": [
                "Scan/report-only knobs exist in analysis tools (e.g. tau_reio and likelihood scan ranges).",
                "This audit verdict is for the core assembled lane used in CI gates, not for exploratory scan binaries."
            ]
        },
        "warnings": warnings
    });

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[parameter_degeneracy_ci_gate]").ok();
    writeln!(txt, "overall_pass = {}", overall_pass).ok();
    writeln!(txt, "verdict = {}", verdict).ok();
    writeln!(txt, "free_parameters = {}", free_count).ok();
    writeln!(txt, "rank_all = {}", rank_all).ok();
    writeln!(txt, "rank_tunable = {}", rank_tunable).ok();
    writeln!(txt, "condition_all = {:.9e}", cond_all).ok();
    writeln!(txt, "condition_tunable = {:.9e}", cond_tunable).ok();
    writeln!(
        txt,
        "tunable_to_transfer_max_abs_sensitivity = {:.9e}",
        tunable_to_transfer_max
    )
    .ok();
    writeln!(
        txt,
        "tunable_to_baryo_max_abs_sensitivity = {:.9e}",
        tunable_to_baryo_max
    )
    .ok();

    if !warnings.is_empty() {
        writeln!(txt, "warnings:").ok();
        for w in &warnings {
            writeln!(txt, "- {}", w).ok();
        }
    }

    let mut jf = File::create(&json_path).expect("create json");
    writeln!(
        jf,
        "{}",
        serde_json::to_string_pretty(&report).expect("serialize json")
    )
    .ok();

    println!(
        "degeneracy gate: pass={} (free={}, rank_tunable={}, transfer_coupling={:.3e})",
        overall_pass, free_count, rank_tunable, tunable_to_transfer_max
    );
    println!("wrote {txt_path}");
    println!("wrote {json_path}");

    if !overall_pass {
        std::process::exit(2);
    }
}
