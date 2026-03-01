//! Inverse-design bridge for CTC/FTL simulation lanes (steps 1-4).
//!
//! This tool turns an existing campaign output into:
//! 1) Minimal door invariants.
//! 2) Sensitivity + ablation results.
//! 3) Dimensionless target windows.
//! 4) Candidate physical analog mappings from those windows.
//!
//! Scope:
//! - Systems/inverse-design analysis on simulation outputs.
//! - Not a physical engine claim.

use serde_json::{json, Value};
use std::fs;
use std::path::{Path, PathBuf};

const FC_VOID_DEFAULT: f64 = 3.0 / 16.0;
const BRANCHING_DEFAULT: f64 = 3.0;
const ETA_DEFAULT: f64 = 4.0 / 6.0;
const INFRA_DEFAULT: f64 = 16.0 / 6.0;
const PDG_ALPHA_INV_2025: f64 = 137.035_999_177;
const PDG_SIN2_THETAW_MZ_MS_2025: f64 = 0.23122;
const PDG_N_NU_2025: f64 = 2.9963;
const PDG_MW_GEV_2025: f64 = 80.3692;
const PDG_MZ_GEV_2025: f64 = 91.1880;

#[derive(Clone, Copy, Debug)]
struct Params {
    beta: f64,
    s: f64,
    kappa: f64,
    f_void: f64,
    radius_m: f64,
    period_s: f64,
    budget_per_door_j: f64,
    n_loops: f64,
    branching: f64,
    eta: f64,
    infra: f64,
}

#[derive(Clone, Copy, Debug)]
struct Metrics {
    threshold_j: f64,
    q_eff: f64,
    coordinate_speed_over_c: f64,
    t_eff_norm: f64,
    predeparture_margin: f64,
    predeparture: bool,
    local_timelike: bool,
    coordinate_superluminal: bool,
    gate_open: bool,
    topology_gain: f64,
    topology_residual: f64,
    stability_margin: f64,
}

#[derive(Clone, Copy, Debug)]
struct ParamSweepRow {
    beta: f64,
    s: f64,
    q_eff: f64,
    n_loops: f64,
    stability_margin: f64,
    objective: f64,
}

fn clamp_positive(x: f64, floor: f64) -> f64 {
    x.abs().max(floor)
}

fn env_f64(name: &str, default: f64) -> f64 {
    std::env::var(name)
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(default)
}

fn path_get_f64(v: &Value, path: &[&str]) -> Option<f64> {
    let mut cur = v;
    for key in path {
        cur = cur.get(*key)?;
    }
    cur.as_f64()
}

fn path_get_bool(v: &Value, path: &[&str]) -> Option<bool> {
    let mut cur = v;
    for key in path {
        cur = cur.get(*key)?;
    }
    cur.as_bool()
}

fn threshold_j(p: &Params) -> f64 {
    clamp_positive(
        p.kappa * p.f_void * p.radius_m.abs() * p.period_s.abs(),
        1e-30,
    )
}

fn metrics(p: &Params) -> Metrics {
    let th = threshold_j(p);
    let q_eff = p.budget_per_door_j / th;
    let coordinate_speed_over_c = 1.0 / p.s;
    let t_eff_norm = 1.0 + p.n_loops * (1.0 - q_eff);
    let predeparture_margin = -t_eff_norm;
    let local_timelike = p.beta < 1.0;
    let coordinate_superluminal = coordinate_speed_over_c > 1.0;
    let gate_open = q_eff >= 1.0;
    let topology_gain = p.branching * p.f_void * p.eta * p.infra;
    let topology_residual = (topology_gain - 1.0).abs();
    let stability_margin = (1.0 - p.beta).min(1.0 - p.s).min(q_eff - 1.0);
    Metrics {
        threshold_j: th,
        q_eff,
        coordinate_speed_over_c,
        t_eff_norm,
        predeparture_margin,
        predeparture: t_eff_norm < 0.0,
        local_timelike,
        coordinate_superluminal,
        gate_open,
        topology_gain,
        topology_residual,
        stability_margin,
    }
}

fn read_campaign_json(path: &Path) -> Option<Value> {
    let raw = fs::read_to_string(path).ok()?;
    serde_json::from_str::<Value>(&raw).ok()
}

fn choose_row(v: &Value) -> Option<&Value> {
    let rows = v.get("rows")?.as_array()?;
    if let Some(row) = rows
        .iter()
        .find(|r| path_get_bool(r, &["predeparture_enabled"]) == Some(true))
    {
        return Some(row);
    }
    rows.last()
}

fn pct_delta(base: f64, shifted: f64) -> f64 {
    if base.abs() < 1e-20 {
        0.0
    } else {
        (shifted / base - 1.0) * 100.0
    }
}

fn perturb_params(mut p: Params, key: &str, delta: f64) -> Params {
    match key {
        "beta" => p.beta = (p.beta * (1.0 + delta)).clamp(1e-6, 0.999_999),
        "s" => p.s = (p.s * (1.0 + delta)).clamp(1e-6, 0.999_999),
        "kappa" => p.kappa = clamp_positive(p.kappa * (1.0 + delta), 1e-30),
        "f_void" => p.f_void = clamp_positive(p.f_void * (1.0 + delta), 1e-9),
        "radius_m" => p.radius_m = clamp_positive(p.radius_m * (1.0 + delta), 1e-30),
        "period_s" => p.period_s = clamp_positive(p.period_s * (1.0 + delta), 1e-30),
        "budget_per_door_j" => {
            p.budget_per_door_j = clamp_positive(p.budget_per_door_j * (1.0 + delta), 1e-30)
        }
        "n_loops" => p.n_loops = clamp_positive(p.n_loops * (1.0 + delta), 1.0),
        _ => {}
    }
    p
}

fn ablation_case(name: &str, p: Params) -> Value {
    let m = metrics(&p);
    let viable = m.local_timelike && m.coordinate_superluminal && m.gate_open && m.predeparture;
    json!({
      "name": name,
      "beta": p.beta,
      "s": p.s,
      "q_eff": m.q_eff,
      "n_loops": p.n_loops,
      "t_eff_norm": m.t_eff_norm,
      "predeparture_margin": m.predeparture_margin,
      "local_timelike": m.local_timelike,
      "coordinate_superluminal": m.coordinate_superluminal,
      "gate_open": m.gate_open,
      "viable": viable,
      "stability_margin": m.stability_margin
    })
}

fn min_max(values: &[f64]) -> Option<(f64, f64)> {
    if values.is_empty() {
        return None;
    }
    let min_v = values.iter().copied().fold(f64::INFINITY, f64::min);
    let max_v = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    Some((min_v, max_v))
}

fn rel_pct(model: f64, measured: f64) -> f64 {
    if measured.abs() < 1e-20 {
        0.0
    } else {
        100.0 * (model / measured - 1.0)
    }
}

fn main() {
    let out_dir = std::env::var("GUTOE_INVERSE_DESIGN_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ctc_inverse_design_steps_1_4".to_string());
    let campaign_path = std::env::var("GUTOE_INVERSE_DESIGN_CAMPAIGN_JSON")
        .unwrap_or_else(|_| "/tmp/bh_renders/ctc_50y_campaign_fast/ctc_50y_campaign_fast.json".to_string());
    let delta = env_f64("GUTOE_INVERSE_DESIGN_SENS_DELTA", 0.05).clamp(0.001, 0.25);
    let n_loops_default = env_f64("GUTOE_INVERSE_DESIGN_N_LOOPS", 120.0).max(1.0);

    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let campaign_json = read_campaign_json(Path::new(&campaign_path)).unwrap_or_else(|| json!({}));
    let row = choose_row(&campaign_json).cloned().unwrap_or_else(|| json!({}));

    let beta = path_get_f64(&campaign_json, &["inputs", "local_beta"])
        .or_else(|| path_get_f64(&row, &["beta"]))
        .unwrap_or(0.8)
        .clamp(1e-6, 0.999_999);
    let s = path_get_f64(&campaign_json, &["inputs", "shortcut_factor_s"])
        .unwrap_or(0.1)
        .clamp(1e-6, 0.999_999);
    let kappa = path_get_f64(&campaign_json, &["inputs", "kappa_j_per_m_s"])
        .unwrap_or(5.645_474_097_135_454e37)
        .max(1e-30);
    let f_void = path_get_f64(&campaign_json, &["inputs", "void_fraction"])
        .unwrap_or(FC_VOID_DEFAULT)
        .max(1e-9);
    let radius_m = path_get_f64(&row, &["radius_m"]).unwrap_or(8.0e-19).max(1e-30);
    let period_s = path_get_f64(&row, &["period_s"]).unwrap_or(5.0e-27).max(1e-30);
    let doors_active = path_get_f64(&row, &["doors_active"]).unwrap_or(1.0).max(1e-12);
    let program_budget_j = path_get_f64(&row, &["program_budget_j"]).unwrap_or(1.0).max(1e-30);
    let budget_per_door_j = (program_budget_j / doors_active).max(1e-30);
    let n_loops = n_loops_default;

    let p0 = Params {
        beta,
        s,
        kappa,
        f_void,
        radius_m,
        period_s,
        budget_per_door_j,
        n_loops,
        branching: BRANCHING_DEFAULT,
        eta: ETA_DEFAULT,
        infra: INFRA_DEFAULT,
    };
    let m0 = metrics(&p0);

    // Step 1: minimal invariants
    let invariants = json!({
      "local_timelike_required": "beta < 1",
      "coordinate_shortcut_required": "s < 1 and chi = 1/s > 1",
      "budget_gate_required": "q_eff = E_door / (kappa*f_void*R*T) >= 1",
      "predeparture_required": "mu = n*(q_eff-1) > 1",
      "topology_criticality_required": "g = branching*f_void*eta*infra ≈ 1",
      "baseline_values": {
        "beta": p0.beta,
        "s": p0.s,
        "chi_coordinate_speed_over_c": m0.coordinate_speed_over_c,
        "q_eff": m0.q_eff,
        "mu_n_q_minus_1": p0.n_loops * (m0.q_eff - 1.0),
        "t_eff_norm": m0.t_eff_norm,
        "topology_gain_g": m0.topology_gain,
        "topology_residual_abs_g_minus_1": m0.topology_residual,
        "stability_margin": m0.stability_margin,
        "threshold_j": m0.threshold_j,
        "budget_per_door_j": p0.budget_per_door_j
      },
      "baseline_gate_checks": {
        "local_timelike": m0.local_timelike,
        "coordinate_superluminal": m0.coordinate_superluminal,
        "budget_gate_open": m0.gate_open,
        "predeparture": m0.predeparture
      }
    });

    // Step 2: sensitivity and ablation
    let knobs = [
        "beta",
        "s",
        "kappa",
        "f_void",
        "radius_m",
        "period_s",
        "budget_per_door_j",
        "n_loops",
    ];
    let mut sensitivities = Vec::new();
    for key in knobs {
        let p_minus = perturb_params(p0, key, -delta);
        let p_plus = perturb_params(p0, key, delta);
        let m_minus = metrics(&p_minus);
        let m_plus = metrics(&p_plus);
        sensitivities.push(json!({
          "parameter": key,
          "delta_fraction": delta,
          "q_eff_pct_on_plus": pct_delta(m0.q_eff, m_plus.q_eff),
          "q_eff_pct_on_minus": pct_delta(m0.q_eff, m_minus.q_eff),
          "predeparture_margin_pct_on_plus": pct_delta(m0.predeparture_margin, m_plus.predeparture_margin),
          "predeparture_margin_pct_on_minus": pct_delta(m0.predeparture_margin, m_minus.predeparture_margin),
          "coord_speed_over_c_pct_on_plus": pct_delta(m0.coordinate_speed_over_c, m_plus.coordinate_speed_over_c),
          "coord_speed_over_c_pct_on_minus": pct_delta(m0.coordinate_speed_over_c, m_minus.coordinate_speed_over_c),
          "stability_margin_pct_on_plus": pct_delta(m0.stability_margin, m_plus.stability_margin),
          "stability_margin_pct_on_minus": pct_delta(m0.stability_margin, m_minus.stability_margin)
        }));
    }

    let th0 = m0.threshold_j;
    let ablations = vec![
        ablation_case("baseline", p0),
        ablation_case(
            "disable_shortcut_channel_s_to_1",
            Params {
                s: 1.0,
                ..p0
            },
        ),
        ablation_case(
            "break_local_timelike_beta_to_1",
            Params {
                beta: 1.0,
                ..p0
            },
        ),
        ablation_case(
            "subgate_budget_90pct",
            Params {
                budget_per_door_j: 0.9 * th0,
                ..p0
            },
        ),
        ablation_case(
            "break_even_budget",
            Params {
                budget_per_door_j: 1.0 * th0,
                ..p0
            },
        ),
        ablation_case("no_loop_accumulation_n0", Params { n_loops: 0.0, ..p0 }),
        ablation_case(
            "visible_weight_f_void_13_over_16",
            Params {
                f_void: 13.0 / 16.0,
                ..p0
            },
        ),
        ablation_case(
            "double_geometry_scale_R_and_T",
            Params {
                radius_m: 2.0 * p0.radius_m,
                period_s: 2.0 * p0.period_s,
                ..p0
            },
        ),
    ];

    // Step 3: dimensionless target windows from grid search.
    let beta_grid = [0.60, 0.70, 0.80, 0.90, 0.95];
    let s_grid = [0.05, 0.10, 0.20, 0.30, 0.40];
    let q_grid = [0.80, 1.00, 1.10, 1.20, 1.50, 2.00];
    let n_grid = [20.0, 60.0, 120.0, 200.0];

    let mut passing = Vec::<ParamSweepRow>::new();
    let mut robust = Vec::<ParamSweepRow>::new();

    for beta_i in beta_grid {
        for s_i in s_grid {
            for q_i in q_grid {
                for n_i in n_grid {
                    let mut p = p0;
                    p.beta = beta_i;
                    p.s = s_i;
                    p.n_loops = n_i;
                    p.budget_per_door_j = q_i * threshold_j(&p);
                    let m = metrics(&p);
                    let objective = m.predeparture_margin.max(0.0)
                        * m.coordinate_speed_over_c
                        * m.stability_margin.max(0.0);
                    let row = ParamSweepRow {
                        beta: beta_i,
                        s: s_i,
                        q_eff: m.q_eff,
                        n_loops: n_i,
                        stability_margin: m.stability_margin,
                        objective,
                    };
                    let pass =
                        m.local_timelike && m.coordinate_superluminal && m.gate_open && m.predeparture;
                    if pass {
                        passing.push(row);
                    }
                    if pass && m.stability_margin > 0.05 && m.predeparture_margin > 0.5 {
                        robust.push(row);
                    }
                }
            }
        }
    }

    let total_grid = (beta_grid.len() * s_grid.len() * q_grid.len() * n_grid.len()) as f64;
    let pass_rate = if total_grid > 0.0 {
        passing.len() as f64 / total_grid
    } else {
        0.0
    };
    let robust_rate = if total_grid > 0.0 {
        robust.len() as f64 / total_grid
    } else {
        0.0
    };

    let robust_betas: Vec<f64> = robust.iter().map(|r| r.beta).collect();
    let robust_ss: Vec<f64> = robust.iter().map(|r| r.s).collect();
    let robust_qs: Vec<f64> = robust.iter().map(|r| r.q_eff).collect();
    let robust_ns: Vec<f64> = robust.iter().map(|r| r.n_loops).collect();

    let robust_windows = json!({
      "beta_local_window": min_max(&robust_betas).map(|(mn,mx)| json!({"min": mn, "max": mx})).unwrap_or(json!(null)),
      "shortcut_factor_s_window": min_max(&robust_ss).map(|(mn,mx)| json!({"min": mn, "max": mx})).unwrap_or(json!(null)),
      "q_eff_window": min_max(&robust_qs).map(|(mn,mx)| json!({"min": mn, "max": mx})).unwrap_or(json!(null)),
      "n_loops_window": min_max(&robust_ns).map(|(mn,mx)| json!({"min": mn, "max": mx})).unwrap_or(json!(null))
    });

    let best = robust
        .iter()
        .max_by(|a, b| a.objective.total_cmp(&b.objective))
        .copied()
        .or_else(|| {
            passing
                .iter()
                .max_by(|a, b| a.objective.total_cmp(&b.objective))
                .copied()
        });

    let target_dimensionless = if let Some(best_row) = best {
        json!({
          "Pi1_q_eff": best_row.q_eff,
          "Pi2_chi_coord_over_c": 1.0 / best_row.s,
          "Pi3_beta_local": best_row.beta,
          "Pi4_mu_loop_margin_n_q_minus_1": best_row.n_loops * (best_row.q_eff - 1.0),
          "Pi5_g_topology_gain": m0.topology_gain,
          "Pi6_eps_criticality_abs_g_minus_1": m0.topology_residual,
          "Pi7_stability_margin": best_row.stability_margin
        })
    } else {
        json!(null)
    };

    let model_alpha_inv = 137.0_f64;
    let model_sin2 = 3.0_f64 / 13.0_f64;
    let model_n_nu = 3.0_f64;
    let pdg_sin2_on_shell = 1.0 - (PDG_MW_GEV_2025 / PDG_MZ_GEV_2025).powi(2);
    let pdg_anchor_residuals = json!({
      "alpha_inv": {
        "model": model_alpha_inv,
        "pdg": PDG_ALPHA_INV_2025,
        "rel_pct_model_minus_pdg": rel_pct(model_alpha_inv, PDG_ALPHA_INV_2025)
      },
      "sin2_thetaW_MS_MZ": {
        "model": model_sin2,
        "pdg": PDG_SIN2_THETAW_MZ_MS_2025,
        "rel_pct_model_minus_pdg": rel_pct(model_sin2, PDG_SIN2_THETAW_MZ_MS_2025)
      },
      "n_nu": {
        "model": model_n_nu,
        "pdg": PDG_N_NU_2025,
        "abs_model_minus_pdg": model_n_nu - PDG_N_NU_2025,
        "rel_pct_model_minus_pdg": rel_pct(model_n_nu, PDG_N_NU_2025)
      },
      "pdg_sin2_on_shell_from_masses": {
        "value": pdg_sin2_on_shell,
        "mw_gev": PDG_MW_GEV_2025,
        "mz_gev": PDG_MZ_GEV_2025
      }
    });

    // Step 4: physical analog map from those Pi groups.
    let pi = target_dimensionless.clone();
    let analog_map = json!([
      {
        "class": "Phase-Locked Photonic Delay Mesh",
        "mapping": {
          "Pi1_q_eff": "loop_phase_gain / loop_loss",
          "Pi2_chi": "tau_reference / tau_shortcut_path",
          "Pi3_beta": "v_group / v_medium_cap",
          "Pi4_mu": "N_loops * (loop_gain-1)"
        },
        "primary_controls": ["coupler ratio", "delay-line length", "phase-lock bandwidth", "loop count gate"],
        "lab_observables": ["group delay histograms", "phase closure residual", "early-arrival index"],
        "target_pi": pi
      },
      {
        "class": "Superconducting Resonator Loop (Flux/Phase Domain)",
        "mapping": {
          "Pi1_q_eff": "pump_power / effective dissipation",
          "Pi2_chi": "normal transit latency / coupled transit latency",
          "Pi3_beta": "signal velocity ratio in bounded line",
          "Pi5_g": "symmetry-constrained gain product"
        },
        "primary_controls": ["Q-factor tuning", "pump detuning", "coupling inductance", "feedback phase"],
        "lab_observables": ["ring-down time", "phase slip count", "causal-order anomaly statistic"],
        "target_pi": pi
      },
      {
        "class": "RF Ring + Digital Twin Controller",
        "mapping": {
          "Pi1_q_eff": "feedback writeback / per-cycle drain",
          "Pi2_chi": "software-routed latency collapse ratio",
          "Pi4_mu": "retry-depth * (effective gain-1)"
        },
        "primary_controls": ["PLL gains", "buffer schedule", "feedback delay taps", "retry scheduler"],
        "lab_observables": ["timestamp inversion counts", "deterministic replay agreement", "energy ledger closure"],
        "target_pi": pi
      },
      {
        "class": "Analog-Gravity Metamaterial Loop",
        "mapping": {
          "Pi2_chi": "effective metric path compression",
          "Pi3_beta": "local wave speed / local cap",
          "Pi7_stability": "minimum margin to local bound and gate"
        },
        "primary_controls": ["index gradient", "modulation depth", "boundary actuation timing"],
        "lab_observables": ["wavefront causal cone checks", "shortcut-factor inference", "bounded-energy closure"],
        "target_pi": pi
      }
    ]);

    let payload = json!({
      "scope": "inverse-design simulation-to-analog bridge (steps 1-4); not physical engine claim",
      "inputs": {
        "campaign_json_path": campaign_path,
        "sensitivity_delta_fraction": delta,
        "n_loops_for_invariants": n_loops_default
      },
      "step_1_minimal_invariants": invariants,
      "step_2_sensitivity_and_ablation": {
        "local_sensitivity": sensitivities,
        "ablation_cases": ablations
      },
      "step_3_dimensionless_targets": {
        "grid_total": total_grid,
        "passing_count": passing.len(),
        "passing_rate": pass_rate,
        "robust_count": robust.len(),
        "robust_rate": robust_rate,
        "robust_windows": robust_windows,
        "recommended_pi_targets": target_dimensionless,
        "pdg_ratio_anchors_2025": pdg_anchor_residuals
      },
      "step_4_physical_analog_map": analog_map
    });

    let txt_path = out.join("ctc_inverse_design_steps_1_4.txt");
    let json_path = out.join("ctc_inverse_design_steps_1_4.json");

    let mut txt = String::new();
    txt.push_str("[ctc_inverse_design_steps_1_4]\n");
    txt.push_str("simulation-to-analog inverse-design bridge\n\n");
    txt.push_str("step1: minimal invariants extracted\n");
    txt.push_str(&format!("  beta={:.6}, s={:.6}, q_eff={:.6e}\n", p0.beta, p0.s, m0.q_eff));
    txt.push_str(&format!(
        "  chi=u/c={:.6}, mu=n(q-1)={:.6e}, g={:.6e}, |g-1|={:.3e}\n",
        m0.coordinate_speed_over_c,
        p0.n_loops * (m0.q_eff - 1.0),
        m0.topology_gain,
        m0.topology_residual
    ));
    txt.push_str(&format!(
        "  gate checks: timelike={} superluminal={} gate_open={} predeparture={}\n",
        m0.local_timelike, m0.coordinate_superluminal, m0.gate_open, m0.predeparture
    ));
    txt.push_str("step2: sensitivity+ablation complete\n");
    txt.push_str("step3: dimensionless windows derived from grid\n");
    txt.push_str(&format!(
        "  passing={}/{} ({:.3}) robust={}/{} ({:.3})\n",
        passing.len(),
        total_grid as usize,
        pass_rate,
        robust.len(),
        total_grid as usize,
        robust_rate
    ));
    txt.push_str("  PDG anchors (2025):\n");
    txt.push_str(&format!(
        "    alpha^-1 model/pdg = {:.6}/{:.9} (rel {:.6}%)\n",
        model_alpha_inv,
        PDG_ALPHA_INV_2025,
        rel_pct(model_alpha_inv, PDG_ALPHA_INV_2025)
    ));
    txt.push_str(&format!(
        "    sin^2(theta_W)_MS model/pdg = {:.9}/{:.5} (rel {:.6}%)\n",
        model_sin2,
        PDG_SIN2_THETAW_MZ_MS_2025,
        rel_pct(model_sin2, PDG_SIN2_THETAW_MZ_MS_2025)
    ));
    txt.push_str(&format!(
        "    N_nu model/pdg = {:.4}/{:.4} (abs {:.6})\n",
        model_n_nu,
        PDG_N_NU_2025,
        model_n_nu - PDG_N_NU_2025
    ));
    txt.push_str("step4: physical analog mappings produced\n");

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(
        &json_path,
        serde_json::to_string_pretty(&payload).expect("serialize json"),
    )
    .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "inverse-design complete: pass_rate={:.3}, robust_rate={:.3}",
        pass_rate, robust_rate
    );
}
