//! Rapid 20-test falsification harness for GUTOE live-fire campaign.
//!
//! Produces a machine-readable and human-readable checklist with:
//! - dataset/source lane,
//! - test statistic,
//! - pass threshold,
//! - kill threshold,
//! - current status (from available run artifacts when possible).
//!
//! Scope:
//! - Program management + falsification tracking.
//! - Not a claim that all tests are executed here.

use serde_json::{json, Value};
use std::fs;
use std::path::PathBuf;

const ALPHA_INV_CODATA_2022: f64 = 137.035_999_177;
const MP_ME_CODATA_2022: f64 = 1836.152_673_43;

#[derive(Clone, Copy)]
enum Status {
    Pass,
    Fail,
    Open,
}

impl Status {
    fn as_str(self) -> &'static str {
        match self {
            Status::Pass => "PASS",
            Status::Fail => "FAIL",
            Status::Open => "OPEN",
        }
    }
}

struct TestCase {
    id: &'static str,
    name: &'static str,
    dataset: &'static str,
    metric: &'static str,
    pass_if: &'static str,
    kill_if: &'static str,
    status: Status,
    observed: String,
}

fn read_json(path: &str) -> Option<Value> {
    let raw = fs::read_to_string(path).ok()?;
    serde_json::from_str::<Value>(&raw).ok()
}

fn weak_model_reduced_chi2(v: &Value, model_name: &str) -> Option<f64> {
    v.get("model_fits")?
        .as_array()?
        .iter()
        .find(|m| m.get("model").and_then(Value::as_str) == Some(model_name))
        .and_then(|m| m.get("reduced_chi2"))
        .and_then(Value::as_f64)
}

fn weak_model_max_abs_pull(v: &Value, model_name: &str) -> Option<f64> {
    let pulls = v
        .get("model_fits")?
        .as_array()?
        .iter()
        .find(|m| m.get("model").and_then(Value::as_str) == Some(model_name))?
        .get("pulls")?
        .as_array()?;
    let mut max_abs = 0.0_f64;
    for p in pulls {
        if let Some(z) = p.get("pull_sigma").and_then(Value::as_f64) {
            max_abs = max_abs.max(z.abs());
        }
    }
    Some(max_abs)
}

fn solve_alpha_inv_structural() -> f64 {
    // x = 137 + 5/x - b/x^2 - c/x^3, with b=9+5/32 and c=1/125.
    let b = 9.0 + 5.0 / 32.0;
    let c = 1.0 / 125.0;
    let mut x = 137.036_f64;
    for _ in 0..40 {
        let x2 = x * x;
        let x3 = x2 * x;
        let x4 = x3 * x;
        let f = x - (137.0 + 5.0 / x - b / x2 - c / x3);
        let df = 1.0 - (-5.0 / x2 + 2.0 * b / x3 + 3.0 * c / x4);
        x -= f / df;
    }
    x
}

fn is_core_algebraic(test_id: &str) -> bool {
    matches!(
        test_id,
        "T01" | "T02" | "T03" | "T11" | "T12" | "T13" | "T14" | "T15" | "T16" | "T21" | "T22"
    )
}

fn main() {
    let out_dir = std::env::var("GUTOE_FALSIFICATION_20_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ctc_falsification_20".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let weak_path = std::env::var("GUTOE_WEAK_FIT_JSON").unwrap_or_else(|_| {
        "/tmp/bh_renders/ctc_weak_angle_identity_fit/ctc_weak_angle_identity_fit.json".to_string()
    });
    let door_path = std::env::var("GUTOE_DOOR_FIT_JSON").unwrap_or_else(|_| {
        "/tmp/bh_renders/ctc_public_door_joint_fit/ctc_public_door_joint_fit.json".to_string()
    });
    let next3_path = std::env::var("GUTOE_NEXT3_JSON").unwrap_or_else(|_| {
        "/tmp/bh_renders/ctc_next3_t07_t11_t13/ctc_next3_t07_t11_t13.json".to_string()
    });
    let t01678_path = std::env::var("GUTOE_T01678_JSON").unwrap_or_else(|_| {
        "/tmp/bh_renders/ctc_t01_t06_t07_t08/ctc_t01_t06_t07_t08.json".to_string()
    });
    let remaining12_path = std::env::var("GUTOE_REMAINING12_JSON").unwrap_or_else(|_| {
        "/tmp/bh_renders/ctc_remaining12_runner/ctc_remaining12_runner.json".to_string()
    });

    let weak = read_json(&weak_path).unwrap_or_else(|| json!({}));
    let door = read_json(&door_path).unwrap_or_else(|| json!({}));
    let next3 = read_json(&next3_path).unwrap_or_else(|| json!({}));
    let t01678 = read_json(&t01678_path).unwrap_or_else(|| json!({}));
    let remaining12 = read_json(&remaining12_path).unwrap_or_else(|| json!({}));

    let red_base = weak_model_reduced_chi2(&weak, "base_fixed");
    let red_zero = weak_model_reduced_chi2(&weak, "rational_plus_alpha_ln10_over_4pi_log_fixed");
    let pull_zero = weak_model_max_abs_pull(&weak, "rational_plus_alpha_ln10_over_4pi_log_fixed");
    let z_stouffer = door
        .get("joint_fit")
        .and_then(|j| j.get("z_stouffer"))
        .and_then(Value::as_f64);
    let door_detected = door
        .get("joint_fit")
        .and_then(|j| j.get("door_detected"))
        .and_then(Value::as_bool);

    let mut tests = vec![
        TestCase {
            id: "T01",
            name: "Kill naive weak-angle identity 3/13",
            dataset: "EW 9-point anchor set",
            metric: "reduced_chi2(base_fixed)",
            pass_if: "> 5.0 (naive lane falsified)",
            kill_if: "<= 5.0",
            status: Status::Open,
            observed: "pending".to_string(),
        },
        TestCase {
            id: "T02",
            name: "Zero-free weak-angle identity survival",
            dataset: "EW 9-point anchor set",
            metric: "reduced_chi2(zero_free_formula)",
            pass_if: "<= 1.5",
            kill_if: "> 2.5",
            status: Status::Open,
            observed: "pending".to_string(),
        },
        TestCase {
            id: "T03",
            name: "Zero-free weak-angle max pull guard",
            dataset: "EW 9-point anchor set",
            metric: "max_abs_pull_sigma(zero_free_formula)",
            pass_if: "<= 2.5",
            kill_if: "> 3.0",
            status: Status::Open,
            observed: "pending".to_string(),
        },
        TestCase {
            id: "T04",
            name: "Spontaneous-door global excess",
            dataset: "Public coarse joint lanes (CB + EHT + low-Q EW)",
            metric: "z_stouffer",
            pass_if: ">= 5.0 with lane coherence",
            kill_if: "< 3.0",
            status: Status::Open,
            observed: "pending".to_string(),
        },
        TestCase {
            id: "T05",
            name: "Cosmic birefringence isotropic non-zero",
            dataset: "Planck/WMAP/ACT combined",
            metric: "beta / sigma(beta)",
            pass_if: ">= 5.0",
            kill_if: "< 3.0",
            status: Status::Open,
            observed: "not-run-in-harness".to_string(),
        },
        TestCase {
            id: "T06",
            name: "Cosmic birefringence anisotropic amplitude",
            dataset: "Anisotropic CB maps",
            metric: "A_CB / sigma(A_CB)",
            pass_if: ">= 5.0",
            kill_if: "< 3.0",
            status: Status::Open,
            observed: "not-run-in-harness".to_string(),
        },
        TestCase {
            id: "T07",
            name: "DESI void scalar anomaly lane",
            dataset: "DESI DR1 void catalogs",
            metric: "z_void_scalar",
            pass_if: ">= 5.0",
            kill_if: "< 3.0",
            status: Status::Open,
            observed: "pending_quantitative_scalar".to_string(),
        },
        TestCase {
            id: "T08",
            name: "DESI void cross-check consistency",
            dataset: "Void counts + lensing",
            metric: "delta_param_between_channels",
            pass_if: "|delta| <= 1 sigma",
            kill_if: "|delta| > 3 sigma",
            status: Status::Open,
            observed: "not-run".to_string(),
        },
        TestCase {
            id: "T09",
            name: "EHT Kerr deviation residual",
            dataset: "M87* + Sgr A*",
            metric: "fractional_deviation / sigma",
            pass_if: ">= 5 sigma coherent anomaly",
            kill_if: "< 3 sigma",
            status: Status::Open,
            observed: "coarse-lane-only".to_string(),
        },
        TestCase {
            id: "T10",
            name: "EHT cross-object parameter consistency",
            dataset: "M87* and Sgr A* joint",
            metric: "shared_param_delta_sigma",
            pass_if: "<= 1 sigma",
            kill_if: "> 3 sigma",
            status: Status::Open,
            observed: "not-run".to_string(),
        },
        TestCase {
            id: "T11",
            name: "Weak-angle holdout forecast",
            dataset: "Chronological holdout EW points",
            metric: "holdout_reduced_chi2",
            pass_if: "<= 1.5",
            kill_if: "> 3.0",
            status: Status::Open,
            observed: "not-run".to_string(),
        },
        TestCase {
            id: "T12",
            name: "Weak-angle leave-one-experiment-out stability",
            dataset: "EW anchors",
            metric: "max_delta_pred_sigma_LOO",
            pass_if: "<= 1.0",
            kill_if: "> 3.0",
            status: Status::Open,
            observed: "not-run".to_string(),
        },
        TestCase {
            id: "T13",
            name: "Weak-angle scheme-coherent fit",
            dataset: "MSbar + effective-angle converted set",
            metric: "reduced_chi2_scheme_clean",
            pass_if: "<= 1.5",
            kill_if: "> 3.0",
            status: Status::Open,
            observed: "not-run".to_string(),
        },
        TestCase {
            id: "T14",
            name: "Weak-angle slope lock",
            dataset: "EW anchors with free log slope",
            metric: "|k_fit - alpha ln10/(4pi)| / sigma_k",
            pass_if: "<= 2 sigma",
            kill_if: "> 3 sigma",
            status: Status::Open,
            observed: "not-run".to_string(),
        },
        TestCase {
            id: "T15",
            name: "Topological-gain lock",
            dataset: "Topology gate",
            metric: "|G - 1|",
            pass_if: "<= 1e-12",
            kill_if: "> 1e-9",
            status: Status::Open,
            observed: "not-run-in-this-harness".to_string(),
        },
        TestCase {
            id: "T16",
            name: "No-door null replication",
            dataset: "Independent rerun of joint lanes",
            metric: "door_detected flag",
            pass_if: "false",
            kill_if: "true without reproducibility",
            status: Status::Open,
            observed: "pending".to_string(),
        },
        TestCase {
            id: "T17",
            name: "Analog photonic Pi-target lock",
            dataset: "Lab analog (photonic delay mesh)",
            metric: "distance_to_Pi_target",
            pass_if: "<= tolerance band",
            kill_if: "> 3x tolerance",
            status: Status::Open,
            observed: "not-run".to_string(),
        },
        TestCase {
            id: "T18",
            name: "Analog superconducting lock threshold",
            dataset: "Resonator loop",
            metric: "lock/no-lock transition sharpness",
            pass_if: "sharp threshold with predicted location",
            kill_if: "no threshold or displaced threshold",
            status: Status::Open,
            observed: "not-run".to_string(),
        },
        TestCase {
            id: "T19",
            name: "Analog RF-PLL causal-order anomaly test",
            dataset: "RF loop + timestamp challenge",
            metric: "predeparture anomaly rate",
            pass_if: "above prereg threshold",
            kill_if: "consistent with null",
            status: Status::Open,
            observed: "not-run".to_string(),
        },
        TestCase {
            id: "T20",
            name: "Global 4-lane consistency",
            dataset: "EW + CB + Voids + EHT",
            metric: "shared-param posterior overlap",
            pass_if: "non-empty overlap at 95% credible level",
            kill_if: "disjoint posteriors",
            status: Status::Open,
            observed: "not-run".to_string(),
        },
        TestCase {
            id: "T21",
            name: "Alpha structural cubic closure",
            dataset: "CODATA 2022",
            metric: "|alpha_inv_pred - alpha_inv_phys| (ppb)",
            pass_if: "<= 0.1 ppb",
            kill_if: "> 1.0 ppb",
            status: Status::Open,
            observed: "not-run".to_string(),
        },
        TestCase {
            id: "T22",
            name: "mp/me shared-series closure",
            dataset: "CODATA 2022",
            metric: "|mp/me_pred - mp/me_phys| (ppb)",
            pass_if: "<= 1.0 ppb",
            kill_if: "> 10.0 ppb",
            status: Status::Open,
            observed: "not-run".to_string(),
        },
    ];

    // Auto-populate the tests we already ran.
    if let Some(v) = red_base {
        tests[0].observed = format!("{v:.6}");
        tests[0].status = if v > 5.0 { Status::Pass } else { Status::Fail };
    }
    if let Some(v) = red_zero {
        tests[1].observed = format!("{v:.6}");
        tests[1].status = if v <= 1.5 { Status::Pass } else { Status::Fail };
    }
    if let Some(v) = pull_zero {
        tests[2].observed = format!("{v:.6}");
        tests[2].status = if v <= 2.5 { Status::Pass } else { Status::Fail };
    }
    if let (Some(z), Some(det)) = (z_stouffer, door_detected) {
        tests[3].observed = format!("z_stouffer={z:.6}, door_detected={det}");
        tests[3].status = if z >= 5.0 && det {
            Status::Pass
        } else if z < 3.0 && !det {
            Status::Fail
        } else {
            Status::Open
        };
    }

    // T07/T11/T13 imports from next3 run.
    if let Some(t13) = next3.get("tests").and_then(|t| t.get("T13_scheme_coherent_fit")) {
        if let (Some(v), Some(s)) = (
            t13.get("value").and_then(Value::as_f64),
            t13.get("status").and_then(Value::as_str),
        ) {
            tests[12].observed = format!("{v:.6}");
            tests[12].status = match s {
                "PASS" => Status::Pass,
                "FAIL" => Status::Fail,
                _ => Status::Open,
            };
        }
    }
    if let Some(t11) = next3
        .get("tests")
        .and_then(|t| t.get("T11_chronological_holdout"))
    {
        if let (Some(v), Some(s)) = (
            t11.get("value").and_then(Value::as_f64),
            t11.get("status").and_then(Value::as_str),
        ) {
            tests[10].observed = format!("{v:.6}");
            tests[10].status = match s {
                "PASS" => Status::Pass,
                "FAIL" => Status::Fail,
                _ => Status::Open,
            };
        }
    }
    if let Some(t07) = next3.get("tests").and_then(|t| t.get("T07_void_scalar_lane")) {
        if let (Some(v), Some(s)) = (
            t07.get("value").and_then(Value::as_f64),
            t07.get("status").and_then(Value::as_str),
        ) {
            tests[6].observed = format!("{v:.6}");
            tests[6].status = match s {
                "PASS" => Status::Pass,
                "FAIL" => Status::Fail,
                _ => Status::Open,
            };
        }
    }

    // T01/T06/T07/T08 imports from immediate-run artifact.
    if let Some(t01) = t01678.get("tests").and_then(|t| t.get("T01")) {
        if let Some(s) = t01.get("status").and_then(Value::as_str) {
            tests[0].status = match s {
                "PASS" => Status::Pass,
                "FAIL" => Status::Fail,
                _ => Status::Open,
            };
        }
        if let Some(obs) = t01.get("observed").and_then(Value::as_str) {
            tests[0].observed = obs.to_string();
        }
    }
    if let Some(t06) = t01678.get("tests").and_then(|t| t.get("T06")) {
        if let Some(s) = t06.get("status").and_then(Value::as_str) {
            tests[5].status = match s {
                "PASS" => Status::Pass,
                "FAIL" => Status::Fail,
                _ => Status::Open,
            };
        }
        if let Some(obs) = t06.get("observed").and_then(Value::as_str) {
            tests[5].observed = obs.to_string();
        }
    }
    if let Some(t07) = t01678.get("tests").and_then(|t| t.get("T07")) {
        if let Some(s) = t07.get("status").and_then(Value::as_str) {
            tests[6].status = match s {
                "PASS" => Status::Pass,
                "FAIL" => Status::Fail,
                _ => Status::Open,
            };
        }
        if let Some(obs) = t07.get("observed").and_then(Value::as_str) {
            tests[6].observed = obs.to_string();
        }
    }
    if let Some(t08) = t01678.get("tests").and_then(|t| t.get("T08")) {
        if let Some(s) = t08.get("status").and_then(Value::as_str) {
            tests[7].status = match s {
                "PASS" => Status::Pass,
                "FAIL" => Status::Fail,
                _ => Status::Open,
            };
        }
        if let Some(obs) = t08.get("observed").and_then(Value::as_str) {
            tests[7].observed = obs.to_string();
        }
    }

    // Remaining-open tests import.
    if let Some(tests_obj) = remaining12.get("tests") {
        let import_test = |slot: &mut TestCase, key: &str| {
            if let Some(t) = tests_obj.get(key) {
                if let Some(s) = t.get("status").and_then(Value::as_str) {
                    slot.status = match s {
                        "PASS" => Status::Pass,
                        "FAIL" => Status::Fail,
                        _ => Status::Open,
                    };
                }
                if let Some(obs) = t.get("observed").and_then(Value::as_str) {
                    slot.observed = obs.to_string();
                }
            }
        };
        import_test(&mut tests[4], "T05");
        import_test(&mut tests[7], "T08");
        import_test(&mut tests[8], "T09");
        import_test(&mut tests[9], "T10");
        import_test(&mut tests[11], "T12");
        import_test(&mut tests[13], "T14");
        import_test(&mut tests[14], "T15");
        import_test(&mut tests[15], "T16");
        import_test(&mut tests[16], "T17");
        import_test(&mut tests[17], "T18");
        import_test(&mut tests[18], "T19");
        import_test(&mut tests[19], "T20");
    }

    // T21/T22: corrected alpha and mp/me structural lanes.
    let alpha_inv_pred = solve_alpha_inv_structural();
    let alpha_ppb = ((alpha_inv_pred - ALPHA_INV_CODATA_2022).abs() / ALPHA_INV_CODATA_2022) * 1.0e9;
    tests[20].observed = format!(
        "alpha_inv_pred={alpha_inv_pred:.12}, alpha_inv_phys={ALPHA_INV_CODATA_2022:.12}, abs_ppb={alpha_ppb:.9}"
    );
    tests[20].status = if alpha_ppb <= 0.1 {
        Status::Pass
    } else if alpha_ppb > 1.0 {
        Status::Fail
    } else {
        Status::Open
    };

    let alpha_pred = 1.0 / alpha_inv_pred;
    let b = 9.0 + 5.0 / 32.0;
    let c = 1.0 / 125.0;
    let g = 4.0 - 8.0 * alpha_pred;
    let mp_me_pred = 6.0 * std::f64::consts::PI.powi(5)
        + 5.0 * alpha_pred
        - g * b * alpha_pred * alpha_pred
        - g * c * alpha_pred * alpha_pred * alpha_pred;
    let mp_ppb = ((mp_me_pred - MP_ME_CODATA_2022).abs() / MP_ME_CODATA_2022) * 1.0e9;
    tests[21].observed = format!(
        "mp_me_pred={mp_me_pred:.12}, mp_me_phys={MP_ME_CODATA_2022:.11}, abs_ppb={mp_ppb:.9}"
    );
    tests[21].status = if mp_ppb <= 1.0 {
        Status::Pass
    } else if mp_ppb > 10.0 {
        Status::Fail
    } else {
        Status::Open
    };

    let pass_n = tests.iter().filter(|t| matches!(t.status, Status::Pass)).count();
    let fail_n = tests.iter().filter(|t| matches!(t.status, Status::Fail)).count();
    let open_n = tests.iter().filter(|t| matches!(t.status, Status::Open)).count();

    let core = tests
        .iter()
        .filter(|t| is_core_algebraic(t.id))
        .collect::<Vec<_>>();
    let speculative = tests
        .iter()
        .filter(|t| !is_core_algebraic(t.id))
        .collect::<Vec<_>>();

    let core_pass = core.iter().filter(|t| matches!(t.status, Status::Pass)).count();
    let core_fail = core.iter().filter(|t| matches!(t.status, Status::Fail)).count();
    let core_open = core.iter().filter(|t| matches!(t.status, Status::Open)).count();

    let spec_pass = speculative
        .iter()
        .filter(|t| matches!(t.status, Status::Pass))
        .count();
    let spec_fail = speculative
        .iter()
        .filter(|t| matches!(t.status, Status::Fail))
        .count();
    let spec_open = speculative
        .iter()
        .filter(|t| matches!(t.status, Status::Open))
        .count();

    let tests_json = tests
        .iter()
        .map(|t| {
            json!({
              "id": t.id,
              "name": t.name,
              "dataset": t.dataset,
              "metric": t.metric,
              "pass_if": t.pass_if,
              "kill_if": t.kill_if,
              "category": if is_core_algebraic(t.id) { "core_algebraic" } else { "speculative_detection" },
              "status": t.status.as_str(),
              "observed": t.observed
            })
        })
        .collect::<Vec<_>>();

    let payload = json!({
      "scope": "rapid live-fire 20-test falsification harness",
      "inputs": {
        "weak_fit_json": weak_path,
        "door_fit_json": door_path,
        "next3_json": next3_path,
        "t01678_json": t01678_path,
        "remaining12_json": remaining12_path
      },
      "summary": {
        "pass": pass_n,
        "fail": fail_n,
        "open": open_n,
        "total": tests.len()
      },
      "core_algebraic_score": {
        "pass": core_pass,
        "fail": core_fail,
        "open": core_open,
        "total": core.len()
      },
      "speculative_detection_score": {
        "pass": spec_pass,
        "fail": spec_fail,
        "open": spec_open,
        "total": speculative.len()
      },
      "tests": tests_json
    });

    let txt_path = out.join("ctc_falsification_20.txt");
    let json_path = out.join("ctc_falsification_20.json");

    let mut txt = String::new();
    txt.push_str("[ctc_falsification_20]\n");
    txt.push_str("rapid live-fire falsification checklist (20 tests)\n\n");
    txt.push_str(&format!(
        "summary: PASS={} FAIL={} OPEN={} TOTAL={}\n\n",
        pass_n,
        fail_n,
        open_n,
        tests.len()
    ));
    txt.push_str(&format!(
        "core_algebraic_score: PASS={} FAIL={} OPEN={} TOTAL={}\n",
        core_pass, core_fail, core_open, core.len()
    ));
    txt.push_str(&format!(
        "speculative_detection_score: PASS={} FAIL={} OPEN={} TOTAL={}\n\n",
        spec_pass, spec_fail, spec_open, speculative.len()
    ));
    for t in &tests {
        txt.push_str(&format!(
            "{} [{}] {}\n  category={}\n  dataset={}\n  metric={}\n  pass_if={}\n  kill_if={}\n  observed={}\n\n",
            t.id,
            t.status.as_str(),
            t.name,
            if is_core_algebraic(t.id) {
                "core_algebraic"
            } else {
                "speculative_detection"
            },
            t.dataset,
            t.metric,
            t.pass_if,
            t.kill_if,
            t.observed
        ));
    }

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(
        &json_path,
        serde_json::to_string_pretty(&payload).expect("serialize json"),
    )
    .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "falsification_20 summary: PASS={} FAIL={} OPEN={}",
        pass_n, fail_n, open_n
    );
}
