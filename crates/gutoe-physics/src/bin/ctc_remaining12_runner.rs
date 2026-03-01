//! Execute remaining-open falsification tests where current artifacts permit closure:
//! T05, T08, T09, T10, T12, T14, T15, T16, T17, T18, T19, T20.

use serde_json::{json, Value};
use std::fs;
use std::path::PathBuf;

const MZ_GEV: f64 = 91.1876;
const ALPHA_INV: f64 = 137.035_999_177;
const A_RATIONAL: f64 = 3.0 / 13.0 + 1.0 / (13.0 * 13.0 * 13.0);

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

fn read_json(path: &str) -> Option<Value> {
    let raw = fs::read_to_string(path).ok()?;
    serde_json::from_str::<Value>(&raw).ok()
}

fn find_lane<'a>(v: &'a Value, name: &str) -> Option<&'a Value> {
    let lanes = v.get("lanes")?.as_array()?;
    lanes.iter().find(|l| l.get("name").and_then(Value::as_str) == Some(name))
}

fn has_quantified_void_lensing_pair(v: &Value) -> bool {
    let Some(lanes) = v.get("lanes").and_then(Value::as_array) else {
        return false;
    };
    let mut void_numeric = false;
    let mut lensing_numeric = false;
    for lane in lanes {
        let name = lane
            .get("name")
            .and_then(Value::as_str)
            .unwrap_or_default()
            .to_ascii_lowercase();
        let has_numeric = lane.get("z").and_then(Value::as_f64).is_some();
        if has_numeric && name.contains("void") {
            void_numeric = true;
        }
        if has_numeric && (name.contains("lensing") || name.contains("kappa")) {
            lensing_numeric = true;
        }
    }
    void_numeric && lensing_numeric
}

fn main() {
    let out_dir = std::env::var("GUTOE_REMAINING12_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ctc_remaining12_runner".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let weak_path = std::env::var("GUTOE_WEAK_FIT_JSON").unwrap_or_else(|_| {
        "/tmp/bh_renders/ctc_weak_angle_identity_fit/ctc_weak_angle_identity_fit.json".to_string()
    });
    let door_path = std::env::var("GUTOE_DOOR_FIT_JSON").unwrap_or_else(|_| {
        "/tmp/bh_renders/ctc_public_door_joint_fit/ctc_public_door_joint_fit.json".to_string()
    });
    let door_rerun_path = std::env::var("GUTOE_DOOR_FIT_RERUN_JSON").unwrap_or_else(|_| {
        "/tmp/bh_renders/ctc_public_door_joint_fit_rerun/ctc_public_door_joint_fit.json".to_string()
    });
    let topo_path = std::env::var("GUTOE_TOPO_JSON").unwrap_or_else(|_| {
        "/tmp/bh_renders/topology_overdetermination_probe/topology_overdetermination_probe.json".to_string()
    });

    let weak = read_json(&weak_path).unwrap_or_else(|| json!({}));
    let door = read_json(&door_path).unwrap_or_else(|| json!({}));
    let door_rerun = read_json(&door_rerun_path).unwrap_or_else(|| json!({}));
    let topo = read_json(&topo_path).unwrap_or_else(|| json!({}));

    // T05: isotropic CB non-zero
    let t05 = if let Some(lane) = find_lane(&door, "cosmic_birefringence_isotropic_beta") {
        if let Some(z) = lane.get("z").and_then(Value::as_f64) {
            let status = if z >= 5.0 {
                Status::Pass
            } else if z < 3.0 {
                Status::Fail
            } else {
                Status::Open
            };
            (status, format!("beta/sigma={z:.6}"))
        } else {
            (Status::Open, "isotropic lane present but no numeric z".to_string())
        }
    } else {
        (
            Status::Open,
            "blocked: isotropic beta lane not quantified in current joint artifact".to_string(),
        )
    };

    // T08: void+lensing cross-check
    let t08 = if has_quantified_void_lensing_pair(&door) {
        (
            Status::Open,
            "quantified pair present; delta metric not computed in this coarse runner".to_string(),
        )
    } else {
        (
            Status::Open,
            "blocked: no quantified void+lensing scalar pair in current artifact".to_string(),
        )
    };

    // T09: EHT residual
    let t09 = if let Some(lane) = find_lane(&door, "eht_shadow_fractional_deviation_from_kerr") {
        if let Some(z) = lane.get("z").and_then(Value::as_f64) {
            let status = if z >= 5.0 {
                Status::Pass
            } else if z < 3.0 {
                Status::Fail
            } else {
                Status::Open
            };
            (status, format!("fractional_deviation/sigma={z:.6}"))
        } else {
            (Status::Open, "EHT lane present but no numeric z".to_string())
        }
    } else {
        (Status::Open, "missing EHT lane in joint artifact".to_string())
    };

    // T10: cross-object consistency (needs separate M87/SgrA inferred parameter rows)
    let t10 = (
        Status::Open,
        "blocked: no separate M87*/SgrA* parameter rows in current coarse artifact".to_string(),
    );

    // T12/T14 from weak artifact points.
    let t12_t14 = if let Some(points) = weak.get("points").and_then(Value::as_array) {
        let mut rows = Vec::new();
        for p in points {
            let q = p.get("Q_GeV").and_then(Value::as_f64).unwrap_or(MZ_GEV);
            let y = p.get("value").and_then(Value::as_f64).unwrap_or(0.0);
            let s = p.get("sigma").and_then(Value::as_f64).unwrap_or(1.0);
            let name = p
                .get("name")
                .and_then(Value::as_str)
                .unwrap_or("unknown")
                .to_string();
            let x = (MZ_GEV / q).ln();
            rows.push((name, x, y, s));
        }

        let den: f64 = rows.iter().map(|(_, x, _, s)| x * x / (s * s)).sum();
        if den <= 0.0 {
            (
                (Status::Open, "blocked: insufficient nonzero-Q leverage for LOO".to_string()),
                (Status::Open, "blocked: insufficient nonzero-Q leverage for slope".to_string()),
            )
        } else {
            let num: f64 = rows
                .iter()
                .map(|(_, x, y, s)| x * (y - A_RATIONAL) / (s * s))
                .sum();
            let k_fit = num / den;
            let sigma_k = (1.0 / den).sqrt();
            let k_exp = (1.0 / ALPHA_INV) * 10.0_f64.ln() / (4.0 * std::f64::consts::PI);
            let z_k = ((k_fit - k_exp) / sigma_k).abs();
            let t14 = (
                if z_k <= 2.0 {
                    Status::Pass
                } else if z_k > 3.0 {
                    Status::Fail
                } else {
                    Status::Open
                },
                format!("k_fit={k_fit:.12e}, k_expected={k_exp:.12e}, sigma_k={sigma_k:.12e}, z={z_k:.6}"),
            );

            // LOO stability: max sigma shift in predictions across all well-defined folds.
            let full_preds: Vec<f64> = rows.iter().map(|(_, x, _, _)| A_RATIONAL + k_fit * x).collect();
            let mut max_delta_sigma = 0.0_f64;
            let mut usable_folds = 0_usize;

            for j in 0..rows.len() {
                let den2: f64 = rows
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| *i != j)
                    .map(|(_, (_, x, _, s))| x * x / (s * s))
                    .sum();
                if den2 <= 0.0 {
                    continue;
                }
                usable_folds += 1;
                let num2: f64 = rows
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| *i != j)
                    .map(|(_, (_, x, y, s))| x * (y - A_RATIONAL) / (s * s))
                    .sum();
                let k2 = num2 / den2;

                for (i, (_, x, _, s)) in rows.iter().enumerate() {
                    let pred2 = A_RATIONAL + k2 * x;
                    let delta_sig = ((pred2 - full_preds[i]) / s).abs();
                    if delta_sig > max_delta_sigma {
                        max_delta_sigma = delta_sig;
                    }
                }
            }

            let t12 = if usable_folds == 0 {
                (
                    Status::Open,
                    "blocked: all LOO folds singular under current Q-leverage".to_string(),
                )
            } else {
                (
                    if max_delta_sigma <= 1.0 {
                        Status::Pass
                    } else if max_delta_sigma > 3.0 {
                        Status::Fail
                    } else {
                        Status::Open
                    },
                    format!(
                        "max_delta_pred_sigma_LOO={max_delta_sigma:.6}, usable_folds={usable_folds}"
                    ),
                )
            };

            (t12, t14)
        }
    } else {
        (
            (Status::Open, "missing weak points in artifact".to_string()),
            (Status::Open, "missing weak points in artifact".to_string()),
        )
    };

    let (t12, t14) = t12_t14;

    // T15: topology gain lock.
    let t15 = if let Some(abs_resid) = topo
        .get("closure")
        .and_then(|c| c.get("abs_residual"))
        .and_then(Value::as_f64)
    {
        (
            if abs_resid <= 1e-12 {
                Status::Pass
            } else if abs_resid > 1e-9 {
                Status::Fail
            } else {
                Status::Open
            },
            format!("|G-1|={abs_resid:.12e}"),
        )
    } else {
        (Status::Open, "missing topology closure residual".to_string())
    };

    // T16: no-door null replication.
    let t16 = {
        let d0 = door
            .get("joint_fit")
            .and_then(|j| j.get("door_detected"))
            .and_then(Value::as_bool);
        let d1 = door_rerun
            .get("joint_fit")
            .and_then(|j| j.get("door_detected"))
            .and_then(Value::as_bool);
        match (d0, d1) {
            (Some(a), Some(b)) => {
                let status = if !a && !b {
                    Status::Pass
                } else if a || b {
                    Status::Fail
                } else {
                    Status::Open
                };
                (status, format!("door_detected run0={a}, run1={b}"))
            }
            _ => (
                Status::Open,
                "blocked: missing one or both door-detection artifacts for replication".to_string(),
            ),
        }
    };

    // T17/T18/T19/T20 remain blocked by absent lab/global posterior artifacts in current local set.
    let t17 = (Status::Open, "blocked: no photonic Pi-target lab artifact in current runset".to_string());
    let t18 = (Status::Open, "blocked: no superconducting lock-threshold artifact in current runset".to_string());
    let t19 = (Status::Open, "blocked: no RF-PLL causal-order artifact in current runset".to_string());
    let t20 = (
        Status::Open,
        "blocked: missing quantified posterior overlap across EW+CB+Voids+EHT lanes".to_string(),
    );

    let payload = json!({
      "scope": "remaining-open tests runner",
      "tests": {
        "T05": {"status": t05.0.as_str(), "observed": t05.1},
        "T08": {"status": t08.0.as_str(), "observed": t08.1},
        "T09": {"status": t09.0.as_str(), "observed": t09.1},
        "T10": {"status": t10.0.as_str(), "observed": t10.1},
        "T12": {"status": t12.0.as_str(), "observed": t12.1},
        "T14": {"status": t14.0.as_str(), "observed": t14.1},
        "T15": {"status": t15.0.as_str(), "observed": t15.1},
        "T16": {"status": t16.0.as_str(), "observed": t16.1},
        "T17": {"status": t17.0.as_str(), "observed": t17.1},
        "T18": {"status": t18.0.as_str(), "observed": t18.1},
        "T19": {"status": t19.0.as_str(), "observed": t19.1},
        "T20": {"status": t20.0.as_str(), "observed": t20.1}
      }
    });

    let txt_path = out.join("ctc_remaining12_runner.txt");
    let json_path = out.join("ctc_remaining12_runner.json");

    let mut txt = String::new();
    txt.push_str("[ctc_remaining12_runner]\n");
    txt.push_str("remaining-open tests execution\n\n");
    for key in ["T05", "T08", "T09", "T10", "T12", "T14", "T15", "T16", "T17", "T18", "T19", "T20"] {
        if let Some(t) = payload.get("tests").and_then(|x| x.get(key)) {
            let s = t.get("status").and_then(Value::as_str).unwrap_or("OPEN");
            let o = t.get("observed").and_then(Value::as_str).unwrap_or("-");
            txt.push_str(&format!("{key} status={s} observed={o}\n"));
        }
    }

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(
        &json_path,
        serde_json::to_string_pretty(&payload).expect("serialize"),
    )
    .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
}
