//! Execute immediate falsification checks requested for tests T01, T06, T07, T08.
//!
//! Inputs (existing artifacts):
//! - weak-angle fit JSON
//! - coarse public door joint-fit JSON
//! - next3 JSON (for T07 if already quantified)
//!
//! Output:
//! - /tmp/bh_renders/ctc_t01_t06_t07_t08/ctc_t01_t06_t07_t08.{txt,json}

use serde_json::{json, Value};
use std::fs;
use std::path::PathBuf;

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

fn weak_model_reduced_chi2(v: &Value, model_name: &str) -> Option<f64> {
    v.get("model_fits")?
        .as_array()?
        .iter()
        .find(|m| m.get("model").and_then(Value::as_str) == Some(model_name))
        .and_then(|m| m.get("reduced_chi2"))
        .and_then(Value::as_f64)
}

fn find_lane_z(v: &Value, lane_name: &str) -> Option<f64> {
    let lanes = v.get("lanes")?.as_array()?;
    for lane in lanes {
        if lane.get("name").and_then(Value::as_str) == Some(lane_name) {
            return lane.get("z").and_then(Value::as_f64);
        }
    }
    None
}

fn has_quantified_void_lensing_pair(v: &Value) -> bool {
    // Current coarse artifact generally carries one DESI void lane without scalar z.
    // We require two quantitative channels for T08: void + lensing (both numeric).
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
    let out_dir = std::env::var("GUTOE_T01678_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ctc_t01_t06_t07_t08".to_string());
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

    let weak = read_json(&weak_path).unwrap_or_else(|| json!({}));
    let door = read_json(&door_path).unwrap_or_else(|| json!({}));
    let next3 = read_json(&next3_path).unwrap_or_else(|| json!({}));

    // T01
    let t01_val = weak_model_reduced_chi2(&weak, "base_fixed");
    let (t01_status, t01_obs) = if let Some(v) = t01_val {
        let s = if v > 5.0 { Status::Pass } else { Status::Fail };
        (s, format!("reduced_chi2(base_fixed)={v:.6}"))
    } else {
        (Status::Open, "missing weak-fit base_fixed model".to_string())
    };

    // T06
    let t06_z = find_lane_z(&door, "cosmic_birefringence_anisotropic_A_CB");
    let (t06_status, t06_obs) = if let Some(z) = t06_z {
        let s = if z >= 5.0 {
            Status::Pass
        } else if z < 3.0 {
            Status::Fail
        } else {
            Status::Open
        };
        (s, format!("A_CB/sigma(A_CB)={z:.6}"))
    } else {
        (
            Status::Open,
            "missing quantified anisotropic birefringence lane".to_string(),
        )
    };

    // T07
    let t07 = next3
        .get("tests")
        .and_then(|t| t.get("T07_void_scalar_lane"));
    let (t07_status, t07_obs) = if let Some(v) = t07 {
        let z = v.get("value").and_then(Value::as_f64);
        let s = v.get("status").and_then(Value::as_str);
        if let (Some(z), Some(s)) = (z, s) {
            let status = match s {
                "PASS" => Status::Pass,
                "FAIL" => Status::Fail,
                _ => Status::Open,
            };
            (status, format!("z_void_scalar={z:.6}"))
        } else {
            (
                Status::Open,
                "next3 T07 present but malformed value/status".to_string(),
            )
        }
    } else {
        (
            Status::Open,
            "missing quantified T07 lane in next3 artifact".to_string(),
        )
    };

    // T08
    let has_pair = has_quantified_void_lensing_pair(&door);
    let (t08_status, t08_obs) = if has_pair {
        (
            Status::Open,
            "quantified void+lensing pair present; delta metric not yet computed in this runner"
                .to_string(),
        )
    } else {
        (
            Status::Open,
            "blocked: no quantified void+lensing scalar pair in current public artifact".to_string(),
        )
    };

    let payload = json!({
      "scope": "immediate execution for tests T01, T06, T07, T08",
      "inputs": {
        "weak_fit_json": weak_path,
        "door_fit_json": door_path,
        "next3_json": next3_path
      },
      "tests": {
        "T01": {
          "metric": "reduced_chi2(base_fixed)",
          "pass_if": "> 5.0",
          "kill_if": "<= 5.0",
          "status": t01_status.as_str(),
          "observed": t01_obs
        },
        "T06": {
          "metric": "A_CB / sigma(A_CB)",
          "pass_if": ">= 5.0",
          "kill_if": "< 3.0",
          "status": t06_status.as_str(),
          "observed": t06_obs
        },
        "T07": {
          "metric": "z_void_scalar",
          "pass_if": ">= 5.0",
          "kill_if": "< 3.0",
          "status": t07_status.as_str(),
          "observed": t07_obs
        },
        "T08": {
          "metric": "delta_param_between_channels",
          "pass_if": "|delta| <= 1 sigma",
          "kill_if": "|delta| > 3 sigma",
          "status": t08_status.as_str(),
          "observed": t08_obs
        }
      }
    });

    let txt_path = out.join("ctc_t01_t06_t07_t08.txt");
    let json_path = out.join("ctc_t01_t06_t07_t08.json");

    let mut txt = String::new();
    txt.push_str("[ctc_t01_t06_t07_t08]\n");
    txt.push_str("immediate execution for tests T01, T06, T07, T08\n\n");
    txt.push_str(&format!("T01 status={} observed={}\n", t01_status.as_str(), t01_obs));
    txt.push_str(&format!("T06 status={} observed={}\n", t06_status.as_str(), t06_obs));
    txt.push_str(&format!("T07 status={} observed={}\n", t07_status.as_str(), t07_obs));
    txt.push_str(&format!("T08 status={} observed={}\n", t08_status.as_str(), t08_obs));

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(
        &json_path,
        serde_json::to_string_pretty(&payload).expect("serialize json"),
    )
    .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
}
