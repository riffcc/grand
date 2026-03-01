//! Dynamic topology creation probe (Path B gate scaffold).
//!
//! Gate model:
//!   threshold = (3/16) * |R| * |T|
//!   pass iff period T > 0 and budget >= threshold
//!
//! Operational toy metric:
//!   local timelike travel to/from patch + n loops of period T,
//!   effective coordinate-time shift = -n T.

use serde_json::json;
use std::fs;
use std::path::PathBuf;

const VOID_FRACTION: f64 = 3.0 / 16.0;

#[derive(Clone)]
struct Case {
    name: &'static str,
    start_x: f64,
    goal_x: f64,
    radius_r: f64,
    period_t: f64,
    budget: f64,
    n_max: usize,
}

fn abs(x: f64) -> f64 {
    x.abs()
}

fn threshold(radius_r: f64, period_t: f64) -> f64 {
    VOID_FRACTION * abs(radius_r) * abs(period_t)
}

fn main() {
    let out_dir = std::env::var("GUTOE_DYNAMIC_TOPOLOGY_CREATION_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/dynamic_topology_creation_probe".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let cases = vec![
        Case {
            name: "under_budget",
            start_x: -500.0,
            goal_x: 500.0,
            radius_r: 20.0,
            period_t: 50.0,
            budget: 100.0,
            n_max: 40,
        },
        Case {
            name: "at_budget",
            start_x: -500.0,
            goal_x: 500.0,
            radius_r: 20.0,
            period_t: 50.0,
            budget: threshold(20.0, 50.0),
            n_max: 40,
        },
        Case {
            name: "over_budget",
            start_x: -1200.0,
            goal_x: 1500.0,
            radius_r: 40.0,
            period_t: 120.0,
            budget: 20_000.0,
            n_max: 80,
        },
    ];

    let mut rows = Vec::new();

    for c in cases {
        let thr = threshold(c.radius_r, c.period_t);
        let gate_pass = c.period_t > 0.0 && c.budget >= thr;

        let baseline_light_time = abs(c.goal_x - c.start_x);

        // Entry/exit via nearest patch boundary (toy local-timelike routing).
        let x_entry = if c.start_x <= 0.0 { -c.radius_r } else { c.radius_r };
        let dt_in = abs(c.start_x - x_entry);
        let dt_out = abs(c.goal_x - x_entry);

        let mut best_eff = f64::INFINITY;
        let mut best_prop = f64::INFINITY;
        let mut best_n = 0usize;

        for n in 0..=c.n_max {
            let loops = n as f64;
            let proper = dt_in + loops * c.period_t + dt_out;
            let eff = dt_in - loops * c.period_t + dt_out;
            if eff < best_eff {
                best_eff = eff;
                best_prop = proper;
                best_n = n;
            }
        }

        let effective_superluminal = gate_pass && best_eff < baseline_light_time;
        let pre_departure = gate_pass && best_eff < 0.0;

        rows.push(json!({
            "case": c.name,
            "start_x": c.start_x,
            "goal_x": c.goal_x,
            "radius_r": c.radius_r,
            "period_t": c.period_t,
            "budget": c.budget,
            "threshold": thr,
            "gate_pass": gate_pass,
            "baseline_light_time": baseline_light_time,
            "best_n": best_n,
            "best_effective_arrival": best_eff,
            "best_local_proper_time": best_prop,
            "effective_superluminal": effective_superluminal,
            "pre_departure": pre_departure
        }));
    }

    let payload = json!({
        "model": {
            "void_fraction": VOID_FRACTION,
            "gate_threshold": "(3/16)*|R|*|T|",
            "operational_assumptions": [
                "local travel timelike",
                "loop identification gives effective coordinate shift -T",
                "no claim of physical implementation mechanism"
            ]
        },
        "results": rows
    });

    let txt_path = out.join("dynamic_topology_creation_probe.txt");
    let json_path = out.join("dynamic_topology_creation_probe.json");

    let mut txt = String::new();
    txt.push_str("[dynamic_topology_creation_probe]\n");
    txt.push_str("threshold = (3/16)*|R|*|T|\n\n");
    for r in payload["results"].as_array().expect("array") {
        txt.push_str(&format!(
            "{}: gate={}, threshold={:.6e}, budget={:.6e}, baseline={:.6e}, best_eff={:.6e}, best_n={}, superluminal={}, pre_departure={}\n",
            r["case"].as_str().unwrap_or("case"),
            r["gate_pass"].as_bool().unwrap_or(false),
            r["threshold"].as_f64().unwrap_or(f64::NAN),
            r["budget"].as_f64().unwrap_or(f64::NAN),
            r["baseline_light_time"].as_f64().unwrap_or(f64::NAN),
            r["best_effective_arrival"].as_f64().unwrap_or(f64::NAN),
            r["best_n"].as_u64().unwrap_or(0),
            r["effective_superluminal"].as_bool().unwrap_or(false),
            r["pre_departure"].as_bool().unwrap_or(false)
        ));
    }

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(&json_path, serde_json::to_string_pretty(&payload).expect("json")).expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
}
