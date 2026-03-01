//! Compact local-topology creation probe.
//!
//! Toy operational model:
//! - Local patch support: x in [-R, R]
//! - Timelike travel to patch entry and out to destination (|dx| <= dt, c=1)
//! - Inside patch: n identified loops, each giving effective coordinate shift -T
//!   while consuming positive local proper-time.
//!
//! Reports whether finite-support creation can yield effective superluminal or
//! pre-departure coordinate arrivals in this toy quotient model.

use serde_json::json;
use std::fs;
use std::path::PathBuf;

#[derive(Clone)]
struct Case {
    name: &'static str,
    x_start: f64,
    x_goal: f64,
    patch_r: f64,
    period_t: f64,
    n_max: usize,
    entry_steps: usize,
}

fn abs(x: f64) -> f64 {
    x.abs()
}

fn main() {
    let out_dir = std::env::var("GUTOE_CTC_LOCAL_PATCH_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ctc_local_patch_creation_probe".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let cases = vec![
        Case {
            name: "through_patch_large_span",
            x_start: -1000.0,
            x_goal: 1000.0,
            patch_r: 20.0,
            period_t: 100.0,
            n_max: 60,
            entry_steps: 401,
        },
        Case {
            name: "same_side_far_patch",
            x_start: 500.0,
            x_goal: 700.0,
            patch_r: 20.0,
            period_t: 100.0,
            n_max: 60,
            entry_steps: 401,
        },
        Case {
            name: "moderate_patch_high_period",
            x_start: -300.0,
            x_goal: 450.0,
            patch_r: 50.0,
            period_t: 250.0,
            n_max: 50,
            entry_steps: 401,
        },
    ];

    let mut rows = Vec::new();

    for c in cases {
        let baseline_light_time = abs(c.x_goal - c.x_start);

        let mut best_effective_arrival = f64::INFINITY;
        let mut best_proper_time = f64::INFINITY;
        let mut best_entry = 0.0;
        let mut best_n = 0usize;

        for k in 0..c.entry_steps {
            let alpha = if c.entry_steps <= 1 {
                0.0
            } else {
                k as f64 / (c.entry_steps - 1) as f64
            };
            let x_entry = -c.patch_r + 2.0 * c.patch_r * alpha;

            let dt_in = abs(x_entry - c.x_start);
            let dt_out = abs(c.x_goal - x_entry);

            for n in 0..=c.n_max {
                let loops = n as f64;

                // Local-proper-time accounting: all positive pieces.
                let proper_time = dt_in + loops * c.period_t + dt_out;

                // Effective coordinate-time arrival in quotient view.
                let eff_arrival = dt_in - loops * c.period_t + dt_out;

                if eff_arrival < best_effective_arrival {
                    best_effective_arrival = eff_arrival;
                    best_proper_time = proper_time;
                    best_entry = x_entry;
                    best_n = n;
                }
            }
        }

        let effective_superluminal = best_effective_arrival < baseline_light_time;
        let pre_departure_arrival = best_effective_arrival < 0.0;

        rows.push(json!({
            "case": c.name,
            "x_start": c.x_start,
            "x_goal": c.x_goal,
            "patch_r": c.patch_r,
            "period_t": c.period_t,
            "n_max": c.n_max,
            "baseline_light_time": baseline_light_time,
            "best_entry_x": best_entry,
            "best_n": best_n,
            "best_effective_arrival_time": best_effective_arrival,
            "best_local_proper_time": best_proper_time,
            "effective_superluminal": effective_superluminal,
            "pre_departure_arrival": pre_departure_arrival
        }));
    }

    let payload = json!({
        "model": {
            "type": "compact_local_identification_patch",
            "assumptions": [
                "local travel timelike (c=1)",
                "identification active only for |x|<=R",
                "each patch loop contributes effective coordinate shift -T"
            ]
        },
        "results": rows
    });

    let txt_path = out.join("ctc_local_patch_creation_probe.txt");
    let json_path = out.join("ctc_local_patch_creation_probe.json");

    let mut txt = String::new();
    txt.push_str("[ctc_local_patch_creation_probe]\n");
    txt.push_str("model = compact_local_identification_patch\n");
    txt.push_str("\n[cases]\n");
    for r in payload["results"].as_array().expect("array") {
        txt.push_str(&format!(
            "{}: baseline={:.12e}, best_eff={:.12e}, best_proper={:.12e}, best_n={}, superluminal={}, pre_departure={}\n",
            r["case"].as_str().unwrap_or("case"),
            r["baseline_light_time"].as_f64().unwrap_or(f64::NAN),
            r["best_effective_arrival_time"].as_f64().unwrap_or(f64::NAN),
            r["best_local_proper_time"].as_f64().unwrap_or(f64::NAN),
            r["best_n"].as_u64().unwrap_or(0),
            r["effective_superluminal"].as_bool().unwrap_or(false),
            r["pre_departure_arrival"].as_bool().unwrap_or(false)
        ));
    }

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(&json_path, serde_json::to_string_pretty(&payload).expect("json")).expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
}
