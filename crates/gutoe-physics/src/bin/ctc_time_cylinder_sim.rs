//! CTC simulation on a 1+1 Minkowski time-cylinder.
//!
//! Worldline parameterization:
//!   t(λ) = t0 + T λ
//!   x(λ) = x0 + A sin(2π λ)
//! with λ ∈ [0, 1], periodic identification t ~ t + nT.
//!
//! In the covering space, local segment interval is:
//!   ds² = -(dt)² + (dx)²
//! and we require ds² < 0 (timelike).
//! Endpoint closure is evaluated both in coordinate space and under
//! time-cylinder identification.

use serde_json::json;
use std::f64::consts::PI;
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

#[derive(Debug, Clone, Copy)]
struct Event {
    t: f64,
    x: f64,
}

fn interval_sq(a: Event, b: Event) -> f64 {
    let dt = b.t - a.t;
    let dx = b.x - a.x;
    -(dt * dt) + (dx * dx)
}

fn env_f64(name: &str, default: f64) -> f64 {
    std::env::var(name)
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(default)
}

fn env_usize(name: &str, default: usize) -> usize {
    std::env::var(name)
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(default)
}

fn main() {
    let out_dir = std::env::var("GUTOE_CTC_SIM_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ctc_time_cylinder".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let t0 = env_f64("GUTOE_CTC_T0", 0.0);
    let x0 = env_f64("GUTOE_CTC_X0", 0.0);
    let period_t = env_f64("GUTOE_CTC_PERIOD_T", 1.0);
    let amp_x = env_f64("GUTOE_CTC_AMP_X", 0.0);
    let steps = env_usize("GUTOE_CTC_STEPS", 200).max(2);

    let mut worldline = Vec::with_capacity(steps + 1);
    for i in 0..=steps {
        let lambda = i as f64 / steps as f64;
        let t = t0 + period_t * lambda;
        let x = x0 + amp_x * (2.0 * PI * lambda).sin();
        worldline.push(Event { t, x });
    }

    let mut proper_time = 0.0_f64;
    let mut timelike_all = true;
    let mut ds2_max = f64::NEG_INFINITY;
    let mut ds2_min = f64::INFINITY;
    let mut violating_segments = 0usize;

    for i in 0..steps {
        let a = worldline[i];
        let b = worldline[i + 1];
        let ds2 = interval_sq(a, b);
        if ds2 >= 0.0 {
            timelike_all = false;
            violating_segments += 1;
        } else {
            proper_time += (-ds2).sqrt();
        }
        ds2_max = ds2_max.max(ds2);
        ds2_min = ds2_min.min(ds2);
    }

    let start = worldline.first().copied().expect("start");
    let end = worldline.last().copied().expect("end");

    let dx_total = end.x - start.x;
    let dt_total = end.t - start.t;
    let coord_closed = dx_total.abs() < 1e-12 && dt_total.abs() < 1e-12;

    let winding_est = dt_total / period_t;
    let winding_rounded = winding_est.round();
    let winding_integer = (winding_est - winding_rounded).abs() < 1e-9;
    let identified_closed = dx_total.abs() < 1e-9 && winding_integer && winding_rounded.abs() >= 1.0;

    let csv_path = out.join("ctc_worldline.csv");
    let txt_path = out.join("ctc_sim_report.txt");
    let json_path = out.join("ctc_sim_report.json");

    let mut csv = File::create(&csv_path).expect("create csv");
    writeln!(csv, "i,lambda,t,x").expect("write csv header");
    for (i, e) in worldline.iter().enumerate() {
        let lambda = i as f64 / steps as f64;
        writeln!(csv, "{i},{lambda:.12},{:.12e},{:.12e}", e.t, e.x).expect("write csv row");
    }

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[inputs]").expect("write");
    writeln!(txt, "t0 = {:.12e}", t0).expect("write");
    writeln!(txt, "x0 = {:.12e}", x0).expect("write");
    writeln!(txt, "period_t = {:.12e}", period_t).expect("write");
    writeln!(txt, "amp_x = {:.12e}", amp_x).expect("write");
    writeln!(txt, "steps = {}", steps).expect("write");
    writeln!(txt).expect("write");

    writeln!(txt, "[local_causality]").expect("write");
    writeln!(txt, "timelike_all_segments = {}", timelike_all).expect("write");
    writeln!(txt, "violating_segments = {}", violating_segments).expect("write");
    writeln!(txt, "segment_ds2_min = {:.12e}", ds2_min).expect("write");
    writeln!(txt, "segment_ds2_max = {:.12e}", ds2_max).expect("write");
    writeln!(txt, "proper_time_sum = {:.12e}", proper_time).expect("write");
    writeln!(txt).expect("write");

    writeln!(txt, "[closure]").expect("write");
    writeln!(txt, "dx_total = {:.12e}", dx_total).expect("write");
    writeln!(txt, "dt_total = {:.12e}", dt_total).expect("write");
    writeln!(txt, "coordinate_closed = {}", coord_closed).expect("write");
    writeln!(txt, "winding_estimate = {:.12e}", winding_est).expect("write");
    writeln!(txt, "winding_integer = {}", winding_integer).expect("write");
    writeln!(txt, "identified_closed = {}", identified_closed).expect("write");

    let payload = json!({
        "inputs": {
            "t0": t0,
            "x0": x0,
            "period_t": period_t,
            "amp_x": amp_x,
            "steps": steps
        },
        "local_causality": {
            "timelike_all_segments": timelike_all,
            "violating_segments": violating_segments,
            "segment_ds2_min": ds2_min,
            "segment_ds2_max": ds2_max,
            "proper_time_sum": proper_time
        },
        "closure": {
            "dx_total": dx_total,
            "dt_total": dt_total,
            "coordinate_closed": coord_closed,
            "winding_estimate": winding_est,
            "winding_integer": winding_integer,
            "identified_closed": identified_closed
        }
    });
    fs::write(
        &json_path,
        serde_json::to_string_pretty(&payload).expect("encode json"),
    )
    .expect("write json");

    println!("wrote {}", csv_path.display());
    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "timelike_all={} identified_closed={} winding≈{:.3}",
        timelike_all, identified_closed, winding_est
    );
}

