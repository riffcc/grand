//! CTC retrocompute speedup probe.
//!
//! Goal:
//! - Quantify apparent wall-clock acceleration when computation results are
//!   returned predeparture through a CTC-like channel.
//!
//! Model:
//! - Finite hardware performs `task_flops / hardware_flops_per_s` seconds of work.
//! - A predeparture return fraction `r` (0..1) shifts result arrival backward by
//!   `r * compute_time`.
//! - Observed external latency is clamped to a finite instrumentation floor.
//! - Retry closure uses `p_eventual = 1 - (1-p)^n`.
//!
//! Interpretation:
//! - Internal compute cost remains finite and positive.
//! - External latency can be arbitrarily small as `r -> 1`.
//! - Apparent speedup can therefore be arbitrarily large on finite hardware.

use serde_json::json;
use std::fs;
use std::path::PathBuf;

fn env_f64(name: &str, default: f64) -> f64 {
    std::env::var(name)
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(default)
}

fn env_u64(name: &str, default: u64) -> u64 {
    std::env::var(name)
        .ok()
        .and_then(|s| s.parse::<u64>().ok())
        .unwrap_or(default)
}

fn clamp01(x: f64) -> f64 {
    x.clamp(0.0, 1.0)
}

fn eventual_success_prob(p: f64, n: u64) -> f64 {
    1.0 - (1.0 - p).powf(n as f64)
}

fn main() {
    let out_dir = std::env::var("GUTOE_RETROCOMPUTE_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ctc_retrocompute_speedup_probe".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    // Workload and hardware
    let task_flops = env_f64("GUTOE_RETROCOMPUTE_TASK_FLOPS", 1.0e20).max(1.0);
    let hardware_flops_per_s = env_f64("GUTOE_RETROCOMPUTE_HW_FLOPS_PER_S", 1.0e15).max(1.0);
    let compute_time_s = task_flops / hardware_flops_per_s;

    // Retro channel
    let predeparture_fraction = clamp01(env_f64("GUTOE_RETROCOMPUTE_PREDEP_FRAC", 0.999_999));
    let instrumentation_floor_s = env_f64("GUTOE_RETROCOMPUTE_FLOOR_S", 1.0e-9).max(1.0e-15);
    let shifted_time_s = compute_time_s * (1.0 - predeparture_fraction);
    let observed_latency_s = shifted_time_s.max(instrumentation_floor_s);
    let apparent_speedup = compute_time_s / observed_latency_s;

    // Retry closure
    let p_single_pass = clamp01(env_f64("GUTOE_RETROCOMPUTE_P_SINGLE", 0.12));
    let retry_depth = env_u64("GUTOE_RETROCOMPUTE_RETRY_DEPTH", 100).max(1);
    let p_eventual = eventual_success_prob(p_single_pass, retry_depth);
    let expected_attempts = if p_single_pass > 0.0 {
        (1.0 / p_single_pass).min(retry_depth as f64)
    } else {
        retry_depth as f64
    };
    let expected_internal_time_s = compute_time_s * expected_attempts;
    let effective_speedup_with_retries = expected_internal_time_s / observed_latency_s;

    // Sweep near r->1 to show unbounded trend (finite hardware, shrinking external latency).
    let sweep = [0.9, 0.99, 0.999, 0.9999, 0.99999, 0.999999]
        .iter()
        .map(|&r| {
            let lat = (compute_time_s * (1.0 - r)).max(instrumentation_floor_s);
            let su = compute_time_s / lat;
            json!({
              "predeparture_fraction": r,
              "observed_latency_s": lat,
              "apparent_speedup": su
            })
        })
        .collect::<Vec<_>>();

    let payload = json!({
      "inputs": {
        "task_flops": task_flops,
        "hardware_flops_per_s": hardware_flops_per_s,
        "predeparture_fraction": predeparture_fraction,
        "instrumentation_floor_s": instrumentation_floor_s,
        "p_single_pass": p_single_pass,
        "retry_depth": retry_depth
      },
      "core": {
        "compute_time_s_internal": compute_time_s,
        "observed_latency_s_external": observed_latency_s,
        "apparent_speedup": apparent_speedup
      },
      "retry_closure": {
        "eventual_success_prob": p_eventual,
        "expected_attempts_truncated": expected_attempts,
        "expected_internal_time_s": expected_internal_time_s,
        "effective_speedup_with_retries": effective_speedup_with_retries
      },
      "sweep_near_unity_predeparture": sweep,
      "scope": "simulation of retrocompute observability on finite hardware; not a physical engine claim"
    });

    let txt_path = out.join("ctc_retrocompute_speedup_probe.txt");
    let json_path = out.join("ctc_retrocompute_speedup_probe.json");

    let mut txt = String::new();
    txt.push_str("[ctc_retrocompute_speedup_probe]\n");
    txt.push_str("finite hardware + predeparture return channel\n\n");
    txt.push_str(&format!("compute_time_s_internal = {:.6e}\n", compute_time_s));
    txt.push_str(&format!("observed_latency_s_external = {:.6e}\n", observed_latency_s));
    txt.push_str(&format!("apparent_speedup = {:.6e}\n", apparent_speedup));
    txt.push_str(&format!("p_single_pass = {:.6}\n", p_single_pass));
    txt.push_str(&format!("retry_depth = {}\n", retry_depth));
    txt.push_str(&format!("eventual_success_prob = {:.12}\n", p_eventual));
    txt.push_str(&format!(
        "effective_speedup_with_retries = {:.6e}\n",
        effective_speedup_with_retries
    ));

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(
        &json_path,
        serde_json::to_string_pretty(&payload).expect("json"),
    )
    .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "retrocompute: speedup={:.3e}, eventual_p={:.9}",
        apparent_speedup, p_eventual
    );
}
