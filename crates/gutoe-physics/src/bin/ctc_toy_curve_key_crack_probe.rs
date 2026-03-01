//! Toy ECC 8-bit key crack probe with retrocompute latency model.
//!
//! What this does:
//! - Defines a tiny pedagogical elliptic curve over a small prime field.
//! - Samples many 8-bit private keys and recovers them by brute-force discrete log.
//! - Measures normal wall-clock cracking time.
//! - Applies a retrocompute observability model to compare apparent external time.
//!
//! Scope:
//! - Educational toy model only (8-bit key space).
//! - Not a claim about practical cryptanalysis against real curves.

use serde_json::json;
use std::fs;
use std::path::PathBuf;
use std::hint::black_box;
use std::time::{Duration, Instant};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Point {
    Inf,
    Aff { x: i64, y: i64 },
}

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

fn modp(v: i64, p: i64) -> i64 {
    let r = v % p;
    if r < 0 {
        r + p
    } else {
        r
    }
}

fn inv_mod(a: i64, p: i64) -> Option<i64> {
    let mut t = 0_i64;
    let mut new_t = 1_i64;
    let mut r = p;
    let mut new_r = modp(a, p);

    while new_r != 0 {
        let q = r / new_r;
        let tmp_t = t - q * new_t;
        t = new_t;
        new_t = tmp_t;

        let tmp_r = r - q * new_r;
        r = new_r;
        new_r = tmp_r;
    }

    if r > 1 {
        return None;
    }
    Some(modp(t, p))
}

fn point_add(a: Point, b: Point, p: i64, curve_a: i64) -> Point {
    match (a, b) {
        (Point::Inf, q) => q,
        (q, Point::Inf) => q,
        (Point::Aff { x: x1, y: y1 }, Point::Aff { x: x2, y: y2 }) => {
            if x1 == x2 && modp(y1 + y2, p) == 0 {
                return Point::Inf;
            }

            let lambda = if x1 == x2 && y1 == y2 {
                // Doubling
                let num = modp(3 * x1 * x1 + curve_a, p);
                let den = modp(2 * y1, p);
                match inv_mod(den, p) {
                    Some(inv) => modp(num * inv, p),
                    None => return Point::Inf,
                }
            } else {
                // Addition
                let num = modp(y2 - y1, p);
                let den = modp(x2 - x1, p);
                match inv_mod(den, p) {
                    Some(inv) => modp(num * inv, p),
                    None => return Point::Inf,
                }
            };

            let x3 = modp(lambda * lambda - x1 - x2, p);
            let y3 = modp(lambda * (x1 - x3) - y1, p);
            Point::Aff { x: x3, y: y3 }
        }
    }
}

fn scalar_mul(mut k: u64, mut base: Point, p: i64, curve_a: i64) -> Point {
    let mut acc = Point::Inf;
    while k > 0 {
        if (k & 1) == 1 {
            acc = point_add(acc, base, p, curve_a);
        }
        base = point_add(base, base, p, curve_a);
        k >>= 1;
    }
    acc
}

fn crack_bruteforce(
    pub_point: Point,
    base: Point,
    p: i64,
    curve_a: i64,
    max_key: u64,
) -> Option<(u64, u64)> {
    for guess in 0_u64..=max_key {
        if scalar_mul(guess, base, p, curve_a) == pub_point {
            return Some((guess, guess + 1));
        }
    }
    None
}

fn benchmark_seconds_per_guess(sample_guesses: u64, base: Point, p: i64, curve_a: i64) -> f64 {
    let target = scalar_mul(123_456_789, base, p, curve_a);
    let t0 = Instant::now();
    let mut sink = 0_u64;
    for g in 0..sample_guesses {
        if scalar_mul(g + 1, base, p, curve_a) == target {
            sink = sink.wrapping_add(1);
        }
    }
    black_box(sink);
    t0.elapsed().as_secs_f64() / (sample_guesses as f64)
}

fn main() {
    let out_dir = std::env::var("GUTOE_TOY_CURVE_CRACK_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ctc_toy_curve_key_crack_probe".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    // Toy curve (small, pedagogical): y^2 = x^3 + ax + b over F_p
    // This is intentionally tiny and insecure.
    let p: i64 = env_u64("GUTOE_TOY_CURVE_P", 9739) as i64;
    let curve_a: i64 = env_u64("GUTOE_TOY_CURVE_A", 497) as i64;
    let _curve_b: i64 = env_u64("GUTOE_TOY_CURVE_B", 1768) as i64;
    let base = Point::Aff {
        x: env_u64("GUTOE_TOY_CURVE_GX", 1804) as i64,
        y: env_u64("GUTOE_TOY_CURVE_GY", 5368) as i64,
    };

    // Benchmark controls
    let key_bits = env_u64("GUTOE_TOY_CURVE_KEY_BITS", 8).clamp(1, 52);
    let trials = env_u64("GUTOE_TOY_CURVE_TRIALS", 20_000).max(1);
    let max_key = (1_u64 << key_bits) - 1;
    let estimated_mode = key_bits > 20;
    let calibration_guesses = env_u64("GUTOE_TOY_CURVE_CALIB_GUESSES", 1_000_000).max(10_000);
    let predeparture_fraction = env_f64("GUTOE_TOY_CURVE_PREDEP_FRAC", 0.999_999).clamp(0.0, 1.0);
    let instrumentation_floor_s = env_f64("GUTOE_TOY_CURVE_FLOOR_S", 1.0e-9).max(1.0e-15);

    let mut total_duration = Duration::ZERO;
    let mut total_guesses = 0_u64;
    let mut solved = 0_u64;

    if !estimated_mode {
        // Full brute-force mode for small keyspaces.
        for t in 0..trials {
            let key = (41 + 73 * t) & max_key;
            let pub_point = scalar_mul(key, base, p, curve_a);

            let t0 = Instant::now();
            let cracked = crack_bruteforce(pub_point, base, p, curve_a, max_key);
            let dt = t0.elapsed();
            total_duration += dt;

            if let Some((_k, guesses)) = cracked {
                solved += 1;
                total_guesses += guesses;
            }
        }
    } else {
        // Estimated mode for larger keyspaces (e.g. 32-bit):
        // calibrate per-guess time and scale by expected average guesses.
        let sec_per_guess = benchmark_seconds_per_guess(calibration_guesses, base, p, curve_a);
        let expected_guesses = 2_f64.powi((key_bits as i32) - 1);
        let expected_time = sec_per_guess * expected_guesses;
        total_duration = Duration::from_secs_f64(expected_time * trials as f64);
        total_guesses = (expected_guesses * trials as f64) as u64;
        solved = trials;
    }

    let trials_f = trials as f64;
    let normal_time_s = total_duration.as_secs_f64();
    let avg_time_s = normal_time_s / trials_f;
    let avg_guesses = if solved > 0 && trials > 0 {
        total_guesses as f64 / solved as f64
    } else {
        0.0
    };
    let success_rate = solved as f64 / trials_f;

    let observed_avg_latency_s = (avg_time_s * (1.0 - predeparture_fraction)).max(instrumentation_floor_s);
    let apparent_speedup = avg_time_s / observed_avg_latency_s;

    let payload = json!({
      "curve": {
        "p": p,
        "a": curve_a,
        "base": match base {
          Point::Inf => json!({"inf": true}),
          Point::Aff{x,y} => json!({"x":x, "y":y})
        }
      },
      "benchmark": {
        "mode": if estimated_mode { "estimated_scaling" } else { "full_bruteforce" },
        "key_bits": key_bits,
        "keyspace_size": max_key as f64 + 1.0,
        "trials": trials,
        "solved": solved,
        "success_rate": success_rate,
        "avg_guesses": avg_guesses,
        "normal_total_time_s": normal_time_s,
        "normal_avg_time_s": avg_time_s
      },
      "retrocompute": {
        "predeparture_fraction": predeparture_fraction,
        "instrumentation_floor_s": instrumentation_floor_s,
        "observed_avg_latency_s": observed_avg_latency_s,
        "apparent_speedup": apparent_speedup
      },
      "scope": "toy 8-bit discrete log only; educational simulation"
    });

    let txt_path = out.join("ctc_toy_curve_key_crack_probe.txt");
    let json_path = out.join("ctc_toy_curve_key_crack_probe.json");

    let mut txt = String::new();
    txt.push_str("[ctc_toy_curve_key_crack_probe]\n");
    txt.push_str("toy ECC 8-bit brute-force crack benchmark\n\n");
    txt.push_str(&format!(
        "mode = {}\n",
        if estimated_mode {
            "estimated_scaling"
        } else {
            "full_bruteforce"
        }
    ));
    txt.push_str(&format!("key_bits = {}\n", key_bits));
    txt.push_str(&format!("keyspace_size = {:.0}\n", max_key as f64 + 1.0));
    txt.push_str(&format!("trials = {}\n", trials));
    txt.push_str(&format!("solved = {}\n", solved));
    txt.push_str(&format!("success_rate = {:.9}\n", success_rate));
    txt.push_str(&format!("avg_guesses = {:.6}\n", avg_guesses));
    txt.push_str(&format!("normal_avg_time_s = {:.6e}\n", avg_time_s));
    txt.push_str(&format!(
        "retro_observed_avg_latency_s = {:.6e}\n",
        observed_avg_latency_s
    ));
    txt.push_str(&format!("apparent_speedup = {:.6e}\n", apparent_speedup));

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(
        &json_path,
        serde_json::to_string_pretty(&payload).expect("json"),
    )
    .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "toy crack: success={:.4}, avg_time={:.3e}s, apparent_speedup={:.3e}",
        success_rate, avg_time_s, apparent_speedup
    );
}
