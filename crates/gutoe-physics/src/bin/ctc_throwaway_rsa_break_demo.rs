//! Throwaway RSA break demo in controlled simulation lane.
//!
//! Safety scope:
//! - Generates a fresh intentionally weak RSA key on host (toy strength).
//! - Cracks only that generated key.
//! - Reports normal host compute time vs retro-observed simulated latency.
//! - No targeting of external/real-world keys.

use rand::rngs::OsRng;
use rand::RngCore;
use serde_json::json;
use std::fs;
use std::path::PathBuf;
use std::time::Instant;

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

fn gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

fn mod_pow(mut base: u64, mut exp: u64, modu: u64) -> u64 {
    if modu == 1 {
        return 0;
    }
    let mut out: u128 = 1;
    let mut b: u128 = (base % modu) as u128;
    let m = modu as u128;
    while exp > 0 {
        if (exp & 1) == 1 {
            out = (out * b) % m;
        }
        b = (b * b) % m;
        exp >>= 1;
    }
    out as u64
}

fn mod_inverse_u64(a: u64, m: u64) -> Option<u64> {
    let mut t: i128 = 0;
    let mut new_t: i128 = 1;
    let mut r: i128 = m as i128;
    let mut new_r: i128 = a as i128;

    while new_r != 0 {
        let q = r / new_r;
        let tmp_t = t - q * new_t;
        t = new_t;
        new_t = tmp_t;
        let tmp_r = r - q * new_r;
        r = new_r;
        new_r = tmp_r;
    }

    if r != 1 {
        return None;
    }
    if t < 0 {
        t += m as i128;
    }
    Some(t as u64)
}

fn is_prime_trial(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    if n % 2 == 0 {
        return n == 2;
    }
    let mut d = 3_u64;
    while (d as u128) * (d as u128) <= n as u128 {
        if n % d == 0 {
            return false;
        }
        d += 2;
    }
    true
}

fn random_odd_in_bits(rng: &mut OsRng, bits: u64) -> u64 {
    let b = bits.clamp(4, 31);
    let mask = (1_u64 << b) - 1;
    let top = 1_u64 << (b - 1);
    let mut x = rng.next_u64() & mask;
    x |= 1;
    x |= top;
    x
}

fn random_prime_bits(rng: &mut OsRng, bits: u64) -> u64 {
    loop {
        let cand = random_odd_in_bits(rng, bits);
        if is_prime_trial(cand) {
            return cand;
        }
    }
}

fn factor_by_trial(n: u64) -> Option<(u64, u64, u64)> {
    if n % 2 == 0 {
        return Some((2, n / 2, 1));
    }
    let mut f = 3_u64;
    let mut tested = 0_u64;
    while (f as u128) * (f as u128) <= n as u128 {
        tested += 1;
        if n % f == 0 {
            return Some((f, n / f, tested));
        }
        f += 2;
    }
    None
}

fn main() {
    let out_dir = std::env::var("GUTOE_RSA_THROWAWAY_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ctc_throwaway_rsa_break_demo".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    // Intentionally weak toy key size for controlled demonstration.
    let n_bits = env_u64("GUTOE_RSA_THROWAWAY_BITS", 48).clamp(16, 52);
    let prime_bits = (n_bits / 2).clamp(8, 31);

    let sim_t_inject = env_f64("GUTOE_RSA_SIM_T_INJECT", 100.0);
    let sim_t_crack_start = env_f64("GUTOE_RSA_SIM_T_CRACK_START", 1_000.0);
    let sim_seconds_per_host_second = env_f64("GUTOE_RSA_SIM_SCALE", 1.0).max(1e-12);
    let retro_shift_factor = env_f64("GUTOE_RSA_RETRO_SHIFT_FACTOR", 1.20).max(0.0);
    let observer_floor_s = env_f64("GUTOE_RSA_OBSERVER_FLOOR_S", 1.0e-9).max(1e-15);

    // 1) Fresh throwaway host key generation.
    let mut rng = OsRng;
    let p = random_prime_bits(&mut rng, prime_bits);
    let mut q = random_prime_bits(&mut rng, prime_bits);
    while q == p {
        q = random_prime_bits(&mut rng, prime_bits);
    }
    let n = p.saturating_mul(q);
    let phi = (p - 1).saturating_mul(q - 1);

    let mut e = 65_537_u64;
    if gcd(e, phi) != 1 {
        e = 3;
        while e < phi && gcd(e, phi) != 1 {
            e += 2;
        }
    }
    let d = mod_inverse_u64(e, phi).expect("e and phi are coprime");

    let mut msg = rng.next_u64() % n.max(3);
    if msg < 2 {
        msg = 2;
    }
    let ciphertext = mod_pow(msg, e, n);
    let public_commit = blake3::hash(format!("n={n};e={e};c={ciphertext}").as_bytes())
        .to_hex()
        .to_string();

    // 2) Simulated forward-time crack (factor n, recover d, decrypt).
    let t0 = Instant::now();
    let (pf, qf, divisors_tested) = factor_by_trial(n).expect("composite n should factor");
    let crack_time_s = t0.elapsed().as_secs_f64();
    let phi_f = (pf - 1).saturating_mul(qf - 1);
    let d_f = mod_inverse_u64(e, phi_f).expect("recovered factors should permit inverse");
    let msg_f = mod_pow(ciphertext, d_f, n);

    let factor_ok = (pf * qf == n) || (qf * pf == n);
    let decrypt_ok = msg_f == msg;
    let d_ok = d_f == d;

    // 3) Simulated retro-observed timing.
    let sim_cover_duration_s = crack_time_s * sim_seconds_per_host_second;
    let sim_t_cover_end = sim_t_crack_start + sim_cover_duration_s;
    let sim_t_result_available = sim_t_cover_end - retro_shift_factor * sim_cover_duration_s;
    let predeparture = sim_t_result_available < sim_t_crack_start;
    let predeparture_margin_s = (sim_t_crack_start - sim_t_result_available).max(0.0);
    let observer_response_latency_s = if predeparture {
        observer_floor_s
    } else {
        (sim_t_result_available - sim_t_crack_start).max(observer_floor_s)
    };
    let normal_latency_s = sim_cover_duration_s.max(observer_floor_s);
    let apparent_speedup = normal_latency_s / observer_response_latency_s;

    let payload = json!({
      "host_generation": {
        "throwaway": true,
        "n_bits_target": n_bits,
        "prime_bits": prime_bits,
        "p": p,
        "q": q,
        "n": n,
        "phi": phi,
        "e": e,
        "d": d,
        "message": msg,
        "ciphertext": ciphertext,
        "public_commit_blake3": public_commit
      },
      "crack_result": {
        "factored_p": pf,
        "factored_q": qf,
        "divisors_tested": divisors_tested,
        "recovered_d": d_f,
        "recovered_message": msg_f,
        "factor_ok": factor_ok,
        "decrypt_ok": decrypt_ok,
        "d_match": d_ok,
        "host_crack_time_s": crack_time_s
      },
      "simulated_reality": {
        "sim_t_inject": sim_t_inject,
        "sim_t_crack_start": sim_t_crack_start,
        "sim_cover_duration_s": sim_cover_duration_s,
        "sim_t_cover_end": sim_t_cover_end,
        "retro_shift_factor": retro_shift_factor,
        "sim_t_result_available": sim_t_result_available,
        "predeparture": predeparture,
        "predeparture_margin_s": predeparture_margin_s,
        "observer_response_latency_s": observer_response_latency_s,
        "normal_latency_s": normal_latency_s,
        "apparent_speedup": apparent_speedup
      },
      "scope": "controlled throwaway RSA only; no external key targeting"
    });

    let txt_path = out.join("ctc_throwaway_rsa_break_demo.txt");
    let json_path = out.join("ctc_throwaway_rsa_break_demo.json");

    let mut txt = String::new();
    txt.push_str("[ctc_throwaway_rsa_break_demo]\n");
    txt.push_str("fresh host-generated weak throwaway RSA, cracked in controlled simulation lane\n\n");
    txt.push_str(&format!("n_bits_target = {}\n", n_bits));
    txt.push_str(&format!("n = {}\n", n));
    txt.push_str(&format!("e = {}\n", e));
    txt.push_str(&format!("ciphertext = {}\n", ciphertext));
    txt.push_str(&format!("public_commit_blake3 = {}\n", public_commit));
    txt.push_str(&format!("factor_ok = {}\n", factor_ok));
    txt.push_str(&format!("d_match = {}\n", d_ok));
    txt.push_str(&format!("decrypt_ok = {}\n", decrypt_ok));
    txt.push_str(&format!("host_crack_time_s = {:.6e}\n", crack_time_s));
    txt.push_str(&format!("predeparture = {}\n", predeparture));
    txt.push_str(&format!("predeparture_margin_s = {:.6e}\n", predeparture_margin_s));
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
        "throwaway_rsa_demo: bits={}, factor_ok={}, decrypt_ok={}, predeparture={}, speedup={:.3e}",
        n_bits, factor_ok, decrypt_ok, predeparture, apparent_speedup
    );
}

