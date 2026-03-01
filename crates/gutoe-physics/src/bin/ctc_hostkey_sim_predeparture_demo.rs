//! Host-key -> simulated-reality predeparture crack demo.
//!
//! This executable demonstrates the exact flow:
//! 1) generate a real host key from OS entropy,
//! 2) commit to it on host,
//! 3) inject public challenge material into simulation,
//! 4) crack in forward compute,
//! 5) make result available at predeparture simulated time,
//! 6) verify commitment and timeline ordering.
//!
//! The challenge is hash-preimage based to guarantee exact key recovery
//! (avoids scalar-equivalence ambiguity of tiny toy EC groups).

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

fn mask_for_bits(bits: u64) -> u64 {
    if bits >= 64 {
        u64::MAX
    } else {
        (1_u64 << bits) - 1
    }
}

fn host_commit_hex(key_bits: u64, key: u64, nonce: u64) -> String {
    let msg = format!("host_key_bits={key_bits};host_key={key};nonce={nonce}");
    blake3::hash(msg.as_bytes()).to_hex().to_string()
}

fn challenge_digest_hex(key_bits: u64, key: u64, nonce: u64) -> String {
    let msg = format!("sim_challenge:key_bits={key_bits};key={key};nonce={nonce}");
    blake3::hash(msg.as_bytes()).to_hex().to_string()
}

fn crack_hash_preimage(
    key_bits: u64,
    nonce: u64,
    target_digest: &str,
    max_key: u64,
) -> Option<(u64, u64)> {
    for guess in 0_u64..=max_key {
        let d = challenge_digest_hex(key_bits, guess, nonce);
        if d == target_digest {
            return Some((guess, guess + 1));
        }
    }
    None
}

fn main() {
    let out_dir = std::env::var("GUTOE_HOSTKEY_SIM_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ctc_hostkey_sim_predeparture_demo".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    // Demo controls
    let key_bits = env_u64("GUTOE_HOSTKEY_BITS", 20).clamp(4, 32);
    let max_key = mask_for_bits(key_bits);
    let sim_t_inject = env_f64("GUTOE_HOSTKEY_SIM_T_INJECT", 100.0);
    let sim_t_crack_start = env_f64("GUTOE_HOSTKEY_SIM_T_CRACK_START", 1_000.0);
    let sim_seconds_per_host_second = env_f64("GUTOE_HOSTKEY_SIM_SCALE", 1.0).max(1e-12);
    let retro_shift_factor = env_f64("GUTOE_HOSTKEY_RETRO_SHIFT_FACTOR", 1.20).max(0.0);
    let observer_floor_s = env_f64("GUTOE_HOSTKEY_OBSERVER_FLOOR_S", 1.0e-9).max(1e-15);

    // 1) Host generates a real key from OS entropy and commits.
    let mut rng = OsRng;
    let mut raw = rng.next_u64() & max_key;
    if raw == 0 {
        raw = 1;
    }
    let host_key = raw;
    let host_nonce = rng.next_u64();
    let host_commit = host_commit_hex(key_bits, host_key, host_nonce);

    // 2) Inject host key challenge material into simulated world.
    let public_challenge_digest = challenge_digest_hex(key_bits, host_key, host_nonce);

    // 3) Crack in forward host compute (represents forward cover time in sim).
    let t0 = Instant::now();
    let cracked = crack_hash_preimage(key_bits, host_nonce, &public_challenge_digest, max_key);
    let host_crack_time_s = t0.elapsed().as_secs_f64();
    let (recovered_key, guesses) = cracked.expect("hash-preimage crack must recover key");

    // 4) Reveal/verify commitment.
    let recovered_commit = host_commit_hex(key_bits, recovered_key, host_nonce);
    let recovered_digest = challenge_digest_hex(key_bits, recovered_key, host_nonce);
    let challenge_verified = recovered_digest == public_challenge_digest;
    let commitment_verified = recovered_commit == host_commit;
    let key_match = recovered_key == host_key;

    // 5) Simulated timeline with retro return.
    let sim_cover_duration_s = host_crack_time_s * sim_seconds_per_host_second;
    let sim_t_cover_end = sim_t_crack_start + sim_cover_duration_s;
    let sim_t_result_available = sim_t_cover_end - retro_shift_factor * sim_cover_duration_s;
    let apparent_latency_from_start = sim_t_result_available - sim_t_crack_start;
    let predeparture = sim_t_result_available < sim_t_crack_start;
    let predeparture_margin_s = (sim_t_crack_start - sim_t_result_available).max(0.0);
    let observer_response_latency_s = if predeparture {
        observer_floor_s
    } else {
        (sim_t_result_available - sim_t_crack_start).max(observer_floor_s)
    };
    let normal_latency_s = sim_cover_duration_s.max(observer_floor_s);
    let apparent_speedup = normal_latency_s / observer_response_latency_s;

    let timeline = vec![
        json!({
          "t_sim": sim_t_inject,
          "event": "inject_host_public_challenge",
          "details": {
            "key_bits": key_bits,
            "nonce": host_nonce,
            "challenge_digest_blake3": public_challenge_digest,
            "host_commit_blake3": host_commit
          }
        }),
        json!({"t_sim": sim_t_crack_start, "event": "observer_query_key"}),
        json!({"t_sim": sim_t_crack_start, "event": "crack_start_cover"}),
        json!({"t_sim": sim_t_cover_end, "event": "crack_end_cover"}),
        json!({"t_sim": sim_t_result_available, "event": "recovered_key_available"}),
    ];

    let payload = json!({
      "host_reality": {
        "challenge_type": "blake3_preimage",
        "key_bits": key_bits,
        "host_key": host_key,
        "host_nonce": host_nonce,
        "host_commit_blake3": host_commit,
        "public_challenge_digest_blake3": public_challenge_digest,
        "host_crack_time_s": host_crack_time_s,
        "key_recovered": recovered_key,
        "key_match": key_match,
        "challenge_verified": challenge_verified,
        "commitment_verified": commitment_verified,
        "guesses_used": guesses
      },
      "simulated_reality": {
        "sim_t_inject": sim_t_inject,
        "sim_t_crack_start": sim_t_crack_start,
        "sim_cover_duration_s": sim_cover_duration_s,
        "sim_t_cover_end": sim_t_cover_end,
        "retro_shift_factor": retro_shift_factor,
        "sim_t_result_available": sim_t_result_available,
        "apparent_latency_from_start_s": apparent_latency_from_start,
        "predeparture": predeparture,
        "predeparture_margin_s": predeparture_margin_s,
        "observer_response_latency_s": observer_response_latency_s,
        "normal_latency_s": normal_latency_s,
        "apparent_speedup": apparent_speedup
      },
      "timeline": timeline,
      "scope": "host-key provenance + in-sim crack + predeparture observer semantics"
    });

    let txt_path = out.join("ctc_hostkey_sim_predeparture_demo.txt");
    let json_path = out.join("ctc_hostkey_sim_predeparture_demo.json");

    let mut txt = String::new();
    txt.push_str("[ctc_hostkey_sim_predeparture_demo]\n");
    txt.push_str(
        "host-generated key injected into sim as public challenge, cracked in sim cover-time, observed predeparture\n\n",
    );
    txt.push_str("challenge_type = blake3_preimage\n");
    txt.push_str(&format!("key_bits = {}\n", key_bits));
    txt.push_str(&format!("host_commit_blake3 = {}\n", host_commit));
    txt.push_str(&format!(
        "public_challenge_digest_blake3 = {}\n",
        public_challenge_digest
    ));
    txt.push_str(&format!("key_match = {}\n", key_match));
    txt.push_str(&format!("challenge_verified = {}\n", challenge_verified));
    txt.push_str(&format!("commitment_verified = {}\n", commitment_verified));
    txt.push_str(&format!("guesses_used = {}\n", guesses));
    txt.push_str(&format!("host_crack_time_s = {:.6e}\n", host_crack_time_s));
    txt.push_str(&format!("sim_cover_duration_s = {:.6e}\n", sim_cover_duration_s));
    txt.push_str(&format!(
        "sim_result_available_t = {:.6e} (start {:.6e})\n",
        sim_t_result_available, sim_t_crack_start
    ));
    txt.push_str(&format!("predeparture = {}\n", predeparture));
    txt.push_str(&format!("predeparture_margin_s = {:.6e}\n", predeparture_margin_s));
    txt.push_str(&format!("normal_latency_s = {:.6e}\n", normal_latency_s));
    txt.push_str(&format!(
        "observer_response_latency_s = {:.6e}\n",
        observer_response_latency_s
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
        "demo: key_bits={}, key_match={}, predeparture={}, speedup={:.3e}",
        key_bits, key_match, predeparture, apparent_speedup
    );
}
