//! Fiber channel capacity / compression probe.
//!
//! Structural lane:
//! - visible dims = 4
//! - hidden fiber dims = 12
//! - ratio hidden:visible = 3:1
//!
//! Reports:
//! - quantized bits/event for chosen bits-per-dimension
//! - simple AWGN MIMO capacity estimate for hidden 12 real channels
//! - passive decode limitation flag (projection-only)

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

fn log2(x: f64) -> f64 {
    x.ln() / 2.0_f64.ln()
}

fn main() {
    let out_dir = std::env::var("GUTOE_FIBER_CAPACITY_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/fiber_channel_capacity_probe".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let visible_dims = 4.0_f64;
    let hidden_dims = 12.0_f64;
    let total_dims = visible_dims + hidden_dims;

    let bits_per_dim = env_u64("GUTOE_FIBER_BITS_PER_DIM", 16);
    let snr_hidden = env_f64("GUTOE_FIBER_HIDDEN_SNR", 100.0).max(0.0);
    let bandwidth_hz = env_f64("GUTOE_FIBER_BW_HZ", 1.0).max(0.0);

    // Quantized structural payloads
    let visible_bits = visible_dims * bits_per_dim as f64;
    let hidden_bits = hidden_dims * bits_per_dim as f64;
    let total_bits = total_dims * bits_per_dim as f64;

    // Real-valued AWGN MIMO estimate for hidden 12 channels, identity coupling:
    // C = (n/2) * log2(1 + SNR) bits/use, and bits/s = BW * C.
    let c_hidden_bits_per_use = 0.5 * hidden_dims * log2(1.0 + snr_hidden);
    let c_hidden_bits_per_s = bandwidth_hz * c_hidden_bits_per_use;

    let payload = json!({
      "structure": {
        "visible_dims": visible_dims,
        "hidden_dims": hidden_dims,
        "total_dims": total_dims,
        "hidden_visible_ratio": hidden_dims / visible_dims
      },
      "quantized_payload": {
        "bits_per_dim": bits_per_dim,
        "visible_bits_per_event": visible_bits,
        "hidden_bits_per_event": hidden_bits,
        "total_bits_per_event": total_bits,
        "hidden_is_three_times_visible": (hidden_bits - 3.0*visible_bits).abs() < 1e-12
      },
      "awgn_hidden_capacity_estimate": {
        "snr_hidden_linear": snr_hidden,
        "bandwidth_hz": bandwidth_hz,
        "hidden_capacity_bits_per_use": c_hidden_bits_per_use,
        "hidden_capacity_bits_per_second": c_hidden_bits_per_s
      },
      "decode_boundary": {
        "passive_projection_only_decode_hidden": false,
        "keyed_decode_hidden_possible": true,
        "note": "projection-only is non-injective; hidden recovery requires fiber key"
      }
    });

    let txt_path = out.join("fiber_channel_capacity_probe.txt");
    let json_path = out.join("fiber_channel_capacity_probe.json");

    let mut txt = String::new();
    txt.push_str("[fiber_channel_capacity_probe]\n");
    txt.push_str(&format!(
        "dims: visible={:.0}, hidden={:.0}, ratio={:.3}\n",
        visible_dims,
        hidden_dims,
        hidden_dims / visible_dims
    ));
    txt.push_str(&format!(
        "quantized: bits_per_dim={}, visible={:.0} bits/event, hidden={:.0} bits/event, total={:.0} bits/event\n",
        bits_per_dim, visible_bits, hidden_bits, total_bits
    ));
    txt.push_str(&format!(
        "awgn_hidden: snr={:.6e}, bw={:.6e}Hz, cap={:.6e} bits/use, {:.6e} bits/s\n",
        snr_hidden, bandwidth_hz, c_hidden_bits_per_use, c_hidden_bits_per_s
    ));
    txt.push_str("decode: passive_hidden=false, keyed_hidden=true\n");

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(&json_path, serde_json::to_string_pretty(&payload).expect("json")).expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
}

