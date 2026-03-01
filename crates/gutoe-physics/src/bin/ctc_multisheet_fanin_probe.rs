//! Multi-sheet fan-in probe for temporal-credit compounding.
//!
//! This probe distinguishes:
//! - split routing (one lineage reaches tracked origin branch),
//! - fan-in routing (many independent sheets target same origin branch).
//!
//! Per-sender send capability recurrence:
//!   s_{k+1} = max(0, eta * infra_gain * s_k + base_inflow - loss)
//!
//! Branch population:
//!   N_k = branching^k
//!
//! Contribution to tracked origin branch:
//! - split mode: c_k = s_k
//! - fan-in mode: c_k = merge_fraction * N_k * s_k
//!
//! Conservation per generation:
//!   c_k <= total_sent_k = N_k * s_k   (for merge_fraction <= 1)

use serde_json::json;
use std::fs;
use std::path::PathBuf;

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

fn simulate(
    generations: usize,
    branching: f64,
    merge_fraction: f64,
    eta: f64,
    infra_gain: f64,
    base_inflow: f64,
    loss: f64,
    s0: f64,
    mode_fanin: bool,
) -> serde_json::Value {
    let mut s = s0.max(0.0);
    let mut cumulative = 0.0_f64;
    let mut rows = Vec::with_capacity(generations + 1);
    let mut max_violation = 0.0_f64;
    let mut prev_c = None::<f64>;
    let mut late_ratios = Vec::new();

    for k in 0..=generations {
        let n_k = branching.powi(k as i32).max(0.0);
        let total_sent = n_k * s;
        let c_k = if mode_fanin {
            merge_fraction * total_sent
        } else {
            s
        };
        cumulative += c_k;

        let violation = (c_k - total_sent).max(0.0);
        max_violation = max_violation.max(violation);

        if let Some(pc) = prev_c {
            if pc > 0.0 {
                let r = c_k / pc;
                if k > generations.saturating_sub(8) {
                    late_ratios.push(r);
                }
            }
        }
        prev_c = Some(c_k);

        rows.push(json!({
            "k": k,
            "send_per_sender_j": s,
            "sheet_count": n_k,
            "total_sent_j": total_sent,
            "contrib_to_origin_j": c_k,
            "cumulative_origin_j": cumulative
        }));

        s = (eta * infra_gain * s + base_inflow - loss).max(0.0);
    }

    let empirical_late_ratio = if late_ratios.is_empty() {
        None
    } else {
        Some(late_ratios.iter().sum::<f64>() / late_ratios.len() as f64)
    };

    json!({
        "mode": if mode_fanin { "fanin" } else { "split" },
        "max_generation_conservation_violation_j": max_violation,
        "empirical_late_ratio": empirical_late_ratio,
        "rows": rows
    })
}

fn main() {
    let out_dir = std::env::var("GUTOE_CTC_MULTISHEET_FANIN_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ctc_multisheet_fanin_probe".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let generations = env_usize("GUTOE_CTC_FANIN_GENERATIONS", 60).max(2);
    let branching = env_f64("GUTOE_CTC_FANIN_BRANCHING", 2.0).max(1.0);
    let merge_fraction = env_f64("GUTOE_CTC_FANIN_MERGE_FRACTION", 1.0).clamp(0.0, 1.0);
    let eta = env_f64("GUTOE_CTC_FANIN_ETA", 0.98).max(0.0);
    let infra_gain = env_f64("GUTOE_CTC_FANIN_INFRA_GAIN", 1.02).max(0.0);
    let base_inflow = env_f64("GUTOE_CTC_FANIN_BASE_INFLOW_J", 1e-9).max(0.0);
    let loss = env_f64("GUTOE_CTC_FANIN_LOSS_J", 1e-10).max(0.0);
    let s0 = env_f64("GUTOE_CTC_FANIN_S0_J", 1e-9).max(0.0);

    let per_sender_multiplier = eta * infra_gain;
    let fanin_seed_multiplier = branching * merge_fraction * per_sender_multiplier;

    let split = simulate(
        generations,
        branching,
        merge_fraction,
        eta,
        infra_gain,
        base_inflow,
        loss,
        s0,
        false,
    );
    let fanin = simulate(
        generations,
        branching,
        merge_fraction,
        eta,
        infra_gain,
        base_inflow,
        loss,
        s0,
        true,
    );

    let split_final = split["rows"]
        .as_array()
        .and_then(|v| v.last())
        .and_then(|r| r["cumulative_origin_j"].as_f64())
        .unwrap_or(f64::NAN);
    let fanin_final = fanin["rows"]
        .as_array()
        .and_then(|v| v.last())
        .and_then(|r| r["cumulative_origin_j"].as_f64())
        .unwrap_or(f64::NAN);

    let payload = json!({
        "inputs": {
            "generations": generations,
            "branching": branching,
            "merge_fraction": merge_fraction,
            "eta": eta,
            "infra_gain": infra_gain,
            "base_inflow_j": base_inflow,
            "loss_j": loss,
            "s0_j": s0
        },
        "derived": {
            "per_sender_multiplier": per_sender_multiplier,
            "fanin_seed_multiplier": fanin_seed_multiplier,
            "criterion_hint": "fanin growth strengthens when branching*merge_fraction*eta*infra_gain > 1"
        },
        "split": split,
        "fanin": fanin,
        "summary": {
            "cumulative_origin_split_j": split_final,
            "cumulative_origin_fanin_j": fanin_final,
            "fanin_over_split_ratio": if split_final > 0.0 { fanin_final / split_final } else { f64::NAN }
        }
    });

    let txt_path = out.join("ctc_multisheet_fanin_probe.txt");
    let json_path = out.join("ctc_multisheet_fanin_probe.json");

    let mut txt = String::new();
    txt.push_str("[ctc_multisheet_fanin_probe]\n");
    txt.push_str(&format!(
        "generations={}, branching={:.6e}, merge_fraction={:.6e}, eta={:.6e}, infra_gain={:.6e}, base_inflow={:.6e}J, loss={:.6e}J\n",
        generations, branching, merge_fraction, eta, infra_gain, base_inflow, loss
    ));
    txt.push_str(&format!(
        "per_sender_multiplier={:.12e}, fanin_seed_multiplier={:.12e}\n",
        per_sender_multiplier, fanin_seed_multiplier
    ));
    txt.push_str(&format!(
        "cumulative_origin_split_j={:.12e}, cumulative_origin_fanin_j={:.12e}, fanin_over_split_ratio={:.12e}\n",
        split_final,
        fanin_final,
        if split_final > 0.0 { fanin_final / split_final } else { f64::NAN }
    ));
    txt.push_str("\n[notes]\n");
    txt.push_str("split mode = one-lineage receiver (no branch merge)\n");
    txt.push_str("fanin mode = many independent sheets contribute to one receiver branch\n");
    txt.push_str("conservation checked per generation: contrib_to_origin <= total_sent\n");

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(&json_path, serde_json::to_string_pretty(&payload).expect("json")).expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
}
