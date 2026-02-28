//! Real-matter CTC paradox simulator.
//!
//! This lane treats paradox logic with explicit matter packets:
//! - ancestor packet mass m_A (always present as matter, alive/dead state)
//! - traveler packet mass m_T (either present on local slice or in loop channel)
//!
//! We evaluate:
//! 1) strict single-history contradiction (forced paradox unsat),
//! 2) branch-split consistency,
//! 3) Deutsch NOT fixed point p*=0.5,
//! 4) Monte Carlo branch frequencies + mass bookkeeping invariants.

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::json;
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

#[derive(Debug, Clone, Copy)]
struct Single {
    a0: u8,
    t: u8,
    k: u8,
    a1: u8,
}

#[derive(Debug, Clone, Copy)]
struct Branch {
    o_a0: u8,
    o_t: u8,
    o_a1: u8,
    t_a0: u8,
    t_k: u8,
    t_a1: u8,
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

fn env_usize(name: &str, default: usize) -> usize {
    std::env::var(name)
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(default)
}

fn main() {
    let out_dir = std::env::var("GUTOE_CTC_REAL_MATTER_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ctc_real_matter_paradox".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let m_ancestor_kg = env_f64("GUTOE_CTC_M_ANCESTOR_KG", 80.0);
    let m_traveler_kg = env_f64("GUTOE_CTC_M_TRAVELER_KG", 80.0);
    let p_input = env_f64("GUTOE_CTC_P_INPUT", 0.37);
    let samples = env_usize("GUTOE_CTC_SAMPLES", 200_000);
    let seed = env_u64("GUTOE_CTC_SEED", 137);

    // ── Test 1: strict single-history unsat for forced paradox branch ───────
    let mut single_valid = Vec::new();
    for a0 in [0u8, 1u8] {
        for t in [0u8, 1u8] {
            for k in [0u8, 1u8] {
                for a1 in [0u8, 1u8] {
                    let c1 = t == a0;
                    let c2 = k == t;
                    let c3 = a1 == (a0 * (1 - k));
                    let c4 = a1 == a0;
                    if c1 && c2 && c3 && c4 {
                        single_valid.push(Single { a0, t, k, a1 });
                    }
                }
            }
        }
    }
    let single_forced_paradox: Vec<_> = single_valid.iter().copied().filter(|x| x.t == 1).collect();
    let single_forced_consistent = !single_forced_paradox.is_empty();

    // ── Test 2: branch-split consistency lane ────────────────────────────────
    let mut branch_valid = Vec::new();
    for o_a0 in [0u8, 1u8] {
        for o_t in [0u8, 1u8] {
            for o_a1 in [0u8, 1u8] {
                for t_a0 in [0u8, 1u8] {
                    for t_k in [0u8, 1u8] {
                        for t_a1 in [0u8, 1u8] {
                            let o1 = o_t == o_a0;
                            let o2 = o_a1 == o_a0;
                            let t1 = t_k == o_t;
                            let t2 = t_a1 == (t_a0 * (1 - t_k));
                            if o1 && o2 && t1 && t2 {
                                branch_valid.push(Branch {
                                    o_a0,
                                    o_t,
                                    o_a1,
                                    t_a0,
                                    t_k,
                                    t_a1,
                                });
                            }
                        }
                    }
                }
            }
        }
    }
    let branch_paradox_style: Vec<_> = branch_valid
        .iter()
        .copied()
        .filter(|a| a.o_a0 == 1 && a.o_t == 1 && a.t_a0 == 1 && a.t_k == 1 && a.t_a1 == 0)
        .collect();
    let branch_paradox_consistent = !branch_paradox_style.is_empty();

    // ── Test 3: Deutsch fixed point for NOT map ──────────────────────────────
    let p_next_input = 1.0 - p_input;
    let input_residual = (p_input - p_next_input).abs();
    let p_star = 0.5_f64;
    let p_star_residual = (p_star - (1.0 - p_star)).abs();

    // ── Test 4: Monte Carlo at p*=0.5 with real-matter bookkeeping ─────────
    let mut rng = StdRng::seed_from_u64(seed);
    let mut count_traveler = 0usize;
    let mut count_target_killed = 0usize;
    let mut max_mass_err = 0.0_f64;
    let total_mass_expected = m_ancestor_kg + m_traveler_kg;
    let mut total_mass_sum = 0.0_f64;

    for _ in 0..samples {
        let traveler_present = rng.gen_bool(p_star);
        if traveler_present {
            count_traveler += 1;
            count_target_killed += 1; // in this lane, traveler-present branch performs kill
        }

        // Real-matter bookkeeping:
        // - ancestor packet mass always present (alive or dead arrangement)
        // - traveler packet is either local (present branch) or in loop channel (complement branch)
        let m_local = m_ancestor_kg + if traveler_present { m_traveler_kg } else { 0.0 };
        let m_channel = if traveler_present { 0.0 } else { m_traveler_kg };
        let m_total = m_local + m_channel;
        let err = (m_total - total_mass_expected).abs();
        if err > max_mass_err {
            max_mass_err = err;
        }
        total_mass_sum += m_total;
    }

    let freq_traveler = count_traveler as f64 / samples as f64;
    let freq_killed = count_target_killed as f64 / samples as f64;
    let monte_residual = (freq_traveler - 0.5).abs();
    let avg_total_mass = total_mass_sum / samples as f64;

    // Expected-value readout at p*=0.5.
    let expected_local_mass = m_ancestor_kg + p_star * m_traveler_kg;
    let expected_channel_mass = (1.0 - p_star) * m_traveler_kg;

    let txt_path = out.join("ctc_real_matter_paradox_report.txt");
    let json_path = out.join("ctc_real_matter_paradox_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[inputs]").expect("write");
    writeln!(txt, "m_ancestor_kg = {:.6}", m_ancestor_kg).expect("write");
    writeln!(txt, "m_traveler_kg = {:.6}", m_traveler_kg).expect("write");
    writeln!(txt, "p_input = {:.12}", p_input).expect("write");
    writeln!(txt, "samples = {}", samples).expect("write");
    writeln!(txt, "seed = {}", seed).expect("write");
    writeln!(txt).expect("write");

    writeln!(txt, "[test_1_single_history]").expect("write");
    writeln!(txt, "valid_assignments = {}", single_valid.len()).expect("write");
    writeln!(
        txt,
        "forced_paradox_assignments = {}",
        single_forced_paradox.len()
    )
    .expect("write");
    writeln!(
        txt,
        "forced_paradox_consistent = {}",
        single_forced_consistent
    )
    .expect("write");
    for (i, a) in single_valid.iter().enumerate() {
        writeln!(txt, "  valid[{}]: a0={} t={} k={} a1={}", i, a.a0, a.t, a.k, a.a1)
            .expect("write");
    }
    writeln!(txt).expect("write");

    writeln!(txt, "[test_2_branch_split]").expect("write");
    writeln!(txt, "valid_assignments = {}", branch_valid.len()).expect("write");
    writeln!(
        txt,
        "paradox_style_assignments = {}",
        branch_paradox_style.len()
    )
    .expect("write");
    writeln!(
        txt,
        "paradox_style_consistent = {}",
        branch_paradox_consistent
    )
    .expect("write");
    for (i, a) in branch_paradox_style.iter().take(8).enumerate() {
        writeln!(
            txt,
            "  paradox[{}]: O(a0={},t={},a1={}) | T(a0={},k={},a1={})",
            i, a.o_a0, a.o_t, a.o_a1, a.t_a0, a.t_k, a.t_a1
        )
        .expect("write");
    }
    writeln!(txt).expect("write");

    writeln!(txt, "[test_3_deutsch_fixed_point]").expect("write");
    writeln!(txt, "p_input_next = {:.12}", p_next_input).expect("write");
    writeln!(txt, "input_not_residual = {:.12e}", input_residual).expect("write");
    writeln!(txt, "p_star = {:.12}", p_star).expect("write");
    writeln!(txt, "p_star_residual = {:.12e}", p_star_residual).expect("write");
    writeln!(txt).expect("write");

    writeln!(txt, "[test_4_monte_carlo_real_matter]").expect("write");
    writeln!(txt, "freq_traveler_present = {:.12}", freq_traveler).expect("write");
    writeln!(txt, "freq_target_killed = {:.12}", freq_killed).expect("write");
    writeln!(txt, "monte_half_residual = {:.12e}", monte_residual).expect("write");
    writeln!(txt, "expected_total_mass_kg = {:.12}", total_mass_expected).expect("write");
    writeln!(txt, "avg_total_mass_kg = {:.12}", avg_total_mass).expect("write");
    writeln!(txt, "max_mass_error_kg = {:.12e}", max_mass_err).expect("write");
    writeln!(txt, "expected_local_mass_kg = {:.12}", expected_local_mass).expect("write");
    writeln!(txt, "expected_channel_mass_kg = {:.12}", expected_channel_mass).expect("write");

    let payload = json!({
        "inputs": {
            "m_ancestor_kg": m_ancestor_kg,
            "m_traveler_kg": m_traveler_kg,
            "p_input": p_input,
            "samples": samples,
            "seed": seed
        },
        "test_1_single_history": {
            "valid_assignment_count": single_valid.len(),
            "forced_paradox_assignment_count": single_forced_paradox.len(),
            "forced_paradox_consistent": single_forced_consistent,
            "valid_assignments": single_valid.iter().map(|a| json!({
                "a0": a.a0, "t": a.t, "k": a.k, "a1": a.a1
            })).collect::<Vec<_>>()
        },
        "test_2_branch_split": {
            "valid_assignment_count": branch_valid.len(),
            "paradox_style_assignment_count": branch_paradox_style.len(),
            "paradox_style_consistent": branch_paradox_consistent,
            "paradox_style_samples": branch_paradox_style.iter().map(|a| json!({
                "origin": {"a0": a.o_a0, "t": a.o_t, "a1": a.o_a1},
                "target": {"a0": a.t_a0, "k": a.t_k, "a1": a.t_a1}
            })).collect::<Vec<_>>()
        },
        "test_3_deutsch_fixed_point": {
            "p_input_next": p_next_input,
            "input_not_residual": input_residual,
            "p_star": p_star,
            "p_star_residual": p_star_residual
        },
        "test_4_monte_carlo_real_matter": {
            "freq_traveler_present": freq_traveler,
            "freq_target_killed": freq_killed,
            "monte_half_residual": monte_residual,
            "expected_total_mass_kg": total_mass_expected,
            "avg_total_mass_kg": avg_total_mass,
            "max_mass_error_kg": max_mass_err,
            "expected_local_mass_kg": expected_local_mass,
            "expected_channel_mass_kg": expected_channel_mass
        }
    });
    fs::write(
        &json_path,
        serde_json::to_string_pretty(&payload).expect("encode json"),
    )
    .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "single_forced_consistent={} branch_paradox_consistent={} p_star={:.3} freq={:.6}",
        single_forced_consistent, branch_paradox_consistent, p_star, freq_traveler
    );
}
