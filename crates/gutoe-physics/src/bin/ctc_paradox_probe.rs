//! CTC paradox probe: deterministic inconsistency vs fixed-point consistency.
//!
//! We test three lanes:
//! 1) Deterministic grandfather map on a closed loop (expected: no fixed point).
//! 2) Brute-force event-logic consistency check on a tiny CTC graph.
//! 3) Deutsch mixed-state fixed point for NOT (expected: p = 0.5).

use serde_json::json;
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

#[derive(Debug, Clone, Copy)]
struct Assignment {
    a0: u8, // ancestor alive before intervention
    t: u8,  // traveler exists
    k: u8,  // kill signal sent
    a1: u8, // ancestor alive after intervention
}

fn main() {
    let out_dir = std::env::var("GUTOE_CTC_PARADOX_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ctc_paradox_probe".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    // ── Lane 1: Direct deterministic map p -> 1-p with loop closure p' = p ──
    let deterministic_fixed_points: Vec<u8> = [0u8, 1u8]
        .into_iter()
        .filter(|&p| (1 - p) == p)
        .collect();
    let deterministic_consistent = !deterministic_fixed_points.is_empty();

    // ── Lane 2: Event-logic brute force on a tiny CTC graph ──────────────────
    // Constraints:
    //   C1: t = a0                      (traveler exists iff ancestor line exists)
    //   C2: k = t                       (if traveler exists, they send kill signal)
    //   C3: a1 = a0 * (1-k)             (kill signal removes ancestor)
    //   C4: a1 = a0                     (CTC loop consistency / identification)
    let mut valid = Vec::new();
    for a0 in [0u8, 1u8] {
        for t in [0u8, 1u8] {
            for k in [0u8, 1u8] {
                for a1 in [0u8, 1u8] {
                    let c1 = t == a0;
                    let c2 = k == t;
                    let c3 = a1 == (a0 * (1 - k));
                    let c4 = a1 == a0;
                    if c1 && c2 && c3 && c4 {
                        valid.push(Assignment { a0, t, k, a1 });
                    }
                }
            }
        }
    }
    let graph_consistent = !valid.is_empty();
    let forced_traveler_valid: Vec<_> = valid.iter().copied().filter(|a| a.t == 1).collect();
    let forced_traveler_consistent = !forced_traveler_valid.is_empty();

    // ── Lane 3: Deutsch mixed-state fixed point for NOT gate ────────────────
    // Classical probability on CTC bit: p' = 1-p.
    // Fixed-point equation: p = 1-p -> p = 0.5.
    let p_star = 0.5_f64;
    let deutsch_residual = (p_star - (1.0 - p_star)).abs();

    let txt_path = out.join("ctc_paradox_probe.txt");
    let json_path = out.join("ctc_paradox_probe.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[deterministic_not_loop]").expect("write");
    writeln!(
        txt,
        "fixed_points = {:?}",
        deterministic_fixed_points
    )
    .expect("write");
    writeln!(txt, "consistent = {}", deterministic_consistent).expect("write");
    writeln!(txt).expect("write");

    writeln!(txt, "[event_graph]").expect("write");
    writeln!(txt, "valid_assignments = {}", valid.len()).expect("write");
    for (i, a) in valid.iter().enumerate() {
        writeln!(
            txt,
            "  {}: a0={} t={} k={} a1={}",
            i, a.a0, a.t, a.k, a.a1
        )
        .expect("write");
    }
    writeln!(txt, "graph_consistent = {}", graph_consistent).expect("write");
    writeln!(
        txt,
        "forced_traveler_valid_assignments = {}",
        forced_traveler_valid.len()
    )
    .expect("write");
    writeln!(
        txt,
        "forced_traveler_consistent = {}",
        forced_traveler_consistent
    )
    .expect("write");
    writeln!(txt).expect("write");

    writeln!(txt, "[deutsch_mixed_state]").expect("write");
    writeln!(txt, "p_star = {:.12}", p_star).expect("write");
    writeln!(txt, "residual = {:.12e}", deutsch_residual).expect("write");
    writeln!(txt, "fixed_point_exists = {}", deutsch_residual < 1e-12).expect("write");

    let payload = json!({
        "deterministic_not_loop": {
            "fixed_points": deterministic_fixed_points,
            "consistent": deterministic_consistent
        },
        "event_graph": {
            "valid_assignment_count": valid.len(),
            "valid_assignments": valid.iter().map(|a| json!({
                "a0": a.a0, "t": a.t, "k": a.k, "a1": a.a1
            })).collect::<Vec<_>>(),
            "graph_consistent": graph_consistent,
            "forced_traveler_valid_assignment_count": forced_traveler_valid.len(),
            "forced_traveler_consistent": forced_traveler_consistent
        },
        "deutsch_mixed_state": {
            "p_star": p_star,
            "residual": deutsch_residual,
            "fixed_point_exists": deutsch_residual < 1e-12
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
        "deterministic_consistent={} graph_consistent={} deutsch_fixed_point={}",
        deterministic_consistent,
        graph_consistent,
        deutsch_residual < 1e-12
    );
}
