//! "Prove otherwise" paradox probe.
//!
//! Compares:
//! 1) strict single-history CTC logic (paradox branch unsat),
//! 2) branch-split CTC logic (paradox-style event can be sat).
//!
//! The second lane is a structural "otherwise": it relaxes the single-history
//! closure axiom and only enforces closure on the origin branch.

use serde_json::json;
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

#[derive(Debug, Clone, Copy)]
struct SingleAssignment {
    a0: u8, // ancestor alive before intervention
    t: u8,  // traveler exists
    k: u8,  // kill signal sent
    a1: u8, // ancestor alive after intervention
}

#[derive(Debug, Clone, Copy)]
struct BranchAssignment {
    // Origin branch O (traveler's home branch)
    o_a0: u8,
    o_t: u8,
    o_a1: u8,
    // Target branch T (where intervention happens)
    t_a0: u8,
    t_k: u8,
    t_a1: u8,
}

fn main() {
    let out_dir = std::env::var("GUTOE_CTC_PARADOX_OTHERWISE_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ctc_paradox_otherwise_probe".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    // ── Lane 1: strict single-history (forced traveler branch) ──────────────
    let mut single_valid = Vec::new();
    for a0 in [0u8, 1u8] {
        for t in [0u8, 1u8] {
            for k in [0u8, 1u8] {
                for a1 in [0u8, 1u8] {
                    let c1 = t == a0; // traveler iff ancestor line exists
                    let c2 = k == t; // traveler always attempts kill
                    let c3 = a1 == (a0 * (1 - k)); // kill flips survival
                    let c4 = a1 == a0; // single-history closure
                    if c1 && c2 && c3 && c4 {
                        single_valid.push(SingleAssignment { a0, t, k, a1 });
                    }
                }
            }
        }
    }
    let single_forced_traveler: Vec<_> = single_valid.iter().copied().filter(|a| a.t == 1).collect();

    // ── Lane 2: branch-split "otherwise" model ──────────────────────────────
    // Constraints:
    //   O1: o_t = o_a0                     (traveler exists on origin branch)
    //   O2: o_a1 = o_a0                    (closure enforced only on origin branch)
    //   T1: t_k = o_t                      (traveler acts in target branch)
    //   T2: t_a1 = t_a0 * (1 - t_k)        (kill can remove target ancestor)
    // Optional paradox-style target event:
    //   t_a0 = 1, o_t = 1, t_k = 1, t_a1 = 0
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
                                branch_valid.push(BranchAssignment {
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

    let paradox_style_branch: Vec<_> = branch_valid
        .iter()
        .copied()
        .filter(|a| a.o_t == 1 && a.t_a0 == 1 && a.t_k == 1 && a.t_a1 == 0 && a.o_a0 == 1)
        .collect();

    let txt_path = out.join("ctc_paradox_otherwise_probe.txt");
    let json_path = out.join("ctc_paradox_otherwise_probe.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[single_history]").expect("write");
    writeln!(txt, "valid_assignments = {}", single_valid.len()).expect("write");
    writeln!(
        txt,
        "forced_traveler_valid_assignments = {}",
        single_forced_traveler.len()
    )
    .expect("write");
    writeln!(
        txt,
        "forced_traveler_consistent = {}",
        !single_forced_traveler.is_empty()
    )
    .expect("write");
    for (i, a) in single_forced_traveler.iter().enumerate() {
        writeln!(txt, "  {}: a0={} t={} k={} a1={}", i, a.a0, a.t, a.k, a.a1).expect("write");
    }
    writeln!(txt).expect("write");

    writeln!(txt, "[branch_split_otherwise]").expect("write");
    writeln!(txt, "valid_assignments = {}", branch_valid.len()).expect("write");
    writeln!(
        txt,
        "paradox_style_assignments = {}",
        paradox_style_branch.len()
    )
    .expect("write");
    writeln!(
        txt,
        "paradox_style_consistent = {}",
        !paradox_style_branch.is_empty()
    )
    .expect("write");
    for (i, a) in paradox_style_branch.iter().take(8).enumerate() {
        writeln!(
            txt,
            "  {}: O(a0={},t={},a1={}) | T(a0={},k={},a1={})",
            i, a.o_a0, a.o_t, a.o_a1, a.t_a0, a.t_k, a.t_a1
        )
        .expect("write");
    }

    let payload = json!({
        "single_history": {
            "valid_assignment_count": single_valid.len(),
            "forced_traveler_valid_assignment_count": single_forced_traveler.len(),
            "forced_traveler_consistent": !single_forced_traveler.is_empty(),
        },
        "branch_split_otherwise": {
            "valid_assignment_count": branch_valid.len(),
            "paradox_style_assignment_count": paradox_style_branch.len(),
            "paradox_style_consistent": !paradox_style_branch.is_empty(),
            "sample": paradox_style_branch.iter().take(8).map(|a| json!({
                "origin": {"a0": a.o_a0, "t": a.o_t, "a1": a.o_a1},
                "target": {"a0": a.t_a0, "k": a.t_k, "a1": a.t_a1}
            })).collect::<Vec<_>>()
        }
    });
    fs::write(&json_path, serde_json::to_string_pretty(&payload).expect("encode json"))
        .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "single_forced_traveler_consistent={} branch_paradox_consistent={}",
        !single_forced_traveler.is_empty(),
        !paradox_style_branch.is_empty()
    );
}

