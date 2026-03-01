//! Even-subalgebra suppression probe for CTC fan-in gain.
//!
//! Tests the hypothesis:
//! - Cl(1,3) even-grade filter gives suppression exactly 1/2 (8 of 16 basis states),
//! - canonical gain lane `branching=2, merge=1` then gives `G_eff = 1`,
//! - measured uncapped gain (`1.9992`) needs suppression near 1/2.

use serde_json::json;
use std::fs;
use std::path::PathBuf;

fn env_f64(name: &str, default: f64) -> f64 {
    std::env::var(name)
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(default)
}

fn main() {
    let out_dir = std::env::var("GUTOE_CTC_EVEN_SUPPRESSION_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ctc_even_subalgebra_suppression_probe".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    // Structural counts in Cl(1,3): grades 0..4 => 1,4,6,4,1.
    let even_count = 1.0 + 6.0 + 1.0; // grades 0,2,4
    let odd_count = 4.0 + 4.0; // grades 1,3
    let total_count = even_count + odd_count;
    let suppression_even = even_count / total_count; // 1/2

    // Canonical branch/merge lane from the user hypothesis.
    let branching_canonical = env_f64("GUTOE_CTC_EVEN_BRANCHING", 2.0);
    let merge_canonical = env_f64("GUTOE_CTC_EVEN_MERGE", 1.0);
    let geff_canonical = branching_canonical * merge_canonical * suppression_even;

    // Measured uncapped lane (from origin closure discussion).
    let g_uncapped = env_f64("GUTOE_CTC_EVEN_G_UNCAPPED", 1.9992);
    let g_target = env_f64("GUTOE_CTC_EVEN_G_TARGET", 1.000_000_000_000_010_2);

    let suppression_for_unit = if g_uncapped > 0.0 {
        1.0 / g_uncapped
    } else {
        f64::NAN
    };
    let suppression_for_target = if g_uncapped > 0.0 {
        g_target / g_uncapped
    } else {
        f64::NAN
    };
    let half = 0.5_f64;
    let unit_offset_from_half = suppression_for_unit - half;
    let target_offset_from_half = suppression_for_target - half;

    let payload = json!({
        "structural_counts": {
            "even_count": even_count,
            "odd_count": odd_count,
            "total_count": total_count,
            "suppression_even": suppression_even
        },
        "canonical_lane": {
            "branching": branching_canonical,
            "merge": merge_canonical,
            "geff": geff_canonical
        },
        "uncapped_lane": {
            "g_uncapped": g_uncapped,
            "g_target": g_target,
            "suppression_for_unit_gain": suppression_for_unit,
            "suppression_for_target_gain": suppression_for_target,
            "half": half,
            "unit_offset_from_half": unit_offset_from_half,
            "target_offset_from_half": target_offset_from_half
        }
    });

    let txt_path = out.join("ctc_even_subalgebra_suppression_probe.txt");
    let json_path = out.join("ctc_even_subalgebra_suppression_probe.json");

    let mut txt = String::new();
    txt.push_str("[ctc_even_subalgebra_suppression_probe]\n");
    txt.push_str(&format!(
        "even_count={:.0}, odd_count={:.0}, total_count={:.0}, suppression_even={:.12e}\n",
        even_count, odd_count, total_count, suppression_even
    ));
    txt.push_str(&format!(
        "canonical_geff = branching*merge*suppression = {:.12e}*{:.12e}*{:.12e} = {:.12e}\n",
        branching_canonical, merge_canonical, suppression_even, geff_canonical
    ));
    txt.push_str(&format!(
        "uncapped: g_uncapped={:.12e}, g_target={:.12e}\n",
        g_uncapped, g_target
    ));
    txt.push_str(&format!(
        "suppression_for_unit_gain = 1/g_uncapped = {:.12e}\n",
        suppression_for_unit
    ));
    txt.push_str(&format!(
        "suppression_for_target_gain = g_target/g_uncapped = {:.12e}\n",
        suppression_for_target
    ));
    txt.push_str(&format!(
        "offset_from_half: unit={:+.12e}, target={:+.12e}\n",
        unit_offset_from_half, target_offset_from_half
    ));

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(&json_path, serde_json::to_string_pretty(&payload).expect("json")).expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
}

