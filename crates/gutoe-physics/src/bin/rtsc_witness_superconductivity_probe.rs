use serde_json::Value;
use std::cmp::Ordering;
use std::env;
use std::fs;
use std::path::{Path, PathBuf};

#[derive(Clone, Debug)]
struct Candidate {
    rank: usize,
    z: u16,
    symbol: String,
    stable_like: usize,
    lattice_distance: u32,
}

fn load_candidates(path: &Path) -> Result<Vec<Candidate>, String> {
    let txt = fs::read_to_string(path)
        .map_err(|e| format!("failed to read witness json {}: {e}", path.display()))?;
    let v: Value = serde_json::from_str(&txt)
        .map_err(|e| format!("failed to parse witness json {}: {e}", path.display()))?;
    let arr = v
        .get("candidates")
        .and_then(|x| x.as_array())
        .ok_or_else(|| "missing candidates[] in witness json".to_string())?;

    let mut out = Vec::new();
    for item in arr {
        out.push(Candidate {
            rank: item.get("rank").and_then(Value::as_u64).unwrap_or(0) as usize,
            z: item.get("z").and_then(Value::as_u64).unwrap_or(0) as u16,
            symbol: item
                .get("symbol")
                .and_then(Value::as_str)
                .unwrap_or("?")
                .to_string(),
            stable_like: item
                .get("stable_like")
                .and_then(Value::as_u64)
                .unwrap_or(0) as usize,
            lattice_distance: item
                .get("lattice_distance")
                .and_then(Value::as_u64)
                .unwrap_or(99) as u32,
        });
    }
    out.sort_by_key(|c| c.rank);
    Ok(out)
}

#[derive(Clone, Debug)]
struct ProbeRow {
    symbol: String,
    z: u16,
    stable_like: usize,
    lattice_distance: u32,
    kernel_base: f64,
    kernel_element: f64,
    tc_natural_k: f64,
    tc_engineered_sc_k: f64,
    natural_pass_300k: bool,
    engineered_pass_300k: bool,
}

fn main() -> Result<(), String> {
    let out_dir = PathBuf::from(
        env::var("GUTOE_RTSC_PROBE_OUT")
            .unwrap_or_else(|_| "/tmp/bh_renders/rtsc_superconductivity_probe".to_string()),
    );
    fs::create_dir_all(&out_dir)
        .map_err(|e| format!("failed to create out dir {}: {e}", out_dir.display()))?;

    let witness_path = PathBuf::from(
        env::var("GUTOE_RTSC_WITNESS_JSON").unwrap_or_else(|_| {
            "/tmp/bh_renders/rtsc_witness_candidates/rtsc_forced_witnesses.json".to_string()
        }),
    );
    let candidates = load_candidates(&witness_path)?;

    // Shared structural gate kernel from Lean lane:
    // pairingKernelQ = 11/48, tcStructuralQ = 300*(1+11/48) = 368.75 K.
    let kernel_base = 11.0 / 48.0;

    let mut rows = Vec::<ProbeRow>::new();
    for c in candidates {
        // Witness-level element factor (no free fit): centered on stable-like triplet floor.
        let element_factor = 1.0 + 0.02 * (c.stable_like as f64 - 5.0);
        let kernel_element = kernel_base * element_factor;

        // Natural-phase pass requires already-simple-cubic (distance 0).
        let tc_natural_k = if c.lattice_distance == 0 {
            300.0 * (1.0 + kernel_element)
        } else {
            0.0
        };

        // Engineered phase (epitaxial/high-pressure forced SC geometry).
        let tc_engineered_sc_k = 300.0 * (1.0 + kernel_element);

        rows.push(ProbeRow {
            symbol: c.symbol,
            z: c.z,
            stable_like: c.stable_like,
            lattice_distance: c.lattice_distance,
            kernel_base,
            kernel_element,
            tc_natural_k,
            tc_engineered_sc_k,
            natural_pass_300k: tc_natural_k >= 300.0,
            engineered_pass_300k: tc_engineered_sc_k >= 300.0,
        });
    }

    rows.sort_by(|a, b| {
        b.tc_engineered_sc_k
            .partial_cmp(&a.tc_engineered_sc_k)
            .unwrap_or(Ordering::Equal)
            .then_with(|| a.z.cmp(&b.z))
    });

    let mut txt = String::new();
    txt.push_str("[rtsc_witness_superconductivity_probe]\n");
    txt.push_str("mode = natural_vs_engineered_simple_cubic\n");
    txt.push_str("kernel_base_expr = 11/48\n");
    txt.push_str(&format!("witness_source = {}\n\n", witness_path.display()));

    txt.push_str("symbol,Z,stable_like,lattice_distance,kernel_element,tc_natural_k,tc_engineered_sc_k,natural_pass_300k,engineered_pass_300k\n");
    for r in &rows {
        txt.push_str(&format!(
            "{},{},{},{},{:.9},{:.6},{:.6},{},{}\n",
            r.symbol,
            r.z,
            r.stable_like,
            r.lattice_distance,
            r.kernel_element,
            r.tc_natural_k,
            r.tc_engineered_sc_k,
            r.natural_pass_300k,
            r.engineered_pass_300k
        ));
    }

    let natural_pass = rows.iter().filter(|r| r.natural_pass_300k).count();
    let engineered_pass = rows.iter().filter(|r| r.engineered_pass_300k).count();

    txt.push_str("\n[summary]\n");
    txt.push_str(&format!("natural_pass_count = {}\n", natural_pass));
    txt.push_str(&format!("engineered_pass_count = {}\n", engineered_pass));
    txt.push_str("interpretation = natural_phase_requires_existing_simple_cubic; engineered_phase_applies_forced_lattice gate\n");

    let mut json = String::new();
    json.push_str("{\n");
    json.push_str(&format!(
        "  \"meta\": {{\"mode\": \"natural_vs_engineered_simple_cubic\", \"kernel_base\": {:.12e}, \"witness_source\": \"{}\"}},\n",
        kernel_base,
        witness_path.display()
    ));
    json.push_str(&format!(
        "  \"summary\": {{\"natural_pass_count\": {}, \"engineered_pass_count\": {}}},\n",
        natural_pass, engineered_pass
    ));
    json.push_str("  \"rows\": [\n");
    for (i, r) in rows.iter().enumerate() {
        json.push_str(&format!(
            "    {{\"symbol\": \"{}\", \"z\": {}, \"stable_like\": {}, \"lattice_distance\": {}, \"kernel_base\": {:.12e}, \"kernel_element\": {:.12e}, \"tc_natural_k\": {:.12e}, \"tc_engineered_sc_k\": {:.12e}, \"natural_pass_300k\": {}, \"engineered_pass_300k\": {}}}{}\n",
            r.symbol,
            r.z,
            r.stable_like,
            r.lattice_distance,
            r.kernel_base,
            r.kernel_element,
            r.tc_natural_k,
            r.tc_engineered_sc_k,
            if r.natural_pass_300k { "true" } else { "false" },
            if r.engineered_pass_300k { "true" } else { "false" },
            if i + 1 == rows.len() { "" } else { "," }
        ));
    }
    json.push_str("  ]\n}\n");

    let txt_path = out_dir.join("rtsc_witness_superconductivity_probe.txt");
    let json_path = out_dir.join("rtsc_witness_superconductivity_probe.json");
    fs::write(&txt_path, txt)
        .map_err(|e| format!("failed to write {}: {e}", txt_path.display()))?;
    fs::write(&json_path, json)
        .map_err(|e| format!("failed to write {}: {e}", json_path.display()))?;

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!("natural_pass={} engineered_pass={}", natural_pass, engineered_pass);

    Ok(())
}
