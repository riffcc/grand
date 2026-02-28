use serde_json::Value;
use std::env;
use std::fs;
use std::path::{Path, PathBuf};

#[derive(Clone, Debug)]
struct Candidate {
    rank: usize,
    z: u16,
    symbol: String,
    stable_like: usize,
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
        });
    }
    out.sort_by_key(|c| c.rank);
    Ok(out)
}

#[derive(Clone, Debug)]
struct Row {
    symbol: String,
    z: u16,
    stable_like: usize,
    tc_ideal_k: f64,
    max_uniform_penalty: f64,
    max_strain_at_defect_10: f64,
    max_defect_at_strain_10: f64,
    pass_frac_grid: f64,
}

fn main() -> Result<(), String> {
    let out_dir = PathBuf::from(
        env::var("GUTOE_RTSC_ROBUST_OUT")
            .unwrap_or_else(|_| "/tmp/bh_renders/rtsc_robustness_sweep".to_string()),
    );
    fs::create_dir_all(&out_dir)
        .map_err(|e| format!("failed to create out dir {}: {e}", out_dir.display()))?;

    let witness_path = PathBuf::from(
        env::var("GUTOE_RTSC_WITNESS_JSON").unwrap_or_else(|_| {
            "/tmp/bh_renders/rtsc_witness_candidates/rtsc_forced_witnesses.json".to_string()
        }),
    );
    let candidates = load_candidates(&witness_path)?;

    let kernel_base = 11.0 / 48.0;

    let strains: Vec<f64> = (0..=50).map(|i| i as f64 / 100.0).collect();
    let defects: Vec<f64> = (0..=50).map(|i| i as f64 / 100.0).collect();

    let mut rows = Vec::<Row>::new();
    for c in candidates {
        let element_factor = 1.0 + 0.02 * (c.stable_like as f64 - 5.0);
        let kernel = kernel_base * element_factor;
        let tc_ideal = 300.0 * (1.0 + kernel);

        // Uniform degradation p with Tc_eff = Tc_ideal * (1 - p)
        let max_uniform_penalty = (1.0 - 300.0 / tc_ideal).clamp(0.0, 1.0);

        // Tc_eff under separable strain/defect penalties.
        let tc_eff = |s: f64, d: f64| tc_ideal * (1.0 - s).max(0.0) * (1.0 - d).max(0.0);

        let d_fixed = 0.10;
        let s_fixed = 0.10;
        let max_strain_at_defect_10 = (1.0 - 300.0 / (tc_ideal * (1.0 - d_fixed))).clamp(0.0, 1.0);
        let max_defect_at_strain_10 = (1.0 - 300.0 / (tc_ideal * (1.0 - s_fixed))).clamp(0.0, 1.0);

        let mut pass = 0usize;
        let mut total = 0usize;
        for &s in &strains {
            for &d in &defects {
                total += 1;
                if tc_eff(s, d) >= 300.0 {
                    pass += 1;
                }
            }
        }

        rows.push(Row {
            symbol: c.symbol,
            z: c.z,
            stable_like: c.stable_like,
            tc_ideal_k: tc_ideal,
            max_uniform_penalty,
            max_strain_at_defect_10,
            max_defect_at_strain_10,
            pass_frac_grid: pass as f64 / total as f64,
        });
    }

    rows.sort_by(|a, b| b.tc_ideal_k.partial_cmp(&a.tc_ideal_k).unwrap());

    let mut txt = String::new();
    txt.push_str("[rtsc_witness_robustness_sweep]\n");
    txt.push_str(&format!("witness_source = {}\n", witness_path.display()));
    txt.push_str("kernel_base_expr = 11/48\n");
    txt.push_str("model = Tc_eff = Tc_ideal * (1-strain) * (1-defect)\n");
    txt.push_str("grid = strain,defect in [0,0.50], step=0.01\n\n");
    txt.push_str("symbol,Z,stable_like,tc_ideal_k,max_uniform_penalty,max_strain_at_defect_10,max_defect_at_strain_10,pass_fraction_grid\n");
    for r in &rows {
        txt.push_str(&format!(
            "{},{},{},{:.6},{:.6},{:.6},{:.6},{:.6}\n",
            r.symbol,
            r.z,
            r.stable_like,
            r.tc_ideal_k,
            r.max_uniform_penalty,
            r.max_strain_at_defect_10,
            r.max_defect_at_strain_10,
            r.pass_frac_grid
        ));
    }

    let avg_pass_frac = rows.iter().map(|r| r.pass_frac_grid).sum::<f64>() / rows.len().max(1) as f64;
    txt.push_str("\n[summary]\n");
    txt.push_str(&format!("candidate_count = {}\n", rows.len()));
    txt.push_str(&format!("avg_pass_fraction = {:.6}\n", avg_pass_frac));

    let mut json = String::new();
    json.push_str("{\n");
    json.push_str(&format!(
        "  \"meta\": {{\"kernel_base\": {:.12e}, \"witness_source\": \"{}\", \"model\": \"Tc_eff = Tc_ideal*(1-strain)*(1-defect)\"}},\n",
        kernel_base,
        witness_path.display()
    ));
    json.push_str(&format!(
        "  \"summary\": {{\"candidate_count\": {}, \"avg_pass_fraction\": {:.12e}}},\n",
        rows.len(),
        avg_pass_frac
    ));
    json.push_str("  \"rows\": [\n");
    for (i, r) in rows.iter().enumerate() {
        json.push_str(&format!(
            "    {{\"symbol\":\"{}\",\"z\":{},\"stable_like\":{},\"tc_ideal_k\":{:.12e},\"max_uniform_penalty\":{:.12e},\"max_strain_at_defect_10\":{:.12e},\"max_defect_at_strain_10\":{:.12e},\"pass_fraction_grid\":{:.12e}}}{}\n",
            r.symbol,
            r.z,
            r.stable_like,
            r.tc_ideal_k,
            r.max_uniform_penalty,
            r.max_strain_at_defect_10,
            r.max_defect_at_strain_10,
            r.pass_frac_grid,
            if i + 1 == rows.len() { "" } else { "," }
        ));
    }
    json.push_str("  ]\n}\n");

    let txt_path = out_dir.join("rtsc_witness_robustness_sweep.txt");
    let json_path = out_dir.join("rtsc_witness_robustness_sweep.json");
    fs::write(&txt_path, txt)
        .map_err(|e| format!("failed to write {}: {e}", txt_path.display()))?;
    fs::write(&json_path, json)
        .map_err(|e| format!("failed to write {}: {e}", json_path.display()))?;

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!("sweep_complete candidates={}", rows.len());

    Ok(())
}
