use serde_json::Value;
use std::cmp::Ordering;
use std::env;
use std::fs;
use std::path::{Path, PathBuf};

#[derive(Clone, Debug)]
struct RobustRow {
    symbol: String,
    z: u16,
    tc_ideal_k: f64,
    max_uniform_penalty: f64,
    pass_fraction_grid: f64,
}

fn load_robust(path: &Path) -> Result<Vec<RobustRow>, String> {
    let txt = fs::read_to_string(path)
        .map_err(|e| format!("failed to read robustness json {}: {e}", path.display()))?;
    let v: Value = serde_json::from_str(&txt)
        .map_err(|e| format!("failed to parse robustness json {}: {e}", path.display()))?;
    let arr = v
        .get("rows")
        .and_then(Value::as_array)
        .ok_or_else(|| "missing rows[]".to_string())?;

    let mut out = Vec::new();
    for r in arr {
        out.push(RobustRow {
            symbol: r
                .get("symbol")
                .and_then(Value::as_str)
                .unwrap_or("?")
                .to_string(),
            z: r.get("z").and_then(Value::as_u64).unwrap_or(0) as u16,
            tc_ideal_k: r.get("tc_ideal_k").and_then(Value::as_f64).unwrap_or(0.0),
            max_uniform_penalty: r
                .get("max_uniform_penalty")
                .and_then(Value::as_f64)
                .unwrap_or(0.0),
            pass_fraction_grid: r
                .get("pass_fraction_grid")
                .and_then(Value::as_f64)
                .unwrap_or(0.0),
        });
    }
    Ok(out)
}

#[derive(Clone, Copy)]
struct ProcessMeta {
    native_phase: &'static str,
    engineering_difficulty: f64, // lower is easier
    hazard: f64,                 // lower is safer
    cap: &'static str,
    first_route: &'static str,
}

fn meta(symbol: &str) -> ProcessMeta {
    match symbol {
        // bcc natives are the easiest to push toward cubic symmetry in thin-film lane.
        "Cr" => ProcessMeta {
            native_phase: "bcc",
            engineering_difficulty: 0.25,
            hazard: 0.20,
            cap: "Al2O3",
            first_route: "MBE/sputter epitaxy on MgO(001) with low-T lock + rapid cap",
        },
        "Mo" => ProcessMeta {
            native_phase: "bcc",
            engineering_difficulty: 0.30,
            hazard: 0.15,
            cap: "SiN",
            first_route: "MBE/sputter epitaxy on MgO(001) or SrTiO3(001), thickness sweep 2-20 nm",
        },
        "Pt" => ProcessMeta {
            native_phase: "fcc",
            engineering_difficulty: 0.55,
            hazard: 0.10,
            cap: "Al2O3",
            first_route: "seed-layer templating + epitaxial compression on perovskite(001)",
        },
        "Hf" => ProcessMeta {
            native_phase: "hcp",
            engineering_difficulty: 0.70,
            hazard: 0.35,
            cap: "HfO2",
            first_route: "high-pressure assisted epitaxy + immediate cap to trap metastable phase",
        },
        "Zn" => ProcessMeta {
            native_phase: "hcp",
            engineering_difficulty: 0.75,
            hazard: 0.40,
            cap: "Al2O3",
            first_route: "low-T epitaxy with strain-lock superlattice on cubic template",
        },
        "Cd" => ProcessMeta {
            native_phase: "hcp",
            engineering_difficulty: 0.85,
            hazard: 0.95,
            cap: "Al2O3",
            first_route: "only in sealed UHV lane; use as late validation due toxicity",
        },
        _ => ProcessMeta {
            native_phase: "unknown",
            engineering_difficulty: 0.8,
            hazard: 0.5,
            cap: "Al2O3",
            first_route: "epitaxy",
        },
    }
}

#[derive(Clone)]
struct RouteRow {
    symbol: String,
    z: u16,
    native_phase: &'static str,
    tc_ideal_k: f64,
    margin_k: f64,
    max_uniform_penalty: f64,
    pass_fraction_grid: f64,
    engineering_difficulty: f64,
    hazard: f64,
    score: f64,
    cap: &'static str,
    first_route: &'static str,
}

fn main() -> Result<(), String> {
    let out_dir = PathBuf::from(
        env::var("GUTOE_RTSC_SYNTH_OUT")
            .unwrap_or_else(|_| "/tmp/bh_renders/rtsc_synthesis_routes".to_string()),
    );
    fs::create_dir_all(&out_dir)
        .map_err(|e| format!("failed to create out dir {}: {e}", out_dir.display()))?;

    let robust_path = PathBuf::from(
        env::var("GUTOE_RTSC_ROBUST_JSON").unwrap_or_else(|_| {
            "/tmp/bh_renders/rtsc_robustness_sweep/rtsc_witness_robustness_sweep.json".to_string()
        }),
    );

    let robust_rows = load_robust(&robust_path)?;
    let mut plan = Vec::<RouteRow>::new();

    for r in robust_rows {
        let m = meta(&r.symbol);
        let margin_k = (r.tc_ideal_k - 300.0).max(0.0);
        // Composite score: high margin/robustness, low difficulty/hazard.
        let score = 1.4 * (margin_k / 100.0)
            + 1.2 * r.max_uniform_penalty
            + 0.9 * r.pass_fraction_grid
            - 0.9 * m.engineering_difficulty
            - 0.8 * m.hazard;
        plan.push(RouteRow {
            symbol: r.symbol,
            z: r.z,
            native_phase: m.native_phase,
            tc_ideal_k: r.tc_ideal_k,
            margin_k,
            max_uniform_penalty: r.max_uniform_penalty,
            pass_fraction_grid: r.pass_fraction_grid,
            engineering_difficulty: m.engineering_difficulty,
            hazard: m.hazard,
            score,
            cap: m.cap,
            first_route: m.first_route,
        });
    }

    plan.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap_or(Ordering::Equal));

    let mut txt = String::new();
    txt.push_str("[rtsc_synthesis_routes]\n");
    txt.push_str(&format!("robustness_source = {}\n", robust_path.display()));
    txt.push_str("goal = prioritize first synthesis shots for engineered simple-cubic RTSC tests\n\n");

    txt.push_str("rank,symbol,Z,native_phase,tc_ideal_k,margin_k,max_uniform_penalty,pass_fraction_grid,engineering_difficulty,hazard,score,cap,first_route\n");
    for (i, p) in plan.iter().enumerate() {
        txt.push_str(&format!(
            "{},{},{},{},{:.3},{:.3},{:.6},{:.6},{:.3},{:.3},{:.6},{},{}\n",
            i + 1,
            p.symbol,
            p.z,
            p.native_phase,
            p.tc_ideal_k,
            p.margin_k,
            p.max_uniform_penalty,
            p.pass_fraction_grid,
            p.engineering_difficulty,
            p.hazard,
            p.score,
            p.cap,
            p.first_route
        ));
    }

    let top3: Vec<&RouteRow> = plan.iter().take(3).collect();
    txt.push_str("\n[top3_execution_order]\n");
    for (i, p) in top3.iter().enumerate() {
        txt.push_str(&format!("{}. {} (Z={}) -> {}\n", i + 1, p.symbol, p.z, p.first_route));
    }

    txt.push_str("\n[measurement_protocol]\n");
    txt.push_str("1) In-situ structure: RHEED + LEED during growth, post-growth XRD reciprocal maps\n");
    txt.push_str("2) Local structure: cross-sectional STEM + FFT indexing to confirm simple-cubic ordering\n");
    txt.push_str("3) Transport: 4-probe R(T) from 400K->2K, zero resistance + I-V critical current\n");
    txt.push_str("4) Magnetic: susceptibility/Meissner check to exclude filamentary artifacts\n");

    let mut json = String::new();
    json.push_str("{\n");
    json.push_str(&format!(
        "  \"meta\": {{\"robustness_source\": \"{}\", \"goal\": \"RTSC synthesis prioritization\"}},\n",
        robust_path.display()
    ));
    json.push_str("  \"plan\": [\n");
    for (i, p) in plan.iter().enumerate() {
        json.push_str(&format!(
            "    {{\"rank\": {}, \"symbol\": \"{}\", \"z\": {}, \"native_phase\": \"{}\", \"tc_ideal_k\": {:.12e}, \"margin_k\": {:.12e}, \"max_uniform_penalty\": {:.12e}, \"pass_fraction_grid\": {:.12e}, \"engineering_difficulty\": {:.12e}, \"hazard\": {:.12e}, \"score\": {:.12e}, \"cap\": \"{}\", \"first_route\": \"{}\"}}{}\n",
            i + 1,
            p.symbol,
            p.z,
            p.native_phase,
            p.tc_ideal_k,
            p.margin_k,
            p.max_uniform_penalty,
            p.pass_fraction_grid,
            p.engineering_difficulty,
            p.hazard,
            p.score,
            p.cap,
            p.first_route,
            if i + 1 == plan.len() { "" } else { "," }
        ));
    }
    json.push_str("  ]\n}\n");

    let txt_path = out_dir.join("rtsc_synthesis_routes.txt");
    let json_path = out_dir.join("rtsc_synthesis_routes.json");
    fs::write(&txt_path, txt)
        .map_err(|e| format!("failed to write {}: {e}", txt_path.display()))?;
    fs::write(&json_path, json)
        .map_err(|e| format!("failed to write {}: {e}", json_path.display()))?;

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!("routes={}", plan.len());

    Ok(())
}
