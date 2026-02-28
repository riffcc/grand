use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::Value;
use std::cmp::Ordering;
use std::env;
use std::fs;
use std::path::{Path, PathBuf};

#[derive(Clone, Debug)]
struct PlanRow {
    rank: usize,
    symbol: String,
    z: u16,
    native_phase: String,
    tc_ideal_k: f64,
    max_uniform_penalty: f64,
    engineering_difficulty: f64,
    hazard: f64,
}

fn load_plan(path: &Path) -> Result<Vec<PlanRow>, String> {
    let txt = fs::read_to_string(path)
        .map_err(|e| format!("failed to read routes json {}: {e}", path.display()))?;
    let v: Value = serde_json::from_str(&txt)
        .map_err(|e| format!("failed to parse routes json {}: {e}", path.display()))?;
    let arr = v
        .get("plan")
        .and_then(Value::as_array)
        .ok_or_else(|| "missing plan[]".to_string())?;
    let mut out = Vec::new();
    for p in arr {
        out.push(PlanRow {
            rank: p.get("rank").and_then(Value::as_u64).unwrap_or(0) as usize,
            symbol: p
                .get("symbol")
                .and_then(Value::as_str)
                .unwrap_or("?")
                .to_string(),
            z: p.get("z").and_then(Value::as_u64).unwrap_or(0) as u16,
            native_phase: p
                .get("native_phase")
                .and_then(Value::as_str)
                .unwrap_or("unknown")
                .to_string(),
            tc_ideal_k: p.get("tc_ideal_k").and_then(Value::as_f64).unwrap_or(0.0),
            max_uniform_penalty: p
                .get("max_uniform_penalty")
                .and_then(Value::as_f64)
                .unwrap_or(0.0),
            engineering_difficulty: p
                .get("engineering_difficulty")
                .and_then(Value::as_f64)
                .unwrap_or(1.0),
            hazard: p.get("hazard").and_then(Value::as_f64).unwrap_or(0.5),
        });
    }
    out.sort_by_key(|r| r.rank);
    Ok(out)
}

#[derive(Clone, Copy)]
struct Regime {
    name: &'static str,
    control: f64, // 0..1
}

#[derive(Clone, Debug)]
struct RowOut {
    regime: String,
    symbol: String,
    z: u16,
    p_simple_cubic: f64,
    p_tc_ge_300: f64,
    p_validated: f64,
    mean_tc_sc: f64,
    expected_validated_per_100: f64,
}

fn tri01(rng: &mut StdRng) -> f64 {
    (rng.gen::<f64>() + rng.gen::<f64>() + rng.gen::<f64>()) / 3.0
}

fn clamp01(x: f64) -> f64 {
    x.clamp(0.0, 1.0)
}

fn sigmoid(x: f64) -> f64 {
    1.0 / (1.0 + (-x).exp())
}

fn base_mismatch(native_phase: &str) -> f64 {
    match native_phase {
        "bcc" => 0.10,
        "fcc" => 0.16,
        "hcp" => 0.28,
        _ => 0.25,
    }
}

fn simulate_row(rng: &mut StdRng, p: &PlanRow, regime: Regime, trials: usize) -> RowOut {
    let d = p.engineering_difficulty;
    let h = p.hazard;
    let c = regime.control;

    let mismatch0 = base_mismatch(&p.native_phase);
    let mismatch = (mismatch0 * (1.0 - 0.30 * c)).max(0.02);

    let mut n_sc = 0usize;
    let mut n_tc300 = 0usize;
    let mut n_valid = 0usize;
    let mut tc_sum_sc = 0.0f64;

    for _ in 0..trials {
        let q_epi = clamp01(0.25 + 0.75 * tri01(rng) * c + 0.15 * (1.0 - d));
        let strain_lock = clamp01(0.20 + 0.80 * tri01(rng) * c - 0.10 * d);
        let oxygen = clamp01(0.18 * (1.0 - c) + 0.25 * (1.0 - tri01(rng)) + 0.10 * h);
        let thickness_dev = clamp01(0.20 * (1.0 - c) + 0.20 * (1.0 - tri01(rng)) + 0.10 * d);
        let cap_integrity = clamp01(0.35 + 0.65 * tri01(rng) * c - 0.05 * h);
        let pressure_assist = if p.native_phase == "hcp" {
            clamp01(0.10 + 0.90 * tri01(rng) * c)
        } else {
            clamp01(0.50 + 0.50 * tri01(rng) * c)
        };

        let score = 1.60 * q_epi
            + 1.25 * strain_lock
            + 0.90 * cap_integrity
            + 0.60 * pressure_assist
            - 1.15 * d
            - 0.95 * mismatch
            - 0.80 * oxygen
            - 0.55 * thickness_dev
            + 0.15 * p.max_uniform_penalty;

        let p_sc = sigmoid(2.0 * (score - 1.6));
        let is_sc = rng.gen::<f64>() < p_sc;
        if !is_sc {
            continue;
        }
        n_sc += 1;

        let strain_pen = (0.01
            + 0.10 * d
            + 0.05 * (1.0 - strain_lock)
            + 0.05 * mismatch
            + 0.03 * thickness_dev
            + 0.03 * oxygen
            - 0.06 * pressure_assist)
            .clamp(0.0, 0.75);

        let defect_pen = (0.02
            + 0.09 * d
            + 0.06 * oxygen
            + 0.05 * (1.0 - cap_integrity)
            + 0.03 * (1.0 - q_epi)
            + 0.03 * thickness_dev)
            .clamp(0.0, 0.75);

        let tc_eff = p.tc_ideal_k * (1.0 - strain_pen) * (1.0 - defect_pen);
        tc_sum_sc += tc_eff;

        if tc_eff >= 300.0 {
            n_tc300 += 1;
        }

        let p_zero_r = sigmoid((tc_eff - 300.0) / 7.0) * (1.0 - 0.35 * defect_pen);
        let p_meissner = sigmoid((tc_eff - 300.0) / 9.0) * (1.0 - 0.25 * defect_pen);
        let validated = tc_eff >= 300.0
            && rng.gen::<f64>() < p_zero_r.clamp(0.0, 1.0)
            && rng.gen::<f64>() < p_meissner.clamp(0.0, 1.0);
        if validated {
            n_valid += 1;
        }
    }

    let p_simple_cubic = n_sc as f64 / trials as f64;
    let p_tc_ge_300 = n_tc300 as f64 / trials as f64;
    let p_validated = n_valid as f64 / trials as f64;
    let mean_tc_sc = if n_sc > 0 { tc_sum_sc / n_sc as f64 } else { 0.0 };

    RowOut {
        regime: regime.name.to_string(),
        symbol: p.symbol.clone(),
        z: p.z,
        p_simple_cubic,
        p_tc_ge_300,
        p_validated,
        mean_tc_sc,
        expected_validated_per_100: 100.0 * p_validated,
    }
}

fn main() -> Result<(), String> {
    let out_dir = PathBuf::from(
        env::var("GUTOE_RTSC_END2END_OUT")
            .unwrap_or_else(|_| "/tmp/bh_renders/rtsc_end_to_end_campaign".to_string()),
    );
    fs::create_dir_all(&out_dir)
        .map_err(|e| format!("failed to create out dir {}: {e}", out_dir.display()))?;

    let plan_path = PathBuf::from(
        env::var("GUTOE_RTSC_ROUTES_JSON").unwrap_or_else(|_| {
            "/tmp/bh_renders/rtsc_synthesis_routes/rtsc_synthesis_routes.json".to_string()
        }),
    );

    let trials = env::var("GUTOE_RTSC_TRIALS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(30_000)
        .max(1000);
    let seed = env::var("GUTOE_RTSC_SEED")
        .ok()
        .and_then(|s| s.parse::<u64>().ok())
        .unwrap_or(137_106_107);

    let plan = load_plan(&plan_path)?;
    let regimes = [
        Regime {
            name: "baseline_lab",
            control: 0.55,
        },
        Regime {
            name: "tight_process",
            control: 0.75,
        },
        Regime {
            name: "heroic_process",
            control: 0.90,
        },
    ];

    let mut rows = Vec::<RowOut>::new();
    let mut rng = StdRng::seed_from_u64(seed);
    for reg in regimes {
        for p in &plan {
            rows.push(simulate_row(&mut rng, p, reg, trials));
        }
    }

    let mut txt = String::new();
    txt.push_str("[rtsc_end_to_end_campaign]\n");
    txt.push_str(&format!("routes_source = {}\n", plan_path.display()));
    txt.push_str(&format!("trials_per_candidate = {}\n", trials));
    txt.push_str(&format!("seed = {}\n", seed));
    txt.push_str("pipeline = growth -> simple_cubic_capture -> Tc_degradation -> zeroR+Meissner validation\n\n");
    txt.push_str("regime,symbol,Z,p_simple_cubic,p_tc_ge_300,p_validated,mean_tc_sc,expected_validated_per_100\n");
    for r in &rows {
        txt.push_str(&format!(
            "{},{},{},{:.6},{:.6},{:.6},{:.3},{:.3}\n",
            r.regime,
            r.symbol,
            r.z,
            r.p_simple_cubic,
            r.p_tc_ge_300,
            r.p_validated,
            r.mean_tc_sc,
            r.expected_validated_per_100
        ));
    }

    txt.push_str("\n[top_by_regime]\n");
    for reg in ["baseline_lab", "tight_process", "heroic_process"] {
        let mut subset: Vec<&RowOut> = rows.iter().filter(|r| r.regime == reg).collect();
        subset.sort_by(|a, b| {
            b.p_validated
                .partial_cmp(&a.p_validated)
                .unwrap_or(Ordering::Equal)
                .then_with(|| a.z.cmp(&b.z))
        });
        if let Some(best) = subset.first() {
            txt.push_str(&format!(
                "{}: {} (Z={}) p_validated={:.4} expected/100={:.2}\n",
                reg,
                best.symbol,
                best.z,
                best.p_validated,
                best.expected_validated_per_100
            ));
        }
    }

    let txt_path = out_dir.join("rtsc_end_to_end_campaign.txt");
    fs::write(&txt_path, txt)
        .map_err(|e| format!("failed to write {}: {e}", txt_path.display()))?;

    let mut json = String::new();
    json.push_str("{\n");
    json.push_str(&format!(
        "  \"meta\": {{\"routes_source\": \"{}\", \"trials_per_candidate\": {}, \"seed\": {}}},\n",
        plan_path.display(),
        trials,
        seed
    ));
    json.push_str("  \"rows\": [\n");
    for (i, r) in rows.iter().enumerate() {
        json.push_str(&format!(
            "    {{\"regime\":\"{}\",\"symbol\":\"{}\",\"z\":{},\"p_simple_cubic\":{:.12e},\"p_tc_ge_300\":{:.12e},\"p_validated\":{:.12e},\"mean_tc_sc\":{:.12e},\"expected_validated_per_100\":{:.12e}}}{}\n",
            r.regime,
            r.symbol,
            r.z,
            r.p_simple_cubic,
            r.p_tc_ge_300,
            r.p_validated,
            r.mean_tc_sc,
            r.expected_validated_per_100,
            if i + 1 == rows.len() { "" } else { "," }
        ));
    }
    json.push_str("  ]\n}\n");

    let json_path = out_dir.join("rtsc_end_to_end_campaign.json");
    fs::write(&json_path, json)
        .map_err(|e| format!("failed to write {}: {e}", json_path.display()))?;

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!("rows={}", rows.len());

    Ok(())
}
