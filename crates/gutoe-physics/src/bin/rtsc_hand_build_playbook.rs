use serde_json::Value;
use std::cmp::Ordering;
use std::env;
use std::fs;
use std::path::{Path, PathBuf};

#[derive(Clone, Debug)]
struct PlanRow {
    rank: usize,
    symbol: String,
    native_phase: String,
    engineering_difficulty: f64,
    hazard: f64,
    cap: String,
    first_route: String,
}

#[derive(Clone, Debug)]
struct CampRow {
    regime: String,
    symbol: String,
    p_validated: f64,
}

#[derive(Clone, Debug)]
struct BuildLane {
    symbol: String,
    native_phase: String,
    p_validated: f64,
    engineering_difficulty: f64,
    hazard: f64,
    cap: String,
    first_route: String,
    hand_build_score: f64,
}

fn load_plan(path: &Path) -> Result<Vec<PlanRow>, String> {
    let txt = fs::read_to_string(path)
        .map_err(|e| format!("failed to read {}: {e}", path.display()))?;
    let v: Value = serde_json::from_str(&txt)
        .map_err(|e| format!("failed to parse {}: {e}", path.display()))?;
    let arr = v
        .get("plan")
        .and_then(Value::as_array)
        .ok_or_else(|| "missing plan[]".to_string())?;

    let mut out = Vec::new();
    for row in arr {
        out.push(PlanRow {
            rank: row.get("rank").and_then(Value::as_u64).unwrap_or(0) as usize,
            symbol: row
                .get("symbol")
                .and_then(Value::as_str)
                .unwrap_or("?")
                .to_string(),
            native_phase: row
                .get("native_phase")
                .and_then(Value::as_str)
                .unwrap_or("unknown")
                .to_string(),
            engineering_difficulty: row
                .get("engineering_difficulty")
                .and_then(Value::as_f64)
                .unwrap_or(1.0),
            hazard: row.get("hazard").and_then(Value::as_f64).unwrap_or(1.0),
            cap: row
                .get("cap")
                .and_then(Value::as_str)
                .unwrap_or("Al2O3")
                .to_string(),
            first_route: row
                .get("first_route")
                .and_then(Value::as_str)
                .unwrap_or("")
                .to_string(),
        });
    }
    out.sort_by_key(|r| r.rank);
    Ok(out)
}

fn load_campaign(path: &Path) -> Result<Vec<CampRow>, String> {
    let txt = fs::read_to_string(path)
        .map_err(|e| format!("failed to read {}: {e}", path.display()))?;
    let v: Value = serde_json::from_str(&txt)
        .map_err(|e| format!("failed to parse {}: {e}", path.display()))?;
    let arr = v
        .get("rows")
        .and_then(Value::as_array)
        .ok_or_else(|| "missing rows[]".to_string())?;

    let mut out = Vec::new();
    for row in arr {
        out.push(CampRow {
            regime: row
                .get("regime")
                .and_then(Value::as_str)
                .unwrap_or("")
                .to_string(),
            symbol: row
                .get("symbol")
                .and_then(Value::as_str)
                .unwrap_or("")
                .to_string(),
            p_validated: row
                .get("p_validated")
                .and_then(Value::as_f64)
                .unwrap_or(0.0),
        });
    }
    Ok(out)
}

fn pick_p(rows: &[CampRow], regime: &str, symbol: &str) -> f64 {
    rows.iter()
        .find(|r| r.regime == regime && r.symbol == symbol)
        .map(|r| r.p_validated)
        .unwrap_or(0.0)
}

fn material_recipe(symbol: &str) -> (&'static str, &'static str, &'static str, &'static str, &'static str) {
    match symbol {
        "Cr" => (
            "MgO(001)",
            "DC magnetron sputter (preferred hand-build lane)",
            "150-260 C",
            "0.05-0.12 nm/s",
            "8-18 nm",
        ),
        "Mo" => (
            "MgO(001) or SrTiO3(001)",
            "RF/DC sputter or MBE",
            "180-320 C",
            "0.04-0.10 nm/s",
            "10-22 nm",
        ),
        "Pt" => (
            "Perovskite(001) + seed templating layer",
            "RF sputter + epitaxial compression",
            "120-220 C",
            "0.03-0.08 nm/s",
            "6-16 nm",
        ),
        _ => (
            "MgO(001)",
            "Sputter epitaxy",
            "150-280 C",
            "0.03-0.10 nm/s",
            "8-20 nm",
        ),
    }
}

fn log_choose(n: usize, k: usize) -> f64 {
    if k > n {
        return f64::NEG_INFINITY;
    }
    let k = k.min(n - k);
    let mut s = 0.0;
    for i in 0..k {
        s += ((n - i) as f64).ln() - ((i + 1) as f64).ln();
    }
    s
}

fn binom_pmf(n: usize, k: usize, p: f64) -> f64 {
    if !(0.0..=1.0).contains(&p) {
        return 0.0;
    }
    if p == 0.0 {
        return if k == 0 { 1.0 } else { 0.0 };
    }
    if p == 1.0 {
        return if k == n { 1.0 } else { 0.0 };
    }
    let l = log_choose(n, k) + (k as f64) * p.ln() + ((n - k) as f64) * (1.0 - p).ln();
    l.exp()
}

fn binom_tail_geq(n: usize, k: usize, p: f64) -> f64 {
    if k == 0 {
        return 1.0;
    }
    if k > n {
        return 0.0;
    }
    (k..=n).map(|i| binom_pmf(n, i, p)).sum()
}

fn n_for_target_successes(p: f64, successes: usize, conf: f64) -> usize {
    let mut n = successes.max(1);
    while n < 10_000 {
        if binom_tail_geq(n, successes, p) >= conf {
            return n;
        }
        n += 1;
    }
    n
}

fn main() -> Result<(), String> {
    let out_dir = PathBuf::from(
        env::var("GUTOE_RTSC_HANDBUILD_OUT")
            .unwrap_or_else(|_| "/tmp/bh_renders/rtsc_hand_build_playbook".to_string()),
    );
    fs::create_dir_all(&out_dir)
        .map_err(|e| format!("failed to create out dir {}: {e}", out_dir.display()))?;

    let plan_json = PathBuf::from(
        env::var("GUTOE_RTSC_ROUTES_JSON").unwrap_or_else(|_| {
            "/tmp/bh_renders/rtsc_synthesis_routes/rtsc_synthesis_routes.json".to_string()
        }),
    );
    let campaign_json = PathBuf::from(
        env::var("GUTOE_RTSC_CAMPAIGN_JSON").unwrap_or_else(|_| {
            "/tmp/bh_renders/rtsc_end_to_end_campaign/rtsc_end_to_end_campaign.json".to_string()
        }),
    );
    let regime = env::var("GUTOE_RTSC_REGIME").unwrap_or_else(|_| "heroic_process".to_string());

    let plan = load_plan(&plan_json)?;
    let campaign = load_campaign(&campaign_json)?;

    let mut lanes: Vec<BuildLane> = plan
        .iter()
        .map(|p| {
            let pv = pick_p(&campaign, &regime, &p.symbol);
            // High p_validated, low hazard, low difficulty.
            let score = pv * (1.0 - 0.45 * p.hazard) * (1.0 - 0.30 * p.engineering_difficulty);
            BuildLane {
                symbol: p.symbol.clone(),
                native_phase: p.native_phase.clone(),
                p_validated: pv,
                engineering_difficulty: p.engineering_difficulty,
                hazard: p.hazard,
                cap: p.cap.clone(),
                first_route: p.first_route.clone(),
                hand_build_score: score,
            }
        })
        .filter(|l| l.p_validated > 0.001)
        .collect();

    lanes.sort_by(|a, b| {
        b.hand_build_score
            .partial_cmp(&a.hand_build_score)
            .unwrap_or(Ordering::Equal)
    });

    let top_n = lanes.len().min(3);
    let selected = &lanes[..top_n];

    let mut txt = String::new();
    txt.push_str("[rtsc_hand_build_playbook]\n");
    txt.push_str("mode = practical_hand_build_lanes\n");
    txt.push_str(&format!("regime = {}\n", regime));
    txt.push_str(&format!("routes_source = {}\n", plan_json.display()));
    txt.push_str(&format!("campaign_source = {}\n\n", campaign_json.display()));

    txt.push_str("[minimum_equipment]\n");
    txt.push_str("- UHV sputter system (base <= 1e-8 mbar) with load-lock\n");
    txt.push_str("- Substrate heater (up to 400 C) + controlled cooldown\n");
    txt.push_str("- RGA residual gas monitor + quartz crystal thickness monitor\n");
    txt.push_str("- Inert transfer or immediate cap chamber\n");
    txt.push_str("- XRD (phase gate), 4-probe transport setup (280-380 K), Meissner check\n\n");

    txt.push_str("[global_process_guardrails]\n");
    txt.push_str("- no open-air transfer before capping\n");
    txt.push_str("- keep O2/H2O partial pressure minimal during growth\n");
    txt.push_str("- start with 12-sample DOE per lane (3 temps x 2 rates x 2 thicknesses)\n");
    txt.push_str("- fail-fast: if simple-cubic signature absent in first 12, retune substrate/buffer only\n\n");

    txt.push_str("[selected_lanes]\n");

    let mut json_rows = String::new();
    json_rows.push('[');

    for (idx, lane) in selected.iter().enumerate() {
        let (substrate, method, temp_c, rate, thickness) = material_recipe(&lane.symbol);
        let n95_1 = n_for_target_successes(lane.p_validated, 1, 0.95);
        let n95_5 = n_for_target_successes(lane.p_validated, 5, 0.95);
        let n95_10 = n_for_target_successes(lane.p_validated, 10, 0.95);

        txt.push_str(&format!("lane_{} = {}\n", idx + 1, lane.symbol));
        txt.push_str(&format!("  p_validated = {:.4}\n", lane.p_validated));
        txt.push_str(&format!("  hazard = {:.2}, engineering_difficulty = {:.2}\n", lane.hazard, lane.engineering_difficulty));
        txt.push_str(&format!("  substrate = {}\n", substrate));
        txt.push_str(&format!("  method = {}\n", method));
        txt.push_str(&format!("  deposition_temp = {}\n", temp_c));
        txt.push_str(&format!("  deposition_rate = {}\n", rate));
        txt.push_str(&format!("  film_thickness = {}\n", thickness));
        txt.push_str(&format!("  cap = {}\n", lane.cap));
        txt.push_str(&format!("  route_hint = {}\n", lane.first_route));
        txt.push_str("  QC gates:\n");
        txt.push_str("    1) XRD shows forced cubic lane (reject otherwise)\n");
        txt.push_str("    2) 4-probe resistivity drop/onset near >=300 K\n");
        txt.push_str("    3) repeatable near-zero-R under current ramp\n");
        txt.push_str("    4) Meissner response check (screening signal)\n");
        txt.push_str(&format!("  planning_counts_95pct: >=1 hit: {}, >=5 hits: {}, >=10 hits: {}\n\n", n95_1, n95_5, n95_10));

        json_rows.push_str(&format!(
            concat!(
                "{{\"lane_rank\":{},\"symbol\":\"{}\",\"native_phase\":\"{}\",\"p_validated\":{:.6},",
                "\"hazard\":{:.6},\"engineering_difficulty\":{:.6},\"substrate\":\"{}\",\"method\":\"{}\",",
                "\"deposition_temp_c\":\"{}\",\"deposition_rate_nm_s\":\"{}\",\"thickness_nm\":\"{}\",\"cap\":\"{}\",",
                "\"n95_ge1\":{},\"n95_ge5\":{},\"n95_ge10\":{}}}{}"
            ),
            idx + 1,
            lane.symbol,
            lane.native_phase,
            lane.p_validated,
            lane.hazard,
            lane.engineering_difficulty,
            substrate,
            method,
            temp_c,
            rate,
            thickness,
            lane.cap,
            n95_1,
            n95_5,
            n95_10,
            if idx + 1 == selected.len() { "" } else { "," }
        ));
    }
    json_rows.push(']');

    txt.push_str("[safety_boundary]\n");
    txt.push_str("- Cd lane is excluded from hand-build shortlist due toxicity burden\n");
    txt.push_str("- high-vacuum + high-temperature operations require trained lab workflow\n");
    txt.push_str("- treat this as R&D process planning, not unsupervised home fabrication\n");

    let json = format!(
        concat!(
            "{{\n",
            "  \"meta\": {{\"mode\": \"practical_hand_build_lanes\", \"regime\": \"{}\", \"routes_source\": \"{}\", \"campaign_source\": \"{}\"}},\n",
            "  \"lanes\": {}\n",
            "}}\n"
        ),
        regime,
        plan_json.display(),
        campaign_json.display(),
        json_rows
    );

    let txt_path = out_dir.join("rtsc_hand_build_playbook.txt");
    let json_path = out_dir.join("rtsc_hand_build_playbook.json");

    fs::write(&txt_path, txt)
        .map_err(|e| format!("failed to write {}: {e}", txt_path.display()))?;
    fs::write(&json_path, json)
        .map_err(|e| format!("failed to write {}: {e}", json_path.display()))?;

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());

    Ok(())
}
