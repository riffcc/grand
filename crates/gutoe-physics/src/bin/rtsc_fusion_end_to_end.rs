use serde_json::{json, Value};
use std::env;
use std::fs;
use std::path::{Path, PathBuf};

const C_LIGHT: f64 = 299_792_458.0;
const C2: f64 = C_LIGHT * C_LIGHT;
const SEC_PER_YEAR: f64 = 365.25 * 24.0 * 3600.0;

const DEFAULT_ANTIMATTER_BOUNDS_JSON: &str =
    "/tmp/bh_renders/antimatter_theoretical_bounds/antimatter_theoretical_yield_bounds.json";
const DEFAULT_ETA_MODELED: f64 = 3.065_432_385_279_946_7e-3;
const DEFAULT_ETA_ROUTE_PP: f64 = 1.667_271_686_855_778_8e-1;
const DEFAULT_ETA_PAIR_ABS: f64 = 5.0e-1;

#[derive(Clone, Debug)]
struct ReactorBest {
    material: String,
    shape: String,
    mesh_quality: f64,
    r_major_m: f64,
    a_minor_m: f64,
    t_kev: f64,
    q_engineering: f64,
}

#[derive(Clone, Debug)]
struct ModulePower {
    p_gross_w: f64,
    p_recirc_w: f64,
    p_net_w: f64,
    volume_m3: f64,
}

#[derive(Clone, Debug)]
struct CampaignRow {
    regime: String,
    symbol: String,
    p_validated: f64,
}

#[derive(Clone, Debug)]
struct TimelineRow {
    t_s: usize,
    phase: &'static str,
    temp_frac: f64,
    field_frac: f64,
    confinement_frac: f64,
    availability: f64,
    p_gross_w: f64,
    p_recirc_w: f64,
    p_net_w: f64,
    e_net_j: f64,
}

#[derive(Clone, Debug)]
struct AntimatterBounds {
    source: String,
    eta_modeled_current: f64,
    eta_route_pp_floor: f64,
    eta_pair_absolute: f64,
}

fn load_reactor_best(path: &Path) -> Result<ReactorBest, String> {
    let txt = fs::read_to_string(path)
        .map_err(|e| format!("failed to read reactor json {}: {e}", path.display()))?;
    let v: Value = serde_json::from_str(&txt)
        .map_err(|e| format!("failed to parse reactor json {}: {e}", path.display()))?;
    let b = v
        .get("best_overall")
        .and_then(Value::as_object)
        .ok_or_else(|| "missing best_overall in reactor json".to_string())?;

    Ok(ReactorBest {
        material: b
            .get("material")
            .and_then(Value::as_str)
            .unwrap_or("Cr")
            .to_string(),
        shape: b
            .get("shape")
            .and_then(Value::as_str)
            .unwrap_or("toroidal_honeycomb")
            .to_string(),
        mesh_quality: b.get("mesh_quality").and_then(Value::as_f64).unwrap_or(1.0),
        r_major_m: b.get("r_major_m").and_then(Value::as_f64).unwrap_or(1.5),
        a_minor_m: b.get("a_minor_m").and_then(Value::as_f64).unwrap_or(0.6),
        t_kev: b.get("t_kev").and_then(Value::as_f64).unwrap_or(150.0),
        q_engineering: b
            .get("q_engineering")
            .and_then(Value::as_f64)
            .unwrap_or(1.0),
    })
}

fn parse_f64_field(headers: &[&str], fields: &[&str], key: &str) -> f64 {
    if let Some(i) = headers.iter().position(|h| *h == key) {
        return fields
            .get(i)
            .and_then(|s| s.parse::<f64>().ok())
            .unwrap_or(0.0);
    }
    0.0
}

fn parse_str_field<'a>(headers: &[&str], fields: &'a [&'a str], key: &str) -> &'a str {
    if let Some(i) = headers.iter().position(|h| *h == key) {
        return fields.get(i).copied().unwrap_or("");
    }
    ""
}

fn load_module_power(best: &ReactorBest, all_points_csv: &Path) -> Result<ModulePower, String> {
    let txt = fs::read_to_string(all_points_csv).map_err(|e| {
        format!(
            "failed to read reactor points csv {}: {e}",
            all_points_csv.display()
        )
    })?;
    let mut lines = txt.lines();
    let header = lines
        .next()
        .ok_or_else(|| "empty reactor points csv".to_string())?;
    let headers: Vec<&str> = header.split(',').collect();

    let mut best_match: Option<(f64, ModulePower)> = None;
    for line in lines {
        let fields: Vec<&str> = line.split(',').collect();
        if fields.len() < headers.len() {
            continue;
        }

        let material = parse_str_field(&headers, &fields, "material");
        if material != best.material {
            continue;
        }

        let shape = parse_str_field(&headers, &fields, "shape");
        if shape != best.shape {
            continue;
        }

        let mesh_quality = parse_f64_field(&headers, &fields, "mesh_quality");
        let r_major = parse_f64_field(&headers, &fields, "r_major_m");
        let a_minor = parse_f64_field(&headers, &fields, "a_minor_m");
        let t_kev = parse_f64_field(&headers, &fields, "t_kev");

        let err = (mesh_quality - best.mesh_quality).abs()
            + (r_major - best.r_major_m).abs()
            + (a_minor - best.a_minor_m).abs()
            + 0.02 * (t_kev - best.t_kev).abs();

        let p_gross = parse_f64_field(&headers, &fields, "p_electric_gross_w");
        let p_recirc = parse_f64_field(&headers, &fields, "p_recirc_w");
        let p_net = parse_f64_field(&headers, &fields, "p_net_w");
        let volume = parse_f64_field(&headers, &fields, "volume_m3");

        let cand = ModulePower {
            p_gross_w: p_gross,
            p_recirc_w: p_recirc,
            p_net_w: p_net,
            volume_m3: volume,
        };

        match &best_match {
            Some((be, _)) if err >= *be => {}
            _ => best_match = Some((err, cand)),
        }
    }

    best_match
        .map(|(_, p)| p)
        .ok_or_else(|| "no matching module power row found in all_points csv".to_string())
}

fn load_campaign(path: &Path) -> Result<Vec<CampaignRow>, String> {
    let txt = fs::read_to_string(path)
        .map_err(|e| format!("failed to read campaign json {}: {e}", path.display()))?;
    let v: Value = serde_json::from_str(&txt)
        .map_err(|e| format!("failed to parse campaign json {}: {e}", path.display()))?;
    let arr = v
        .get("rows")
        .and_then(Value::as_array)
        .ok_or_else(|| "missing rows[] in campaign json".to_string())?;

    let mut out = Vec::new();
    for r in arr {
        out.push(CampaignRow {
            regime: r
                .get("regime")
                .and_then(Value::as_str)
                .unwrap_or("")
                .to_string(),
            symbol: r
                .get("symbol")
                .and_then(Value::as_str)
                .unwrap_or("")
                .to_string(),
            p_validated: r.get("p_validated").and_then(Value::as_f64).unwrap_or(0.0),
        });
    }
    Ok(out)
}

fn choose_p_validated(rows: &[CampaignRow], regime: &str, material: &str) -> f64 {
    rows.iter()
        .find(|r| r.regime == regime && r.symbol == material)
        .map(|r| r.p_validated)
        .unwrap_or(0.0)
}

fn json_f64_at(v: &Value, path: &[&str]) -> Option<f64> {
    let mut cur = v;
    for key in path {
        cur = cur.get(*key)?;
    }
    cur.as_f64()
}

fn json_string_at(v: &Value, path: &[&str]) -> Option<String> {
    let mut cur = v;
    for key in path {
        cur = cur.get(*key)?;
    }
    cur.as_str().map(|s| s.to_string())
}

fn load_antimatter_bounds(path: &Path) -> Option<AntimatterBounds> {
    let txt = fs::read_to_string(path).ok()?;
    let v: Value = serde_json::from_str(&txt).ok()?;

    Some(AntimatterBounds {
        source: json_string_at(&v, &["fusion_reference", "source"])
            .unwrap_or_else(|| "antimatter_bounds_file".to_string()),
        eta_modeled_current: json_f64_at(&v, &["efficiencies", "modeled_current"])?,
        eta_route_pp_floor: json_f64_at(&v, &["efficiencies", "route_pp_floor"])?,
        eta_pair_absolute: json_f64_at(&v, &["efficiencies", "pair_absolute"])?,
    })
}

fn default_antimatter_bounds() -> AntimatterBounds {
    AntimatterBounds {
        source: "embedded_defaults".to_string(),
        eta_modeled_current: DEFAULT_ETA_MODELED,
        eta_route_pp_floor: DEFAULT_ETA_ROUTE_PP,
        eta_pair_absolute: DEFAULT_ETA_PAIR_ABS,
    }
}

fn antimatter_stream_ng_per_year_for_annihilation_equivalent(power_w: f64) -> f64 {
    if power_w <= 0.0 {
        return 0.0;
    }
    (power_w / (2.0 * C2)) * SEC_PER_YEAR * 1.0e12
}

fn beam_mw_for_annihilation_parity(fusion_power_mw: f64, eta: f64) -> f64 {
    if eta <= 0.0 {
        return f64::INFINITY;
    }
    fusion_power_mw / (2.0 * eta)
}

fn annihilation_density_ratio_to_fusion(beam_power_mw: f64, eta: f64, fusion_power_mw: f64) -> f64 {
    if fusion_power_mw <= 0.0 {
        return f64::INFINITY;
    }
    2.0 * eta * beam_power_mw / fusion_power_mw
}

fn antimatter_annihilation_output_mw(beam_power_mw: f64, eta: f64) -> f64 {
    2.0 * beam_power_mw.max(0.0) * eta.max(0.0)
}

fn main() -> Result<(), String> {
    let out_dir = PathBuf::from(
        env::var("GUTOE_FUSION_END2END_OUT")
            .unwrap_or_else(|_| "/tmp/bh_renders/rtsc_fusion_end_to_end".to_string()),
    );
    fs::create_dir_all(&out_dir)
        .map_err(|e| format!("failed to create output dir {}: {e}", out_dir.display()))?;

    let reactor_json = PathBuf::from(env::var("GUTOE_FUSION_REACTOR_JSON").unwrap_or_else(|_| {
        "/tmp/bh_renders/rtsc_honeycomb_fusion_reactor/rtsc_honeycomb_fusion_reactor.json"
            .to_string()
    }));
    let reactor_points = PathBuf::from(
        env::var("GUTOE_FUSION_REACTOR_POINTS").unwrap_or_else(|_| {
            "/tmp/bh_renders/rtsc_honeycomb_fusion_reactor/rtsc_honeycomb_fusion_reactor_all_points.csv"
                .to_string()
        }),
    );
    let campaign_json = PathBuf::from(env::var("GUTOE_RTSC_CAMPAIGN_JSON").unwrap_or_else(|_| {
        "/tmp/bh_renders/rtsc_end_to_end_campaign/rtsc_end_to_end_campaign.json".to_string()
    }));
    let antimatter_bounds_json = PathBuf::from(
        env::var("GUTOE_ANTIMATTER_BOUNDS_JSON")
            .unwrap_or_else(|_| DEFAULT_ANTIMATTER_BOUNDS_JSON.to_string()),
    );

    let regime = env::var("GUTOE_RTSC_REGIME").unwrap_or_else(|_| "heroic_process".to_string());
    let target_grid_mw = env::var("GUTOE_GRID_TARGET_MW")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(1000.0)
        .max(50.0);
    let fab_attempts = env::var("GUTOE_FAB_ATTEMPTS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(120)
        .max(10);
    let batch_days = env::var("GUTOE_BATCH_DAYS")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(14.0)
        .max(1.0);
    let anti_beam_budget_mw = env::var("GUTOE_ANTI_BEAM_BUDGET_MW")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(100.0)
        .max(0.0);

    let best = load_reactor_best(&reactor_json)?;
    let module = load_module_power(&best, &reactor_points)?;
    let campaign_rows = load_campaign(&campaign_json)?;
    let p_validated = choose_p_validated(&campaign_rows, &regime, &best.material);

    if p_validated <= 0.0 {
        return Err(format!(
            "no validated synthesis probability for material={} regime={} in {}",
            best.material,
            regime,
            campaign_json.display()
        ));
    }

    let target_grid_w = target_grid_mw * 1.0e6;
    let modules_required = (target_grid_w / module.p_net_w).ceil().max(1.0) as usize;
    let expected_modules_built = (fab_attempts as f64 * p_validated).floor() as usize;
    let modules_online = expected_modules_built.min(modules_required);

    let fabrication_days = if p_validated > 0.0 {
        (modules_required as f64 / (fab_attempts as f64 * p_validated)).ceil() * batch_days
    } else {
        f64::INFINITY
    };

    // Operational timeline (24h): startup -> ramp -> burn -> downramp.
    let dt_s = 1usize;
    let t_startup = 30 * 60;
    let t_ramp = 45 * 60;
    let t_burn = 18 * 3600;
    let t_down = 30 * 60;
    let total_t = t_startup + t_ramp + t_burn + t_down;

    let module_availability_nominal =
        (0.96 + 0.03 * (best.q_engineering / 2.0).clamp(0.0, 1.0)).clamp(0.90, 0.995);
    let startup_overhead_factor = 1.45;

    let mut timeline = Vec::<TimelineRow>::new();
    let mut e_net_j = 0.0f64;

    for t in (0..=total_t).step_by(dt_s) {
        let (phase, temp_frac, field_frac, conf_frac, phase_recirc_boost) = if t <= t_startup {
            let u = t as f64 / t_startup as f64;
            (
                "startup",
                0.20 * u,
                0.50 * u,
                0.10 * u,
                startup_overhead_factor,
            )
        } else if t <= t_startup + t_ramp {
            let u = (t - t_startup) as f64 / t_ramp as f64;
            (
                "ramp",
                0.20 + 0.80 * u,
                0.50 + 0.50 * u,
                0.10 + 0.90 * u,
                1.15,
            )
        } else if t <= t_startup + t_ramp + t_burn {
            let burn_t = t - t_startup - t_ramp;
            let burn_u = burn_t as f64 / t_burn as f64;
            let slow_decay = 1.0 - 0.05 * burn_u;
            ("burn", 1.0, 1.0, slow_decay, 1.00)
        } else {
            let u = (t - t_startup - t_ramp - t_burn) as f64 / t_down as f64;
            (
                "downramp",
                (1.0 - u).max(0.0),
                (1.0 - 0.8 * u).max(0.0),
                (1.0 - u).max(0.0),
                1.10,
            )
        };

        let availability = module_availability_nominal * (0.995 - 0.010 * (1.0 - conf_frac));
        let on_modules = modules_online as f64 * availability;

        let gross = on_modules * module.p_gross_w * temp_frac * field_frac * conf_frac;
        let recirc = on_modules
            * module.p_recirc_w
            * (0.65 + 0.35 * field_frac)
            * (0.70 + 0.30 * temp_frac)
            * phase_recirc_boost;
        let net = gross - recirc;

        e_net_j += net * dt_s as f64;

        timeline.push(TimelineRow {
            t_s: t,
            phase,
            temp_frac,
            field_frac,
            confinement_frac: conf_frac,
            availability,
            p_gross_w: gross,
            p_recirc_w: recirc,
            p_net_w: net,
            e_net_j,
        });
    }

    let p_net_peak_w = timeline
        .iter()
        .map(|r| r.p_net_w)
        .fold(f64::NEG_INFINITY, f64::max);
    let p_net_avg_w = timeline.iter().map(|r| r.p_net_w).sum::<f64>() / timeline.len() as f64;
    let p_recirc_avg_w = timeline.iter().map(|r| r.p_recirc_w).sum::<f64>() / timeline.len() as f64;
    let p_recirc_peak_w = timeline
        .iter()
        .map(|r| r.p_recirc_w)
        .fold(f64::NEG_INFINITY, f64::max);
    let e_net_mwh = e_net_j / 3.6e9;
    let module_net_mw = module.p_net_w / 1.0e6;
    let module_density_w_m3 = if module.volume_m3 > 0.0 {
        module.p_net_w / module.volume_m3
    } else {
        f64::INFINITY
    };

    let antimatter_bounds =
        load_antimatter_bounds(&antimatter_bounds_json).unwrap_or_else(default_antimatter_bounds);
    let antimatter_scenarios = vec![
        ("modeled_current", antimatter_bounds.eta_modeled_current),
        ("route_pp_floor", antimatter_bounds.eta_route_pp_floor),
        ("pair_absolute", antimatter_bounds.eta_pair_absolute),
    ];

    let module_ann_stream_ng_per_year =
        antimatter_stream_ng_per_year_for_annihilation_equivalent(module.p_net_w);
    let target_grid_ann_stream_ng_per_year =
        antimatter_stream_ng_per_year_for_annihilation_equivalent(target_grid_w);

    let mut antimatter_rows_txt = String::new();
    let mut antimatter_rows_json = Vec::<Value>::new();
    let mut hybrid_rows_txt = String::new();
    let mut hybrid_rows_json = Vec::<Value>::new();
    for (name, eta) in antimatter_scenarios {
        let ann_ratio_module_beam =
            annihilation_density_ratio_to_fusion(module_net_mw, eta, module_net_mw);
        let ann_ratio_grid_beam =
            annihilation_density_ratio_to_fusion(target_grid_mw, eta, module_net_mw);
        let beam_parity_mw = beam_mw_for_annihilation_parity(module_net_mw, eta);
        let beam_parity_target_ratio = beam_parity_mw / target_grid_mw.max(1.0e-30);

        antimatter_rows_txt.push_str(&format!(
            "{},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}\n",
            name,
            eta,
            ann_ratio_module_beam,
            ann_ratio_grid_beam,
            beam_parity_mw,
            beam_parity_target_ratio
        ));

        antimatter_rows_json.push(json!({
            "scenario": name,
            "eta_net": eta,
            "annihilation_density_ratio_at_module_net_beam": ann_ratio_module_beam,
            "annihilation_density_ratio_at_target_grid_beam": ann_ratio_grid_beam,
            "beam_mw_for_annihilation_density_parity_with_module": beam_parity_mw,
            "beam_mw_parity_over_target_grid": beam_parity_target_ratio
        }));

        let ann_output_budget_mw = antimatter_annihilation_output_mw(anti_beam_budget_mw, eta);
        let ann_output_budget_w = ann_output_budget_mw * 1.0e6;
        let recirc_offset_frac_avg = if p_recirc_avg_w > 0.0 {
            ann_output_budget_w / p_recirc_avg_w
        } else {
            f64::INFINITY
        };
        let recirc_offset_frac_peak = if p_recirc_peak_w > 0.0 {
            ann_output_budget_w / p_recirc_peak_w
        } else {
            f64::INFINITY
        };
        let beam_for_10pct_recirc =
            beam_mw_for_annihilation_parity((0.10 * p_recirc_avg_w) / 1.0e6, eta);
        let beam_for_50pct_recirc =
            beam_mw_for_annihilation_parity((0.50 * p_recirc_avg_w) / 1.0e6, eta);
        let beam_for_100pct_recirc = beam_mw_for_annihilation_parity(p_recirc_avg_w / 1.0e6, eta);
        let stream_ng_per_year_budget =
            antimatter_stream_ng_per_year_for_annihilation_equivalent(ann_output_budget_w);

        hybrid_rows_txt.push_str(&format!(
            "{},{:.12e},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6e}\n",
            name,
            eta,
            ann_output_budget_mw,
            recirc_offset_frac_avg,
            recirc_offset_frac_peak,
            beam_for_10pct_recirc,
            beam_for_50pct_recirc,
            beam_for_100pct_recirc,
            stream_ng_per_year_budget
        ));

        hybrid_rows_json.push(json!({
            "scenario": name,
            "eta_net": eta,
            "beam_budget_mw": anti_beam_budget_mw,
            "annihilation_output_at_budget_mw": ann_output_budget_mw,
            "recirc_offset_fraction_avg": recirc_offset_frac_avg,
            "recirc_offset_fraction_peak": recirc_offset_frac_peak,
            "beam_mw_for_10pct_avg_recirc_offset": beam_for_10pct_recirc,
            "beam_mw_for_50pct_avg_recirc_offset": beam_for_50pct_recirc,
            "beam_mw_for_100pct_avg_recirc_offset": beam_for_100pct_recirc,
            "antimatter_stream_ng_per_year_at_budget": stream_ng_per_year_budget
        }));
    }

    let mut csv = String::from(
        "t_s,phase,temp_frac,field_frac,confinement_frac,availability,p_gross_w,p_recirc_w,p_net_w,e_net_j\n",
    );
    for r in &timeline {
        csv.push_str(&format!(
            "{},{},{:.6},{:.6},{:.6},{:.6},{:.6e},{:.6e},{:.6e},{:.6e}\n",
            r.t_s,
            r.phase,
            r.temp_frac,
            r.field_frac,
            r.confinement_frac,
            r.availability,
            r.p_gross_w,
            r.p_recirc_w,
            r.p_net_w,
            r.e_net_j
        ));
    }

    let mut txt = String::new();
    txt.push_str("[rtsc_fusion_end_to_end]\n");
    txt.push_str(&format!("regime = {}\n", regime));
    txt.push_str(&format!("material = {}\n", best.material));
    txt.push_str(&format!("shape = {}\n", best.shape));
    txt.push_str(&format!("target_grid_mw = {:.3}\n", target_grid_mw));
    txt.push_str(&format!("module_net_mw = {:.3}\n", module_net_mw));
    txt.push_str(&format!(
        "module_gross_mw = {:.3}\n",
        module.p_gross_w / 1.0e6
    ));
    txt.push_str(&format!(
        "module_recirc_mw = {:.3}\n",
        module.p_recirc_w / 1.0e6
    ));
    txt.push_str(&format!(
        "q_engineering_module = {:.3}\n",
        best.q_engineering
    ));
    txt.push_str(&format!("p_validated = {:.6}\n", p_validated));
    txt.push_str(&format!("fab_attempts_per_batch = {}\n", fab_attempts));
    txt.push_str(&format!("batch_days = {:.3}\n", batch_days));
    txt.push_str(&format!("modules_required = {}\n", modules_required));
    txt.push_str(&format!(
        "expected_modules_built = {}\n",
        expected_modules_built
    ));
    txt.push_str(&format!("modules_online = {}\n", modules_online));
    txt.push_str(&format!(
        "expected_fabrication_days_for_target = {:.3}\n",
        fabrication_days
    ));
    txt.push('\n');

    txt.push_str("[operation_24h]\n");
    txt.push_str(&format!("timeline_seconds = {}\n", total_t));
    txt.push_str(&format!(
        "net_power_peak_mw = {:.3}\n",
        p_net_peak_w / 1.0e6
    ));
    txt.push_str(&format!("net_power_avg_mw = {:.3}\n", p_net_avg_w / 1.0e6));
    txt.push_str(&format!("net_energy_mwh = {:.3}\n", e_net_mwh));
    txt.push('\n');

    txt.push_str("[antimatter_density_comparison]\n");
    txt.push_str(&format!("bounds_source = {}\n", antimatter_bounds.source));
    txt.push_str(&format!(
        "bounds_json = {}\n",
        antimatter_bounds_json.display()
    ));
    txt.push_str(&format!(
        "fusion_module_volume_m3 = {:.12e}\n",
        module.volume_m3
    ));
    txt.push_str(&format!(
        "fusion_module_net_density_w_m3 = {:.12e}\n",
        module_density_w_m3
    ));
    txt.push_str(&format!(
        "annihilation_equivalent_stream_ng_per_year_for_module_net = {:.12e}\n",
        module_ann_stream_ng_per_year
    ));
    txt.push_str(&format!(
        "annihilation_equivalent_stream_ng_per_year_for_target_grid = {:.12e}\n",
        target_grid_ann_stream_ng_per_year
    ));
    txt.push_str("scenario,eta_net,ann_ratio_at_module_net_beam,ann_ratio_at_target_grid_beam,beam_mw_for_ann_density_parity_with_module,beam_mw_parity_over_target_grid\n");
    txt.push_str(&antimatter_rows_txt);
    txt.push('\n');

    txt.push_str("[antimatter_fusion_hybrid_assist]\n");
    txt.push_str(
        "assumption = antimatter_used_as_auxiliary_recirculation_offset_not_primary_generation\n",
    );
    txt.push_str(&format!(
        "antimatter_beam_budget_mw = {:.6}\n",
        anti_beam_budget_mw
    ));
    txt.push_str(&format!(
        "fusion_recirc_avg_mw = {:.6}\n",
        p_recirc_avg_w / 1.0e6
    ));
    txt.push_str(&format!(
        "fusion_recirc_peak_mw = {:.6}\n",
        p_recirc_peak_w / 1.0e6
    ));
    txt.push_str("scenario,eta_net,annihilation_output_at_budget_mw,recirc_offset_fraction_avg,recirc_offset_fraction_peak,beam_mw_for_10pct_avg_recirc_offset,beam_mw_for_50pct_avg_recirc_offset,beam_mw_for_100pct_avg_recirc_offset,antimatter_stream_ng_per_year_at_budget\n");
    txt.push_str(&hybrid_rows_txt);
    txt.push('\n');

    txt.push_str("artifacts:\n");
    txt.push_str(&format!(
        "timeline_csv = {}/rtsc_fusion_end_to_end_timeline.csv\n",
        out_dir.display()
    ));
    txt.push_str(&format!(
        "summary_json = {}/rtsc_fusion_end_to_end_summary.json\n",
        out_dir.display()
    ));

    let json = json!({
        "meta": {
            "regime": regime,
            "material": best.material,
            "shape": best.shape
        },
        "module": {
            "p_net_w": module.p_net_w,
            "p_gross_w": module.p_gross_w,
            "p_recirc_w": module.p_recirc_w,
            "q_engineering": best.q_engineering,
            "volume_m3": module.volume_m3,
            "net_density_w_m3": module_density_w_m3
        },
        "synthesis": {
            "p_validated": p_validated,
            "fab_attempts_per_batch": fab_attempts,
            "batch_days": batch_days,
            "modules_required": modules_required,
            "expected_modules_built": expected_modules_built,
            "modules_online": modules_online,
            "expected_fabrication_days_for_target": fabrication_days
        },
        "operation_24h": {
            "timeline_seconds": total_t,
            "net_power_peak_w": p_net_peak_w,
            "net_power_avg_w": p_net_avg_w,
            "recirc_power_avg_w": p_recirc_avg_w,
            "recirc_power_peak_w": p_recirc_peak_w,
            "net_energy_mwh": e_net_mwh
        },
        "antimatter_comparison": {
            "bounds_source": antimatter_bounds.source,
            "bounds_json": antimatter_bounds_json.display().to_string(),
            "annihilation_equivalent_stream_ng_per_year_for_module_net": module_ann_stream_ng_per_year,
            "annihilation_equivalent_stream_ng_per_year_for_target_grid": target_grid_ann_stream_ng_per_year,
            "scenarios": antimatter_rows_json
        },
        "antimatter_hybrid_assist": {
            "assumption": "antimatter_used_as_auxiliary_recirculation_offset_not_primary_generation",
            "beam_budget_mw": anti_beam_budget_mw,
            "scenarios": hybrid_rows_json
        }
    })
    .to_string();

    let txt_path = out_dir.join("rtsc_fusion_end_to_end_summary.txt");
    let json_path = out_dir.join("rtsc_fusion_end_to_end_summary.json");
    let csv_path = out_dir.join("rtsc_fusion_end_to_end_timeline.csv");

    fs::write(&txt_path, txt)
        .map_err(|e| format!("failed to write {}: {e}", txt_path.display()))?;
    fs::write(&json_path, json)
        .map_err(|e| format!("failed to write {}: {e}", json_path.display()))?;
    fs::write(&csv_path, csv)
        .map_err(|e| format!("failed to write {}: {e}", csv_path.display()))?;

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!("wrote {}", csv_path.display());

    Ok(())
}
