use serde_json::{json, Value};
use std::env;
use std::f64::consts::PI;
use std::fs;
use std::path::Path;

const C_LIGHT: f64 = 299_792_458.0;
const C2: f64 = C_LIGHT * C_LIGHT;
const SEC_PER_YEAR: f64 = 365.25 * 24.0 * 3600.0;
const J_PER_KILOTON_TNT: f64 = 4.184e12;
const J_PER_TWH: f64 = 3.6e15;

const M_PROTON_GEV: f64 = 0.938_272_088_16;
const M_ELECTRON_GEV: f64 = 0.000_510_998_95;

// Fallback reference from the current RTSC fusion lane best-overall point.
const DEFAULT_FUSION_REF_NET_MW: f64 = 29.372_22;
const DEFAULT_FUSION_REF_VOLUME_M3: f64 = 12.448_7;
const DEFAULT_FUSION_REACTOR_JSON: &str =
    "/tmp/bh_renders/rtsc_honeycomb_fusion_reactor/rtsc_honeycomb_fusion_reactor.json";
const DEFAULT_FUSION_REACTOR_POINTS: &str =
    "/tmp/bh_renders/rtsc_honeycomb_fusion_reactor/rtsc_honeycomb_fusion_reactor_all_points.csv";

#[derive(Clone, Debug)]
struct FusionBestPoint {
    material: String,
    shape: String,
    mesh_quality: f64,
    r_major_m: f64,
    a_minor_m: f64,
    t_kev: f64,
    p_net_w: f64,
}

#[derive(Clone, Debug)]
struct FusionReference {
    source: String,
    net_power_w: f64,
    volume_m3: f64,
}

fn env_f64(name: &str, default: f64) -> f64 {
    env::var(name)
        .ok()
        .and_then(|v| v.parse::<f64>().ok())
        .unwrap_or(default)
}

fn env_f64_opt(name: &str) -> Option<f64> {
    env::var(name).ok().and_then(|v| v.parse::<f64>().ok())
}

fn antiproton_threshold_lab_kinetic_gev() -> f64 {
    // p + p(rest) -> p + p + p + pbar
    // Lab threshold kinetic energy for incident proton: 6 m_p.
    6.0 * M_PROTON_GEV
}

fn positron_pair_floor_input_gev_per_pos() -> f64 {
    // Minimum input to create one e+ via pair production is 2 m_e.
    2.0 * M_ELECTRON_GEV
}

fn antihydrogen_rest_gev() -> f64 {
    M_PROTON_GEV + M_ELECTRON_GEV
}

fn antihydrogen_route_floor_input_gev() -> f64 {
    antiproton_threshold_lab_kinetic_gev() + positron_pair_floor_input_gev_per_pos()
}

fn eta_route_pp_floor() -> f64 {
    antihydrogen_rest_gev() / antihydrogen_route_floor_input_gev()
}

fn g_per_mw_year_at_efficiency(eta: f64) -> f64 {
    let kg_per_mw_year_100 = 1.0e6 * SEC_PER_YEAR / C2;
    kg_per_mw_year_100 * 1000.0 * eta
}

fn g_per_year(power_mw: f64, eta: f64) -> f64 {
    power_mw.max(0.0) * g_per_mw_year_at_efficiency(eta.max(0.0))
}

fn power_mw_for_target_g_per_year(target_g_per_year: f64, eta: f64) -> f64 {
    if eta <= 0.0 {
        return f64::INFINITY;
    }
    target_g_per_year.max(0.0) / g_per_mw_year_at_efficiency(eta)
}

fn om_gap(from: f64, to: f64) -> f64 {
    if from <= 0.0 || to <= 0.0 {
        return f64::INFINITY;
    }
    (to / from).log10()
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

fn load_fusion_best(json_path: &Path) -> Option<FusionBestPoint> {
    let txt = fs::read_to_string(json_path).ok()?;
    let v: Value = serde_json::from_str(&txt).ok()?;
    let b = v.get("best_overall")?.as_object()?;

    Some(FusionBestPoint {
        material: b
            .get("material")
            .and_then(Value::as_str)
            .unwrap_or("")
            .to_string(),
        shape: b
            .get("shape")
            .and_then(Value::as_str)
            .unwrap_or("")
            .to_string(),
        mesh_quality: b.get("mesh_quality").and_then(Value::as_f64).unwrap_or(0.0),
        r_major_m: b.get("r_major_m").and_then(Value::as_f64).unwrap_or(0.0),
        a_minor_m: b.get("a_minor_m").and_then(Value::as_f64).unwrap_or(0.0),
        t_kev: b.get("t_kev").and_then(Value::as_f64).unwrap_or(0.0),
        p_net_w: b.get("p_net_w").and_then(Value::as_f64).unwrap_or(0.0),
    })
}

fn load_fusion_volume_from_points_csv(best: &FusionBestPoint, points_csv: &Path) -> Option<f64> {
    let txt = fs::read_to_string(points_csv).ok()?;
    let mut lines = txt.lines();
    let header = lines.next()?;
    let headers: Vec<&str> = header.split(',').collect();

    let mut best_match: Option<(f64, f64)> = None;
    for line in lines {
        let fields: Vec<&str> = line.split(',').collect();
        if fields.len() < headers.len() {
            continue;
        }

        if parse_str_field(&headers, &fields, "material") != best.material {
            continue;
        }
        if parse_str_field(&headers, &fields, "shape") != best.shape {
            continue;
        }

        let mesh_quality = parse_f64_field(&headers, &fields, "mesh_quality");
        let r_major = parse_f64_field(&headers, &fields, "r_major_m");
        let a_minor = parse_f64_field(&headers, &fields, "a_minor_m");
        let t_kev = parse_f64_field(&headers, &fields, "t_kev");
        let volume_m3 = parse_f64_field(&headers, &fields, "volume_m3");

        if volume_m3 <= 0.0 || !volume_m3.is_finite() {
            continue;
        }

        let err = (mesh_quality - best.mesh_quality).abs()
            + (r_major - best.r_major_m).abs()
            + (a_minor - best.a_minor_m).abs()
            + 0.02 * (t_kev - best.t_kev).abs();

        match best_match {
            Some((best_err, _)) if err >= best_err => {}
            _ => best_match = Some((err, volume_m3)),
        }
    }

    best_match.map(|(_, volume)| volume)
}

fn inferred_fusion_volume_m3(best: &FusionBestPoint) -> f64 {
    let torus_like = 2.0 * PI * PI * best.r_major_m.max(0.0) * best.a_minor_m.max(0.0).powi(2);
    match best.shape.as_str() {
        "spherical_honeycomb" => {
            let radius = 1.35 * best.a_minor_m.max(0.0);
            (4.0 / 3.0) * PI * radius.powi(3)
        }
        "toroidal_honeycomb" | "cylindrical_honeycomb" => torus_like,
        _ => torus_like,
    }
}

fn load_fusion_reference() -> FusionReference {
    let net_mw_override = env_f64_opt("GUTOE_ANTI_FUSION_REF_NET_MW")
        .filter(|v| *v > 0.0)
        .map(|v| v * 1.0e6);
    let volume_override = env_f64_opt("GUTOE_ANTI_FUSION_REF_VOLUME_M3").filter(|v| *v > 0.0);

    if net_mw_override.is_some() && volume_override.is_some() {
        return FusionReference {
            source: "env_override".to_string(),
            net_power_w: net_mw_override.unwrap_or(DEFAULT_FUSION_REF_NET_MW * 1.0e6),
            volume_m3: volume_override.unwrap_or(DEFAULT_FUSION_REF_VOLUME_M3),
        };
    }

    let reactor_json = env::var("GUTOE_ANTI_FUSION_REACTOR_JSON")
        .unwrap_or_else(|_| DEFAULT_FUSION_REACTOR_JSON.to_string());
    let reactor_points = env::var("GUTOE_ANTI_FUSION_REACTOR_POINTS")
        .unwrap_or_else(|_| DEFAULT_FUSION_REACTOR_POINTS.to_string());

    if let Some(best) = load_fusion_best(Path::new(&reactor_json)) {
        let net_power_w = net_mw_override.unwrap_or(best.p_net_w.max(0.0));
        let volume_m3 = volume_override.unwrap_or_else(|| {
            load_fusion_volume_from_points_csv(&best, Path::new(&reactor_points))
                .unwrap_or_else(|| inferred_fusion_volume_m3(&best))
        });

        let source = if net_mw_override.is_some() || volume_override.is_some() {
            "mixed_env_plus_rtsc_artifacts"
        } else {
            "rtsc_artifacts"
        };

        return FusionReference {
            source: source.to_string(),
            net_power_w,
            volume_m3,
        };
    }

    let source = if net_mw_override.is_some() || volume_override.is_some() {
        "env_partial_fallback"
    } else {
        "embedded_default"
    };

    FusionReference {
        source: source.to_string(),
        net_power_w: net_mw_override.unwrap_or(DEFAULT_FUSION_REF_NET_MW * 1.0e6),
        volume_m3: volume_override.unwrap_or(DEFAULT_FUSION_REF_VOLUME_M3),
    }
}

fn rest_output_power_mw(beam_power_mw: f64, eta: f64) -> f64 {
    beam_power_mw.max(0.0) * eta.max(0.0)
}

fn annihilation_output_power_mw(beam_power_mw: f64, eta: f64) -> f64 {
    2.0 * rest_output_power_mw(beam_power_mw, eta)
}

fn power_density_w_m3(power_mw: f64, volume_m3: f64) -> f64 {
    if volume_m3 <= 0.0 {
        return f64::INFINITY;
    }
    power_mw.max(0.0) * 1.0e6 / volume_m3
}

fn ratio(num: f64, den: f64) -> f64 {
    if den <= 0.0 {
        return f64::INFINITY;
    }
    num / den
}

fn beam_mw_for_rest_output_mw(target_rest_output_mw: f64, eta: f64) -> f64 {
    if eta <= 0.0 {
        return f64::INFINITY;
    }
    target_rest_output_mw.max(0.0) / eta
}

fn beam_mw_for_annihilation_output_mw(target_ann_output_mw: f64, eta: f64) -> f64 {
    if eta <= 0.0 {
        return f64::INFINITY;
    }
    target_ann_output_mw.max(0.0) / (2.0 * eta)
}

fn ng_per_year_for_rest_output_power_w(power_w: f64) -> f64 {
    if power_w <= 0.0 {
        return 0.0;
    }
    (power_w / C2) * SEC_PER_YEAR * 1.0e12
}

fn ng_per_year_for_annihilation_output_power_w(power_w: f64) -> f64 {
    if power_w <= 0.0 {
        return 0.0;
    }
    (power_w / (2.0 * C2)) * SEC_PER_YEAR * 1.0e12
}

fn roundtrip_output_over_input(eta: f64) -> f64 {
    2.0 * eta.max(0.0)
}

fn roundtrip_penalty_input_over_output(eta: f64) -> f64 {
    let gain = roundtrip_output_over_input(eta);
    if gain <= 0.0 {
        return f64::INFINITY;
    }
    1.0 / gain
}

fn thermodynamic_verdict(eta: f64) -> &'static str {
    let gain = roundtrip_output_over_input(eta);
    if (gain - 1.0).abs() < 1.0e-12 {
        "break_even_ceiling"
    } else if gain < 1.0 {
        "net_energy_sink"
    } else {
        "net_positive_assumption"
    }
}

fn main() {
    let out_dir = env::var("GUTOE_ANTIMATTER_THEORY_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/antimatter_theoretical_bounds".to_string());
    fs::create_dir_all(&out_dir).expect("create output dir");

    // Absolute symmetry ceiling (energy-only): at most half pair rest-energy ends up as antimatter.
    let eta_pair_absolute = env_f64("GUTOE_ANTI_ETA_PAIR_ABS", 0.5).clamp(0.0, 1.0);

    // Route-constrained kinematic floor from p+p antiproton production + positron pair floor.
    let eta_pp_route = env_f64("GUTOE_ANTI_ETA_ROUTE_PP", eta_route_pp_floor()).clamp(0.0, 1.0);

    // Current best modeled net efficiency from existing lane (eta_floor * eta_chain_top).
    let eta_modeled = env_f64("GUTOE_ANTI_ETA_MODELED", 3.065_432_385_279_946_7e-3).clamp(0.0, 1.0);

    let fusion_ref = load_fusion_reference();
    let fusion_net_mw = fusion_ref.net_power_w / 1.0e6;
    let fusion_density_w_m3 = fusion_ref.net_power_w / fusion_ref.volume_m3.max(1.0e-30);
    let fusion_year_energy_j = fusion_ref.net_power_w * SEC_PER_YEAR;
    let fusion_year_energy_twh = fusion_year_energy_j / J_PER_TWH;
    let fusion_year_energy_kt_tnt = fusion_year_energy_j / J_PER_KILOTON_TNT;

    let power_grid_mw = [0.1, 1.0, 5.0, 20.0, 100.0, 500.0, 1_000.0, 5_000.0];
    let targets_g_per_year = [1.0e-6, 1.0e-3, 1.0e-2, 1.0e-1, 1.0, 5.0, 10.0];

    let scenarios = [
        ("modeled_current", eta_modeled),
        ("route_pp_floor", eta_pp_route),
        ("pair_absolute", eta_pair_absolute),
    ];

    let mut txt = String::new();
    txt.push_str("[antimatter_theoretical_yield_bounds]\n");
    txt.push_str("scope = theoretical_output_bounds_only (no operational methodology)\n\n");

    txt.push_str("[physics_inputs]\n");
    txt.push_str(&format!(
        "antiproton_threshold_kinetic_gev = {:.12e}\n",
        antiproton_threshold_lab_kinetic_gev()
    ));
    txt.push_str(&format!(
        "positron_pair_floor_input_gev = {:.12e}\n",
        positron_pair_floor_input_gev_per_pos()
    ));
    txt.push_str(&format!(
        "antihydrogen_rest_gev = {:.12e}\n",
        antihydrogen_rest_gev()
    ));
    txt.push_str(&format!(
        "antihydrogen_route_floor_input_gev = {:.12e}\n",
        antihydrogen_route_floor_input_gev()
    ));
    txt.push_str(&format!("eta_route_pp_floor = {:.12e}\n", eta_pp_route));
    txt.push_str(&format!("eta_pair_absolute = {:.12e}\n", eta_pair_absolute));
    txt.push_str(&format!("eta_modeled_current = {:.12e}\n\n", eta_modeled));

    txt.push_str("[fusion_reference_for_density_comparison]\n");
    txt.push_str(&format!("source = {}\n", fusion_ref.source));
    txt.push_str(&format!("fusion_net_power_mw = {:.12e}\n", fusion_net_mw));
    txt.push_str(&format!(
        "fusion_reference_volume_m3 = {:.12e}\n",
        fusion_ref.volume_m3
    ));
    txt.push_str(&format!(
        "fusion_net_power_density_w_m3 = {:.12e}\n\n",
        fusion_density_w_m3
    ));

    txt.push_str("[efficiency_gaps]\n");
    txt.push_str(&format!(
        "om_gap_modeled_to_route_pp = {:.6}\n",
        om_gap(eta_modeled, eta_pp_route)
    ));
    txt.push_str(&format!(
        "om_gap_modeled_to_pair_absolute = {:.6}\n\n",
        om_gap(eta_modeled, eta_pair_absolute)
    ));

    txt.push_str("[yield_formula]\n");
    txt.push_str("g_per_year = 0.351125654089 * eta_net * P_MW\n\n");

    txt.push_str("[comparison_formulas]\n");
    txt.push_str("rest_output_mw = eta_net * P_beam_mw\n");
    txt.push_str("annihilation_output_mw = 2 * eta_net * P_beam_mw\n");
    txt.push_str("annihilation_density_w_m3 = annihilation_output_mw * 1e6 / V_ref_m3\n");
    txt.push_str(
        "ratio_to_fusion_density = annihilation_density_w_m3 / fusion_net_power_density_w_m3\n\n",
    );
    txt.push_str("[thermodynamic_boundary]\n");
    txt.push_str("roundtrip_output_over_input = 2 * eta_net\n");
    txt.push_str(
        "eta_net <= 0.5 implies roundtrip_output_over_input <= 1.0 (no net-positive power)\n",
    );
    txt.push_str("eta_net < 0.5 implies net energy sink for power-generation use\n\n");

    let mut csv_yield =
        String::from("scenario,eta_net,power_mw,g_per_year,mg_per_year,mg_per_day\n");
    let mut csv_density = String::from(
        "scenario,eta_net,power_mw,ng_per_year,rest_output_mw,annihilation_output_mw,rest_density_w_m3,annihilation_density_w_m3,ratio_rest_to_fusion_density,ratio_annihilation_to_fusion_density\n",
    );
    let mut csv_thermo = String::from(
        "scenario,eta_net,roundtrip_output_over_input,roundtrip_net_fraction,input_over_output_penalty,beam_mw_for_annihilation_output_equal_to_fusion_net,beam_penalty_vs_fusion_net,thermodynamic_verdict\n",
    );
    let mut thermo_rows_txt = String::new();
    let mut thermo_rows_json = Vec::<Value>::new();

    for (name, eta) in scenarios {
        txt.push_str(&format!("[yield_vs_power:{}]\n", name));
        txt.push_str(&format!("eta_net = {:.12e}\n", eta));
        txt.push_str("power_mw,g_per_year,mg_per_year,ng_per_year,mg_per_day\n");
        for p in power_grid_mw {
            let g = g_per_year(p, eta);
            let mg = g * 1000.0;
            let ng = g * 1.0e9;
            let mg_day = mg / 365.25;
            txt.push_str(&format!(
                "{:.3},{:.12e},{:.12e},{:.12e},{:.12e}\n",
                p, g, mg, ng, mg_day
            ));
            csv_yield.push_str(&format!(
                "{},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}\n",
                name, eta, p, g, mg, mg_day
            ));
        }
        txt.push('\n');

        txt.push_str(&format!("[fusion_density_vs_power:{}]\n", name));
        txt.push_str("power_mw,rest_output_mw,annihilation_output_mw,rest_density_w_m3,annihilation_density_w_m3,ratio_rest_to_fusion_density,ratio_annihilation_to_fusion_density\n");
        for p in power_grid_mw {
            let ng_per_year = g_per_year(p, eta) * 1.0e9;
            let rest_output_mw = rest_output_power_mw(p, eta);
            let ann_output_mw = annihilation_output_power_mw(p, eta);
            let rest_density = power_density_w_m3(rest_output_mw, fusion_ref.volume_m3);
            let ann_density = power_density_w_m3(ann_output_mw, fusion_ref.volume_m3);
            let rest_ratio = ratio(rest_density, fusion_density_w_m3);
            let ann_ratio = ratio(ann_density, fusion_density_w_m3);
            txt.push_str(&format!(
                "{:.3},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}\n",
                p, rest_output_mw, ann_output_mw, rest_density, ann_density, rest_ratio, ann_ratio
            ));
            csv_density.push_str(&format!(
                "{},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}\n",
                name,
                eta,
                p,
                ng_per_year,
                rest_output_mw,
                ann_output_mw,
                rest_density,
                ann_density,
                rest_ratio,
                ann_ratio
            ));
        }
        let p_req_ann = beam_mw_for_annihilation_output_mw(fusion_net_mw, eta);
        let p_req_rest = beam_mw_for_rest_output_mw(fusion_net_mw, eta);
        let roundtrip = roundtrip_output_over_input(eta);
        let roundtrip_net = roundtrip - 1.0;
        let penalty = roundtrip_penalty_input_over_output(eta);
        let beam_penalty_vs_fusion = ratio(p_req_ann, fusion_net_mw);
        let verdict = thermodynamic_verdict(eta);
        txt.push_str(&format!(
            "beam_mw_for_annihilation_output_equal_to_fusion_net = {:.12e}\n",
            p_req_ann
        ));
        txt.push_str(&format!(
            "beam_mw_for_rest_output_equal_to_fusion_net = {:.12e}\n",
            p_req_rest
        ));
        txt.push_str(&format!(
            "annihilation_match_ng_per_year = {:.12e}\n\n",
            g_per_year(p_req_ann, eta) * 1.0e9
        ));

        csv_thermo.push_str(&format!(
            "{},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{}\n",
            name,
            eta,
            roundtrip,
            roundtrip_net,
            penalty,
            p_req_ann,
            beam_penalty_vs_fusion,
            verdict
        ));
        thermo_rows_txt.push_str(&format!(
            "{},{:.6},{:.6},{:.6},{:.3},{:.3},{}\n",
            name, eta, roundtrip, roundtrip_net, p_req_ann, beam_penalty_vs_fusion, verdict
        ));
        thermo_rows_json.push(json!({
            "scenario": name,
            "eta_net": eta,
            "roundtrip_output_over_input": roundtrip,
            "roundtrip_net_fraction": roundtrip_net,
            "input_over_output_penalty": penalty,
            "beam_mw_for_annihilation_output_equal_to_fusion_net": p_req_ann,
            "beam_penalty_vs_fusion_net": beam_penalty_vs_fusion,
            "thermodynamic_verdict": verdict
        }));
    }

    let mut csv_power = String::from("scenario,eta_net,target_g_per_year,power_mw_required\n");
    txt.push_str("[power_required_for_targets]\n");
    for (name, eta) in scenarios {
        txt.push_str(&format!("scenario = {} (eta={:.12e})\n", name, eta));
        txt.push_str("target_g_per_year,power_mw_required\n");
        for t in targets_g_per_year {
            let p_req = power_mw_for_target_g_per_year(t, eta);
            txt.push_str(&format!("{:.12e},{:.12e}\n", t, p_req));
            csv_power.push_str(&format!(
                "{},{:.12e},{:.12e},{:.12e}\n",
                name, eta, t, p_req
            ));
        }
        txt.push('\n');
    }

    let target_five_g = 5.0;
    let p_modeled_5g = power_mw_for_target_g_per_year(target_five_g, eta_modeled);
    let p_pp_5g = power_mw_for_target_g_per_year(target_five_g, eta_pp_route);
    let p_abs_5g = power_mw_for_target_g_per_year(target_five_g, eta_pair_absolute);

    txt.push_str("[five_gram_reference]\n");
    txt.push_str(&format!("target_g_per_year = {:.6}\n", target_five_g));
    txt.push_str(&format!(
        "power_mw_modeled_current = {:.12e}\n",
        p_modeled_5g
    ));
    txt.push_str(&format!("power_mw_route_pp_floor = {:.12e}\n", p_pp_5g));
    txt.push_str(&format!("power_mw_pair_absolute = {:.12e}\n", p_abs_5g));
    txt.push_str(&format!(
        "target_ng_per_year = {:.12e}\n\n",
        target_five_g * 1.0e9
    ));

    txt.push_str("[plain_language_summary]\n");
    txt.push_str("interpretation = antimatter_as_primary_power_is_bounded_by_roundtrip_2eta\n");
    txt.push_str(&format!("fusion_reference_net_mw = {:.3}\n", fusion_net_mw));
    txt.push_str(&format!(
        "fusion_reference_year_energy_twh = {:.6}\n",
        fusion_year_energy_twh
    ));
    txt.push_str(&format!(
        "fusion_reference_year_energy_kilotons_tnt = {:.3}\n",
        fusion_year_energy_kt_tnt
    ));
    txt.push_str("scenario,eta_net,roundtrip_output_over_input,roundtrip_net_fraction,beam_mw_for_fusion_parity,beam_penalty_vs_fusion_net,thermodynamic_verdict\n");
    txt.push_str(&thermo_rows_txt);
    txt.push('\n');

    txt.push_str("[fusion_match_requirements]\n");
    txt.push_str(&format!("fusion_net_power_mw = {:.12e}\n", fusion_net_mw));
    txt.push_str(&format!(
        "antimatter_stream_ng_per_year_for_rest_equivalent = {:.12e}\n",
        ng_per_year_for_rest_output_power_w(fusion_ref.net_power_w)
    ));
    txt.push_str(&format!(
        "antimatter_stream_ng_per_year_for_annihilation_equivalent = {:.12e}\n",
        ng_per_year_for_annihilation_output_power_w(fusion_ref.net_power_w)
    ));

    let p_match_modeled_ann = beam_mw_for_annihilation_output_mw(fusion_net_mw, eta_modeled);
    let p_match_route_ann = beam_mw_for_annihilation_output_mw(fusion_net_mw, eta_pp_route);
    let p_match_abs_ann = beam_mw_for_annihilation_output_mw(fusion_net_mw, eta_pair_absolute);

    let json = json!({
        "meta": {
            "scope": "theoretical_output_bounds_only"
        },
        "efficiencies": {
            "modeled_current": eta_modeled,
            "route_pp_floor": eta_pp_route,
            "pair_absolute": eta_pair_absolute
        },
        "om_gaps": {
            "modeled_to_route_pp": om_gap(eta_modeled, eta_pp_route),
            "modeled_to_pair_absolute": om_gap(eta_modeled, eta_pair_absolute)
        },
        "five_gram_power_mw": {
            "modeled_current": p_modeled_5g,
            "route_pp_floor": p_pp_5g,
            "pair_absolute": p_abs_5g
        },
        "fusion_reference": {
            "source": fusion_ref.source,
            "net_power_mw": fusion_net_mw,
            "volume_m3": fusion_ref.volume_m3,
            "net_density_w_m3": fusion_density_w_m3
        },
        "fusion_match_requirements": {
            "antimatter_stream_ng_per_year_rest":
                ng_per_year_for_rest_output_power_w(fusion_ref.net_power_w),
            "antimatter_stream_ng_per_year_annihilation":
                ng_per_year_for_annihilation_output_power_w(fusion_ref.net_power_w),
            "beam_mw_for_annihilation_match": {
                "modeled_current": p_match_modeled_ann,
                "route_pp_floor": p_match_route_ann,
                "pair_absolute": p_match_abs_ann
            }
        },
        "thermodynamic_boundary": {
            "max_roundtrip_output_over_input_for_eta_le_half": 1.0,
            "statement": "eta_net <= 0.5 implies antimatter production cannot be net-positive as an energy source"
        },
        "fusion_reference_year_energy": {
            "energy_j": fusion_year_energy_j,
            "energy_twh": fusion_year_energy_twh,
            "energy_kilotons_tnt": fusion_year_energy_kt_tnt
        },
        "scenario_thermodynamic_verdicts": thermo_rows_json
    })
    .to_string();

    let txt_path = format!("{out_dir}/antimatter_theoretical_yield_bounds.txt");
    let csv_yield_path = format!("{out_dir}/antimatter_theoretical_yield_vs_power.csv");
    let csv_density_path = format!("{out_dir}/antimatter_theoretical_density_vs_power.csv");
    let csv_thermo_path = format!("{out_dir}/antimatter_theoretical_thermodynamic_verdicts.csv");
    let csv_power_path = format!("{out_dir}/antimatter_theoretical_power_for_targets.csv");
    let json_path = format!("{out_dir}/antimatter_theoretical_yield_bounds.json");

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(&csv_yield_path, csv_yield).expect("write yield csv");
    fs::write(&csv_density_path, csv_density).expect("write density csv");
    fs::write(&csv_thermo_path, csv_thermo).expect("write thermo csv");
    fs::write(&csv_power_path, csv_power).expect("write target csv");
    fs::write(&json_path, json).expect("write json");

    println!("wrote {}", txt_path);
    println!("wrote {}", csv_yield_path);
    println!("wrote {}", csv_density_path);
    println!("wrote {}", csv_thermo_path);
    println!("wrote {}", csv_power_path);
    println!("wrote {}", json_path);
}
