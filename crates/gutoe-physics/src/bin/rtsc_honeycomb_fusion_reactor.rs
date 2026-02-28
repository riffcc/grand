use gutoe_physics::{gamow_factor, ALPHA_LEADING_ORDER};
use serde_json::Value;
use std::env;
use std::f64::consts::PI;
use std::fs;
use std::path::Path;

const MU0: f64 = 4.0 * PI * 1.0e-7;
const KB: f64 = 1.380_649e-23;
const KEV_TO_K: f64 = 1.160_451_812e7;
const MEV_TO_J: f64 = 1.602_176_634e-13;
const E_FUSION_DHE3_MEV: f64 = 18.353_055;
const E_FUSION_DHE3_J: f64 = E_FUSION_DHE3_MEV * MEV_TO_J;
const PARTICLE_PRESSURE_FACTOR: f64 = 2.5; // ions + electrons in D-He3 mix
const CALIBRATED_SIGMAV_AT_100KEV: f64 = 1.0e-22; // m^3/s, reference proxy anchor
const MAGNET_LOAD_COEFF: f64 = 60.0; // W/(m^2*T^2) RTSC coil overhead proxy

#[derive(Clone, Debug)]
struct Material {
    symbol: String,
    z: u16,
    tc_ideal_k: f64,
    margin_k: f64,
    pass_fraction_grid: f64,
    engineering_difficulty: f64,
    hazard: f64,
}

#[derive(Clone, Copy, Debug)]
struct ShapeCfg {
    name: &'static str,
    confinement_base: f64,
    beta_factor: f64,
    field_utilization: f64,
    electric_eff_factor: f64,
}

#[derive(Clone, Debug)]
struct ReactorPoint {
    material: String,
    shape: &'static str,
    mesh_quality: f64,
    r_major_m: f64,
    aspect: f64,
    a_minor_m: f64,
    t_kev: f64,
    b_operating_t: f64,
    n_fuel_m3: f64,
    volume_m3: f64,
    surface_m2: f64,
    confinement: f64,
    p_fusion_w: f64,
    p_electric_gross_w: f64,
    p_recirc_w: f64,
    p_net_w: f64,
    q_engineering: f64,
    hoop_stress_pa: f64,
    allowable_stress_pa: f64,
    score: f64,
}

fn reduced_mass_mev(a1: f64, a2: f64) -> f64 {
    (a1 * a2 / (a1 + a2)) * 931.494
}

fn gamow_dhe3(t_kev: f64) -> f64 {
    if t_kev <= 0.0 {
        return 0.0;
    }
    let alpha_eff = ALPHA_LEADING_ORDER * 1.0 * 2.0;
    let m_reduced = reduced_mass_mev(2.0, 3.0);
    let e_cm_mev = t_kev / 1000.0;
    gamow_factor(alpha_eff, m_reduced, e_cm_mev).unwrap_or(0.0)
}

fn sigma_v_dhe3(t_kev: f64) -> f64 {
    let g_ref = gamow_dhe3(100.0).max(1.0e-30);
    let g = gamow_dhe3(t_kev);
    let thermal = (t_kev / 100.0).max(0.05).sqrt();
    let v = CALIBRATED_SIGMAV_AT_100KEV * (g / g_ref) * thermal;
    v.clamp(1.0e-26, 5.0e-22)
}

fn fallback_materials() -> Vec<Material> {
    vec![
        Material {
            symbol: "Mo".to_string(),
            z: 42,
            tc_ideal_k: 370.125,
            margin_k: 70.125,
            pass_fraction_grid: 0.079969,
            engineering_difficulty: 0.30,
            hazard: 0.15,
        },
        Material {
            symbol: "Cr".to_string(),
            z: 24,
            tc_ideal_k: 366.0,
            margin_k: 66.0,
            pass_fraction_grid: 0.073049,
            engineering_difficulty: 0.25,
            hazard: 0.20,
        },
        Material {
            symbol: "Pt".to_string(),
            z: 78,
            tc_ideal_k: 371.5,
            margin_k: 71.5,
            pass_fraction_grid: 0.084198,
            engineering_difficulty: 0.55,
            hazard: 0.10,
        },
        Material {
            symbol: "Hf".to_string(),
            z: 72,
            tc_ideal_k: 372.875,
            margin_k: 72.875,
            pass_fraction_grid: 0.086505,
            engineering_difficulty: 0.70,
            hazard: 0.35,
        },
        Material {
            symbol: "Zn".to_string(),
            z: 30,
            tc_ideal_k: 370.125,
            margin_k: 70.125,
            pass_fraction_grid: 0.079969,
            engineering_difficulty: 0.75,
            hazard: 0.40,
        },
        Material {
            symbol: "Cd".to_string(),
            z: 48,
            tc_ideal_k: 368.75,
            margin_k: 68.75,
            pass_fraction_grid: 0.078431,
            engineering_difficulty: 0.85,
            hazard: 0.95,
        },
    ]
}

fn load_materials(plan_json_path: &Path) -> Vec<Material> {
    let parsed = fs::read_to_string(plan_json_path)
        .ok()
        .and_then(|s| serde_json::from_str::<Value>(&s).ok());

    let mut out = Vec::new();
    if let Some(v) = parsed {
        if let Some(arr) = v.get("plan").and_then(|x| x.as_array()) {
            for row in arr {
                let symbol = row
                    .get("symbol")
                    .and_then(|x| x.as_str())
                    .unwrap_or("")
                    .to_string();
                if symbol.is_empty() {
                    continue;
                }
                let z = row.get("z").and_then(|x| x.as_u64()).unwrap_or(0) as u16;
                let tc_ideal_k = row
                    .get("tc_ideal_k")
                    .and_then(|x| x.as_f64())
                    .unwrap_or(0.0);
                let margin_k = row.get("margin_k").and_then(|x| x.as_f64()).unwrap_or(0.0);
                let pass_fraction_grid = row
                    .get("pass_fraction_grid")
                    .and_then(|x| x.as_f64())
                    .unwrap_or(0.0);
                let engineering_difficulty = row
                    .get("engineering_difficulty")
                    .and_then(|x| x.as_f64())
                    .unwrap_or(0.5);
                let hazard = row.get("hazard").and_then(|x| x.as_f64()).unwrap_or(0.5);

                out.push(Material {
                    symbol,
                    z,
                    tc_ideal_k,
                    margin_k,
                    pass_fraction_grid,
                    engineering_difficulty,
                    hazard,
                });
            }
        }
    }

    if out.is_empty() {
        fallback_materials()
    } else {
        out
    }
}

fn shape_catalog() -> [ShapeCfg; 3] {
    [
        ShapeCfg {
            name: "toroidal_honeycomb",
            confinement_base: 1.00,
            beta_factor: 1.00,
            field_utilization: 0.95,
            electric_eff_factor: 1.00,
        },
        ShapeCfg {
            name: "spherical_honeycomb",
            confinement_base: 0.87,
            beta_factor: 0.90,
            field_utilization: 0.90,
            electric_eff_factor: 0.92,
        },
        ShapeCfg {
            name: "cylindrical_honeycomb",
            confinement_base: 0.64,
            beta_factor: 0.72,
            field_utilization: 0.82,
            electric_eff_factor: 0.86,
        },
    ]
}

fn torus_volume_surface(r_major: f64, a_minor: f64) -> (f64, f64) {
    let volume = 2.0 * PI * PI * r_major * a_minor.powi(2);
    let surface = 4.0 * PI * PI * r_major * a_minor;
    (volume, surface)
}

fn sphere_volume_surface(radius: f64) -> (f64, f64) {
    let volume = (4.0 / 3.0) * PI * radius.powi(3);
    let surface = 4.0 * PI * radius.powi(2);
    (volume, surface)
}

fn cylinder_volume_surface(radius: f64, length: f64) -> (f64, f64) {
    let volume = PI * radius.powi(2) * length;
    let surface = 2.0 * PI * radius * length + 2.0 * PI * radius.powi(2);
    (volume, surface)
}

fn evaluate_point(
    m: &Material,
    s: ShapeCfg,
    mesh_quality: f64,
    r_major_m: f64,
    aspect: f64,
    t_kev: f64,
) -> Option<ReactorPoint> {
    if !(0.8..=220.0).contains(&r_major_m) || !(2.2..=7.5).contains(&aspect) || t_kev <= 0.0 {
        return None;
    }

    let a_minor_m = (r_major_m / aspect).max(0.10);

    let base_b_max = 20.0 + 0.17 * m.margin_k + 80.0 * m.pass_fraction_grid;
    let engineering_factor = (1.0 - 0.22 * m.engineering_difficulty - 0.08 * m.hazard).max(0.35);
    let mesh_field_gain = 0.90 + 0.10 * mesh_quality.clamp(0.0, 1.0);
    let b_max_t = base_b_max * engineering_factor * mesh_field_gain;
    let b_operating_t = 0.82 * b_max_t * s.field_utilization;

    let beta_target = 0.045 * s.beta_factor * mesh_quality.powf(1.2);
    let p_magnetic = b_operating_t * b_operating_t / (2.0 * MU0);
    let p_plasma = beta_target * p_magnetic;

    let t_k = t_kev * KEV_TO_K;
    let n_fuel_m3 = p_plasma / (PARTICLE_PRESSURE_FACTOR * KB * t_k);

    let (volume_m3, surface_m2) = match s.name {
        "toroidal_honeycomb" => torus_volume_surface(r_major_m, a_minor_m),
        "spherical_honeycomb" => sphere_volume_surface(a_minor_m * 1.35),
        "cylindrical_honeycomb" => cylinder_volume_surface(a_minor_m, 2.0 * PI * r_major_m),
        _ => return None,
    };

    let finite_small = 1.0 - (-a_minor_m / 0.8).exp();
    let finite_large = (-(r_major_m / 120.0).powf(1.2)).exp();
    let confinement = (s.confinement_base * mesh_quality.powf(1.4) * finite_small * finite_large)
        .clamp(0.0, 1.0);

    let sigma_v = sigma_v_dhe3(t_kev);
    let reaction_rate = 0.25 * n_fuel_m3.powi(2) * sigma_v;
    let p_fusion_w = reaction_rate * E_FUSION_DHE3_J * volume_m3 * confinement;

    let electric_eff = 0.62 * s.electric_eff_factor;
    let p_electric_gross_w = p_fusion_w * electric_eff;

    let p_magnet_overhead = surface_m2
        * b_operating_t.powi(2)
        * MAGNET_LOAD_COEFF
        * (1.0 + 0.6 * (1.0 - mesh_quality));

    let p_heating = p_fusion_w * (0.12 + 0.28 * (1.0 - confinement));
    let p_aux = p_fusion_w * (0.04 + 0.06 * m.engineering_difficulty + 0.03 * m.hazard);
    let p_recirc_w = p_magnet_overhead + p_heating + p_aux;

    let p_net_w = p_electric_gross_w - p_recirc_w;
    let q_engineering = if p_recirc_w > 0.0 {
        p_electric_gross_w / p_recirc_w
    } else {
        f64::INFINITY
    };

    let coil_thickness_m = 0.35 + 0.12 * a_minor_m * mesh_quality;
    let hoop_stress_pa = p_magnetic * (r_major_m / coil_thickness_m);
    let allowable_stress_pa = 1.15e9 * (1.0 - 0.35 * m.engineering_difficulty).max(0.40);

    if hoop_stress_pa > allowable_stress_pa {
        return None;
    }

    let score = p_net_w.max(0.0) * (q_engineering / 8.0).clamp(0.0, 1.8) * (0.4 + 0.6 * confinement);

    Some(ReactorPoint {
        material: m.symbol.clone(),
        shape: s.name,
        mesh_quality,
        r_major_m,
        aspect,
        a_minor_m,
        t_kev,
        b_operating_t,
        n_fuel_m3,
        volume_m3,
        surface_m2,
        confinement,
        p_fusion_w,
        p_electric_gross_w,
        p_recirc_w,
        p_net_w,
        q_engineering,
        hoop_stress_pa,
        allowable_stress_pa,
        score,
    })
}

fn fmt_mw(w: f64) -> f64 {
    w / 1.0e6
}

fn main() {
    let out_dir = env::var("GUTOE_FUSION_REACTOR_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/rtsc_honeycomb_fusion_reactor".to_string());
    fs::create_dir_all(&out_dir).expect("create output dir");

    let routes_json = env::var("GUTOE_RTSC_ROUTES_JSON").unwrap_or_else(|_| {
        "/tmp/bh_renders/rtsc_synthesis_routes/rtsc_synthesis_routes.json".to_string()
    });
    let materials = load_materials(Path::new(&routes_json));

    let meshes = [0.86, 0.90, 0.94, 0.97, 1.00];
    let aspects = [2.4, 3.0, 3.6, 4.2, 5.0, 6.0, 7.0];
    let temps_kev = [60.0, 80.0, 100.0, 120.0, 150.0, 180.0];

    let mut radii = Vec::new();
    let mut r = 0.8;
    while r <= 220.0 {
        radii.push(r);
        r *= 1.115;
    }

    let mut points = Vec::<ReactorPoint>::new();
    for material in &materials {
        for shape in shape_catalog() {
            for &mesh in &meshes {
                for &aspect in &aspects {
                    for &t in &temps_kev {
                        for &r_major in &radii {
                            if let Some(p) = evaluate_point(material, shape, mesh, r_major, aspect, t) {
                                points.push(p);
                            }
                        }
                    }
                }
            }
        }
    }

    points.sort_by(|a, b| b.score.total_cmp(&a.score));

    let best = points.first().cloned();
    let best_torus_honeycomb = points
        .iter()
        .filter(|p| p.shape == "toroidal_honeycomb" && p.mesh_quality >= 0.99)
        .max_by(|a, b| a.score.total_cmp(&b.score))
        .cloned();

    let smallest_viable = points
        .iter()
        .filter(|p| p.p_net_w >= 10.0e6 && p.q_engineering >= 1.7 && p.confinement >= 0.42)
        .min_by(|a, b| a.r_major_m.total_cmp(&b.r_major_m))
        .cloned();

    let largest_viable = points
        .iter()
        .filter(|p| p.p_net_w >= 1.0e6 && p.q_engineering >= 1.5 && p.confinement >= 0.30)
        .max_by(|a, b| a.r_major_m.total_cmp(&b.r_major_m))
        .cloned();

    let largest_positive = points
        .iter()
        .filter(|p| p.p_net_w >= 0.0 && p.q_engineering >= 1.0 && p.confinement >= 0.20)
        .max_by(|a, b| a.r_major_m.total_cmp(&b.r_major_m))
        .cloned();

    let mut by_material: Vec<(String, usize, f64)> = materials
        .iter()
        .map(|m| {
            let pool: Vec<&ReactorPoint> = points.iter().filter(|p| p.material == m.symbol).collect();
            let max_net = pool
                .iter()
                .map(|p| p.p_net_w)
                .fold(0.0_f64, |a, b| a.max(b));
            (m.symbol.clone(), pool.len(), max_net)
        })
        .collect();
    by_material.sort_by(|a, b| b.2.total_cmp(&a.2));

    let mut csv = String::from(
        "material,shape,mesh_quality,r_major_m,aspect,a_minor_m,t_kev,b_operating_t,n_fuel_m3,volume_m3,surface_m2,confinement,p_fusion_w,p_electric_gross_w,p_recirc_w,p_net_w,q_engineering,hoop_stress_pa,allowable_stress_pa,score\n",
    );
    for p in &points {
        csv.push_str(&format!(
            "{},{},{:.3},{:.6},{:.3},{:.6},{:.3},{:.6},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e}\n",
            p.material,
            p.shape,
            p.mesh_quality,
            p.r_major_m,
            p.aspect,
            p.a_minor_m,
            p.t_kev,
            p.b_operating_t,
            p.n_fuel_m3,
            p.volume_m3,
            p.surface_m2,
            p.confinement,
            p.p_fusion_w,
            p.p_electric_gross_w,
            p.p_recirc_w,
            p.p_net_w,
            p.q_engineering,
            p.hoop_stress_pa,
            p.allowable_stress_pa,
            p.score
        ));
    }

    let mut txt = String::new();
    txt.push_str("[rtsc_honeycomb_fusion_reactor]\n");
    txt.push_str("fuel = D+He3\n");
    txt.push_str("mode = beta_limited_magnetic_confinement_with_rtsc_mesh\n");
    txt.push_str(&format!("materials_source = {}\n", routes_json));
    txt.push_str(&format!("candidate_materials = {}\n", materials.len()));
    txt.push_str(&format!("evaluated_points = {}\n\n", points.len()));

    if let Some(p) = &best {
        txt.push_str("[best_overall]\n");
        txt.push_str(&format!("material = {}\n", p.material));
        txt.push_str(&format!("shape = {}\n", p.shape));
        txt.push_str(&format!("mesh_quality = {:.3}\n", p.mesh_quality));
        txt.push_str(&format!("R_major_m = {:.3}\n", p.r_major_m));
        txt.push_str(&format!("a_minor_m = {:.3}\n", p.a_minor_m));
        txt.push_str(&format!("aspect = {:.3}\n", p.aspect));
        txt.push_str(&format!("T_keV = {:.3}\n", p.t_kev));
        txt.push_str(&format!("B_operating_T = {:.3}\n", p.b_operating_t));
        txt.push_str(&format!("n_fuel_m^-3 = {:.6e}\n", p.n_fuel_m3));
        txt.push_str(&format!("confinement = {:.4}\n", p.confinement));
        txt.push_str(&format!("P_fusion_MW = {:.3}\n", fmt_mw(p.p_fusion_w)));
        txt.push_str(&format!("P_electric_gross_MW = {:.3}\n", fmt_mw(p.p_electric_gross_w)));
        txt.push_str(&format!("P_recirc_MW = {:.3}\n", fmt_mw(p.p_recirc_w)));
        txt.push_str(&format!("P_net_MW = {:.3}\n", fmt_mw(p.p_net_w)));
        txt.push_str(&format!("Q_engineering = {:.3}\n", p.q_engineering));
        txt.push_str(&format!("hoop_stress_GPa = {:.3}\n", p.hoop_stress_pa / 1.0e9));
        txt.push_str(&format!(
            "allowable_stress_GPa = {:.3}\n\n",
            p.allowable_stress_pa / 1.0e9
        ));
    }

    if let Some(p) = &best_torus_honeycomb {
        txt.push_str("[best_toroidal_honeycomb_perfect_mesh]\n");
        txt.push_str(&format!("material = {}\n", p.material));
        txt.push_str(&format!("R_major_m = {:.3}\n", p.r_major_m));
        txt.push_str(&format!("a_minor_m = {:.3}\n", p.a_minor_m));
        txt.push_str(&format!("T_keV = {:.3}\n", p.t_kev));
        txt.push_str(&format!("B_operating_T = {:.3}\n", p.b_operating_t));
        txt.push_str(&format!("P_net_MW = {:.3}\n", fmt_mw(p.p_net_w)));
        txt.push_str(&format!("Q_engineering = {:.3}\n\n", p.q_engineering));
    }

    if let Some(p) = &smallest_viable {
        txt.push_str("[smallest_viable]\n");
        txt.push_str("criteria = P_net>=10MW && Q_engineering>=1.7 && confinement>=0.42\n");
        txt.push_str(&format!("material = {}\n", p.material));
        txt.push_str(&format!("shape = {}\n", p.shape));
        txt.push_str(&format!("mesh_quality = {:.3}\n", p.mesh_quality));
        txt.push_str(&format!("R_major_m = {:.3}\n", p.r_major_m));
        txt.push_str(&format!("a_minor_m = {:.3}\n", p.a_minor_m));
        txt.push_str(&format!("P_net_MW = {:.3}\n", fmt_mw(p.p_net_w)));
        txt.push_str(&format!("Q_engineering = {:.3}\n\n", p.q_engineering));
    } else {
        txt.push_str("[smallest_viable]\nnone\n\n");
    }

    if let Some(p) = &largest_viable {
        txt.push_str("[largest_viable]\n");
        txt.push_str("criteria = P_net>=1MW && Q_engineering>=1.5 && confinement>=0.30\n");
        txt.push_str(&format!("material = {}\n", p.material));
        txt.push_str(&format!("shape = {}\n", p.shape));
        txt.push_str(&format!("mesh_quality = {:.3}\n", p.mesh_quality));
        txt.push_str(&format!("R_major_m = {:.3}\n", p.r_major_m));
        txt.push_str(&format!("a_minor_m = {:.3}\n", p.a_minor_m));
        txt.push_str(&format!("P_net_MW = {:.3}\n", fmt_mw(p.p_net_w)));
        txt.push_str(&format!("Q_engineering = {:.3}\n\n", p.q_engineering));
    } else {
        txt.push_str("[largest_viable]\nnone\n\n");
    }

    if let Some(p) = &largest_positive {
        txt.push_str("[largest_positive]\n");
        txt.push_str("criteria = P_net>=0 && Q_engineering>=1.0 && confinement>=0.20\n");
        txt.push_str(&format!("material = {}\n", p.material));
        txt.push_str(&format!("shape = {}\n", p.shape));
        txt.push_str(&format!("mesh_quality = {:.3}\n", p.mesh_quality));
        txt.push_str(&format!("R_major_m = {:.3}\n", p.r_major_m));
        txt.push_str(&format!("a_minor_m = {:.3}\n", p.a_minor_m));
        txt.push_str(&format!("P_net_MW = {:.3}\n", fmt_mw(p.p_net_w)));
        txt.push_str(&format!("Q_engineering = {:.3}\n\n", p.q_engineering));
    } else {
        txt.push_str("[largest_positive]\nnone\n\n");
    }

    txt.push_str("[material_ranking_by_max_net]\n");
    txt.push_str("rank,material,max_net_mw,evaluated_points\n");
    for (i, (symbol, n, max_net)) in by_material.iter().enumerate() {
        txt.push_str(&format!("{},{},{:.3},{}\n", i + 1, symbol, fmt_mw(*max_net), n));
    }

    let best_json = best.as_ref().map(|p| {
        format!(
            concat!(
                "{{\"material\":\"{}\",\"shape\":\"{}\",\"mesh_quality\":{:.6},\"r_major_m\":{:.6},\"a_minor_m\":{:.6},\"aspect\":{:.6},\"t_kev\":{:.6},",
                "\"b_operating_t\":{:.6},\"n_fuel_m3\":{:.6e},\"confinement\":{:.6},\"p_net_w\":{:.6e},\"q_engineering\":{:.6}}}"
            ),
            p.material,
            p.shape,
            p.mesh_quality,
            p.r_major_m,
            p.a_minor_m,
            p.aspect,
            p.t_kev,
            p.b_operating_t,
            p.n_fuel_m3,
            p.confinement,
            p.p_net_w,
            p.q_engineering
        )
    });

    let smallest_json = smallest_viable.as_ref().map(|p| {
        format!(
            "{{\"material\":\"{}\",\"shape\":\"{}\",\"mesh_quality\":{:.6},\"r_major_m\":{:.6},\"a_minor_m\":{:.6},\"p_net_w\":{:.6e},\"q_engineering\":{:.6}}}",
            p.material, p.shape, p.mesh_quality, p.r_major_m, p.a_minor_m, p.p_net_w, p.q_engineering
        )
    });

    let largest_json = largest_viable.as_ref().map(|p| {
        format!(
            "{{\"material\":\"{}\",\"shape\":\"{}\",\"mesh_quality\":{:.6},\"r_major_m\":{:.6},\"a_minor_m\":{:.6},\"p_net_w\":{:.6e},\"q_engineering\":{:.6}}}",
            p.material, p.shape, p.mesh_quality, p.r_major_m, p.a_minor_m, p.p_net_w, p.q_engineering
        )
    });

    let largest_positive_json = largest_positive.as_ref().map(|p| {
        format!(
            "{{\"material\":\"{}\",\"shape\":\"{}\",\"mesh_quality\":{:.6},\"r_major_m\":{:.6},\"a_minor_m\":{:.6},\"p_net_w\":{:.6e},\"q_engineering\":{:.6}}}",
            p.material, p.shape, p.mesh_quality, p.r_major_m, p.a_minor_m, p.p_net_w, p.q_engineering
        )
    });

    let mut by_mat_json = String::new();
    by_mat_json.push('[');
    for (i, (symbol, n, max_net)) in by_material.iter().enumerate() {
        by_mat_json.push_str(&format!(
            "{{\"rank\":{},\"material\":\"{}\",\"max_net_w\":{:.6e},\"evaluated_points\":{}}}{}",
            i + 1,
            symbol,
            max_net,
            n,
            if i + 1 == by_material.len() { "" } else { "," }
        ));
    }
    by_mat_json.push(']');

    let json = format!(
        concat!(
            "{{\n",
            "  \"meta\": {{\"fuel\": \"D+He3\", \"mode\": \"beta_limited_magnetic_confinement_with_rtsc_mesh\", \"materials_source\": \"{}\", \"candidate_materials\": {}, \"evaluated_points\": {}}},\n",
            "  \"best_overall\": {},\n",
            "  \"smallest_viable\": {},\n",
            "  \"largest_viable\": {},\n",
            "  \"largest_positive\": {},\n",
            "  \"material_ranking_by_max_net\": {}\n",
            "}}\n"
        ),
        routes_json,
        materials.len(),
        points.len(),
        best_json.unwrap_or_else(|| "null".to_string()),
        smallest_json.unwrap_or_else(|| "null".to_string()),
        largest_json.unwrap_or_else(|| "null".to_string()),
        largest_positive_json.unwrap_or_else(|| "null".to_string()),
        by_mat_json
    );

    let txt_path = format!("{out_dir}/rtsc_honeycomb_fusion_reactor.txt");
    let json_path = format!("{out_dir}/rtsc_honeycomb_fusion_reactor.json");
    let csv_path = format!("{out_dir}/rtsc_honeycomb_fusion_reactor_all_points.csv");

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(&json_path, json).expect("write json");
    fs::write(&csv_path, csv).expect("write csv");

    println!("wrote {}", txt_path);
    println!("wrote {}", json_path);
    println!("wrote {}", csv_path);
}
