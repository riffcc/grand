//! One-way geometry-FTL loophole probe.
//!
//! Hypothesis under test:
//! - keep local propagation subluminal (`v_local <= c`)
//! - avoid CTC by forcing one-way corridor topology
//! - avoid negative energy by sourcing only positive-energy dark-sector budget
//! - ask whether useful FTL (v_eff > c) is energetically feasible

use gutoe_physics::constants::{C, G, PLANCK_LENGTH};
use gutoe_physics::dark_sector::{dark_density_particle, vacuum_energy_density_structural};
use gutoe_physics::equations::ModifiedEinstein;
use gutoe_physics::singularity_resolution::lattice_core_radius_m;
use serde_json::json;
use std::f64::consts::PI;
use std::fs;
use std::path::PathBuf;

#[derive(Debug, Clone, Copy)]
struct ProbePoint {
    v_eff_over_c: f64,
    radius_m: f64,
    shortcut_s: f64,
    rho_required_classical_j_m3: f64,
    rho_required_gutoe_scale_radius_j_m3: f64,
    rho_required_gutoe_scale_planck_j_m3: f64,
    ratio_classical_over_source: f64,
    ratio_gutoe_scale_radius_over_source: f64,
    ratio_gutoe_scale_planck_over_source: f64,
    ratio_gutoe_scale_radius_rear10_over_source: f64,
    ratio_gutoe_scale_planck_rear10_over_source: f64,
}

fn required_curvature_energy_density_j_m3_with_g(
    shortcut_s: f64,
    radius_m: f64,
    g_eff: f64,
) -> f64 {
    // Einstein scale: ρ ~ (c^4 / 8πG) * K, with K ~ ((1/s)-1)/R² for a corridor strain proxy.
    let pref = C.powi(4) / (8.0 * std::f64::consts::PI * g_eff);
    let amp = (1.0 / shortcut_s - 1.0).max(0.0);
    pref * amp / radius_m.powi(2)
}

fn main() {
    let out_dir = std::env::var("GUTOE_FTL_ONE_WAY_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ftl_one_way_lentz_probe".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    // Optimistic positive-energy source lane:
    // visible mass density = dense interstellar cloud (~1e-20 kg/m³, optimistic),
    // dark density from structural ratio (5/11), plus structural vacuum density.
    let rho_visible_kg_m3 = 1.0e-20_f64;
    let rho_dark_kg_m3 = dark_density_particle(rho_visible_kg_m3);
    let rho_vac_kg_m3 = vacuum_energy_density_structural();
    let rho_source_kg_m3 = rho_dark_kg_m3 + rho_vac_kg_m3;
    let rho_source_j_m3 = rho_source_kg_m3 * C * C;

    // "Useful" FTL regime scan: 2c..100c and corridor radius 100 m..10,000 km.
    let v_grid = [2.0, 3.0, 5.0, 10.0, 20.0, 50.0, 100.0];
    let r_grid = [1.0e2, 1.0e3, 1.0e4, 1.0e5, 1.0e6, 1.0e7];

    let rear_face_factor = 0.1_f64;
    let einstein = ModifiedEinstein::new();
    let g_eff_planck = einstein.effective_g(PLANCK_LENGTH);

    let mut points = Vec::new();
    let mut best: Option<ProbePoint> = None;
    for &v_eff_over_c in &v_grid {
        let shortcut_s = 1.0 / v_eff_over_c; // coordinate speed = c/s
        for &radius_m in &r_grid {
            let g_eff_radius = einstein.effective_g(radius_m);
            let rho_required_classical_j_m3 =
                required_curvature_energy_density_j_m3_with_g(shortcut_s, radius_m, G);
            let rho_required_gutoe_scale_radius_j_m3 =
                required_curvature_energy_density_j_m3_with_g(shortcut_s, radius_m, g_eff_radius);
            let rho_required_gutoe_scale_planck_j_m3 =
                required_curvature_energy_density_j_m3_with_g(shortcut_s, radius_m, g_eff_planck);

            let ratio_classical_over_source = rho_required_classical_j_m3 / rho_source_j_m3;
            let ratio_gutoe_scale_radius_over_source =
                rho_required_gutoe_scale_radius_j_m3 / rho_source_j_m3;
            let ratio_gutoe_scale_planck_over_source =
                rho_required_gutoe_scale_planck_j_m3 / rho_source_j_m3;
            let ratio_gutoe_scale_radius_rear10_over_source =
                rear_face_factor * ratio_gutoe_scale_radius_over_source;
            let ratio_gutoe_scale_planck_rear10_over_source =
                rear_face_factor * ratio_gutoe_scale_planck_over_source;

            let p = ProbePoint {
                v_eff_over_c,
                radius_m,
                shortcut_s,
                rho_required_classical_j_m3,
                rho_required_gutoe_scale_radius_j_m3,
                rho_required_gutoe_scale_planck_j_m3,
                ratio_classical_over_source,
                ratio_gutoe_scale_radius_over_source,
                ratio_gutoe_scale_planck_over_source,
                ratio_gutoe_scale_radius_rear10_over_source,
                ratio_gutoe_scale_planck_rear10_over_source,
            };
            if best
                .map(|b| {
                    p.ratio_gutoe_scale_radius_over_source < b.ratio_gutoe_scale_radius_over_source
                })
                .unwrap_or(true)
            {
                best = Some(p);
            }
            points.push(p);
        }
    }

    // Large-geometry desperation scan: if anything passes, it will pass only here.
    let mega_r_grid = [1.0e10, 1.0e12, 1.0e14, 1.0e16, 1.0e18, 1.0e20, 1.0e22];
    let mut best_mega: Option<ProbePoint> = None;
    for &radius_m in &mega_r_grid {
        let v_eff_over_c = 2.0;
        let shortcut_s = 1.0 / v_eff_over_c;
        let g_eff_radius = einstein.effective_g(radius_m);
        let rho_required_classical_j_m3 =
            required_curvature_energy_density_j_m3_with_g(shortcut_s, radius_m, G);
        let rho_required_gutoe_scale_radius_j_m3 =
            required_curvature_energy_density_j_m3_with_g(shortcut_s, radius_m, g_eff_radius);
        let rho_required_gutoe_scale_planck_j_m3 =
            required_curvature_energy_density_j_m3_with_g(shortcut_s, radius_m, g_eff_planck);

        let ratio_classical_over_source = rho_required_classical_j_m3 / rho_source_j_m3;
        let ratio_gutoe_scale_radius_over_source =
            rho_required_gutoe_scale_radius_j_m3 / rho_source_j_m3;
        let ratio_gutoe_scale_planck_over_source =
            rho_required_gutoe_scale_planck_j_m3 / rho_source_j_m3;
        let ratio_gutoe_scale_radius_rear10_over_source =
            rear_face_factor * ratio_gutoe_scale_radius_over_source;
        let ratio_gutoe_scale_planck_rear10_over_source =
            rear_face_factor * ratio_gutoe_scale_planck_over_source;

        let p = ProbePoint {
            v_eff_over_c,
            radius_m,
            shortcut_s,
            rho_required_classical_j_m3,
            rho_required_gutoe_scale_radius_j_m3,
            rho_required_gutoe_scale_planck_j_m3,
            ratio_classical_over_source,
            ratio_gutoe_scale_radius_over_source,
            ratio_gutoe_scale_planck_over_source,
            ratio_gutoe_scale_radius_rear10_over_source,
            ratio_gutoe_scale_planck_rear10_over_source,
        };
        if best_mega
            .map(|b| {
                p.ratio_gutoe_scale_radius_over_source < b.ratio_gutoe_scale_radius_over_source
            })
            .unwrap_or(true)
        {
            best_mega = Some(p);
        }
    }

    let best = best.expect("nonempty probe");
    let best_mega = best_mega.expect("nonempty mega probe");
    let practical_feasible = best.ratio_gutoe_scale_radius_over_source <= 1.0;
    let practical_feasible_rear10 = best.ratio_gutoe_scale_radius_rear10_over_source <= 1.0;
    let practical_feasible_planck_rear10 = best.ratio_gutoe_scale_planck_rear10_over_source <= 1.0;
    let mega_feasible = best_mega.ratio_gutoe_scale_radius_over_source <= 1.0;
    let mega_feasible_rear10 = best_mega.ratio_gutoe_scale_radius_rear10_over_source <= 1.0;
    let mega_feasible_planck_rear10 = best_mega.ratio_gutoe_scale_planck_rear10_over_source <= 1.0;

    let qg_radius_gain =
        best.ratio_classical_over_source / best.ratio_gutoe_scale_radius_over_source;
    let qg_planck_gain =
        best.ratio_classical_over_source / best.ratio_gutoe_scale_planck_over_source;

    // ── Tube-geometry zoom (user-requested) ──────────────────────────────────
    let tube_length_m = 1.0e22_f64;
    let tube_radius_m = 100.0_f64;
    let tube_shortcut_s = 0.5_f64; // 2c effective
    let tube_volume_m3 = PI * tube_radius_m.powi(2) * tube_length_m;
    let l_sun_w = 3.828e26_f64;
    let seconds_per_year = 365.25_f64 * 24.0 * 3600.0;

    // Curvature scale assumptions:
    // (A) optimistic longitudinal mode: curvature set by tube length
    // (B) conservative transverse mode: curvature set by tube radius
    let rho_req_long_classical =
        required_curvature_energy_density_j_m3_with_g(tube_shortcut_s, tube_length_m, G);
    let rho_req_trans_classical =
        required_curvature_energy_density_j_m3_with_g(tube_shortcut_s, tube_radius_m, G);

    let rho_req_long_rear10 = rear_face_factor * rho_req_long_classical;
    let rho_req_trans_rear10 = rear_face_factor * rho_req_trans_classical;

    let e_req_long_classical_j = rho_req_long_classical * tube_volume_m3;
    let e_req_trans_classical_j = rho_req_trans_classical * tube_volume_m3;
    let e_req_long_rear10_j = rho_req_long_rear10 * tube_volume_m3;
    let e_req_trans_rear10_j = rho_req_trans_rear10 * tube_volume_m3;

    let sun_seconds_long_classical = e_req_long_classical_j / l_sun_w;
    let sun_seconds_trans_classical = e_req_trans_classical_j / l_sun_w;
    let sun_seconds_long_rear10 = e_req_long_rear10_j / l_sun_w;
    let sun_seconds_trans_rear10 = e_req_trans_rear10_j / l_sun_w;

    let sun_years_long_classical = sun_seconds_long_classical / seconds_per_year;
    let sun_years_trans_classical = sun_seconds_trans_classical / seconds_per_year;
    let sun_years_long_rear10 = sun_seconds_long_rear10 / seconds_per_year;
    let sun_years_trans_rear10 = sun_seconds_trans_rear10 / seconds_per_year;

    // ── Compact self-contained shell sweep ───────────────────────────────────
    let compact_radius_grid_m = [10.0_f64, 50.0, 100.0];
    let compact_thickness_grid_m = [0.1_f64, 1.0, 5.0];
    let compact_v_eff_grid = [2.0_f64, 10.0];

    let mut compact_rows_txt = String::new();
    compact_rows_txt.push_str(
        "columns=radius_m,thickness_m,v_eff_over_c,shortcut_s,curvature_scale_model,shell_volume_m3,rho_req_j_m3,e_req_j,sun_years,e_req_rear10_j,sun_years_rear10\n",
    );
    let mut compact_rows_json = Vec::new();
    let mut compact_best_radius_model: Option<(f64, f64, f64, f64, f64, f64, f64)> = None;
    let mut compact_best_thickness_model: Option<(f64, f64, f64, f64, f64, f64, f64)> = None;
    let mut compact_best_core_model: Option<(f64, f64, f64, f64, f64, f64, f64)> = None;
    let r_core_m = lattice_core_radius_m(PLANCK_LENGTH);
    for &r_bubble in &compact_radius_grid_m {
        for &thickness in &compact_thickness_grid_m {
            for &v_eff in &compact_v_eff_grid {
                let shortcut_s = 1.0 / v_eff;
                let shell_volume_m3 =
                    (4.0 / 3.0) * PI * ((r_bubble + thickness).powi(3) - r_bubble.powi(3));

                // Model A (optimistic for compact): curvature set by bubble radius.
                let rho_req_r =
                    required_curvature_energy_density_j_m3_with_g(shortcut_s, r_bubble, G);
                let e_req_r = rho_req_r * shell_volume_m3;
                let e_req_r_rear10 = rear_face_factor * e_req_r;
                let sun_years_r = (e_req_r / l_sun_w) / seconds_per_year;
                let sun_years_r_rear10 = (e_req_r_rear10 / l_sun_w) / seconds_per_year;

                compact_rows_txt.push_str(&format!(
                    "{:.6e},{:.6e},{:.3},{:.6e},radius,{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e}\n",
                    r_bubble,
                    thickness,
                    v_eff,
                    shortcut_s,
                    shell_volume_m3,
                    rho_req_r,
                    e_req_r,
                    sun_years_r,
                    e_req_r_rear10,
                    sun_years_r_rear10
                ));
                compact_rows_json.push(json!({
                    "radius_m": r_bubble,
                    "thickness_m": thickness,
                    "v_eff_over_c": v_eff,
                    "shortcut_s": shortcut_s,
                    "curvature_scale_model": "radius",
                    "shell_volume_m3": shell_volume_m3,
                    "rho_req_j_m3": rho_req_r,
                    "e_req_j": e_req_r,
                    "sun_years": sun_years_r,
                    "e_req_rear10_j": e_req_r_rear10,
                    "sun_years_rear10": sun_years_r_rear10
                }));
                if compact_best_radius_model
                    .map(|b| e_req_r_rear10 < b.4)
                    .unwrap_or(true)
                {
                    compact_best_radius_model = Some((
                        r_bubble,
                        thickness,
                        v_eff,
                        shell_volume_m3,
                        e_req_r_rear10,
                        sun_years_r_rear10,
                        rho_req_r,
                    ));
                }

                // Model B (conservative): curvature set by wall thickness.
                let rho_req_t =
                    required_curvature_energy_density_j_m3_with_g(shortcut_s, thickness, G);
                let e_req_t = rho_req_t * shell_volume_m3;
                let e_req_t_rear10 = rear_face_factor * e_req_t;
                let sun_years_t = (e_req_t / l_sun_w) / seconds_per_year;
                let sun_years_t_rear10 = (e_req_t_rear10 / l_sun_w) / seconds_per_year;

                compact_rows_txt.push_str(&format!(
                    "{:.6e},{:.6e},{:.3},{:.6e},thickness,{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e}\n",
                    r_bubble,
                    thickness,
                    v_eff,
                    shortcut_s,
                    shell_volume_m3,
                    rho_req_t,
                    e_req_t,
                    sun_years_t,
                    e_req_t_rear10,
                    sun_years_t_rear10
                ));
                compact_rows_json.push(json!({
                    "radius_m": r_bubble,
                    "thickness_m": thickness,
                    "v_eff_over_c": v_eff,
                    "shortcut_s": shortcut_s,
                    "curvature_scale_model": "thickness",
                    "shell_volume_m3": shell_volume_m3,
                    "rho_req_j_m3": rho_req_t,
                    "e_req_j": e_req_t,
                    "sun_years": sun_years_t,
                    "e_req_rear10_j": e_req_t_rear10,
                    "sun_years_rear10": sun_years_t_rear10
                }));
                if compact_best_thickness_model
                    .map(|b| e_req_t_rear10 < b.4)
                    .unwrap_or(true)
                {
                    compact_best_thickness_model = Some((
                        r_bubble,
                        thickness,
                        v_eff,
                        shell_volume_m3,
                        e_req_t_rear10,
                        sun_years_t_rear10,
                        rho_req_t,
                    ));
                }

                // Model C (near-core exotic geometry): curvature pinned at lattice core scale.
                let rho_req_core =
                    required_curvature_energy_density_j_m3_with_g(shortcut_s, r_core_m, G);
                let e_req_core = rho_req_core * shell_volume_m3;
                let e_req_core_rear10 = rear_face_factor * e_req_core;
                let sun_years_core = (e_req_core / l_sun_w) / seconds_per_year;
                let sun_years_core_rear10 = (e_req_core_rear10 / l_sun_w) / seconds_per_year;

                compact_rows_txt.push_str(&format!(
                    "{:.6e},{:.6e},{:.3},{:.6e},core,{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e}\n",
                    r_bubble,
                    thickness,
                    v_eff,
                    shortcut_s,
                    shell_volume_m3,
                    rho_req_core,
                    e_req_core,
                    sun_years_core,
                    e_req_core_rear10,
                    sun_years_core_rear10
                ));
                compact_rows_json.push(json!({
                    "radius_m": r_bubble,
                    "thickness_m": thickness,
                    "v_eff_over_c": v_eff,
                    "shortcut_s": shortcut_s,
                    "curvature_scale_model": "core",
                    "curvature_scale_m": r_core_m,
                    "shell_volume_m3": shell_volume_m3,
                    "rho_req_j_m3": rho_req_core,
                    "e_req_j": e_req_core,
                    "sun_years": sun_years_core,
                    "e_req_rear10_j": e_req_core_rear10,
                    "sun_years_rear10": sun_years_core_rear10
                }));
                if compact_best_core_model
                    .map(|b| e_req_core_rear10 < b.4)
                    .unwrap_or(true)
                {
                    compact_best_core_model = Some((
                        r_bubble,
                        thickness,
                        v_eff,
                        shell_volume_m3,
                        e_req_core_rear10,
                        sun_years_core_rear10,
                        rho_req_core,
                    ));
                }
            }
        }
    }
    let compact_best_radius_model = compact_best_radius_model.expect("compact rows nonempty");
    let compact_best_thickness_model = compact_best_thickness_model.expect("compact rows nonempty");
    let compact_best_core_model = compact_best_core_model.expect("compact rows nonempty");

    let txt = format!(
        "[source]\n\
rho_visible_kg_m3={:.6e}\n\
rho_dark_kg_m3={:.6e}\n\
rho_vac_kg_m3={:.6e}\n\
rho_source_j_m3={:.6e}\n\
lambda_qg={:.12}\n\
g_eff_planck_over_g={:.12}\n\n\
[best_practical]\n\
v_eff_over_c={:.3}\n\
shortcut_s={:.6e}\n\
radius_m={:.6e}\n\
rho_required_classical_j_m3={:.6e}\n\
rho_required_gutoe_scale_radius_j_m3={:.6e}\n\
rho_required_gutoe_scale_planck_j_m3={:.6e}\n\
ratio_classical_over_source={:.6e}\n\
ratio_gutoe_scale_radius_over_source={:.6e}\n\
ratio_gutoe_scale_planck_over_source={:.6e}\n\
ratio_gutoe_scale_radius_rear10_over_source={:.6e}\n\
ratio_gutoe_scale_planck_rear10_over_source={:.6e}\n\
qg_gain_radius={:.12}\n\
qg_gain_planck={:.12}\n\
practical_feasible={}\n\
practical_feasible_rear10={}\n\
practical_feasible_planck_rear10={}\n\n\
[best_mega]\n\
v_eff_over_c={:.3}\n\
shortcut_s={:.6e}\n\
radius_m={:.6e}\n\
rho_required_classical_j_m3={:.6e}\n\
rho_required_gutoe_scale_radius_j_m3={:.6e}\n\
rho_required_gutoe_scale_planck_j_m3={:.6e}\n\
ratio_classical_over_source={:.6e}\n\
ratio_gutoe_scale_radius_over_source={:.6e}\n\
ratio_gutoe_scale_planck_over_source={:.6e}\n\
ratio_gutoe_scale_radius_rear10_over_source={:.6e}\n\
ratio_gutoe_scale_planck_rear10_over_source={:.6e}\n\
mega_feasible={}\n\
mega_feasible_rear10={}\n\
mega_feasible_planck_rear10={}\n\n\
[tube_geometry_zoom]\n\
tube_length_m={:.6e}\n\
tube_radius_m={:.6e}\n\
tube_volume_m3={:.6e}\n\
shortcut_s={:.6e}\n\
rho_req_longitudinal_classical_j_m3={:.6e}\n\
rho_req_transverse_classical_j_m3={:.6e}\n\
rho_req_longitudinal_rear10_j_m3={:.6e}\n\
rho_req_transverse_rear10_j_m3={:.6e}\n\
e_req_longitudinal_classical_j={:.6e}\n\
e_req_transverse_classical_j={:.6e}\n\
e_req_longitudinal_rear10_j={:.6e}\n\
e_req_transverse_rear10_j={:.6e}\n\
sun_years_longitudinal_classical={:.6e}\n\
sun_years_transverse_classical={:.6e}\n\
sun_years_longitudinal_rear10={:.6e}\n\
sun_years_transverse_rear10={:.6e}\n\n\
[compact_shell_zoom]\n\
r_core_m={:.6e}\n\
best_radius_model_rear10: radius_m={:.6e} thickness_m={:.6e} v_eff_over_c={:.3} shell_volume_m3={:.6e} rho_req_j_m3={:.6e} e_req_rear10_j={:.6e} sun_years_rear10={:.6e}\n\
best_thickness_model_rear10: radius_m={:.6e} thickness_m={:.6e} v_eff_over_c={:.3} shell_volume_m3={:.6e} rho_req_j_m3={:.6e} e_req_rear10_j={:.6e} sun_years_rear10={:.6e}\n\
best_core_model_rear10: radius_m={:.6e} thickness_m={:.6e} v_eff_over_c={:.3} shell_volume_m3={:.6e} rho_req_j_m3={:.6e} e_req_rear10_j={:.6e} sun_years_rear10={:.6e}\n\
{}\n",
        rho_visible_kg_m3,
        rho_dark_kg_m3,
        rho_vac_kg_m3,
        rho_source_j_m3,
        einstein.lambda_qg,
        g_eff_planck / G,
        best.v_eff_over_c,
        best.shortcut_s,
        best.radius_m,
        best.rho_required_classical_j_m3,
        best.rho_required_gutoe_scale_radius_j_m3,
        best.rho_required_gutoe_scale_planck_j_m3,
        best.ratio_classical_over_source,
        best.ratio_gutoe_scale_radius_over_source,
        best.ratio_gutoe_scale_planck_over_source,
        best.ratio_gutoe_scale_radius_rear10_over_source,
        best.ratio_gutoe_scale_planck_rear10_over_source,
        qg_radius_gain,
        qg_planck_gain,
        practical_feasible,
        practical_feasible_rear10,
        practical_feasible_planck_rear10,
        best_mega.v_eff_over_c,
        best_mega.shortcut_s,
        best_mega.radius_m,
        best_mega.rho_required_classical_j_m3,
        best_mega.rho_required_gutoe_scale_radius_j_m3,
        best_mega.rho_required_gutoe_scale_planck_j_m3,
        best_mega.ratio_classical_over_source,
        best_mega.ratio_gutoe_scale_radius_over_source,
        best_mega.ratio_gutoe_scale_planck_over_source,
        best_mega.ratio_gutoe_scale_radius_rear10_over_source,
        best_mega.ratio_gutoe_scale_planck_rear10_over_source,
        mega_feasible,
        mega_feasible_rear10,
        mega_feasible_planck_rear10,
        tube_length_m,
        tube_radius_m,
        tube_volume_m3,
        tube_shortcut_s,
        rho_req_long_classical,
        rho_req_trans_classical,
        rho_req_long_rear10,
        rho_req_trans_rear10,
        e_req_long_classical_j,
        e_req_trans_classical_j,
        e_req_long_rear10_j,
        e_req_trans_rear10_j,
        sun_years_long_classical,
        sun_years_trans_classical,
        sun_years_long_rear10,
        sun_years_trans_rear10,
        r_core_m,
        compact_best_radius_model.0,
        compact_best_radius_model.1,
        compact_best_radius_model.2,
        compact_best_radius_model.3,
        compact_best_radius_model.6,
        compact_best_radius_model.4,
        compact_best_radius_model.5,
        compact_best_thickness_model.0,
        compact_best_thickness_model.1,
        compact_best_thickness_model.2,
        compact_best_thickness_model.3,
        compact_best_thickness_model.6,
        compact_best_thickness_model.4,
        compact_best_thickness_model.5,
        compact_best_core_model.0,
        compact_best_core_model.1,
        compact_best_core_model.2,
        compact_best_core_model.3,
        compact_best_core_model.6,
        compact_best_core_model.4,
        compact_best_core_model.5,
        compact_rows_txt
    );
    fs::write(out.join("ftl_one_way_lentz_probe.txt"), txt).expect("write txt");

    let payload = json!({
        "assumptions": {
            "one_way_corridor": true,
            "local_signal_bound_enforced": true,
            "negative_energy_used": false,
            "ctc_blocked_by_topology": true,
            "rho_visible_kg_m3": rho_visible_kg_m3,
            "rear_face_factor_hypothesis": rear_face_factor
        },
        "source": {
            "rho_dark_kg_m3": rho_dark_kg_m3,
            "rho_vac_kg_m3": rho_vac_kg_m3,
            "rho_source_j_m3": rho_source_j_m3,
            "lambda_qg": einstein.lambda_qg,
            "g_eff_planck_over_g": g_eff_planck / G
        },
        "best_practical": {
            "v_eff_over_c": best.v_eff_over_c,
            "shortcut_s": best.shortcut_s,
            "radius_m": best.radius_m,
            "rho_required_classical_j_m3": best.rho_required_classical_j_m3,
            "rho_required_gutoe_scale_radius_j_m3": best.rho_required_gutoe_scale_radius_j_m3,
            "rho_required_gutoe_scale_planck_j_m3": best.rho_required_gutoe_scale_planck_j_m3,
            "ratio_classical_over_source": best.ratio_classical_over_source,
            "ratio_gutoe_scale_radius_over_source": best.ratio_gutoe_scale_radius_over_source,
            "ratio_gutoe_scale_planck_over_source": best.ratio_gutoe_scale_planck_over_source,
            "ratio_gutoe_scale_radius_rear10_over_source": best.ratio_gutoe_scale_radius_rear10_over_source,
            "ratio_gutoe_scale_planck_rear10_over_source": best.ratio_gutoe_scale_planck_rear10_over_source,
            "qg_gain_radius": qg_radius_gain,
            "qg_gain_planck": qg_planck_gain,
            "feasible": practical_feasible,
            "feasible_rear10": practical_feasible_rear10,
            "feasible_planck_rear10": practical_feasible_planck_rear10
        },
        "best_mega": {
            "v_eff_over_c": best_mega.v_eff_over_c,
            "shortcut_s": best_mega.shortcut_s,
            "radius_m": best_mega.radius_m,
            "rho_required_classical_j_m3": best_mega.rho_required_classical_j_m3,
            "rho_required_gutoe_scale_radius_j_m3": best_mega.rho_required_gutoe_scale_radius_j_m3,
            "rho_required_gutoe_scale_planck_j_m3": best_mega.rho_required_gutoe_scale_planck_j_m3,
            "ratio_classical_over_source": best_mega.ratio_classical_over_source,
            "ratio_gutoe_scale_radius_over_source": best_mega.ratio_gutoe_scale_radius_over_source,
            "ratio_gutoe_scale_planck_over_source": best_mega.ratio_gutoe_scale_planck_over_source,
            "ratio_gutoe_scale_radius_rear10_over_source": best_mega.ratio_gutoe_scale_radius_rear10_over_source,
            "ratio_gutoe_scale_planck_rear10_over_source": best_mega.ratio_gutoe_scale_planck_rear10_over_source,
            "feasible": mega_feasible,
            "feasible_rear10": mega_feasible_rear10,
            "feasible_planck_rear10": mega_feasible_planck_rear10
        },
        "tube_geometry_zoom": {
            "tube_length_m": tube_length_m,
            "tube_radius_m": tube_radius_m,
            "tube_volume_m3": tube_volume_m3,
            "shortcut_s": tube_shortcut_s,
            "rho_req_longitudinal_classical_j_m3": rho_req_long_classical,
            "rho_req_transverse_classical_j_m3": rho_req_trans_classical,
            "rho_req_longitudinal_rear10_j_m3": rho_req_long_rear10,
            "rho_req_transverse_rear10_j_m3": rho_req_trans_rear10,
            "e_req_longitudinal_classical_j": e_req_long_classical_j,
            "e_req_transverse_classical_j": e_req_trans_classical_j,
            "e_req_longitudinal_rear10_j": e_req_long_rear10_j,
            "e_req_transverse_rear10_j": e_req_trans_rear10_j,
            "sun_years_longitudinal_classical": sun_years_long_classical,
            "sun_years_transverse_classical": sun_years_trans_classical,
            "sun_years_longitudinal_rear10": sun_years_long_rear10,
            "sun_years_transverse_rear10": sun_years_trans_rear10
        },
        "compact_shell_zoom": {
            "r_core_m": r_core_m,
            "best_radius_model_rear10": {
                "radius_m": compact_best_radius_model.0,
                "thickness_m": compact_best_radius_model.1,
                "v_eff_over_c": compact_best_radius_model.2,
                "shell_volume_m3": compact_best_radius_model.3,
                "e_req_rear10_j": compact_best_radius_model.4,
                "sun_years_rear10": compact_best_radius_model.5,
                "rho_req_j_m3": compact_best_radius_model.6
            },
            "best_thickness_model_rear10": {
                "radius_m": compact_best_thickness_model.0,
                "thickness_m": compact_best_thickness_model.1,
                "v_eff_over_c": compact_best_thickness_model.2,
                "shell_volume_m3": compact_best_thickness_model.3,
                "e_req_rear10_j": compact_best_thickness_model.4,
                "sun_years_rear10": compact_best_thickness_model.5,
                "rho_req_j_m3": compact_best_thickness_model.6
            },
            "best_core_model_rear10": {
                "radius_m": compact_best_core_model.0,
                "thickness_m": compact_best_core_model.1,
                "v_eff_over_c": compact_best_core_model.2,
                "shell_volume_m3": compact_best_core_model.3,
                "e_req_rear10_j": compact_best_core_model.4,
                "sun_years_rear10": compact_best_core_model.5,
                "rho_req_j_m3": compact_best_core_model.6
            },
            "rows": compact_rows_json
        },
        "point_count": points.len()
    });
    fs::write(
        out.join("ftl_one_way_lentz_probe.json"),
        serde_json::to_string_pretty(&payload).expect("json encode"),
    )
    .expect("write json");

    println!(
        "best practical ratio (gutoe scale=R) required/source = {:.3e} (feasible={})",
        best.ratio_gutoe_scale_radius_over_source, practical_feasible
    );
    println!(
        "best practical ratio with rear10 hypothesis = {:.3e} (feasible={})",
        best.ratio_gutoe_scale_radius_rear10_over_source, practical_feasible_rear10
    );
    println!(
        "best mega ratio (gutoe scale=R) required/source = {:.3e} at R={:.3e} m (feasible={})",
        best_mega.ratio_gutoe_scale_radius_over_source, best_mega.radius_m, mega_feasible
    );
    println!(
        "best mega ratio with rear10 hypothesis = {:.3e} (feasible={})",
        best_mega.ratio_gutoe_scale_radius_rear10_over_source, mega_feasible_rear10
    );
    println!(
        "tube zoom (longitudinal curvature): E={:.3e} J ({:.3e} solar-years), rear10 E={:.3e} J ({:.3e} solar-years)",
        e_req_long_classical_j,
        sun_years_long_classical,
        e_req_long_rear10_j,
        sun_years_long_rear10
    );
    println!(
        "tube zoom (transverse curvature): E={:.3e} J ({:.3e} solar-years), rear10 E={:.3e} J ({:.3e} solar-years)",
        e_req_trans_classical_j,
        sun_years_trans_classical,
        e_req_trans_rear10_j,
        sun_years_trans_rear10
    );
    println!(
        "compact shell best (radius-model, rear10): E={:.3e} J ({:.3e} solar-years) at R={:.3e} m, t={:.3e} m, v_eff/c={:.1}",
        compact_best_radius_model.4,
        compact_best_radius_model.5,
        compact_best_radius_model.0,
        compact_best_radius_model.1,
        compact_best_radius_model.2
    );
    println!(
        "compact shell best (thickness-model, rear10): E={:.3e} J ({:.3e} solar-years) at R={:.3e} m, t={:.3e} m, v_eff/c={:.1}",
        compact_best_thickness_model.4,
        compact_best_thickness_model.5,
        compact_best_thickness_model.0,
        compact_best_thickness_model.1,
        compact_best_thickness_model.2
    );
    println!(
        "compact shell best (core-model, rear10): E={:.3e} J ({:.3e} solar-years) at R={:.3e} m, t={:.3e} m, v_eff/c={:.1}, r_core={:.3e} m",
        compact_best_core_model.4,
        compact_best_core_model.5,
        compact_best_core_model.0,
        compact_best_core_model.1,
        compact_best_core_model.2,
        r_core_m
    );
    println!(
        "wrote {}",
        out.join("ftl_one_way_lentz_probe.txt").display()
    );
    println!(
        "wrote {}",
        out.join("ftl_one_way_lentz_probe.json").display()
    );
}
