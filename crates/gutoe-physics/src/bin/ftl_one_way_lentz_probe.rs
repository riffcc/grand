//! One-way geometry-FTL loophole probe.
//!
//! Hypothesis under test:
//! - keep local propagation subluminal (`v_local <= c`)
//! - avoid CTC by forcing one-way corridor topology
//! - avoid negative energy by sourcing only positive-energy dark-sector budget
//! - ask whether useful FTL (v_eff > c) is energetically feasible

use gutoe_physics::constants::{C, G};
use gutoe_physics::dark_sector::{dark_density_particle, vacuum_energy_density_structural};
use serde_json::json;
use std::fs;
use std::path::PathBuf;

#[derive(Debug, Clone, Copy)]
struct ProbePoint {
    v_eff_over_c: f64,
    radius_m: f64,
    shortcut_s: f64,
    rho_required_j_m3: f64,
    rho_source_j_m3: f64,
    density_ratio_required_over_source: f64,
}

fn required_curvature_energy_density_j_m3(shortcut_s: f64, radius_m: f64) -> f64 {
    // Einstein scale: ρ ~ (c^4 / 8πG) * K, with K ~ ((1/s)-1)/R² for a corridor strain proxy.
    let pref = C.powi(4) / (8.0 * std::f64::consts::PI * G);
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

    let mut points = Vec::new();
    let mut best: Option<ProbePoint> = None;
    for &v_eff_over_c in &v_grid {
        let shortcut_s = 1.0 / v_eff_over_c; // coordinate speed = c/s
        for &radius_m in &r_grid {
            let rho_required_j_m3 = required_curvature_energy_density_j_m3(shortcut_s, radius_m);
            let density_ratio_required_over_source = rho_required_j_m3 / rho_source_j_m3;
            let p = ProbePoint {
                v_eff_over_c,
                radius_m,
                shortcut_s,
                rho_required_j_m3,
                rho_source_j_m3,
                density_ratio_required_over_source,
            };
            if best
                .map(|b| p.density_ratio_required_over_source < b.density_ratio_required_over_source)
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
        let rho_required_j_m3 = required_curvature_energy_density_j_m3(shortcut_s, radius_m);
        let density_ratio_required_over_source = rho_required_j_m3 / rho_source_j_m3;
        let p = ProbePoint {
            v_eff_over_c,
            radius_m,
            shortcut_s,
            rho_required_j_m3,
            rho_source_j_m3,
            density_ratio_required_over_source,
        };
        if best_mega
            .map(|b| p.density_ratio_required_over_source < b.density_ratio_required_over_source)
            .unwrap_or(true)
        {
            best_mega = Some(p);
        }
    }

    let best = best.expect("nonempty probe");
    let best_mega = best_mega.expect("nonempty mega probe");
    let practical_feasible = best.density_ratio_required_over_source <= 1.0;
    let mega_feasible = best_mega.density_ratio_required_over_source <= 1.0;

    let txt = format!(
        "[source]\n\
rho_visible_kg_m3={:.6e}\n\
rho_dark_kg_m3={:.6e}\n\
rho_vac_kg_m3={:.6e}\n\
rho_source_j_m3={:.6e}\n\n\
[best_practical]\n\
v_eff_over_c={:.3}\n\
shortcut_s={:.6e}\n\
radius_m={:.6e}\n\
rho_required_j_m3={:.6e}\n\
ratio_required_over_source={:.6e}\n\
practical_feasible={}\n\n\
[best_mega]\n\
v_eff_over_c={:.3}\n\
shortcut_s={:.6e}\n\
radius_m={:.6e}\n\
rho_required_j_m3={:.6e}\n\
ratio_required_over_source={:.6e}\n\
mega_feasible={}\n",
        rho_visible_kg_m3,
        rho_dark_kg_m3,
        rho_vac_kg_m3,
        rho_source_j_m3,
        best.v_eff_over_c,
        best.shortcut_s,
        best.radius_m,
        best.rho_required_j_m3,
        best.density_ratio_required_over_source,
        practical_feasible,
        best_mega.v_eff_over_c,
        best_mega.shortcut_s,
        best_mega.radius_m,
        best_mega.rho_required_j_m3,
        best_mega.density_ratio_required_over_source,
        mega_feasible
    );
    fs::write(out.join("ftl_one_way_lentz_probe.txt"), txt).expect("write txt");

    let payload = json!({
        "assumptions": {
            "one_way_corridor": true,
            "local_signal_bound_enforced": true,
            "negative_energy_used": false,
            "ctc_blocked_by_topology": true,
            "rho_visible_kg_m3": rho_visible_kg_m3
        },
        "source": {
            "rho_dark_kg_m3": rho_dark_kg_m3,
            "rho_vac_kg_m3": rho_vac_kg_m3,
            "rho_source_j_m3": rho_source_j_m3
        },
        "best_practical": {
            "v_eff_over_c": best.v_eff_over_c,
            "shortcut_s": best.shortcut_s,
            "radius_m": best.radius_m,
            "rho_required_j_m3": best.rho_required_j_m3,
            "ratio_required_over_source": best.density_ratio_required_over_source,
            "feasible": practical_feasible
        },
        "best_mega": {
            "v_eff_over_c": best_mega.v_eff_over_c,
            "shortcut_s": best_mega.shortcut_s,
            "radius_m": best_mega.radius_m,
            "rho_required_j_m3": best_mega.rho_required_j_m3,
            "ratio_required_over_source": best_mega.density_ratio_required_over_source,
            "feasible": mega_feasible
        },
        "point_count": points.len()
    });
    fs::write(
        out.join("ftl_one_way_lentz_probe.json"),
        serde_json::to_string_pretty(&payload).expect("json encode"),
    )
    .expect("write json");

    println!(
        "best practical ratio required/source = {:.3e} (feasible={})",
        best.density_ratio_required_over_source, practical_feasible
    );
    println!(
        "best mega ratio required/source = {:.3e} at R={:.3e} m (feasible={})",
        best_mega.density_ratio_required_over_source, best_mega.radius_m, mega_feasible
    );
    println!("wrote {}", out.join("ftl_one_way_lentz_probe.txt").display());
    println!("wrote {}", out.join("ftl_one_way_lentz_probe.json").display());
}

