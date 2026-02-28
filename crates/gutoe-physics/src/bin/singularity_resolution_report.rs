//! GRAND-128: singularity-resolution report (black-hole core + Big-Bang bounce).

use gutoe_physics::{
    bounce_kernel, bounce_scale_factor, bounce_volume_fraction, hubble_sq_bounce_si,
    kretschmann_classical_m4, kretschmann_regularized_m4, lattice_core_radius_m,
    lattice_critical_density_kg_m3, planck_density_kg_m3, regularized_g_tt, schwarzschild_radius_m,
    C_INF, PLANCK_LENGTH,
};
use serde_json::json;
use std::fs;
use std::path::PathBuf;

const SOLAR_MASS_KG: f64 = 1.988_47e30;

fn env_f64(name: &str, default: f64) -> f64 {
    std::env::var(name)
        .ok()
        .and_then(|v| v.parse::<f64>().ok())
        .unwrap_or(default)
}

fn main() {
    let out_dir = std::env::var("GUTOE_SINGULARITY_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/singularity_resolution".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let bh_mass_solar = env_f64("GUTOE_BH_MASS_SOLAR", 10.0).max(1.0e-9);
    let bh_mass_kg = bh_mass_solar * SOLAR_MASS_KG;

    let l_p = PLANCK_LENGTH;
    let r_core = lattice_core_radius_m(l_p);
    let r_s = schwarzschild_radius_m(bh_mass_kg);
    let g_tt_origin = regularized_g_tt(0.0, r_s, l_p);
    let k_origin_reg = kretschmann_regularized_m4(0.0, r_s, l_p);

    let r_probe = (1.0e-6 * r_core).max(1.0e-300);
    let k_probe_classical = kretschmann_classical_m4(r_probe, r_s).unwrap_or(f64::INFINITY);
    let k_probe_reg = kretschmann_regularized_m4(r_probe, r_s, l_p);
    let k_suppression_ratio = if k_probe_reg > 0.0 {
        k_probe_classical / k_probe_reg
    } else {
        f64::INFINITY
    };

    let rho_p = planck_density_kg_m3();
    let rho_crit = lattice_critical_density_kg_m3();
    let h2_0 = hubble_sq_bounce_si(0.0, rho_crit);
    let h2_half = hubble_sq_bounce_si(0.5 * rho_crit, rho_crit);
    let h2_crit = hubble_sq_bounce_si(rho_crit, rho_crit);
    let h2_over = hubble_sq_bounce_si(1.1 * rho_crit, rho_crit);

    let w = 1.0 / 3.0;
    let rho_a1 = env_f64("GUTOE_BOUNCE_RHO_A1_FRAC", 1.0e-120) * rho_crit;
    let a_bounce = bounce_scale_factor(rho_a1, w, rho_crit).unwrap_or(f64::NAN);
    let v_bounce = bounce_volume_fraction(rho_a1, w, rho_crit).unwrap_or(f64::NAN);

    let bh_origin_finite =
        g_tt_origin.is_finite() && k_origin_reg.is_finite() && k_origin_reg > 0.0;
    let bh_regularization_active = k_suppression_ratio > 1.0;
    let bounce_kernel_anchors_ok = bounce_kernel(0.0, rho_crit) == 0.0
        && bounce_kernel(rho_crit, rho_crit).abs() < 1.0e-12 * rho_crit;
    let bounce_turnaround_ok = h2_half > 0.0 && h2_crit.abs() < 1.0e-18 && h2_over.abs() < 1.0e-18;
    let bounce_scale_ok = a_bounce.is_finite() && a_bounce > 0.0 && a_bounce < 1.0;

    let passes_all = bh_origin_finite
        && bh_regularization_active
        && bounce_kernel_anchors_ok
        && bounce_turnaround_ok
        && bounce_scale_ok;

    let mut bounce_csv = String::from("rho_over_rho_crit,kernel_kg_m3,h2_s2\n");
    for i in 0..=12 {
        let frac = i as f64 / 10.0;
        let rho = frac * rho_crit;
        bounce_csv.push_str(&format!(
            "{:.6},{:.12e},{:.12e}\n",
            frac,
            bounce_kernel(rho, rho_crit),
            hubble_sq_bounce_si(rho, rho_crit)
        ));
    }

    let mut bh_csv =
        String::from("r_over_r_core,r_m,kretschmann_classical_m4,kretschmann_regularized_m4\n");
    for frac in [0.0, 1.0e-9, 1.0e-6, 1.0e-3, 1.0, 10.0] {
        let r = frac * r_core;
        let k_class = kretschmann_classical_m4(r, r_s).unwrap_or(f64::INFINITY);
        let k_reg = kretschmann_regularized_m4(r, r_s, l_p);
        bh_csv.push_str(&format!(
            "{:.12e},{:.12e},{:.12e},{:.12e}\n",
            frac, r, k_class, k_reg
        ));
    }

    let txt = format!(
        "[meta]\nlane = GRAND-128_singularity_resolution\nC_inf = {:.12e}\n\n[black_hole]\nmass_solar = {:.6}\nmass_kg = {:.12e}\nplanck_length_m = {:.12e}\nr_core_m = {:.12e}\nr_s_m = {:.12e}\ng_tt_origin_regularized = {:.12e}\nK_origin_regularized_m^-4 = {:.12e}\nK_probe_classical_m^-4 = {:.12e}\nK_probe_regularized_m^-4 = {:.12e}\nK_probe_classical_over_regularized = {:.12e}\n\n[big_bang_bounce]\nrho_planck_kg_m^-3 = {:.12e}\nrho_crit_kg_m^-3 = {:.12e}\nrho_crit_over_planck = {:.12e}\nH2_rho0_s^-2 = {:.12e}\nH2_rhohalf_s^-2 = {:.12e}\nH2_rhocrit_s^-2 = {:.12e}\nH2_rho1p1crit_s^-2 = {:.12e}\nrho_a1_kg_m^-3 = {:.12e}\nw = {:.12e}\na_bounce = {:.12e}\nvolume_fraction_bounce = {:.12e}\n\n[gate]\nbh_origin_finite = {}\nbh_regularization_active = {}\nbounce_kernel_anchors_ok = {}\nbounce_turnaround_ok = {}\nbounce_scale_ok = {}\npasses_all = {}\n",
        C_INF,
        bh_mass_solar,
        bh_mass_kg,
        l_p,
        r_core,
        r_s,
        g_tt_origin,
        k_origin_reg,
        k_probe_classical,
        k_probe_reg,
        k_suppression_ratio,
        rho_p,
        rho_crit,
        rho_crit / rho_p,
        h2_0,
        h2_half,
        h2_crit,
        h2_over,
        rho_a1,
        w,
        a_bounce,
        v_bounce,
        bh_origin_finite,
        bh_regularization_active,
        bounce_kernel_anchors_ok,
        bounce_turnaround_ok,
        bounce_scale_ok,
        passes_all
    );

    let payload = json!({
        "meta": {
            "lane": "GRAND-128_singularity_resolution",
            "c_inf": C_INF,
        },
        "black_hole": {
            "mass_solar": bh_mass_solar,
            "mass_kg": bh_mass_kg,
            "planck_length_m": l_p,
            "r_core_m": r_core,
            "r_s_m": r_s,
            "g_tt_origin_regularized": g_tt_origin,
            "kretschmann_origin_regularized_m4": k_origin_reg,
            "probe": {
                "r_probe_m": r_probe,
                "kretschmann_classical_m4": k_probe_classical,
                "kretschmann_regularized_m4": k_probe_reg,
                "classical_over_regularized": k_suppression_ratio,
            }
        },
        "big_bang_bounce": {
            "rho_planck_kg_m3": rho_p,
            "rho_crit_kg_m3": rho_crit,
            "rho_crit_over_planck": rho_crit / rho_p,
            "h2": {
                "rho0": h2_0,
                "rho_half": h2_half,
                "rho_crit": h2_crit,
                "rho_1p1_crit": h2_over,
            },
            "reference_eos_w": w,
            "rho_a1_kg_m3": rho_a1,
            "a_bounce": a_bounce,
            "volume_fraction_bounce": v_bounce,
        },
        "gate": {
            "bh_origin_finite": bh_origin_finite,
            "bh_regularization_active": bh_regularization_active,
            "bounce_kernel_anchors_ok": bounce_kernel_anchors_ok,
            "bounce_turnaround_ok": bounce_turnaround_ok,
            "bounce_scale_ok": bounce_scale_ok,
            "passes_all": passes_all,
        }
    });

    let txt_path = out.join("singularity_resolution_report.txt");
    let json_path = out.join("singularity_resolution_report.json");
    let bounce_csv_path = out.join("singularity_resolution_bounce_curve.csv");
    let bh_csv_path = out.join("singularity_resolution_bh_profile.csv");

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(
        &json_path,
        serde_json::to_string_pretty(&payload).expect("json serialize"),
    )
    .expect("write json");
    fs::write(&bounce_csv_path, bounce_csv).expect("write bounce csv");
    fs::write(&bh_csv_path, bh_csv).expect("write bh csv");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!("wrote {}", bounce_csv_path.display());
    println!("wrote {}", bh_csv_path.display());
    println!(
        "singularity_resolution: pass={} (r_core={:.3e} m, rho_crit={:.3e} kg/m^3)",
        passes_all, r_core, rho_crit
    );
}
