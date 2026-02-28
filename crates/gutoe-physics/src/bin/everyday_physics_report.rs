//! Everyday-physics multi-lane report.

use gutoe_physics::{
    evaluate_cat_purr_resonance, evaluate_coffee_flavor_shift, evaluate_default_bird_wing_efficiency,
    evaluate_rayleigh_scattering, evaluate_soap_bubble_optimum, CatPurrResonanceInput,
    CoffeeChemistryInput, RayleighModelInput, SoapBubbleInput,
};
use serde_json::json;
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

fn main() {
    let rayleigh = evaluate_rayleigh_scattering(RayleighModelInput::default());
    let soap = evaluate_soap_bubble_optimum(SoapBubbleInput::default());
    let cat = evaluate_cat_purr_resonance(CatPurrResonanceInput::default());
    let coffee_cfg = CoffeeChemistryInput::default();
    let coffee_rows = [0.0, 1_000.0, 2_000.0, 3_000.0, 4_000.0]
        .into_iter()
        .map(|h| evaluate_coffee_flavor_shift(h, coffee_cfg))
        .collect::<Vec<_>>();
    let wing_rows = evaluate_default_bird_wing_efficiency();
    let wing_winner = wing_rows
        .first()
        .expect("default_bird_wing_efficiency must produce at least one row");

    let out_dir = std::env::var("GUTOE_EVERYDAY_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/everyday_physics".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let txt_path = out.join("everyday_physics_report.txt");
    let json_path = out.join("everyday_physics_report.json");
    let coffee_csv_path = out.join("coffee_altitude_sweep.csv");
    let wing_csv_path = out.join("wing_efficiency_rank.csv");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[everyday_physics]").expect("write");
    writeln!(txt, "topic_count = 5").expect("write");
    writeln!(txt, "").expect("write");

    writeln!(txt, "[rayleigh]").expect("write");
    writeln!(
        txt,
        "blue_to_red_scattering_ratio = {:.6}",
        rayleigh.blue_to_red_scattering_ratio
    )
    .expect("write");
    writeln!(
        txt,
        "midday_blue_share_of_scattered_light = {:.6}",
        rayleigh.midday_blue_share_of_scattered_light
    )
    .expect("write");
    writeln!(
        txt,
        "sunset_red_to_blue_direct_ratio = {:.6}",
        rayleigh.sunset_red_to_blue_direct_ratio
    )
    .expect("write");
    writeln!(txt, "").expect("write");

    writeln!(txt, "[soap_bubble]").expect("write");
    writeln!(
        txt,
        "sphere_double_surface_energy_j = {:.9e}",
        soap.sphere_double_surface_energy_j
    )
    .expect("write");
    writeln!(
        txt,
        "cube_energy_penalty_percent = {:.6}",
        soap.cube_energy_penalty_percent
    )
    .expect("write");
    writeln!(
        txt,
        "prolate_energy_penalty_percent = {:.6}",
        soap.prolate_energy_penalty_percent
    )
    .expect("write");
    writeln!(txt, "").expect("write");

    writeln!(txt, "[cat_purr]").expect("write");
    writeln!(
        txt,
        "predicted_purr_frequency_hz = {:.6}",
        cat.predicted_purr_frequency_hz
    )
    .expect("write");
    writeln!(txt, "in_healing_band = {}", cat.in_healing_band).expect("write");
    writeln!(txt, "healing_overlap_score = {:.6}", cat.healing_overlap_score).expect("write");
    writeln!(txt, "").expect("write");

    writeln!(txt, "[coffee_altitude]").expect("write");
    for row in &coffee_rows {
        writeln!(
            txt,
            "alt={:>4.0}m boil={:.2}C bitter_rel={:.4} acid_rel={:.4} receptor_rel={:.4} acid_to_bitter={:.4}",
            row.altitude_m,
            row.boiling_temperature_c,
            row.bitter_intensity_relative,
            row.acidic_intensity_relative,
            row.receptor_affinity_relative,
            row.acidity_to_bitterness_ratio
        )
        .expect("write");
    }
    writeln!(txt, "").expect("write");

    writeln!(txt, "[wing_efficiency]").expect("write");
    writeln!(txt, "winner = {}", wing_winner.name).expect("write");
    writeln!(txt, "winner_ld_max = {:.6}", wing_winner.ld_max).expect("write");
    for row in &wing_rows {
        writeln!(
            txt,
            "{}: ar={:.4} cl_opt={:.4} ld_max={:.4}",
            row.name, row.aspect_ratio, row.cl_opt_for_ld_max, row.ld_max
        )
        .expect("write");
    }

    let mut coffee_csv = String::from(
        "altitude_m,pressure_pa,boiling_temperature_c,bitter_intensity_relative,acidic_intensity_relative,receptor_affinity_relative,acidity_to_bitterness_ratio\n",
    );
    for row in &coffee_rows {
        coffee_csv.push_str(&format!(
            "{:.0},{:.3},{:.6},{:.9},{:.9},{:.9},{:.9}\n",
            row.altitude_m,
            row.pressure_pa,
            row.boiling_temperature_c,
            row.bitter_intensity_relative,
            row.acidic_intensity_relative,
            row.receptor_affinity_relative,
            row.acidity_to_bitterness_ratio
        ));
    }
    fs::write(&coffee_csv_path, coffee_csv).expect("write coffee csv");

    let mut wing_csv =
        String::from("rank,name,aspect_ratio,cl_opt_for_ld_max,induced_drag_at_opt,ld_max\n");
    for (idx, row) in wing_rows.iter().enumerate() {
        wing_csv.push_str(&format!(
            "{},{},{:.6},{:.6},{:.6},{:.6}\n",
            idx + 1,
            row.name,
            row.aspect_ratio,
            row.cl_opt_for_ld_max,
            row.induced_drag_at_opt,
            row.ld_max
        ));
    }
    fs::write(&wing_csv_path, wing_csv).expect("write wing csv");

    let payload = json!({
        "meta": {
            "lane": "everyday_physics",
            "topics": [
                "rayleigh_sky_sunset",
                "soap_bubble_minimal_surface",
                "cat_purr_resonance_band",
                "coffee_altitude_flavor_chain",
                "bird_wing_lift_drag_efficiency"
            ],
            "note": "Reduced-order transduction report."
        },
        "rayleigh": {
            "blue_to_red_scattering_ratio": rayleigh.blue_to_red_scattering_ratio,
            "midday_blue_share_of_scattered_light": rayleigh.midday_blue_share_of_scattered_light,
            "sunset_red_to_blue_direct_ratio": rayleigh.sunset_red_to_blue_direct_ratio,
            "blue_cross_section_m2": rayleigh.blue_cross_section_m2,
            "red_cross_section_m2": rayleigh.red_cross_section_m2
        },
        "soap_bubble": {
            "sphere_area_m2": soap.sphere_area_m2,
            "cube_area_m2": soap.cube_area_m2,
            "prolate_area_m2": soap.prolate_area_m2,
            "sphere_double_surface_energy_j": soap.sphere_double_surface_energy_j,
            "cube_double_surface_energy_j": soap.cube_double_surface_energy_j,
            "prolate_double_surface_energy_j": soap.prolate_double_surface_energy_j,
            "cube_energy_penalty_percent": soap.cube_energy_penalty_percent,
            "prolate_energy_penalty_percent": soap.prolate_energy_penalty_percent
        },
        "cat_purr": {
            "predicted_purr_frequency_hz": cat.predicted_purr_frequency_hz,
            "in_healing_band": cat.in_healing_band,
            "distance_from_band_center_hz": cat.distance_from_band_center_hz,
            "healing_overlap_score": cat.healing_overlap_score,
            "effective_stiffness_n_per_m": cat.effective_stiffness_n_per_m
        },
        "coffee_altitude": coffee_rows.iter().map(|r| json!({
            "altitude_m": r.altitude_m,
            "pressure_pa": r.pressure_pa,
            "boiling_temperature_c": r.boiling_temperature_c,
            "bitter_intensity_relative": r.bitter_intensity_relative,
            "acidic_intensity_relative": r.acidic_intensity_relative,
            "receptor_affinity_relative": r.receptor_affinity_relative,
            "acidity_to_bitterness_ratio": r.acidity_to_bitterness_ratio
        })).collect::<Vec<_>>(),
        "wing_efficiency": {
            "winner": wing_winner.name,
            "winner_ld_max": wing_winner.ld_max,
            "ranked": wing_rows.iter().map(|r| json!({
                "name": r.name,
                "aspect_ratio": r.aspect_ratio,
                "cl_opt_for_ld_max": r.cl_opt_for_ld_max,
                "ld_max": r.ld_max
            })).collect::<Vec<_>>()
        }
    });
    fs::write(
        &json_path,
        serde_json::to_string_pretty(&payload).expect("serialize"),
    )
    .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!("wrote {}", coffee_csv_path.display());
    println!("wrote {}", wing_csv_path.display());
    println!(
        "everyday_physics: rayleigh_ratio={:.3} cat_purr={:.2}Hz wing_winner={} ld_max={:.2}",
        rayleigh.blue_to_red_scattering_ratio,
        cat.predicted_purr_frequency_hz,
        wing_winner.name,
        wing_winner.ld_max
    );
}
