//! Everyday extremes report lane.

use gutoe_physics::{
    default_ice_temperature_sweep, evaluate_mpemba, evaluate_mpemba_small_sweep,
    evaluate_popcorn_popping, evaluate_raindrop_shape_sweep, MpembaInput, PopcornInput,
};
use serde_json::json;
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

fn main() {
    let ice_rows = default_ice_temperature_sweep();
    let popcorn = evaluate_popcorn_popping(PopcornInput::default());
    let raindrop = evaluate_raindrop_shape_sweep();
    let mpemba = evaluate_mpemba(MpembaInput::default());
    let mpemba_sweep = evaluate_mpemba_small_sweep();

    let out_dir = std::env::var("GUTOE_EVERYDAY_EXTREMES_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/everyday_extremes".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let txt_path = out.join("everyday_extremes_report.txt");
    let json_path = out.join("everyday_extremes_report.json");
    let ice_csv_path = out.join("ice_slipperiness_sweep.csv");
    let rain_csv_path = out.join("raindrop_shape_sweep.csv");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[everyday_extremes]").expect("write");
    writeln!(txt, "topic_count = 4").expect("write");
    writeln!(txt, "").expect("write");

    writeln!(txt, "[ice_slipperiness]").expect("write");
    for (temp_c, row) in &ice_rows {
        writeln!(
            txt,
            "T={:>5.1}C qll={:.3}nm mu={:.4} effT={:.3}C (dT_pressure={:.3}, dT_friction={:.3})",
            temp_c,
            row.quasi_liquid_layer_thickness_nm,
            row.friction_coefficient,
            row.effective_surface_temperature_c,
            row.pressure_melting_shift_c,
            row.frictional_heating_shift_c
        )
        .expect("write");
    }
    writeln!(txt, "").expect("write");

    writeln!(txt, "[popcorn]").expect("write");
    writeln!(txt, "pops = {}", popcorn.pops).expect("write");
    writeln!(txt, "ready_temperature_c = {:.3}", popcorn.ready_temperature_c).expect("write");
    writeln!(txt, "burst_temperature_c = {:.3}", popcorn.burst_temperature_c).expect("write");
    writeln!(txt, "hysteresis_delta_c = {:.3}", popcorn.hysteresis_delta_c).expect("write");
    writeln!(txt, "ready_time_s = {:.3}", popcorn.ready_time_s).expect("write");
    writeln!(txt, "burst_time_s = {:.3}", popcorn.burst_time_s).expect("write");
    writeln!(
        txt,
        "internal_pressure_ready_mpa = {:.6}",
        popcorn.internal_pressure_ready_mpa
    )
    .expect("write");
    writeln!(
        txt,
        "internal_pressure_burst_mpa = {:.6}",
        popcorn.internal_pressure_burst_mpa
    )
    .expect("write");
    writeln!(
        txt,
        "rupture_threshold_ready_mpa = {:.6}",
        popcorn.rupture_threshold_ready_mpa
    )
    .expect("write");
    writeln!(
        txt,
        "rupture_threshold_burst_mpa = {:.6}",
        popcorn.rupture_threshold_burst_mpa
    )
    .expect("write");
    writeln!(
        txt,
        "pressure_margin_ready_mpa = {:.6}",
        popcorn.pressure_margin_ready_mpa
    )
    .expect("write");
    writeln!(
        txt,
        "pressure_margin_burst_mpa = {:.6}",
        popcorn.pressure_margin_burst_mpa
    )
    .expect("write");
    writeln!(
        txt,
        "estimated_expansion_ratio = {:.3}",
        popcorn.estimated_expansion_ratio
    )
    .expect("write");
    writeln!(txt, "").expect("write");

    writeln!(txt, "[raindrop_shape]").expect("write");
    writeln!(
        txt,
        "optimal_diameter_mm = {:.3}",
        raindrop.optimal.diameter_mm
    )
    .expect("write");
    writeln!(
        txt,
        "optimal_aspect_ratio = {:.4}",
        raindrop.optimal.aspect_ratio
    )
    .expect("write");
    writeln!(
        txt,
        "optimal_weber_number = {:.4}",
        raindrop.optimal.weber_number
    )
    .expect("write");
    writeln!(
        txt,
        "optimal_terminal_velocity_m_s = {:.4}",
        raindrop.optimal.terminal_velocity_m_s
    )
    .expect("write");
    writeln!(txt, "optimal_stable = {}", raindrop.optimal.stable).expect("write");
    writeln!(txt, "").expect("write");

    writeln!(txt, "[mpemba]").expect("write");
    writeln!(
        txt,
        "hot_total_freeze_time_min = {:.3}",
        mpemba.hot_total_freeze_time_s / 60.0
    )
    .expect("write");
    writeln!(
        txt,
        "cold_total_freeze_time_min = {:.3}",
        mpemba.cold_total_freeze_time_s / 60.0
    )
    .expect("write");
    writeln!(txt, "hot_faster = {}", mpemba.hot_faster).expect("write");
    writeln!(
        txt,
        "time_advantage_minutes = {:.3}",
        mpemba.time_advantage_minutes
    )
    .expect("write");
    writeln!(
        txt,
        "small_sweep_hot_faster_fraction = {:.4} ({}/{})",
        mpemba_sweep.hot_faster_fraction,
        mpemba_sweep.hot_faster_count,
        mpemba_sweep.sample_count
    )
    .expect("write");

    let mut ice_csv = String::from(
        "ice_temperature_c,effective_surface_temperature_c,quasi_liquid_layer_thickness_nm,friction_coefficient,pressure_melting_shift_c,frictional_heating_shift_c\n",
    );
    for (temp_c, row) in &ice_rows {
        ice_csv.push_str(&format!(
            "{:.3},{:.6},{:.6},{:.9},{:.6},{:.6}\n",
            temp_c,
            row.effective_surface_temperature_c,
            row.quasi_liquid_layer_thickness_nm,
            row.friction_coefficient,
            row.pressure_melting_shift_c,
            row.frictional_heating_shift_c
        ));
    }
    fs::write(&ice_csv_path, ice_csv).expect("write ice csv");

    let mut rain_csv = String::from(
        "diameter_mm,aspect_ratio,bond_number,weber_number,terminal_velocity_m_s,drag_coefficient,transport_score,stable\n",
    );
    for p in &raindrop.points {
        rain_csv.push_str(&format!(
            "{:.3},{:.6},{:.6},{:.6},{:.6},{:.6},{:.9},{}\n",
            p.diameter_mm,
            p.aspect_ratio,
            p.bond_number,
            p.weber_number,
            p.terminal_velocity_m_s,
            p.drag_coefficient,
            p.transport_score,
            p.stable
        ));
    }
    fs::write(&rain_csv_path, rain_csv).expect("write rain csv");

    let payload = json!({
        "meta": {
            "lane": "everyday_extremes",
            "topics": [
                "ice_slipperiness",
                "popcorn_popping",
                "raindrop_shape_optimum",
                "mpemba_regime"
            ],
            "note": "Reduced-order structural transduction lane."
        },
        "ice": ice_rows.iter().map(|(t, r)| json!({
            "ice_temperature_c": t,
            "effective_surface_temperature_c": r.effective_surface_temperature_c,
            "quasi_liquid_layer_thickness_nm": r.quasi_liquid_layer_thickness_nm,
            "friction_coefficient": r.friction_coefficient,
            "pressure_melting_shift_c": r.pressure_melting_shift_c,
            "frictional_heating_shift_c": r.frictional_heating_shift_c
        })).collect::<Vec<_>>(),
        "popcorn": {
            "pops": popcorn.pops,
            "ready_temperature_c": popcorn.ready_temperature_c,
            "burst_temperature_c": popcorn.burst_temperature_c,
            "hysteresis_delta_c": popcorn.hysteresis_delta_c,
            "ready_time_s": popcorn.ready_time_s,
            "burst_time_s": popcorn.burst_time_s,
            "internal_pressure_ready_mpa": popcorn.internal_pressure_ready_mpa,
            "internal_pressure_burst_mpa": popcorn.internal_pressure_burst_mpa,
            "rupture_threshold_ready_mpa": popcorn.rupture_threshold_ready_mpa,
            "rupture_threshold_burst_mpa": popcorn.rupture_threshold_burst_mpa,
            "pressure_margin_ready_mpa": popcorn.pressure_margin_ready_mpa,
            "pressure_margin_burst_mpa": popcorn.pressure_margin_burst_mpa,
            "estimated_expansion_ratio": popcorn.estimated_expansion_ratio
        },
        "raindrop_optimal": {
            "diameter_mm": raindrop.optimal.diameter_mm,
            "aspect_ratio": raindrop.optimal.aspect_ratio,
            "bond_number": raindrop.optimal.bond_number,
            "weber_number": raindrop.optimal.weber_number,
            "terminal_velocity_m_s": raindrop.optimal.terminal_velocity_m_s,
            "stable": raindrop.optimal.stable
        },
        "mpemba": {
            "hot_total_freeze_time_s": mpemba.hot_total_freeze_time_s,
            "cold_total_freeze_time_s": mpemba.cold_total_freeze_time_s,
            "hot_faster": mpemba.hot_faster,
            "time_advantage_minutes": mpemba.time_advantage_minutes,
            "small_sweep_hot_faster_fraction": mpemba_sweep.hot_faster_fraction,
            "small_sweep_hot_faster_count": mpemba_sweep.hot_faster_count,
            "small_sweep_sample_count": mpemba_sweep.sample_count
        }
    });
    fs::write(
        &json_path,
        serde_json::to_string_pretty(&payload).expect("serialize"),
    )
    .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!("wrote {}", ice_csv_path.display());
    println!("wrote {}", rain_csv_path.display());
    println!(
        "everyday_extremes: popcorn_ready={:.1}C burst={:.1}C raindrop_opt={:.2}mm mpemba_hot_faster={} sweep_frac={:.3}",
        popcorn.ready_temperature_c,
        popcorn.burst_temperature_c,
        raindrop.optimal.diameter_mm,
        mpemba.hot_faster,
        mpemba_sweep.hot_faster_fraction
    );
}
