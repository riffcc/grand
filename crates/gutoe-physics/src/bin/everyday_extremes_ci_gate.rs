//! CI gate for everyday extremes lane.

use gutoe_physics::{
    evaluate_ice_slipperiness, evaluate_mpemba, evaluate_mpemba_small_sweep,
    evaluate_popcorn_popping, evaluate_raindrop_shape_sweep, IceSlipperinessInput, MpembaInput,
    PopcornInput,
};
use serde_json::json;
use std::fs::{self, File};
use std::io::Write;
use std::process;

fn env_f64(key: &str, default: f64) -> f64 {
    std::env::var(key)
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(default)
}

fn main() {
    let warm_ice = evaluate_ice_slipperiness(IceSlipperinessInput {
        ice_temperature_c: -2.0,
        ..IceSlipperinessInput::default()
    });
    let cold_ice = evaluate_ice_slipperiness(IceSlipperinessInput {
        ice_temperature_c: -20.0,
        ..IceSlipperinessInput::default()
    });

    let popcorn = evaluate_popcorn_popping(PopcornInput::default());
    let raindrop = evaluate_raindrop_shape_sweep();
    let mpemba_default = evaluate_mpemba(MpembaInput::default());
    let mpemba_control = evaluate_mpemba(MpembaInput {
        evap_fraction_hot: 0.0,
        evap_fraction_cold: 0.0,
        convection_boost_hot: 0.0,
        convection_boost_cold: 0.0,
        supercool_hot_c: 0.0,
        supercool_cold_c: 0.0,
        freezing_flux_boost_hot: 0.0,
        freezing_flux_boost_cold: 0.0,
        ..MpembaInput::default()
    });
    let mpemba_sweep = evaluate_mpemba_small_sweep();

    let min_ice_mu_drop = env_f64("GUTOE_EXTREMES_MIN_ICE_MU_DROP", 0.10);
    let min_pop_ready_temp_c = env_f64("GUTOE_EXTREMES_MIN_POP_READY_TEMP_C", 145.0);
    let max_pop_ready_temp_c = env_f64("GUTOE_EXTREMES_MAX_POP_READY_TEMP_C", 175.0);
    let min_pop_burst_temp_c = env_f64("GUTOE_EXTREMES_MIN_POP_BURST_TEMP_C", 170.0);
    let max_pop_burst_temp_c = env_f64("GUTOE_EXTREMES_MAX_POP_BURST_TEMP_C", 195.0);
    let min_pop_hysteresis_c = env_f64("GUTOE_EXTREMES_MIN_POP_HYSTERESIS_C", 5.0);
    let min_expansion_ratio = env_f64("GUTOE_EXTREMES_MIN_POP_EXPANSION", 10.0);
    let min_rain_opt_d_mm = env_f64("GUTOE_EXTREMES_MIN_RAIN_OPT_D_MM", 2.5);
    let max_rain_opt_d_mm = env_f64("GUTOE_EXTREMES_MAX_RAIN_OPT_D_MM", 6.0);
    let min_mpemba_fraction = env_f64("GUTOE_EXTREMES_MIN_MPEMBA_SWEEP_FRACTION", 0.2);

    let ice_mu_drop = cold_ice.friction_coefficient - warm_ice.friction_coefficient;
    let ice_ok = ice_mu_drop >= min_ice_mu_drop;

    let popcorn_ok = popcorn.pops
        && popcorn.ready_temperature_c >= min_pop_ready_temp_c
        && popcorn.ready_temperature_c <= max_pop_ready_temp_c
        && popcorn.burst_temperature_c >= min_pop_burst_temp_c
        && popcorn.burst_temperature_c <= max_pop_burst_temp_c
        && popcorn.hysteresis_delta_c >= min_pop_hysteresis_c
        && popcorn.burst_temperature_c > popcorn.ready_temperature_c
        && popcorn.estimated_expansion_ratio >= min_expansion_ratio;

    let rain_ok = raindrop.optimal.stable
        && raindrop.optimal.aspect_ratio < 1.0
        && raindrop.optimal.diameter_mm >= min_rain_opt_d_mm
        && raindrop.optimal.diameter_mm <= max_rain_opt_d_mm;

    let mpemba_regime_ok = mpemba_default.hot_faster
        && !mpemba_control.hot_faster
        && mpemba_sweep.hot_faster_fraction >= min_mpemba_fraction;

    let overall_pass = ice_ok && popcorn_ok && rain_ok && mpemba_regime_ok;

    let out_dir =
        std::env::var("GUTOE_EVERYDAY_EXTREMES_GATE_OUT").unwrap_or_else(|_| "/tmp/bh_renders".to_string());
    let _ = fs::create_dir_all(&out_dir);
    let json_path = format!("{out_dir}/everyday_extremes_ci_gate.json");

    let payload = json!({
        "overall_pass": overall_pass,
        "windows": {
            "min_ice_mu_drop": min_ice_mu_drop,
            "min_pop_ready_temp_c": min_pop_ready_temp_c,
            "max_pop_ready_temp_c": max_pop_ready_temp_c,
            "min_pop_burst_temp_c": min_pop_burst_temp_c,
            "max_pop_burst_temp_c": max_pop_burst_temp_c,
            "min_pop_hysteresis_c": min_pop_hysteresis_c,
            "min_expansion_ratio": min_expansion_ratio,
            "min_rain_opt_d_mm": min_rain_opt_d_mm,
            "max_rain_opt_d_mm": max_rain_opt_d_mm,
            "min_mpemba_sweep_fraction": min_mpemba_fraction
        },
        "summary": {
            "ice_mu_cold_minus_warm": ice_mu_drop,
            "warm_mu": warm_ice.friction_coefficient,
            "cold_mu": cold_ice.friction_coefficient,
            "pop_ready_temperature_c": popcorn.ready_temperature_c,
            "pop_burst_temperature_c": popcorn.burst_temperature_c,
            "pop_hysteresis_delta_c": popcorn.hysteresis_delta_c,
            "pop_expansion_ratio": popcorn.estimated_expansion_ratio,
            "raindrop_opt_diameter_mm": raindrop.optimal.diameter_mm,
            "raindrop_opt_aspect_ratio": raindrop.optimal.aspect_ratio,
            "mpemba_default_hot_faster": mpemba_default.hot_faster,
            "mpemba_control_hot_faster": mpemba_control.hot_faster,
            "mpemba_sweep_hot_faster_fraction": mpemba_sweep.hot_faster_fraction
        },
        "gate": {
            "ice_ok": ice_ok,
            "popcorn_ok": popcorn_ok,
            "raindrop_ok": rain_ok,
            "mpemba_regime_ok": mpemba_regime_ok
        }
    });

    let mut json_file = File::create(&json_path).expect("create gate json");
    writeln!(
        json_file,
        "{}",
        serde_json::to_string_pretty(&payload).expect("serialize gate")
    )
    .expect("write gate");

    println!(
        "everyday_extremes_gate: pass={} ice_drop={:.3} pop_ready={:.1} pop_burst={:.1} rain_opt={:.2}mm mpemba_frac={:.3}",
        overall_pass,
        ice_mu_drop,
        popcorn.ready_temperature_c,
        popcorn.burst_temperature_c,
        raindrop.optimal.diameter_mm,
        mpemba_sweep.hot_faster_fraction
    );
    println!("wrote {json_path}");

    if !overall_pass {
        eprintln!(
            "FAIL: ice_ok={} popcorn_ok={} raindrop_ok={} mpemba_regime_ok={}",
            ice_ok, popcorn_ok, rain_ok, mpemba_regime_ok
        );
        process::exit(2);
    }
}
