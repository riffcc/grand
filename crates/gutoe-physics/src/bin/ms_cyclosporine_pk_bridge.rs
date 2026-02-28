//! Cyclosporine PK bridge: effect-site concentration -> blood concentration proxy.
//!
//! This is a reduced-order translational model with explicit uncertainty factors,
//! not clinical dosing guidance.

use gutoe_physics::{
    default_cyclosporine_pk_bridge_input, simulate_cyclosporine_pk_bridge,
    summarize_cyclosporine_pk_bridge,
};
use serde_json::json;
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

fn env_f64(key: &str, default: f64) -> f64 {
    std::env::var(key)
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(default)
}

fn env_u64(key: &str, default: u64) -> u64 {
    std::env::var(key)
        .ok()
        .and_then(|s| s.parse::<u64>().ok())
        .unwrap_or(default)
}

fn env_usize(key: &str, default: usize) -> usize {
    std::env::var(key)
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(default)
}

fn main() {
    let mut input = default_cyclosporine_pk_bridge_input();
    input.site_target_nanomolar = env_f64("GUTOE_MS_SITE_TARGET_NM", input.site_target_nanomolar);
    input.molecular_weight_g_mol = env_f64("GUTOE_CYCLOSPORINE_MW_G_MOL", input.molecular_weight_g_mol);
    input.blood_to_site_gain_median =
        env_f64("GUTOE_MS_PK_GAIN_MEDIAN", input.blood_to_site_gain_median);
    input.blood_to_site_gain_gsd = env_f64("GUTOE_MS_PK_GAIN_GSD", input.blood_to_site_gain_gsd);
    input.samples = env_usize("GUTOE_MS_PK_SAMPLES", input.samples);
    input.seed = env_u64("GUTOE_MS_PK_SEED", input.seed);

    let ensemble = simulate_cyclosporine_pk_bridge(input);
    let s = summarize_cyclosporine_pk_bridge(&ensemble);

    let out_dir = std::env::var("GUTOE_MS_CYCLOSPORINE_PK_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ms_cyclosporine_pk_bridge".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let txt_path = out.join("ms_cyclosporine_pk_bridge_report.txt");
    let json_path = out.join("ms_cyclosporine_pk_bridge_report.json");
    let csv_path = out.join("ms_cyclosporine_pk_bridge_quantiles.csv");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[ms_cyclosporine_pk_bridge]").expect("write");
    writeln!(txt, "site_target_nM = {:.9}", input.site_target_nanomolar).expect("write");
    writeln!(txt, "molecular_weight_g_mol = {:.9}", input.molecular_weight_g_mol).expect("write");
    writeln!(txt, "blood_to_site_gain_median = {:.9}", input.blood_to_site_gain_median).expect("write");
    writeln!(txt, "blood_to_site_gain_gsd = {:.9}", input.blood_to_site_gain_gsd).expect("write");
    writeln!(txt, "samples = {}", input.samples).expect("write");
    writeln!(txt, "seed = {}", input.seed).expect("write");
    writeln!(txt, "p05_nM = {:.9}", s.p05_nanomolar).expect("write");
    writeln!(txt, "p50_nM = {:.9}", s.p50_nanomolar).expect("write");
    writeln!(txt, "p95_nM = {:.9}", s.p95_nanomolar).expect("write");
    writeln!(txt, "p05_ng_mL = {:.9}", s.p05_ng_ml).expect("write");
    writeln!(txt, "p50_ng_mL = {:.9}", s.p50_ng_ml).expect("write");
    writeln!(txt, "p95_ng_mL = {:.9}", s.p95_ng_ml).expect("write");
    writeln!(txt, "mean_ng_mL = {:.9}", s.mean_ng_ml).expect("write");

    let csv = format!(
        "quantile,nanomolar,ng_mL\n0.05,{:.9},{:.9}\n0.25,{:.9},{:.9}\n0.50,{:.9},{:.9}\n0.75,{:.9},{:.9}\n0.95,{:.9},{:.9}\n",
        s.p05_nanomolar,
        s.p05_ng_ml,
        s.p25_nanomolar,
        s.p25_ng_ml,
        s.p50_nanomolar,
        s.p50_ng_ml,
        s.p75_nanomolar,
        s.p75_ng_ml,
        s.p95_nanomolar,
        s.p95_ng_ml,
    );
    fs::write(&csv_path, csv).expect("write csv");

    let payload = json!({
        "meta": {
            "lane": "ms_cyclosporine_pk_bridge",
            "note": "reduced-order translational exposure bridge, not dosing guidance"
        },
        "input": {
            "site_target_nM": input.site_target_nanomolar,
            "molecular_weight_g_mol": input.molecular_weight_g_mol,
            "blood_to_site_gain_median": input.blood_to_site_gain_median,
            "blood_to_site_gain_gsd": input.blood_to_site_gain_gsd,
            "samples": input.samples,
            "seed": input.seed
        },
        "summary": {
            "p05_nM": s.p05_nanomolar,
            "p25_nM": s.p25_nanomolar,
            "p50_nM": s.p50_nanomolar,
            "p75_nM": s.p75_nanomolar,
            "p95_nM": s.p95_nanomolar,
            "p05_ng_mL": s.p05_ng_ml,
            "p25_ng_mL": s.p25_ng_ml,
            "p50_ng_mL": s.p50_ng_ml,
            "p75_ng_mL": s.p75_ng_ml,
            "p95_ng_mL": s.p95_ng_ml,
            "mean_ng_mL": s.mean_ng_ml
        }
    });
    fs::write(&json_path, serde_json::to_string_pretty(&payload).expect("serialize"))
        .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", csv_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "ms_cyclosporine_pk_bridge: p50={:.3} ng/mL p95={:.3} ng/mL (site_target={:.3} nM)",
        s.p50_ng_ml, s.p95_ng_ml, input.site_target_nanomolar
    );
}
