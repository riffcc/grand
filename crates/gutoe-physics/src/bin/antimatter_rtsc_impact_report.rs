use std::env;
use std::fs;

const C_LIGHT: f64 = 299_792_458.0;
const C2: f64 = C_LIGHT * C_LIGHT;

#[derive(Clone, Copy, Debug)]
struct Scenario {
    name: &'static str,
    // Multipliers relative to today's best five-angle stack.
    capture_uplift: f64,
    cooling_uplift: f64,
    recombination_uplift: f64,
    storage_uplift: f64,
}

fn env_f64(name: &str, default: f64) -> f64 {
    env::var(name)
        .ok()
        .and_then(|v| v.parse::<f64>().ok())
        .unwrap_or(default)
}

fn ng_per_year(beam_power_mw: f64, eta_floor: f64, eta_chain: f64) -> f64 {
    let beam_power_w = beam_power_mw * 1.0e6;
    let eta_net = eta_floor * eta_chain;
    let rest_power_w = beam_power_w * eta_net;
    let mass_rate_kg_s = rest_power_w / C2;
    mass_rate_kg_s * 365.25 * 24.0 * 3600.0 * 1.0e12
}

fn main() {
    let out_dir = env::var("GUTOE_ANTIMATTER_RTSC_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/antimatter_rtsc".to_string());
    fs::create_dir_all(&out_dir).expect("create output dir");

    let eta_floor = env_f64("GUTOE_ANTI_ETA_FLOOR", 0.166_727_168_685_6);
    let beam_power_mw = env_f64("GUTOE_ANTI_BEAM_POWER_MW", 5.0);

    // Baseline chain factors from prior lane.
    let eta_chain_baseline = 1.0e-10;

    // Today's best non-RTSC stack from five-angle sweep:
    // production x capture x cooling x recombination x storage
    // = 6 * 10 * 6 * 5 * 6 = 10800.
    let best_non_rtsc_multiplier = 10_800.0;
    let eta_chain_best_non_rtsc = eta_chain_baseline * best_non_rtsc_multiplier;
    let ng_best_non_rtsc = ng_per_year(beam_power_mw, eta_floor, eta_chain_best_non_rtsc);

    // RTSC material assumptions are hypothetical and exploratory.
    let scenarios = [
        Scenario {
            name: "rtsc_conservative",
            capture_uplift: 1.5,
            cooling_uplift: 1.4,
            recombination_uplift: 1.2,
            storage_uplift: 1.6,
        },
        Scenario {
            name: "rtsc_plausible",
            capture_uplift: 2.5,
            cooling_uplift: 2.0,
            recombination_uplift: 1.6,
            storage_uplift: 2.8,
        },
        Scenario {
            name: "rtsc_aggressive",
            capture_uplift: 4.0,
            cooling_uplift: 3.0,
            recombination_uplift: 2.2,
            storage_uplift: 4.5,
        },
    ];

    let mut txt = String::new();
    let mut csv = String::from(
        "scenario,total_uplift_vs_best_non_rtsc,total_multiplier_vs_baseline,orders_of_magnitude_vs_baseline,eta_chain,ng_per_year,ng_per_year_gain_vs_best_non_rtsc\n",
    );
    let mut sweep_csv =
        String::from("capture_uplift,cooling_uplift,recombination_uplift,storage_uplift,total_uplift_vs_best,total_multiplier_vs_baseline,orders_of_magnitude,ng_per_year\n");

    txt.push_str("[antimatter_rtsc_impact]\n");
    txt.push_str("Assumption: RTSC enables higher field duty-cycle + lower parasitic cryo overhead in capture/cooling/storage.\n");
    txt.push_str("This is a hypothetical material-science lane, not a claim of discovered ambient RT superconductivity.\n\n");
    txt.push_str(&format!("beam_power_mw = {:.6}\n", beam_power_mw));
    txt.push_str(&format!("eta_floor = {:.12e}\n", eta_floor));
    txt.push_str(&format!("eta_chain_baseline = {:.12e}\n", eta_chain_baseline));
    txt.push_str(&format!(
        "best_non_rtsc_multiplier = {:.6}\n",
        best_non_rtsc_multiplier
    ));
    txt.push_str(&format!(
        "best_non_rtsc_eta_chain = {:.12e}\n",
        eta_chain_best_non_rtsc
    ));
    txt.push_str(&format!(
        "best_non_rtsc_ng_per_year = {:.12e}\n\n",
        ng_best_non_rtsc
    ));

    for s in scenarios {
        let uplift = s.capture_uplift * s.cooling_uplift * s.recombination_uplift * s.storage_uplift;
        let total_mult = best_non_rtsc_multiplier * uplift;
        let eta_chain = eta_chain_baseline * total_mult;
        let ng = ng_per_year(beam_power_mw, eta_floor, eta_chain);
        let gain = ng / ng_best_non_rtsc.max(1e-300);
        let om = total_mult.log10();

        txt.push_str(&format!("[scenario:{}]\n", s.name));
        txt.push_str(&format!("capture_uplift = {:.6}\n", s.capture_uplift));
        txt.push_str(&format!("cooling_uplift = {:.6}\n", s.cooling_uplift));
        txt.push_str(&format!(
            "recombination_uplift = {:.6}\n",
            s.recombination_uplift
        ));
        txt.push_str(&format!("storage_uplift = {:.6}\n", s.storage_uplift));
        txt.push_str(&format!("total_uplift_vs_best_non_rtsc = {:.6}\n", uplift));
        txt.push_str(&format!(
            "total_multiplier_vs_baseline = {:.12e}\n",
            total_mult
        ));
        txt.push_str(&format!("orders_of_magnitude_vs_baseline = {:.6}\n", om));
        txt.push_str(&format!("eta_chain = {:.12e}\n", eta_chain));
        txt.push_str(&format!("ng_per_year = {:.12e}\n", ng));
        txt.push_str(&format!(
            "ng_per_year_gain_vs_best_non_rtsc = {:.12e}\n\n",
            gain
        ));

        csv.push_str(&format!(
            "{},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}\n",
            s.name, uplift, total_mult, om, eta_chain, ng, gain
        ));
    }

    // Parameter sweep to probe feasible RTSC window.
    let cap_vals = [1.2, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0];
    let cool_vals = [1.1, 1.3, 1.6, 2.0, 2.5, 3.0, 4.0];
    let rec_vals = [1.05, 1.1, 1.2, 1.4, 1.6, 2.0, 2.5];
    let store_vals = [1.2, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0];

    let mut best_ng = 0.0;
    let mut best_tuple = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    for &cap in &cap_vals {
        for &cool in &cool_vals {
            for &rec in &rec_vals {
                for &store in &store_vals {
                    let uplift = cap * cool * rec * store;
                    let total_mult = best_non_rtsc_multiplier * uplift;
                    let eta_chain = eta_chain_baseline * total_mult;
                    let ng = ng_per_year(beam_power_mw, eta_floor, eta_chain);
                    let om = total_mult.log10();
                    if ng > best_ng {
                        best_ng = ng;
                        best_tuple = (cap, cool, rec, store, uplift, total_mult, om);
                    }
                    sweep_csv.push_str(&format!(
                        "{:.6},{:.6},{:.6},{:.6},{:.12e},{:.12e},{:.12e},{:.12e}\n",
                        cap, cool, rec, store, uplift, total_mult, om, ng
                    ));
                }
            }
        }
    }

    txt.push_str("[rtsc_parameter_sweep]\n");
    txt.push_str("capture_uplift ∈ {1.2..5.0}, cooling_uplift ∈ {1.1..4.0}, recombination_uplift ∈ {1.05..2.5}, storage_uplift ∈ {1.2..5.0}\n");
    txt.push_str(&format!(
        "best_capture_uplift = {:.6}\n",
        best_tuple.0
    ));
    txt.push_str(&format!(
        "best_cooling_uplift = {:.6}\n",
        best_tuple.1
    ));
    txt.push_str(&format!(
        "best_recombination_uplift = {:.6}\n",
        best_tuple.2
    ));
    txt.push_str(&format!(
        "best_storage_uplift = {:.6}\n",
        best_tuple.3
    ));
    txt.push_str(&format!(
        "best_total_uplift_vs_best_non_rtsc = {:.12e}\n",
        best_tuple.4
    ));
    txt.push_str(&format!(
        "best_total_multiplier_vs_baseline = {:.12e}\n",
        best_tuple.5
    ));
    txt.push_str(&format!(
        "best_orders_of_magnitude_vs_baseline = {:.12e}\n",
        best_tuple.6
    ));
    txt.push_str(&format!("best_ng_per_year = {:.12e}\n", best_ng));
    txt.push_str(&format!(
        "best_ng_gain_vs_best_non_rtsc = {:.12e}\n",
        best_ng / ng_best_non_rtsc.max(1e-300)
    ));

    let txt_path = format!("{out_dir}/antimatter_rtsc_impact_report.txt");
    let csv_path = format!("{out_dir}/antimatter_rtsc_impact_report.csv");
    let sweep_path = format!("{out_dir}/antimatter_rtsc_impact_sweep.csv");
    fs::write(&txt_path, txt).expect("write txt");
    fs::write(&csv_path, csv).expect("write csv");
    fs::write(&sweep_path, sweep_csv).expect("write sweep csv");
    println!("wrote {}", txt_path);
    println!("wrote {}", csv_path);
    println!("wrote {}", sweep_path);
}

