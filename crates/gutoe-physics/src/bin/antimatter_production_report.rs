use std::env;
use std::fs;

const C_LIGHT: f64 = 299_792_458.0;
const C2: f64 = C_LIGHT * C_LIGHT;
const EV_TO_J: f64 = 1.602_176_634e-19;
const GEV_TO_J: f64 = 1.0e9 * EV_TO_J;

const M_PROTON_GEV: f64 = 0.938_272_088_16;
const M_ELECTRON_GEV: f64 = 0.000_510_998_95;
const M_ANTI_HYDROGEN_KG: f64 = 1.673_532_84e-27; // ~= m_p + m_e

#[derive(Clone, Copy, Debug)]
struct Scenario {
    name: &'static str,
    beam_power_mw: f64,
    // Engineering chain efficiency beyond kinematic floor:
    // production x collection x cooling/deceleration x trapping/recombination.
    eta_chain: f64,
}

fn env_f64(name: &str, default: f64) -> f64 {
    env::var(name)
        .ok()
        .and_then(|v| v.parse::<f64>().ok())
        .unwrap_or(default)
}

fn antiproton_threshold_lab_kinetic_gev() -> f64 {
    // p + p(rest) -> p + p + p + pbar
    // In lab with target at rest: E_beam,total,thr = 7 m_p
    // => K_thr = 6 m_p
    6.0 * M_PROTON_GEV
}

fn positron_pair_floor_input_gev_per_pos() -> f64 {
    // gamma-gamma (or gamma-nucleus) pair production floor:
    // 1.022 MeV input creates e+ and e-; allocate full energy to one produced e+.
    2.0 * M_ELECTRON_GEV
}

fn antihydrogen_floor_input_gev() -> f64 {
    antiproton_threshold_lab_kinetic_gev() + positron_pair_floor_input_gev_per_pos()
}

fn antihydrogen_rest_gev() -> f64 {
    M_PROTON_GEV + M_ELECTRON_GEV
}

fn kinematic_floor_efficiency() -> f64 {
    antihydrogen_rest_gev() / antihydrogen_floor_input_gev()
}

fn main() {
    let out_dir = env::var("GUTOE_ANTIMATTER_PROD_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/antimatter_production".to_string());
    fs::create_dir_all(&out_dir).expect("create output dir");

    let k_floor = kinematic_floor_efficiency();
    let scenarios = vec![
        Scenario {
            name: "today_conservative",
            beam_power_mw: env_f64("GUTOE_ANTI_P_MW_TODAY", 1.0),
            eta_chain: env_f64("GUTOE_ANTI_ETA_CHAIN_TODAY", 1.0e-12),
        },
        Scenario {
            name: "today_optimistic",
            beam_power_mw: env_f64("GUTOE_ANTI_P_MW_OPT", 5.0),
            eta_chain: env_f64("GUTOE_ANTI_ETA_CHAIN_OPT", 1.0e-10),
        },
        Scenario {
            name: "near_term_stretch",
            beam_power_mw: env_f64("GUTOE_ANTI_P_MW_STRETCH", 20.0),
            eta_chain: env_f64("GUTOE_ANTI_ETA_CHAIN_STRETCH", 1.0e-8),
        },
    ];

    let mut csv = String::from(
        "scenario,beam_power_mw,eta_chain,eta_net,rest_power_w,atoms_per_s,kg_per_year,ng_per_year,years_per_1ng,years_per_1ug\n",
    );
    let mut txt = String::new();

    txt.push_str("[antimatter_production_report]\n");
    txt.push_str("physical route modeled:\n");
    txt.push_str("1) antiprotons via p+p -> p+p+p+pbar near threshold\n");
    txt.push_str("2) positrons via pair production (1.022 MeV floor per produced e+)\n");
    txt.push_str("3) trapping/deceleration/recombination into antihydrogen\n\n");

    let k_thr = antiproton_threshold_lab_kinetic_gev();
    let e_pos = positron_pair_floor_input_gev_per_pos();
    let e_floor = antihydrogen_floor_input_gev();
    let e_rest = antihydrogen_rest_gev();
    txt.push_str(&format!(
        "antiproton kinetic threshold (lab) = {:.9} GeV\n",
        k_thr
    ));
    txt.push_str(&format!("positron floor input/particle = {:.9} GeV\n", e_pos));
    txt.push_str(&format!(
        "antihydrogen floor input/atom = {:.9} GeV\n",
        e_floor
    ));
    txt.push_str(&format!("antihydrogen rest/atom = {:.9} GeV\n", e_rest));
    txt.push_str(&format!(
        "kinematic floor efficiency (rest/input) = {:.9}\n\n",
        k_floor
    ));

    txt.push_str("Assumption: eta_net = eta_chain * kinematic_floor_efficiency\n");
    txt.push_str("Caveat: this ignores many second-order machine losses and hadronic multiplicity penalties.\n\n");

    for s in scenarios {
        let beam_power_w = s.beam_power_mw * 1.0e6;
        let eta_net = (s.eta_chain * k_floor).max(0.0);
        let rest_power_w = beam_power_w * eta_net;
        let mass_rate_kg_s = rest_power_w / C2;
        let atoms_per_s = mass_rate_kg_s / M_ANTI_HYDROGEN_KG;
        let kg_per_year = mass_rate_kg_s * 365.25 * 24.0 * 3600.0;
        let ng_per_year = kg_per_year * 1.0e12;
        let years_per_1ng = if ng_per_year > 0.0 {
            1.0 / ng_per_year
        } else {
            f64::INFINITY
        };
        let years_per_1ug = if ng_per_year > 0.0 {
            1000.0 / ng_per_year
        } else {
            f64::INFINITY
        };

        csv.push_str(&format!(
            "{},{:.6},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}\n",
            s.name,
            s.beam_power_mw,
            s.eta_chain,
            eta_net,
            rest_power_w,
            atoms_per_s,
            kg_per_year,
            ng_per_year,
            years_per_1ng,
            years_per_1ug
        ));

        txt.push_str(&format!("[scenario:{}]\n", s.name));
        txt.push_str(&format!("beam_power_mw = {:.6}\n", s.beam_power_mw));
        txt.push_str(&format!("eta_chain = {:.12e}\n", s.eta_chain));
        txt.push_str(&format!("eta_net = {:.12e}\n", eta_net));
        txt.push_str(&format!("rest_power_w = {:.12e}\n", rest_power_w));
        txt.push_str(&format!("atoms_per_s = {:.12e}\n", atoms_per_s));
        txt.push_str(&format!("kg_per_year = {:.12e}\n", kg_per_year));
        txt.push_str(&format!("ng_per_year = {:.12e}\n", ng_per_year));
        txt.push_str(&format!("years_per_1ng = {:.12e}\n", years_per_1ng));
        txt.push_str(&format!("years_per_1ug = {:.12e}\n\n", years_per_1ug));
    }

    // Pure kinematic lower-bound power required for target masses.
    let mut lower_bound_csv = String::from("target_mass,kg,input_energy_j,input_energy_twh\n");
    txt.push_str("[kinematic_lower_bound_only]\n");
    for (name, kg) in [
        ("1_ng", 1.0e-12),
        ("1_ug", 1.0e-9),
        ("1_mg", 1.0e-6),
        ("1_g", 1.0e-3),
    ] {
        let e_rest_j = kg * C2;
        let e_input_j = e_rest_j / k_floor.max(1.0e-30);
        let e_twh = e_input_j / 3.6e15;
        lower_bound_csv.push_str(&format!("{},{:.12e},{:.12e},{:.12e}\n", name, kg, e_input_j, e_twh));
        txt.push_str(&format!(
            "{}: input_energy_j={:.12e} ({:.6} TWh)\n",
            name, e_input_j, e_twh
        ));
    }

    // Cross-check conversion for one anti-H atom.
    let e_floor_j = e_floor * GEV_TO_J;
    let e_rest_j = e_rest * GEV_TO_J;
    txt.push_str(&format!(
        "\nper_antihydrogen_atom: floor_input={:.12e} J, rest={:.12e} J\n",
        e_floor_j, e_rest_j
    ));

    let txt_path = format!("{out_dir}/antimatter_production_report.txt");
    let csv_path = format!("{out_dir}/antimatter_production_report.csv");
    let lb_path = format!("{out_dir}/antimatter_production_lower_bound.csv");
    fs::write(&txt_path, txt).expect("write txt");
    fs::write(&csv_path, csv).expect("write csv");
    fs::write(&lb_path, lower_bound_csv).expect("write lower-bound csv");
    println!("wrote {}", txt_path);
    println!("wrote {}", csv_path);
    println!("wrote {}", lb_path);
}

