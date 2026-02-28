use std::cmp::Ordering;
use std::env;
use std::fs;

const C_LIGHT: f64 = 299_792_458.0;
const C2: f64 = C_LIGHT * C_LIGHT;

#[derive(Clone, Copy, Debug)]
struct Level {
    label: &'static str,
    multiplier: f64,
    effort: u32,
}

#[derive(Clone, Copy, Debug)]
struct Angle {
    name: &'static str,
    baseline: f64,
    levels: &'static [Level],
}

#[derive(Clone, Debug)]
struct Combo {
    indices: Vec<usize>,
    multiplier: f64,
    effort: u32,
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

fn om_gain(multiplier: f64) -> f64 {
    multiplier.max(1.0e-300).log10()
}

fn main() {
    let out_dir = env::var("GUTOE_ANTIMATTER_OM_SWEEP_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/antimatter_om_reduction".to_string());
    fs::create_dir_all(&out_dir).expect("create output dir");

    let eta_floor = env_f64("GUTOE_ANTI_ETA_FLOOR", 0.166_727_168_685_6);
    let beam_power_mw_baseline = env_f64("GUTOE_ANTI_BEAM_POWER_MW_BASE", 5.0);

    static PROD: &[Level] = &[
        Level {
            label: "baseline",
            multiplier: 1.0,
            effort: 0,
        },
        Level {
            label: "w_target_microfocus",
            multiplier: 2.0,
            effort: 1,
        },
        Level {
            label: "iridium_segmented",
            multiplier: 3.5,
            effort: 2,
        },
        Level {
            label: "liquid_metal_jet",
            multiplier: 6.0,
            effort: 4,
        },
    ];
    static CAP: &[Level] = &[
        Level {
            label: "baseline",
            multiplier: 1.0,
            effort: 0,
        },
        Level {
            label: "nbti_lens",
            multiplier: 3.0,
            effort: 1,
        },
        Level {
            label: "nb3sn_high_field",
            multiplier: 5.0,
            effort: 2,
        },
        Level {
            label: "rebco_20t_array",
            multiplier: 10.0,
            effort: 4,
        },
    ];
    static COOL: &[Level] = &[
        Level {
            label: "baseline",
            multiplier: 1.0,
            effort: 0,
        },
        Level {
            label: "sic_multi_foil",
            multiplier: 2.2,
            effort: 1,
        },
        Level {
            label: "diamond_membrane_stack",
            multiplier: 3.8,
            effort: 2,
        },
        Level {
            label: "cryo_rfq_stack",
            multiplier: 6.0,
            effort: 4,
        },
    ];
    static REC: &[Level] = &[
        Level {
            label: "baseline",
            multiplier: 1.0,
            effort: 0,
        },
        Level {
            label: "dense_positron_stack",
            multiplier: 2.0,
            effort: 1,
        },
        Level {
            label: "plasma_lattice_control",
            multiplier: 3.0,
            effort: 2,
        },
        Level {
            label: "rydberg_assisted",
            multiplier: 5.0,
            effort: 3,
        },
    ];
    static STORE: &[Level] = &[
        Level {
            label: "baseline",
            multiplier: 1.0,
            effort: 0,
        },
        Level {
            label: "graphene_lined",
            multiplier: 2.0,
            effort: 1,
        },
        Level {
            label: "cryo_dlc",
            multiplier: 3.5,
            effort: 2,
        },
        Level {
            label: "hybrid_magnetic_graphene",
            multiplier: 6.0,
            effort: 4,
        },
    ];

    // New throughput levers to wipe additional OMs.
    static DUTY: &[Level] = &[
        Level {
            label: "baseline_pulsed",
            multiplier: 1.0,
            effort: 0,
        },
        Level {
            label: "improved_rf_stability",
            multiplier: 1.6,
            effort: 1,
        },
        Level {
            label: "rtsc_high_duty_capture",
            multiplier: 2.6,
            effort: 3,
        },
        Level {
            label: "continuous_near_cw",
            multiplier: 3.8,
            effort: 5,
        },
    ];
    static POWER: &[Level] = &[
        Level {
            label: "5mw_baseline",
            multiplier: 1.0,
            effort: 0,
        },
        Level {
            label: "20mw_modular",
            multiplier: 4.0,
            effort: 2,
        },
        Level {
            label: "100mw_multiline",
            multiplier: 20.0,
            effort: 5,
        },
        Level {
            label: "500mw_facility",
            multiplier: 100.0,
            effort: 9,
        },
    ];
    static PARALLEL: &[Level] = &[
        Level {
            label: "single_line",
            multiplier: 1.0,
            effort: 0,
        },
        Level {
            label: "dual_line",
            multiplier: 2.0,
            effort: 1,
        },
        Level {
            label: "quad_line",
            multiplier: 4.0,
            effort: 2,
        },
        Level {
            label: "hexadec_line",
            multiplier: 16.0,
            effort: 5,
        },
    ];
    static OPS: &[Level] = &[
        Level {
            label: "manual_ops",
            multiplier: 1.0,
            effort: 0,
        },
        Level {
            label: "closed_loop_tuning",
            multiplier: 1.4,
            effort: 1,
        },
        Level {
            label: "adaptive_beam_ai",
            multiplier: 2.0,
            effort: 2,
        },
        Level {
            label: "full_realtime_feedback",
            multiplier: 2.8,
            effort: 4,
        },
    ];

    let angles: [Angle; 9] = [
        Angle {
            name: "production",
            baseline: 1.0e-4,
            levels: PROD,
        },
        Angle {
            name: "capture",
            baseline: 1.0e-3,
            levels: CAP,
        },
        Angle {
            name: "cooling",
            baseline: 1.0e-2,
            levels: COOL,
        },
        Angle {
            name: "recombination",
            baseline: 1.0e-1,
            levels: REC,
        },
        Angle {
            name: "storage",
            baseline: 1.0,
            levels: STORE,
        },
        Angle {
            name: "duty_cycle",
            baseline: 1.0,
            levels: DUTY,
        },
        Angle {
            name: "beam_power",
            baseline: 1.0,
            levels: POWER,
        },
        Angle {
            name: "parallel_lines",
            baseline: 1.0,
            levels: PARALLEL,
        },
        Angle {
            name: "ops_feedback",
            baseline: 1.0,
            levels: OPS,
        },
    ];

    let eta_chain_baseline: f64 = angles.iter().map(|a| a.baseline).product();
    let ng_baseline = ng_per_year(beam_power_mw_baseline, eta_floor, eta_chain_baseline);

    let mut combos: Vec<Combo> = Vec::new();
    for i0 in 0..angles[0].levels.len() {
        for i1 in 0..angles[1].levels.len() {
            for i2 in 0..angles[2].levels.len() {
                for i3 in 0..angles[3].levels.len() {
                    for i4 in 0..angles[4].levels.len() {
                        for i5 in 0..angles[5].levels.len() {
                            for i6 in 0..angles[6].levels.len() {
                                for i7 in 0..angles[7].levels.len() {
                                    for i8 in 0..angles[8].levels.len() {
                                        let idx = vec![i0, i1, i2, i3, i4, i5, i6, i7, i8];
                                        let mut m = 1.0;
                                        let mut e = 0u32;
                                        for (aidx, lidx) in idx.iter().enumerate() {
                                            let lv = angles[aidx].levels[*lidx];
                                            m *= lv.multiplier;
                                            e += lv.effort;
                                        }
                                        combos.push(Combo {
                                            indices: idx,
                                            multiplier: m,
                                            effort: e,
                                        });
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    combos.sort_by(|a, b| {
        b.multiplier
            .partial_cmp(&a.multiplier)
            .unwrap_or(Ordering::Equal)
            .then_with(|| a.effort.cmp(&b.effort))
    });

    let top = combos.first().cloned().expect("non-empty combos");

    let budgets = [8u32, 12u32, 16u32, 20u32, 24u32, 28u32];
    let mut budget_best: Vec<(u32, Combo)> = Vec::new();
    for b in budgets {
        let best = combos
            .iter()
            .filter(|c| c.effort <= b)
            .max_by(|x, y| {
                x.multiplier
                    .partial_cmp(&y.multiplier)
                    .unwrap_or(Ordering::Equal)
            })
            .cloned()
            .expect("at least one combo under budget");
        budget_best.push((b, best));
    }

    let mut txt = String::new();
    txt.push_str("[antimatter_om_reduction_sweep]\n");
    txt.push_str(&format!("eta_floor = {:.12e}\n", eta_floor));
    txt.push_str(&format!("beam_power_mw_baseline = {:.6}\n", beam_power_mw_baseline));
    txt.push_str(&format!("eta_chain_baseline = {:.12e}\n", eta_chain_baseline));
    txt.push_str(&format!("ng_per_year_baseline = {:.12e}\n\n", ng_baseline));

    txt.push_str("[top_unconstrained]\n");
    let eta_top = eta_chain_baseline * top.multiplier;
    let ng_top = ng_per_year(beam_power_mw_baseline, eta_floor, eta_top);
    txt.push_str(&format!("multiplier_vs_baseline = {:.12e}\n", top.multiplier));
    txt.push_str(&format!("om_gain_vs_baseline = {:.6}\n", om_gain(top.multiplier)));
    txt.push_str(&format!("effort = {}\n", top.effort));
    txt.push_str(&format!("eta_chain = {:.12e}\n", eta_top));
    txt.push_str(&format!("ng_per_year = {:.12e}\n", ng_top));
    txt.push_str(&format!("years_per_1ng = {:.12e}\n\n", 1.0 / ng_top.max(1e-300)));

    txt.push_str("[budget_frontier]\n");
    for (b, c) in &budget_best {
        let eta = eta_chain_baseline * c.multiplier;
        let ng = ng_per_year(beam_power_mw_baseline, eta_floor, eta);
        txt.push_str(&format!(
            "budget={} effort={} multiplier={:.12e} om_gain={:.6} ng_per_year={:.12e} years_per_1ng={:.12e}\n",
            b,
            c.effort,
            c.multiplier,
            om_gain(c.multiplier),
            ng,
            1.0 / ng.max(1e-300)
        ));
    }
    txt.push('\n');

    // Print explicit lane for practical budget 16 and stretch budget 24.
    for target_budget in [16u32, 24u32] {
        if let Some((_, c)) = budget_best.iter().find(|(b, _)| *b == target_budget) {
            txt.push_str(&format!("[lane_budget_{}]\n", target_budget));
            for (aidx, lidx) in c.indices.iter().enumerate() {
                let lv = angles[aidx].levels[*lidx];
                txt.push_str(&format!(
                    "{} = {} (x{:.3}, effort={})\n",
                    angles[aidx].name, lv.label, lv.multiplier, lv.effort
                ));
            }
            txt.push('\n');
        }
    }

    let mut csv = String::from(
        "rank,multiplier,om_gain,effort,ng_per_year,years_per_1ng,production,capture,cooling,recombination,storage,duty_cycle,beam_power,parallel_lines,ops_feedback\n",
    );
    for (rank, c) in combos.iter().take(200).enumerate() {
        let eta = eta_chain_baseline * c.multiplier;
        let ng = ng_per_year(beam_power_mw_baseline, eta_floor, eta);
        let labels: Vec<&str> = c
            .indices
            .iter()
            .enumerate()
            .map(|(aidx, lidx)| angles[aidx].levels[*lidx].label)
            .collect();
        csv.push_str(&format!(
            "{},{:.12e},{:.12e},{},{:.12e},{:.12e},{},{},{},{},{},{},{},{},{}\n",
            rank + 1,
            c.multiplier,
            om_gain(c.multiplier),
            c.effort,
            ng,
            1.0 / ng.max(1e-300),
            labels[0],
            labels[1],
            labels[2],
            labels[3],
            labels[4],
            labels[5],
            labels[6],
            labels[7],
            labels[8]
        ));
    }

    let txt_path = format!("{out_dir}/antimatter_om_reduction_sweep.txt");
    let csv_path = format!("{out_dir}/antimatter_om_reduction_sweep_top200.csv");
    fs::write(&txt_path, txt).expect("write txt");
    fs::write(&csv_path, csv).expect("write csv");
    println!("wrote {}", txt_path);
    println!("wrote {}", csv_path);
}
