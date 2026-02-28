use std::env;
use std::fs;

const C_LIGHT: f64 = 299_792_458.0;
const C2: f64 = C_LIGHT * C_LIGHT;

#[derive(Clone, Copy, Debug)]
struct FactorLevel {
    label: &'static str,
    multiplier: f64,
    effort: u32,
}

#[derive(Clone, Copy, Debug)]
struct ChainFactor {
    name: &'static str,
    baseline: f64,
    levels: &'static [FactorLevel],
}

fn env_f64(name: &str, default: f64) -> f64 {
    env::var(name)
        .ok()
        .and_then(|v| v.parse::<f64>().ok())
        .unwrap_or(default)
}

fn main() {
    let out_dir = env::var("GUTOE_ANTIMATTER_TWO_OM_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/antimatter_two_om_plan".to_string());
    fs::create_dir_all(&out_dir).expect("create output dir");

    // From the prior production lane:
    let eta_floor = env_f64("GUTOE_ANTI_ETA_FLOOR", 0.1667271686856);
    let beam_power_mw = env_f64("GUTOE_ANTI_BEAM_POWER_MW", 5.0); // "today_optimistic" anchor

    // Baseline chain decomposition (product = 1e-10).
    // These are engineering factors, separated so we can test upgrade combinations.
    static PROD_LEVELS: &[FactorLevel] = &[
        FactorLevel {
            label: "baseline",
            multiplier: 1.0,
            effort: 0,
        },
        FactorLevel {
            label: "target_optics",
            multiplier: 2.0,
            effort: 1,
        },
        FactorLevel {
            label: "target_optics_plus",
            multiplier: 3.0,
            effort: 2,
        },
    ];
    static CAPTURE_LEVELS: &[FactorLevel] = &[
        FactorLevel {
            label: "baseline",
            multiplier: 1.0,
            effort: 0,
        },
        FactorLevel {
            label: "high_acceptance_lens",
            multiplier: 3.0,
            effort: 1,
        },
        FactorLevel {
            label: "high_acceptance_lens_plus",
            multiplier: 5.0,
            effort: 2,
        },
    ];
    static COOL_LEVELS: &[FactorLevel] = &[
        FactorLevel {
            label: "baseline",
            multiplier: 1.0,
            effort: 0,
        },
        FactorLevel {
            label: "stacked_deceleration",
            multiplier: 3.0,
            effort: 1,
        },
        FactorLevel {
            label: "stacked_deceleration_plus",
            multiplier: 5.0,
            effort: 2,
        },
    ];
    static RECOMB_LEVELS: &[FactorLevel] = &[
        FactorLevel {
            label: "baseline",
            multiplier: 1.0,
            effort: 0,
        },
        FactorLevel {
            label: "nested_trap_recomb",
            multiplier: 2.0,
            effort: 1,
        },
        FactorLevel {
            label: "nested_trap_recomb_plus",
            multiplier: 4.0,
            effort: 2,
        },
    ];
    static STORAGE_LEVELS: &[FactorLevel] = &[
        FactorLevel {
            label: "baseline",
            multiplier: 1.0,
            effort: 0,
        },
        FactorLevel {
            label: "lifetime_control",
            multiplier: 2.0,
            effort: 1,
        },
        FactorLevel {
            label: "lifetime_control_plus",
            multiplier: 3.0,
            effort: 2,
        },
    ];

    let factors = [
        ChainFactor {
            name: "production",
            baseline: 1.0e-4,
            levels: PROD_LEVELS,
        },
        ChainFactor {
            name: "capture",
            baseline: 1.0e-3,
            levels: CAPTURE_LEVELS,
        },
        ChainFactor {
            name: "cooling",
            baseline: 1.0e-2,
            levels: COOL_LEVELS,
        },
        ChainFactor {
            name: "recombination",
            baseline: 1.0e-1,
            levels: RECOMB_LEVELS,
        },
        ChainFactor {
            name: "storage",
            baseline: 1.0,
            levels: STORAGE_LEVELS,
        },
    ];

    let eta_chain_baseline: f64 = factors.iter().map(|f| f.baseline).product();
    let eta_net_baseline = eta_floor * eta_chain_baseline;
    let target_multiplier = 100.0;

    #[derive(Clone, Debug)]
    struct Candidate {
        labels: Vec<(&'static str, &'static str)>,
        upgrade_multiplier: f64,
        effort: u32,
        changed: u32,
    }

    let mut best: Option<Candidate> = None;

    for p in factors[0].levels {
        for c in factors[1].levels {
            for d in factors[2].levels {
                for r in factors[3].levels {
                    for s in factors[4].levels {
                        let levels = [p, c, d, r, s];
                        let upgrade_multiplier: f64 =
                            levels.iter().map(|x| x.multiplier).product();
                        if upgrade_multiplier < target_multiplier {
                            continue;
                        }
                        let effort: u32 = levels.iter().map(|x| x.effort).sum();
                        let changed: u32 = levels.iter().filter(|x| x.multiplier > 1.0).count()
                            as u32;
                        let labels = vec![
                            (factors[0].name, p.label),
                            (factors[1].name, c.label),
                            (factors[2].name, d.label),
                            (factors[3].name, r.label),
                            (factors[4].name, s.label),
                        ];
                        let cand = Candidate {
                            labels,
                            upgrade_multiplier,
                            effort,
                            changed,
                        };
                        let better = match &best {
                            None => true,
                            Some(cur) => {
                                (cand.effort, cand.changed, cand.upgrade_multiplier)
                                    < (cur.effort, cur.changed, cur.upgrade_multiplier)
                            }
                        };
                        if better {
                            best = Some(cand);
                        }
                    }
                }
            }
        }
    }

    let best = best.expect("at least one >=100x candidate");
    let eta_chain_upgraded = eta_chain_baseline * best.upgrade_multiplier;
    let eta_net_upgraded = eta_floor * eta_chain_upgraded;

    let beam_power_w = beam_power_mw * 1.0e6;
    let rest_power_baseline = beam_power_w * eta_net_baseline;
    let rest_power_upgraded = beam_power_w * eta_net_upgraded;

    let mass_rate_baseline = rest_power_baseline / C2;
    let mass_rate_upgraded = rest_power_upgraded / C2;
    let ng_per_year_baseline = mass_rate_baseline * 365.25 * 24.0 * 3600.0 * 1.0e12;
    let ng_per_year_upgraded = mass_rate_upgraded * 365.25 * 24.0 * 3600.0 * 1.0e12;

    let mut txt = String::new();
    txt.push_str("[antimatter_two_om_plan]\n");
    txt.push_str(&format!("eta_floor = {:.12e}\n", eta_floor));
    txt.push_str(&format!("beam_power_mw = {:.6}\n", beam_power_mw));
    txt.push_str(&format!("eta_chain_baseline = {:.12e}\n", eta_chain_baseline));
    txt.push_str(&format!("eta_net_baseline = {:.12e}\n", eta_net_baseline));
    txt.push_str(&format!("target_multiplier = {:.3}\n\n", target_multiplier));

    txt.push_str("[best_candidate_ge_100x]\n");
    txt.push_str(&format!(
        "upgrade_multiplier = {:.6}\n",
        best.upgrade_multiplier
    ));
    txt.push_str(&format!("effort_score = {}\n", best.effort));
    txt.push_str(&format!("changed_factors = {}\n", best.changed));
    for (name, label) in &best.labels {
        txt.push_str(&format!("{} = {}\n", name, label));
    }
    txt.push('\n');
    txt.push_str(&format!("eta_chain_upgraded = {:.12e}\n", eta_chain_upgraded));
    txt.push_str(&format!("eta_net_upgraded = {:.12e}\n", eta_net_upgraded));
    txt.push_str(&format!(
        "rest_power_baseline_w = {:.12e}\n",
        rest_power_baseline
    ));
    txt.push_str(&format!(
        "rest_power_upgraded_w = {:.12e}\n",
        rest_power_upgraded
    ));
    txt.push_str(&format!(
        "ng_per_year_baseline = {:.12e}\n",
        ng_per_year_baseline
    ));
    txt.push_str(&format!(
        "ng_per_year_upgraded = {:.12e}\n",
        ng_per_year_upgraded
    ));
    txt.push_str(&format!(
        "ng_per_year_gain = {:.12e}\n",
        ng_per_year_upgraded / ng_per_year_baseline.max(1e-300)
    ));

    let mut csv = String::from(
        "factor,baseline,selected_label,selected_multiplier,selected_effort\n",
    );
    let selected = &best.labels;
    let selected_mults = [
        selected[0].1,
        selected[1].1,
        selected[2].1,
        selected[3].1,
        selected[4].1,
    ];
    for (idx, f) in factors.iter().enumerate() {
        let level = f
            .levels
            .iter()
            .find(|lv| lv.label == selected_mults[idx])
            .expect("selected level exists");
        csv.push_str(&format!(
            "{},{:.12e},{},{:.6},{}\n",
            f.name, f.baseline, level.label, level.multiplier, level.effort
        ));
    }

    let txt_path = format!("{out_dir}/antimatter_two_om_plan.txt");
    let csv_path = format!("{out_dir}/antimatter_two_om_plan.csv");
    fs::write(&txt_path, txt).expect("write txt");
    fs::write(&csv_path, csv).expect("write csv");
    println!("wrote {}", txt_path);
    println!("wrote {}", csv_path);
}

