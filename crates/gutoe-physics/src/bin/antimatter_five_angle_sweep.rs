use rand::Rng;
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
    material_stack: &'static str,
}

#[derive(Clone, Copy, Debug)]
struct Angle {
    name: &'static str,
    baseline: f64,
    levels: &'static [Level],
}

#[derive(Clone, Debug)]
struct Combo {
    indices: [usize; 5],
    multiplier: f64,
    effort: u32,
}

fn env_f64(name: &str, default: f64) -> f64 {
    env::var(name)
        .ok()
        .and_then(|v| v.parse::<f64>().ok())
        .unwrap_or(default)
}

fn ng_per_year_from_chain(beam_power_mw: f64, eta_floor: f64, eta_chain: f64) -> f64 {
    let beam_power_w = beam_power_mw * 1.0e6;
    let eta_net = eta_floor * eta_chain;
    let rest_power_w = beam_power_w * eta_net;
    let mass_rate_kg_s = rest_power_w / C2;
    mass_rate_kg_s * 365.25 * 24.0 * 3600.0 * 1.0e12
}

fn main() {
    let out_dir = env::var("GUTOE_ANTIMATTER_FIVE_ANGLE_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/antimatter_five_angle".to_string());
    fs::create_dir_all(&out_dir).expect("create output dir");

    let eta_floor = env_f64("GUTOE_ANTI_ETA_FLOOR", 0.166_727_168_685_6);
    let beam_power_mw = env_f64("GUTOE_ANTI_BEAM_POWER_MW", 5.0);

    static PROD_LEVELS: &[Level] = &[
        Level {
            label: "baseline",
            multiplier: 1.0,
            effort: 0,
            material_stack: "existing_target",
        },
        Level {
            label: "w_target_microfocus",
            multiplier: 2.0,
            effort: 1,
            material_stack: "tungsten+microchannel_heat_sink",
        },
        Level {
            label: "iridium_segmented",
            multiplier: 3.5,
            effort: 2,
            material_stack: "iridium+segment_cooling",
        },
        Level {
            label: "liquid_metal_jet",
            multiplier: 6.0,
            effort: 4,
            material_stack: "li_pb_bi_jet+active_mhd_stabilization",
        },
    ];
    static CAPTURE_LEVELS: &[Level] = &[
        Level {
            label: "baseline",
            multiplier: 1.0,
            effort: 0,
            material_stack: "existing_capture_optics",
        },
        Level {
            label: "nbti_lens",
            multiplier: 3.0,
            effort: 1,
            material_stack: "nbti_superconducting_lens",
        },
        Level {
            label: "nb3sn_high_field",
            multiplier: 5.0,
            effort: 2,
            material_stack: "nb3sn_quadrupole_capture",
        },
        Level {
            label: "rebco_20t_array",
            multiplier: 10.0,
            effort: 4,
            material_stack: "rebco_hts_20t_multi_lens",
        },
    ];
    static COOL_LEVELS: &[Level] = &[
        Level {
            label: "baseline",
            multiplier: 1.0,
            effort: 0,
            material_stack: "single_stage_deceleration",
        },
        Level {
            label: "sic_multi_foil",
            multiplier: 2.2,
            effort: 1,
            material_stack: "sic_foil_stack",
        },
        Level {
            label: "diamond_membrane_stack",
            multiplier: 3.8,
            effort: 2,
            material_stack: "cvd_diamond_membranes+rf_phasing",
        },
        Level {
            label: "cryo_rfq_stack",
            multiplier: 6.0,
            effort: 4,
            material_stack: "cryogenic_rfq+diamond_window_train",
        },
    ];
    static RECOMB_LEVELS: &[Level] = &[
        Level {
            label: "baseline",
            multiplier: 1.0,
            effort: 0,
            material_stack: "single_trap_recombination",
        },
        Level {
            label: "dense_positron_stack",
            multiplier: 2.0,
            effort: 1,
            material_stack: "nanoporous_accumulator+cold_buffer_gas",
        },
        Level {
            label: "plasma_lattice_control",
            multiplier: 3.0,
            effort: 2,
            material_stack: "phase_locked_penning_array",
        },
        Level {
            label: "rydberg_assisted",
            multiplier: 5.0,
            effort: 3,
            material_stack: "rydberg_assist+adiabatic_field_ramp",
        },
    ];
    static STORAGE_LEVELS: &[Level] = &[
        Level {
            label: "baseline",
            multiplier: 1.0,
            effort: 0,
            material_stack: "standard_trap_walls",
        },
        Level {
            label: "graphene_lined",
            multiplier: 2.0,
            effort: 1,
            material_stack: "graphene_coated_surfaces",
        },
        Level {
            label: "cryo_dlc",
            multiplier: 3.5,
            effort: 2,
            material_stack: "cryogenic_dlc+ultra_low_outgassing",
        },
        Level {
            label: "hybrid_magnetic_graphene",
            multiplier: 6.0,
            effort: 4,
            material_stack: "superconducting_bottle+graphene_inner_shell",
        },
    ];

    let angles = [
        Angle {
            name: "production",
            baseline: 1.0e-4,
            levels: PROD_LEVELS,
        },
        Angle {
            name: "capture",
            baseline: 1.0e-3,
            levels: CAPTURE_LEVELS,
        },
        Angle {
            name: "cooling",
            baseline: 1.0e-2,
            levels: COOL_LEVELS,
        },
        Angle {
            name: "recombination",
            baseline: 1.0e-1,
            levels: RECOMB_LEVELS,
        },
        Angle {
            name: "storage",
            baseline: 1.0,
            levels: STORAGE_LEVELS,
        },
    ];

    let eta_chain_baseline: f64 = angles.iter().map(|a| a.baseline).product();
    let ng_baseline = ng_per_year_from_chain(beam_power_mw, eta_floor, eta_chain_baseline);

    let mut combos: Vec<Combo> = Vec::new();
    for i0 in 0..angles[0].levels.len() {
        for i1 in 0..angles[1].levels.len() {
            for i2 in 0..angles[2].levels.len() {
                for i3 in 0..angles[3].levels.len() {
                    for i4 in 0..angles[4].levels.len() {
                        let indices = [i0, i1, i2, i3, i4];
                        let mut m = 1.0;
                        let mut e = 0u32;
                        for (aidx, lidx) in indices.iter().enumerate() {
                            let lv = angles[aidx].levels[*lidx];
                            m *= lv.multiplier;
                            e += lv.effort;
                        }
                        combos.push(Combo {
                            indices,
                            multiplier: m,
                            effort: e,
                        });
                    }
                }
            }
        }
    }

    let mut by_best = combos.clone();
    by_best.sort_by(|a, b| {
        b.multiplier
            .partial_cmp(&a.multiplier)
            .unwrap_or(Ordering::Equal)
            .then_with(|| a.effort.cmp(&b.effort))
    });

    let targets = [100.0, 1_000.0, 10_000.0];
    let mut best_for_target: Vec<Option<Combo>> = vec![None; targets.len()];
    for c in &combos {
        for (tidx, t) in targets.iter().enumerate() {
            if c.multiplier < *t {
                continue;
            }
            let improve = match &best_for_target[tidx] {
                None => true,
                Some(cur) => (c.effort, c.multiplier) < (cur.effort, cur.multiplier),
            };
            if improve {
                best_for_target[tidx] = Some(c.clone());
            }
        }
    }

    let top = by_best.first().cloned().expect("non-empty combinations");
    let top_chain = eta_chain_baseline * top.multiplier;
    let ng_top = ng_per_year_from_chain(beam_power_mw, eta_floor, top_chain);

    // Robustness probe for top combo: random 20% multiplicative uncertainty each lane.
    let mut rng = rand::thread_rng();
    let samples = 20_000;
    let mut sampled: Vec<f64> = Vec::with_capacity(samples);
    for _ in 0..samples {
        let mut m = 1.0;
        for (aidx, lidx) in top.indices.iter().enumerate() {
            let base = angles[aidx].levels[*lidx].multiplier;
            let jitter: f64 = rng.gen_range(0.8..=1.2);
            m *= base * jitter;
        }
        sampled.push(m);
    }
    sampled.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
    let p10 = sampled[(samples as f64 * 0.10) as usize];
    let p50 = sampled[(samples as f64 * 0.50) as usize];
    let p90 = sampled[(samples as f64 * 0.90) as usize];

    let mut txt = String::new();
    txt.push_str("[antimatter_five_angle_sweep]\n");
    txt.push_str(&format!("eta_floor = {:.12e}\n", eta_floor));
    txt.push_str(&format!("beam_power_mw = {:.6}\n", beam_power_mw));
    txt.push_str(&format!("eta_chain_baseline = {:.12e}\n", eta_chain_baseline));
    txt.push_str(&format!("ng_per_year_baseline = {:.12e}\n", ng_baseline));
    txt.push('\n');

    txt.push_str("[angle_caps]\n");
    for angle in &angles {
        let best_level = angle
            .levels
            .iter()
            .max_by(|x, y| {
                x.multiplier
                    .partial_cmp(&y.multiplier)
                    .unwrap_or(Ordering::Equal)
            })
            .expect("angle has levels");
        txt.push_str(&format!(
            "{}: best_level={} multiplier={:.3} effort={} material={}\n",
            angle.name,
            best_level.label,
            best_level.multiplier,
            best_level.effort,
            best_level.material_stack
        ));
    }
    txt.push('\n');

    txt.push_str("[target_min_effort]\n");
    for (tidx, target) in targets.iter().enumerate() {
        match &best_for_target[tidx] {
            None => txt.push_str(&format!("target_{:.0}x: unattained\n", target)),
            Some(c) => {
                let eta_chain = eta_chain_baseline * c.multiplier;
                let ng = ng_per_year_from_chain(beam_power_mw, eta_floor, eta_chain);
                txt.push_str(&format!(
                    "target_{:.0}x: effort={} multiplier={:.6} eta_chain={:.12e} ng_per_year={:.12e}\n",
                    target, c.effort, c.multiplier, eta_chain, ng
                ));
                for (aidx, lidx) in c.indices.iter().enumerate() {
                    let lv = angles[aidx].levels[*lidx];
                    txt.push_str(&format!(
                        "  {} = {} (x{:.3}, material={})\n",
                        angles[aidx].name, lv.label, lv.multiplier, lv.material_stack
                    ));
                }
            }
        }
        txt.push('\n');
    }

    txt.push_str("[top_unconstrained_combo]\n");
    txt.push_str(&format!(
        "multiplier={:.6} effort={} eta_chain={:.12e} ng_per_year={:.12e}\n",
        top.multiplier, top.effort, top_chain, ng_top
    ));
    for (aidx, lidx) in top.indices.iter().enumerate() {
        let lv = angles[aidx].levels[*lidx];
        txt.push_str(&format!(
            "{} = {} (x{:.3}, material={})\n",
            angles[aidx].name, lv.label, lv.multiplier, lv.material_stack
        ));
    }
    txt.push('\n');

    txt.push_str("[top_combo_robustness_monte_carlo]\n");
    txt.push_str("assumption: per-angle multiplicative uncertainty = ±20% uniform\n");
    txt.push_str(&format!("samples = {}\n", samples));
    txt.push_str(&format!("p10_multiplier = {:.6}\n", p10));
    txt.push_str(&format!("p50_multiplier = {:.6}\n", p50));
    txt.push_str(&format!("p90_multiplier = {:.6}\n", p90));
    txt.push_str(&format!("p10_om = {:.6}\n", p10.log10()));
    txt.push_str(&format!("p50_om = {:.6}\n", p50.log10()));
    txt.push_str(&format!("p90_om = {:.6}\n", p90.log10()));

    let mut csv = String::from(
        "rank,multiplier,orders_of_magnitude,effort,production,capture,cooling,recombination,storage\n",
    );
    for (rank, c) in by_best.iter().take(120).enumerate() {
        let p = angles[0].levels[c.indices[0]].label;
        let cap = angles[1].levels[c.indices[1]].label;
        let cool = angles[2].levels[c.indices[2]].label;
        let rec = angles[3].levels[c.indices[3]].label;
        let st = angles[4].levels[c.indices[4]].label;
        csv.push_str(&format!(
            "{},{:.12e},{:.12e},{},{},{},{},{},{}\n",
            rank + 1,
            c.multiplier,
            c.multiplier.log10(),
            c.effort,
            p,
            cap,
            cool,
            rec,
            st
        ));
    }

    let txt_path = format!("{out_dir}/antimatter_five_angle_sweep.txt");
    let csv_path = format!("{out_dir}/antimatter_five_angle_sweep_top120.csv");
    fs::write(&txt_path, txt).expect("write txt");
    fs::write(&csv_path, csv).expect("write csv");
    println!("wrote {}", txt_path);
    println!("wrote {}", csv_path);
}

