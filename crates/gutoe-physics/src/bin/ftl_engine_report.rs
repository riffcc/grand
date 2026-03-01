//! FTL engine matrix against current GUTOE causality + geometry gates.
//!
//! This report separates:
//! - local propagation bounds (`v_g <= c`)
//! - geometry/topology shortcuts (warp metrics / wormholes)
//! and evaluates major FTL proposals under one consistent constraint set.

use gutoe_physics::constants::{G, PLANCK_LENGTH};
use gutoe_physics::dark_sector::{dark_density_particle, vacuum_energy_density_structural};
use gutoe_physics::singularity_resolution::lattice_core_radius_m;
use gutoe_physics::C;
use serde_json::json;
use std::f64::consts::PI;
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

fn env_bool(name: &str, default: bool) -> bool {
    match std::env::var(name) {
        Ok(v) => matches!(v.as_str(), "1" | "true" | "TRUE" | "yes" | "YES"),
        Err(_) => default,
    }
}

const REAR_FACE_FACTOR: f64 = 0.1;
const COMPACT_RADIUS_M: f64 = 100.0;
const COMPACT_THICKNESS_M: f64 = 0.1;
const COMPACT_SHORTCUT_S: f64 = 0.5; // 2c effective
const VISIBLE_DENSITY_OPT_KG_M3: f64 = 1.0e-20;

fn required_curvature_energy_density_j_m3(shortcut_s: f64, curvature_scale_m: f64) -> f64 {
    // Einstein scale proxy: ρ ~ (c^4 / 8πG) * K, K ~ ((1/s)-1)/R^2
    let pref = C.powi(4) / (8.0 * PI * G);
    let amp = (1.0 / shortcut_s - 1.0).max(0.0);
    pref * amp / curvature_scale_m.powi(2)
}

fn compact_shell_volume_m3(radius_m: f64, thickness_m: f64) -> f64 {
    (4.0 / 3.0) * PI * ((radius_m + thickness_m).powi(3) - radius_m.powi(3))
}

#[derive(Debug, Clone, Copy)]
enum EngineClass {
    WarpMetric,
    Wormhole,
    Tube,
    QuantumShortcut,
    Tachyonic,
}

impl EngineClass {
    fn as_str(self) -> &'static str {
        match self {
            EngineClass::WarpMetric => "warp_metric",
            EngineClass::Wormhole => "wormhole",
            EngineClass::Tube => "tube",
            EngineClass::QuantumShortcut => "quantum_shortcut",
            EngineClass::Tachyonic => "tachyonic",
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct EngineSpec {
    name: &'static str,
    class: EngineClass,
    local_signal_speed_gt_c_required: bool,
    requires_negative_energy: bool,
    requires_anec_violation: bool,
    ctc_risk: bool,
    traversable_as_stated: bool,
    requires_compact_positive_budget: bool,
}

#[derive(Debug, Clone, Copy)]
struct FrameworkGates {
    local_signal_bound_enforced: bool,
    supports_macroscopic_negative_energy: bool,
    supports_anec_violation: bool,
    allows_ctc: bool,
    lentz_compact_budget_met: bool,
    lentz_compact_floor_j: f64,
    lentz_compact_shortfall_ratio: f64,
    near_core_dead_end_multiplier: f64,
}

#[derive(Debug, Clone)]
struct EngineVerdict {
    spec: EngineSpec,
    pass: bool,
    blocked_by: Vec<&'static str>,
}

fn evaluate_engine(spec: EngineSpec, gates: FrameworkGates) -> EngineVerdict {
    let mut blocked_by = Vec::new();

    if spec.local_signal_speed_gt_c_required && gates.local_signal_bound_enforced {
        blocked_by.push("local_signal_bound_vg_le_c");
    }
    if spec.requires_negative_energy && !gates.supports_macroscopic_negative_energy {
        blocked_by.push("negative_energy_not_available");
    }
    if spec.requires_anec_violation && !gates.supports_anec_violation {
        blocked_by.push("anec_violation_not_available");
    }
    if spec.ctc_risk && !gates.allows_ctc {
        blocked_by.push("chronology_protection_ctc_block");
    }
    if spec.requires_compact_positive_budget && !gates.lentz_compact_budget_met {
        blocked_by.push("compact_shell_energy_floor_not_met");
    }
    if !spec.traversable_as_stated {
        blocked_by.push("not_traversable_as_stated");
    }

    EngineVerdict {
        spec,
        pass: blocked_by.is_empty(),
        blocked_by,
    }
}

fn best_ftl_engines() -> Vec<EngineSpec> {
    vec![
        EngineSpec {
            name: "Alcubierre warp bubble",
            class: EngineClass::WarpMetric,
            local_signal_speed_gt_c_required: false,
            requires_negative_energy: true,
            requires_anec_violation: true,
            ctc_risk: true,
            traversable_as_stated: true,
            requires_compact_positive_budget: false,
        },
        EngineSpec {
            name: "Natario warp metric",
            class: EngineClass::WarpMetric,
            local_signal_speed_gt_c_required: false,
            requires_negative_energy: true,
            requires_anec_violation: true,
            ctc_risk: true,
            traversable_as_stated: true,
            requires_compact_positive_budget: false,
        },
        EngineSpec {
            name: "Van den Broeck warp variant",
            class: EngineClass::WarpMetric,
            local_signal_speed_gt_c_required: false,
            requires_negative_energy: true,
            requires_anec_violation: true,
            ctc_risk: true,
            traversable_as_stated: true,
            requires_compact_positive_budget: false,
        },
        EngineSpec {
            name: "Lentz-style warp soliton",
            class: EngineClass::WarpMetric,
            local_signal_speed_gt_c_required: false,
            requires_negative_energy: false,
            requires_anec_violation: false,
            ctc_risk: true,
            traversable_as_stated: true,
            requires_compact_positive_budget: true,
        },
        EngineSpec {
            name: "Morris-Thorne traversable wormhole",
            class: EngineClass::Wormhole,
            local_signal_speed_gt_c_required: false,
            requires_negative_energy: true,
            requires_anec_violation: true,
            ctc_risk: true,
            traversable_as_stated: true,
            requires_compact_positive_budget: false,
        },
        EngineSpec {
            name: "Einstein-Rosen bridge (classical)",
            class: EngineClass::Wormhole,
            local_signal_speed_gt_c_required: false,
            requires_negative_energy: false,
            requires_anec_violation: false,
            ctc_risk: false,
            traversable_as_stated: false,
            requires_compact_positive_budget: false,
        },
        EngineSpec {
            name: "Krasnikov tube",
            class: EngineClass::Tube,
            local_signal_speed_gt_c_required: false,
            requires_negative_energy: true,
            requires_anec_violation: true,
            ctc_risk: true,
            traversable_as_stated: true,
            requires_compact_positive_budget: false,
        },
        EngineSpec {
            name: "Quantum entanglement messaging",
            class: EngineClass::QuantumShortcut,
            local_signal_speed_gt_c_required: false,
            requires_negative_energy: false,
            requires_anec_violation: false,
            ctc_risk: false,
            traversable_as_stated: false,
            requires_compact_positive_budget: false,
        },
        EngineSpec {
            name: "Tachyonic signaling drive",
            class: EngineClass::Tachyonic,
            local_signal_speed_gt_c_required: true,
            requires_negative_energy: false,
            requires_anec_violation: false,
            ctc_risk: true,
            traversable_as_stated: true,
            requires_compact_positive_budget: false,
        },
    ]
}

fn main() {
    let out_dir = std::env::var("GUTOE_FTL_REPORT_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ftl_engine_report".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    // Shared gate assumptions from current lane state:
    // 1) Lean theorem + runtime tests enforce local no-FTL group velocity.
    // 2) Current cosmology lane yields positive vacuum density (not macroscopic
    //    engineered negative-energy support).
    let rho_vac = vacuum_energy_density_structural();
    let shell_volume_m3 = compact_shell_volume_m3(COMPACT_RADIUS_M, COMPACT_THICKNESS_M);
    let rho_required_compact =
        required_curvature_energy_density_j_m3(COMPACT_SHORTCUT_S, COMPACT_RADIUS_M);
    let lentz_compact_floor_j = rho_required_compact * shell_volume_m3 * REAR_FACE_FACTOR;

    let rho_dark = dark_density_particle(VISIBLE_DENSITY_OPT_KG_M3);
    let rho_source_j_m3 = (rho_dark + rho_vac) * C * C;
    let source_energy_same_shell_j = rho_source_j_m3 * shell_volume_m3;
    let lentz_compact_shortfall_ratio = lentz_compact_floor_j / source_energy_same_shell_j;
    let lentz_compact_budget_met = lentz_compact_shortfall_ratio <= 1.0;

    let r_core = lattice_core_radius_m(PLANCK_LENGTH);
    let rho_required_core = required_curvature_energy_density_j_m3(COMPACT_SHORTCUT_S, r_core);
    let lentz_core_floor_j = rho_required_core * shell_volume_m3 * REAR_FACE_FACTOR;
    let near_core_dead_end_multiplier = lentz_core_floor_j / lentz_compact_floor_j;

    let gates = FrameworkGates {
        local_signal_bound_enforced: true,
        supports_macroscopic_negative_energy: false,
        supports_anec_violation: false,
        allows_ctc: env_bool("GUTOE_ALLOW_CTC", false),
        lentz_compact_budget_met,
        lentz_compact_floor_j,
        lentz_compact_shortfall_ratio,
        near_core_dead_end_multiplier,
    };

    let specs = best_ftl_engines();
    let verdicts: Vec<EngineVerdict> = specs
        .iter()
        .copied()
        .map(|s| evaluate_engine(s, gates))
        .collect();

    let pass_count = verdicts.iter().filter(|v| v.pass).count();
    let fail_count = verdicts.len().saturating_sub(pass_count);

    let txt_path = out.join("ftl_engine_report.txt");
    let json_path = out.join("ftl_engine_report.json");

    let mut txt = File::create(&txt_path).expect("create txt");
    writeln!(txt, "[framework_gates]").expect("write");
    writeln!(txt, "c_m_per_s = {:.6}", C).expect("write");
    writeln!(
        txt,
        "vacuum_energy_density_structural_kg_m3 = {:.12e}",
        rho_vac
    )
    .expect("write");
    writeln!(
        txt,
        "local_signal_bound_enforced = {}",
        gates.local_signal_bound_enforced
    )
    .expect("write");
    writeln!(
        txt,
        "supports_macroscopic_negative_energy = {}",
        gates.supports_macroscopic_negative_energy
    )
    .expect("write");
    writeln!(
        txt,
        "supports_anec_violation = {}",
        gates.supports_anec_violation
    )
    .expect("write");
    writeln!(txt, "allows_ctc = {}", gates.allows_ctc).expect("write");
    writeln!(
        txt,
        "lentz_compact_budget_met = {}",
        gates.lentz_compact_budget_met
    )
    .expect("write");
    writeln!(
        txt,
        "lentz_compact_floor_j = {:.12e}",
        gates.lentz_compact_floor_j
    )
    .expect("write");
    writeln!(
        txt,
        "lentz_compact_shortfall_ratio = {:.12e}",
        gates.lentz_compact_shortfall_ratio
    )
    .expect("write");
    writeln!(
        txt,
        "near_core_dead_end_multiplier = {:.12e}",
        gates.near_core_dead_end_multiplier
    )
    .expect("write");
    writeln!(txt).expect("write");

    writeln!(txt, "[engine_verdicts]").expect("write");
    for v in &verdicts {
        writeln!(txt, "- engine = {}", v.spec.name).expect("write");
        writeln!(txt, "  class = {}", v.spec.class.as_str()).expect("write");
        writeln!(txt, "  pass = {}", v.pass).expect("write");
        if v.blocked_by.is_empty() {
            writeln!(txt, "  blocked_by = []").expect("write");
        } else {
            writeln!(txt, "  blocked_by = {:?}", v.blocked_by).expect("write");
        }
    }
    writeln!(txt).expect("write");
    writeln!(txt, "[summary]").expect("write");
    writeln!(txt, "total = {}", verdicts.len()).expect("write");
    writeln!(txt, "pass = {}", pass_count).expect("write");
    writeln!(txt, "fail = {}", fail_count).expect("write");

    let rows: Vec<_> = verdicts
        .iter()
        .map(|v| {
            json!({
                "name": v.spec.name,
                "class": v.spec.class.as_str(),
                "assumptions": {
                    "local_signal_speed_gt_c_required": v.spec.local_signal_speed_gt_c_required,
                    "requires_negative_energy": v.spec.requires_negative_energy,
                    "requires_anec_violation": v.spec.requires_anec_violation,
                    "ctc_risk": v.spec.ctc_risk,
                    "traversable_as_stated": v.spec.traversable_as_stated,
                    "requires_compact_positive_budget": v.spec.requires_compact_positive_budget,
                },
                "pass": v.pass,
                "blocked_by": v.blocked_by,
            })
        })
        .collect();

    let payload = json!({
        "framework_gates": {
            "c_m_per_s": C,
            "vacuum_energy_density_structural_kg_m3": rho_vac,
            "local_signal_bound_enforced": gates.local_signal_bound_enforced,
            "supports_macroscopic_negative_energy": gates.supports_macroscopic_negative_energy,
            "supports_anec_violation": gates.supports_anec_violation,
            "allows_ctc": gates.allows_ctc,
            "lentz_compact_budget_met": gates.lentz_compact_budget_met,
            "lentz_compact_floor_j": gates.lentz_compact_floor_j,
            "lentz_compact_shortfall_ratio": gates.lentz_compact_shortfall_ratio,
            "near_core_dead_end_multiplier": gates.near_core_dead_end_multiplier
        },
        "summary": {
            "total": verdicts.len(),
            "pass": pass_count,
            "fail": fail_count
        },
        "engines": rows
    });
    fs::write(
        &json_path,
        serde_json::to_string_pretty(&payload).expect("json encode"),
    )
    .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "ftl matrix: pass={} fail={} (total={})",
        pass_count,
        fail_count,
        verdicts.len()
    );
}
