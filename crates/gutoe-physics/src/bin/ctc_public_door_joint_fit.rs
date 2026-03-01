//! Public-data spontaneous-door joint fit (coarse pass).
//!
//! Purpose:
//! - Aggregate current public constraints across the strongest door lanes.
//! - Output a single coarse decision: `door_detected` true/false.
//!
//! Notes:
//! - This is a lightweight synthesis probe, not a publication-grade global fit.
//! - Uses lane-level residuals where quantitative public numbers are available.
//! - Keeps non-quantified lanes explicit as pending/inconclusive.

use serde_json::json;
use std::fs;
use std::path::PathBuf;

#[derive(Clone, Copy)]
struct NumericLane {
    name: &'static str,
    measured_minus_baseline: f64,
    sigma: f64,
    include_in_joint: bool,
    source: &'static str,
    notes: &'static str,
}

fn main() {
    let out_dir = std::env::var("GUTOE_DOOR_JOINT_OUT")
        .unwrap_or_else(|_| "/tmp/bh_renders/ctc_public_door_joint_fit".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    // Lane A: anisotropic cosmic birefringence amplitude A_CB
    // Best fit reported around 0.42e-4 with asymmetric +/- (0.40,0.34)e-4.
    // We use a symmetric sigma proxy 0.37e-4 for this coarse synthesis.
    let lane_biref = NumericLane {
        name: "cosmic_birefringence_anisotropic_A_CB",
        measured_minus_baseline: 0.42e-4,
        sigma: 0.37e-4,
        include_in_joint: true,
        source: "https://www.emergentmind.com/articles/2504.13154",
        notes: "Best-fit is consistent with zero (about <2σ); coarse symmetric sigma proxy used.",
    };

    // Lane B: EHT metric-deviation residual.
    // Public summary: observed size within ~10% of Kerr prediction.
    let lane_eht = NumericLane {
        name: "eht_shadow_fractional_deviation_from_kerr",
        measured_minus_baseline: 0.0,
        sigma: 0.10,
        include_in_joint: true,
        source: "https://eventhorizontelescope.org/publications/first-sagittarius-event-horizon-telescope-results-vi-testing-black-hole-metric",
        notes: "Using a conservative coarse residual model centered at zero with 10% scale.",
    };

    // Lane C: low-energy weak-charge residual (Qweak vs SM).
    let qweak_meas: f64 = 0.0719;
    let qweak_meas_sigma: f64 = 0.0045;
    let qweak_sm: f64 = 0.0708;
    let qweak_sm_sigma: f64 = 0.0003;
    let qweak_resid = qweak_meas - qweak_sm;
    let qweak_sigma = (qweak_meas_sigma.powi(2) + qweak_sm_sigma.powi(2)).sqrt();
    let lane_qweak = NumericLane {
        name: "electroweak_running_lowQ_weak_charge_residual",
        measured_minus_baseline: qweak_resid,
        sigma: qweak_sigma,
        include_in_joint: true,
        source: "https://www.jlab.org/research/qweak",
        notes: "Qweak result is consistent with SM within quoted uncertainty.",
    };

    // Lane D: void topology statistics (DESI DR1 void catalogs).
    // Public source currently provides consistency statements, but no single
    // residual/sigma scalar for a one-parameter door fit in this probe.
    let lane_void_source = "https://doi.org/10.3847/1538-4357/adb559";
    let lane_void_notes =
        "Generally consistent void properties; scalar anomaly residual not uniquely defined here.";

    let lanes = [lane_biref, lane_eht, lane_qweak];

    let mut included = Vec::new();
    let mut z_sum = 0.0_f64;
    let mut chi2_null = 0.0_f64;
    let mut signs = Vec::new();

    for lane in lanes {
        let z = if lane.sigma > 0.0 {
            lane.measured_minus_baseline / lane.sigma
        } else {
            0.0
        };
        let chi2 = z * z;
        if lane.include_in_joint {
            included.push(z);
            z_sum += z;
            chi2_null += chi2;
            signs.push(z.signum());
        }
    }

    let n = included.len() as f64;
    let z_stouffer = if n > 0.0 { z_sum / n.sqrt() } else { 0.0 };
    let mean_abs_z = if n > 0.0 {
        included.iter().map(|z| z.abs()).sum::<f64>() / n
    } else {
        0.0
    };
    let max_abs_z = included
        .iter()
        .map(|z| z.abs())
        .fold(0.0_f64, f64::max);
    let same_sign = signs
        .iter()
        .all(|s| *s >= 0.0)
        || signs.iter().all(|s| *s <= 0.0);

    // Coarse discovery gate (intentionally strict):
    // - strong combined significance
    // - at least one strong single-lane excursion
    // - consistent sign across included lanes
    let door_detected = z_stouffer.abs() >= 5.0 && max_abs_z >= 3.0 && same_sign;

    let lane_rows = vec![
        json!({
          "name": lane_biref.name,
          "measured_minus_baseline": lane_biref.measured_minus_baseline,
          "sigma": lane_biref.sigma,
          "z": lane_biref.measured_minus_baseline / lane_biref.sigma,
          "included_in_joint": lane_biref.include_in_joint,
          "source": lane_biref.source,
          "notes": lane_biref.notes
        }),
        json!({
          "name": lane_eht.name,
          "measured_minus_baseline": lane_eht.measured_minus_baseline,
          "sigma": lane_eht.sigma,
          "z": lane_eht.measured_minus_baseline / lane_eht.sigma,
          "included_in_joint": lane_eht.include_in_joint,
          "source": lane_eht.source,
          "notes": lane_eht.notes
        }),
        json!({
          "name": lane_qweak.name,
          "measured_minus_baseline": lane_qweak.measured_minus_baseline,
          "sigma": lane_qweak.sigma,
          "z": lane_qweak.measured_minus_baseline / lane_qweak.sigma,
          "included_in_joint": lane_qweak.include_in_joint,
          "source": lane_qweak.source,
          "notes": lane_qweak.notes
        }),
        json!({
          "name": "desi_void_topology_statistics",
          "included_in_joint": false,
          "source": lane_void_source,
          "notes": lane_void_notes,
          "status": "pending_quantitative_scalar"
        }),
    ];

    let summary = json!({
      "scope": "coarse public-data door hunt synthesis; not a publication-grade global fit",
      "joint_fit": {
        "n_numeric_lanes_included": n,
        "z_stouffer": z_stouffer,
        "chi2_null": chi2_null,
        "mean_abs_z": mean_abs_z,
        "max_abs_z": max_abs_z,
        "same_sign": same_sign,
        "door_detected": door_detected
      },
      "lanes": lane_rows,
      "interpretation": {
        "result": if door_detected {
          "door-compatible global excess detected (coarse gate)"
        } else {
          "no robust spontaneous-door signal in current coarse public-lane synthesis"
        },
        "next_required": "upgrade to full likelihood-level fit for each lane; include quantified void-stat residuals"
      }
    });

    let txt_path = out.join("ctc_public_door_joint_fit.txt");
    let json_path = out.join("ctc_public_door_joint_fit.json");

    let mut txt = String::new();
    txt.push_str("[ctc_public_door_joint_fit]\n");
    txt.push_str("coarse public-data spontaneous-door synthesis\n\n");
    txt.push_str(&format!("n_numeric_lanes_included = {}\n", n as usize));
    txt.push_str(&format!("z_stouffer = {:.6}\n", z_stouffer));
    txt.push_str(&format!("chi2_null = {:.6}\n", chi2_null));
    txt.push_str(&format!("mean_abs_z = {:.6}\n", mean_abs_z));
    txt.push_str(&format!("max_abs_z = {:.6}\n", max_abs_z));
    txt.push_str(&format!("same_sign = {}\n", same_sign));
    txt.push_str(&format!("door_detected = {}\n", door_detected));

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(
        &json_path,
        serde_json::to_string_pretty(&summary).expect("serialize json"),
    )
    .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "joint door hunt: z_stouffer={:.3}, max|z|={:.3}, detected={}",
        z_stouffer, max_abs_z, door_detected
    );
}
