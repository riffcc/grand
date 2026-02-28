//! GRAND-132: FCNC suppression report from CKM/GIM structure.

use gutoe_em::{
    channel_label, fcnc_gim_from_clifford, fcnc_gim_from_textures, up_flavors,
    FCNC_LOOP_PROXY_EXPECTED,
};
use serde_json::json;
use std::fs;
use std::path::PathBuf;

const TREE_EPS: f64 = 1.0e-12;
const GIM_SUM_EPS: f64 = 1.0e-12;
const MASS_DIFF_EPS: f64 = 1.0e-12;
const LOOP_SUPPRESSION_MAX: f64 = 0.10;
const LOOP_PROXY_TOL: f64 = 1.0e-15;

fn all_channels_below(metrics: &gutoe_em::FcncGimMetrics, max_ratio: f64) -> bool {
    metrics
        .channels
        .iter()
        .all(|ch| ch.gim_suppression_ratio < max_ratio)
}

fn all_mass_diff_residuals_ok(metrics: &gutoe_em::FcncGimMetrics, eps: f64) -> bool {
    metrics
        .channels
        .iter()
        .all(|ch| ch.mass_difference_form_residual_abs < eps)
}

fn channel_rows(metrics: &gutoe_em::FcncGimMetrics) -> Vec<serde_json::Value> {
    metrics
        .channels
        .iter()
        .map(|ch| {
            json!({
                "channel": channel_label(ch),
                "from": ch.from,
                "to": ch.to,
                "lambda_u_abs": ch.lambda_u_abs,
                "lambda_c_abs": ch.lambda_c_abs,
                "lambda_t_abs": ch.lambda_t_abs,
                "lambda_sum_abs": ch.lambda_sum_abs,
                "degenerate_kernel_abs": ch.degenerate_kernel_abs,
                "split_kernel_abs": ch.split_kernel_abs,
                "split_kernel_no_gim_abs": ch.split_kernel_no_gim_abs,
                "gim_suppression_ratio": ch.gim_suppression_ratio,
                "mass_difference_form_residual_abs": ch.mass_difference_form_residual_abs,
            })
        })
        .collect()
}

fn main() {
    let out_dir =
        std::env::var("GUTOE_FCNC_OUT").unwrap_or_else(|_| "/tmp/bh_renders/fcnc_gim".to_string());
    let out = PathBuf::from(out_dir);
    let _ = fs::create_dir_all(&out);

    let direct = fcnc_gim_from_clifford();
    let texture = fcnc_gim_from_textures();

    let tree_level_ok = direct.neutral_current_offdiag_max_abs < TREE_EPS
        && texture.neutral_current_offdiag_max_abs < TREE_EPS
        && direct.neutral_current_diag_drift_max_abs < TREE_EPS
        && texture.neutral_current_diag_drift_max_abs < TREE_EPS;

    let gim_sum_rule_ok = direct.gim_sum_rule_residual_max_abs < GIM_SUM_EPS
        && texture.gim_sum_rule_residual_max_abs < GIM_SUM_EPS;

    let mass_difference_rewrite_ok = all_mass_diff_residuals_ok(&direct, MASS_DIFF_EPS)
        && all_mass_diff_residuals_ok(&texture, MASS_DIFF_EPS);

    let loop_suppression_ok = all_channels_below(&direct, LOOP_SUPPRESSION_MAX)
        && all_channels_below(&texture, LOOP_SUPPRESSION_MAX);

    let structural_proxy_ok =
        (direct.structural_loop_proxy - FCNC_LOOP_PROXY_EXPECTED).abs() < LOOP_PROXY_TOL;

    let passes_all = tree_level_ok
        && gim_sum_rule_ok
        && mass_difference_rewrite_ok
        && loop_suppression_ok
        && structural_proxy_ok;

    let txt_path = out.join("fcnc_gim_report.txt");
    let json_path = out.join("fcnc_gim_report.json");

    let mut txt = String::new();
    txt.push_str("[meta]\n");
    txt.push_str("lane = GRAND-132_fcnc_gim\n");
    txt.push_str(&format!("up_flavors = {:?}\n", up_flavors()));
    txt.push_str("\n");

    txt.push_str("[direct_ckm]\n");
    txt.push_str(&format!(
        "s23 = {:.12e}\ns13 = {:.12e}\nstructural_loop_proxy = {:.12e}\nneutral_current_offdiag_max_abs = {:.12e}\nneutral_current_diag_drift_max_abs = {:.12e}\ngim_sum_rule_residual_max_abs = {:.12e}\nkernel_u = {:.12e}\nkernel_c = {:.12e}\nkernel_t = {:.12e}\n",
        direct.ckm_s23,
        direct.ckm_s13,
        direct.structural_loop_proxy,
        direct.neutral_current_offdiag_max_abs,
        direct.neutral_current_diag_drift_max_abs,
        direct.gim_sum_rule_residual_max_abs,
        direct.loop_kernel_u,
        direct.loop_kernel_c,
        direct.loop_kernel_t
    ));
    for ch in &direct.channels {
        txt.push_str(&format!(
            "channel_{}: lambda_sum_abs={:.12e}, split_abs={:.12e}, no_gim_abs={:.12e}, suppression={:.12e}, mass_diff_residual={:.12e}\n",
            channel_label(ch),
            ch.lambda_sum_abs,
            ch.split_kernel_abs,
            ch.split_kernel_no_gim_abs,
            ch.gim_suppression_ratio,
            ch.mass_difference_form_residual_abs
        ));
    }
    txt.push_str("\n");

    txt.push_str("[texture_ckm]\n");
    txt.push_str(&format!(
        "s23 = {:.12e}\ns13 = {:.12e}\nstructural_loop_proxy = {:.12e}\nneutral_current_offdiag_max_abs = {:.12e}\nneutral_current_diag_drift_max_abs = {:.12e}\ngim_sum_rule_residual_max_abs = {:.12e}\nkernel_u = {:.12e}\nkernel_c = {:.12e}\nkernel_t = {:.12e}\n",
        texture.ckm_s23,
        texture.ckm_s13,
        texture.structural_loop_proxy,
        texture.neutral_current_offdiag_max_abs,
        texture.neutral_current_diag_drift_max_abs,
        texture.gim_sum_rule_residual_max_abs,
        texture.loop_kernel_u,
        texture.loop_kernel_c,
        texture.loop_kernel_t
    ));
    for ch in &texture.channels {
        txt.push_str(&format!(
            "channel_{}: lambda_sum_abs={:.12e}, split_abs={:.12e}, no_gim_abs={:.12e}, suppression={:.12e}, mass_diff_residual={:.12e}\n",
            channel_label(ch),
            ch.lambda_sum_abs,
            ch.split_kernel_abs,
            ch.split_kernel_no_gim_abs,
            ch.gim_suppression_ratio,
            ch.mass_difference_form_residual_abs
        ));
    }
    txt.push_str("\n");

    txt.push_str("[gate]\n");
    txt.push_str(&format!("tree_level_ok = {}\n", tree_level_ok));
    txt.push_str(&format!("gim_sum_rule_ok = {}\n", gim_sum_rule_ok));
    txt.push_str(&format!(
        "mass_difference_rewrite_ok = {}\n",
        mass_difference_rewrite_ok
    ));
    txt.push_str(&format!("loop_suppression_ok = {}\n", loop_suppression_ok));
    txt.push_str(&format!("structural_proxy_ok = {}\n", structural_proxy_ok));
    txt.push_str(&format!("passes_all = {}\n", passes_all));

    let payload = json!({
        "meta": {
            "lane": "GRAND-132_fcnc_gim",
            "up_flavors": up_flavors(),
            "constants": {
                "fcnc_loop_proxy_expected": FCNC_LOOP_PROXY_EXPECTED,
                "tree_eps": TREE_EPS,
                "gim_sum_eps": GIM_SUM_EPS,
                "mass_diff_eps": MASS_DIFF_EPS,
                "loop_suppression_max": LOOP_SUPPRESSION_MAX,
                "loop_proxy_tol": LOOP_PROXY_TOL,
            }
        },
        "direct_ckm": {
            "summary": direct,
            "channels": channel_rows(&direct),
        },
        "texture_ckm": {
            "summary": texture,
            "channels": channel_rows(&texture),
        },
        "gate": {
            "tree_level_ok": tree_level_ok,
            "gim_sum_rule_ok": gim_sum_rule_ok,
            "mass_difference_rewrite_ok": mass_difference_rewrite_ok,
            "loop_suppression_ok": loop_suppression_ok,
            "structural_proxy_ok": structural_proxy_ok,
            "passes_all": passes_all,
        }
    });

    fs::write(&txt_path, txt).expect("write txt");
    fs::write(
        &json_path,
        serde_json::to_string_pretty(&payload).expect("serialize json"),
    )
    .expect("write json");

    println!("wrote {}", txt_path.display());
    println!("wrote {}", json_path.display());
    println!(
        "fcnc_gim: pass={} direct_offdiag={:.3e} texture_offdiag={:.3e}",
        passes_all, direct.neutral_current_offdiag_max_abs, texture.neutral_current_offdiag_max_abs
    );
}
