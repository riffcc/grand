//! Emit a machine-readable GUTOE falsifiability scorecard.
//!
//! Usage:
//!   cargo run -p gutoe-physics --bin falsification_report
//!   GUTOE_SIN2_EW=0.23122 GUTOE_MZ_OVER_MW=1.1345 GUTOE_ALPHA_INV=137.036 \
//!     cargo run -p gutoe-physics --bin falsification_report

use gutoe_physics::{
    evaluate_structural, evaluate_with_corrected, provisional_corrected_observables,
    CorrectedObservables, FalsificationWindows,
};

fn env_f64(key: &str) -> Option<f64> {
    std::env::var(key).ok()?.parse::<f64>().ok()
}

fn main() {
    let structural = evaluate_structural();
    let corrected = if let (Some(sin2), Some(mzmw), Some(alpha_inv)) = (
        env_f64("GUTOE_SIN2_EW"),
        env_f64("GUTOE_MZ_OVER_MW"),
        env_f64("GUTOE_ALPHA_INV"),
    ) {
        CorrectedObservables {
            sin2_theta_w_ew: sin2,
            mz_over_mw: mzmw,
            alpha_inverse: alpha_inv,
            theta_qcd: env_f64("GUTOE_THETA_QCD").unwrap_or(0.0),
            neutron_edm_e_cm: env_f64("GUTOE_NEUTRON_EDM_E_CM").unwrap_or(0.0),
        }
    } else {
        provisional_corrected_observables()
    };
    let windows = FalsificationWindows::default();
    let full = evaluate_with_corrected(corrected, windows);

    println!(
        "{{\n  \"structural_ok\": {},\n  \"gauge_count_ok\": {},\n  \"lambda_qg_ok\": {},\n  \"corrected_ok\": {},\n  \"passes_all\": {},\n  \"corrected\": {{\"sin2_theta_w_ew\": {:.9}, \"mz_over_mw\": {:.9}, \"alpha_inverse\": {:.9}, \"theta_qcd\": {:.3e}, \"neutron_edm_e_cm\": {:.3e}}},\n  \"windows\": {{\"sin2_theta_w_ew\": [{:.5}, {:.5}], \"mz_over_mw\": [{:.4}, {:.4}], \"alpha_inverse\": [{:.3}, {:.3}], \"theta_qcd_abs_max\": {:.3e}, \"neutron_edm_abs_max_e_cm\": {:.3e}}}\n}}",
        structural.structural_ok,
        structural.gauge_count_ok,
        structural.lambda_qg_ok,
        full.corrected_ok,
        full.passes_all(),
        corrected.sin2_theta_w_ew,
        corrected.mz_over_mw,
        corrected.alpha_inverse,
        corrected.theta_qcd,
        corrected.neutron_edm_e_cm,
        windows.sin2_theta_w_ew.min,
        windows.sin2_theta_w_ew.max,
        windows.mz_over_mw.min,
        windows.mz_over_mw.max,
        windows.alpha_inverse.min,
        windows.alpha_inverse.max,
        windows.theta_qcd_abs_max,
        windows.neutron_edm_abs_max_e_cm
    );
}
