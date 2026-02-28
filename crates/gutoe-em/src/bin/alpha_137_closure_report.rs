use gutoe_em::alpha::{
    electron_mass_from_clifford_improved_with_alpha, lepton_masses_from_electron_structural_alpha,
    lepton_masses_from_electron_with_alpha, triangular, ALPHA_INVERSE_PHYSICAL,
    ALPHA_INVERSE_STRUCTURAL, CLIFFORD_DIM,
};

const M_E_OBS: f64 = 0.510_998_950;
const M_MU_OBS: f64 = 105.658_375_5;
const M_TAU_OBS: f64 = 1776.93;

fn main() {
    let alpha_struct = 1.0 / ALPHA_INVERSE_STRUCTURAL;
    let alpha_phys = 1.0 / ALPHA_INVERSE_PHYSICAL;

    // Pure algebra closure statement for alpha^{-1}.
    let t16 = triangular(1 << 4);
    let alpha_inv_from_algebra = t16 + 1;

    // Compare lepton lane under structural vs physical alpha.
    let masses_struct = lepton_masses_from_electron_structural_alpha(M_E_OBS);
    let masses_phys = lepton_masses_from_electron_with_alpha(M_E_OBS, alpha_phys);

    let [me_s, mmu_s, mtau_s] = masses_struct;
    let [me_p, mmu_p, mtau_p] = masses_phys;

    let mu_rel_s = (mmu_s - M_MU_OBS) / M_MU_OBS;
    let tau_rel_s = (mtau_s - M_TAU_OBS) / M_TAU_OBS;
    let mu_rel_p = (mmu_p - M_MU_OBS) / M_MU_OBS;
    let tau_rel_p = (mtau_p - M_TAU_OBS) / M_TAU_OBS;

    // Electron-from-(mu,tau) inversion in both lanes.
    let me_from_mu_tau_struct =
        electron_mass_from_clifford_improved_with_alpha(M_MU_OBS, M_TAU_OBS, alpha_struct);
    let me_from_mu_tau_phys =
        electron_mass_from_clifford_improved_with_alpha(M_MU_OBS, M_TAU_OBS, alpha_phys);

    println!("[alpha_closure]");
    println!(
        "clifford_dim = {}  T(16) = {}  alpha_inv = T(16)+1 = {}",
        CLIFFORD_DIM, t16, alpha_inv_from_algebra
    );
    println!(
        "alpha_inv_structural = {:.12}  alpha_inv_physical = {:.12}",
        ALPHA_INVERSE_STRUCTURAL, ALPHA_INVERSE_PHYSICAL
    );
    println!(
        "alpha_rel_offset = {:+.12e}",
        (ALPHA_INVERSE_PHYSICAL - ALPHA_INVERSE_STRUCTURAL) / ALPHA_INVERSE_STRUCTURAL
    );

    println!();
    println!("[lepton_lane_from_me_obs]");
    println!("structural_alpha: m_e={:.9} m_mu={:.9} m_tau={:.9}", me_s, mmu_s, mtau_s);
    println!(
        "structural_alpha rel_err: mu={:+.9e} tau={:+.9e}",
        mu_rel_s, tau_rel_s
    );
    println!("physical_alpha:   m_e={:.9} m_mu={:.9} m_tau={:.9}", me_p, mmu_p, mtau_p);
    println!(
        "physical_alpha rel_err:   mu={:+.9e} tau={:+.9e}",
        mu_rel_p, tau_rel_p
    );

    println!();
    println!("[electron_from_mu_tau]");
    println!(
        "structural_alpha m_e_pred={:.9} rel={:+.9e}",
        me_from_mu_tau_struct,
        (me_from_mu_tau_struct - M_E_OBS) / M_E_OBS
    );
    println!(
        "physical_alpha   m_e_pred={:.9} rel={:+.9e}",
        me_from_mu_tau_phys,
        (me_from_mu_tau_phys - M_E_OBS) / M_E_OBS
    );
}
