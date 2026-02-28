/*!
 * GUTOE Physics - Field Equations and Predictions
 * Copyright (C) 2026  Riff Labs
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

//! Physics equations from GUTOE framework
//!
//! Key equations to verify:
//! - Veracity wave equation: g^μν∇_μ∇_ν φ + λ_QG l_P² ∇⁴φ = 0
//! - Modified Einstein: G_μν + λ_QG l_P² H_μν = κ T_μν + ξ Λ g_μν
//! - Unification: G = v²/κ, c = v, ħ = l_P² κ
//! - Entropy: S = A/4G + α ln(A/l_P²)
//! - Dispersion: ω² = v² k² - λ_QG l_P² k⁴

pub mod abiogenesis;
pub mod baryogenesis;
pub mod bbn;
pub mod chemical_thermo;
pub mod chiral_symmetry_breaking;
pub mod cmb_class;
pub mod cmb_damping;
pub mod cmb_differential;
pub mod cmb_reionization;
pub mod constants;
pub mod cosmo_transfer;
pub mod dark_matter_falsification;
pub mod dark_sector;
pub mod dynamics_map;
pub mod energy_accounting;
pub mod entropy_progression;
pub mod equations;
pub mod falsifiability;
pub mod few_body_qm;
pub mod galactic_life_map;
pub mod great_filter;
pub mod homochirality;
pub mod inflation;
pub mod lithium7_stellar_depletion;
pub mod mass_gap;
pub mod microphysics;
pub mod multi_zone;
pub mod nuclear_chart;
pub mod nuclear_first_principles;
pub mod reaction_rates;
pub mod runtime_emulator;
pub mod single_zone;
pub mod singularity_resolution;
pub mod spectral_synthesis;
pub mod star_catalog;
pub mod stellar_reactions;
pub mod stiff_integrator;
pub mod uncertainty;
pub mod universe;

pub use abiogenesis::*;
pub use baryogenesis::*;
pub use bbn::*;
pub use chemical_thermo::*;
pub use chiral_symmetry_breaking::*;
pub use cmb_class::*;
pub use cmb_damping::*;
pub use cmb_differential::*;
pub use cmb_reionization::*;
pub use constants::*;
pub use cosmo_transfer::*;
pub use dark_matter_falsification::*;
pub use dark_sector::*;
pub use dynamics_map::*;
pub use energy_accounting::*;
pub use entropy_progression::*;
pub use equations::*;
pub use falsifiability::*;
pub use few_body_qm::*;
pub use galactic_life_map::*;
pub use great_filter::*;
pub use homochirality::*;
pub use inflation::*;
pub use lithium7_stellar_depletion::*;
pub use mass_gap::*;
pub use microphysics::*;
pub use multi_zone::*;
pub use nuclear_chart::*;
pub use nuclear_first_principles::*;
pub use reaction_rates::*;
pub use runtime_emulator::*;
pub use single_zone::*;
pub use singularity_resolution::*;
pub use spectral_synthesis::*;
pub use star_catalog::*;
pub use stellar_reactions::*;
pub use stiff_integrator::*;
pub use uncertainty::*;
pub use universe::*;
