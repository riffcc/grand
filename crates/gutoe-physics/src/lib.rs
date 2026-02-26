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

pub mod equations;
pub mod constants;
pub mod dynamics_map;
pub mod falsifiability;
pub mod runtime_emulator;
pub mod stellar_reactions;
pub mod reaction_rates;
pub mod stiff_integrator;
pub mod energy_accounting;
pub mod single_zone;
pub mod multi_zone;
pub mod spectral_synthesis;
pub mod star_catalog;
pub mod nuclear_chart;
pub mod mass_gap;
pub mod dark_sector;
pub mod dark_matter_falsification;

pub use equations::*;
pub use constants::*;
pub use dynamics_map::*;
pub use falsifiability::*;
pub use runtime_emulator::*;
pub use stellar_reactions::*;
pub use reaction_rates::*;
pub use stiff_integrator::*;
pub use energy_accounting::*;
pub use single_zone::*;
pub use multi_zone::*;
pub use spectral_synthesis::*;
pub use star_catalog::*;
pub use nuclear_chart::*;
pub use mass_gap::*;
pub use dark_sector::*;
pub use dark_matter_falsification::*;
