/*!
 * GUTOE Core - Tripartite Quantum System
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

//! GUTOE Core - Tripartite Quantum System
//!
//! Implements the tripartite quantum system from the Grand Unifying Theory of Everything
//!
//! Core states:
//! - VOID: Pure undefined state (the "null" reference)
//! - SINE: |1⟩ (existence/presence)
//! - COSINE: |0⟩ (absence/orthogonality)
//! - TANGENT: tan = sin/cos — the slope/ratio state (the "relationship" state)
//!
//! Also implements the 12-state hexagonal system from VOID-DIFFERENTIATION.md:
//! "The system evolves into a hexagonal structure with twelve states (six on each face)"
//!
//! See: /mnt/castle/garage/gutoe-research/GUTOE.md

pub mod states;
pub mod gates;
pub mod errors;
pub mod hex_states;

pub use states::*;
pub use gates::*;
pub use errors::*;
pub use hex_states::*;
