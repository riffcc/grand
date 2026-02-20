/*!
 * GUTOE Lattice - Hexagonal Spatial Structure
 * Copyright (C) 2026  Wings
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

//! Hexagonal lattice structure for vector rail propagation
//!
//! From GUTOE.md:
//! "Uses cubic coordinates (q, r, s) where q+r+s=0"
//! "O(1) distance calculation: max(|q1-q2|, |r1-r2|, |s1-s2|)"

pub mod hex_lattice;
pub mod vector_rails;
pub mod lattice_3d;

pub use hex_lattice::*;
pub use vector_rails::*;
pub use lattice_3d::*;
