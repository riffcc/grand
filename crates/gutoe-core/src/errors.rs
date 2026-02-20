/*!
 * GUTOE Core - Error Types
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

/*!
 * GUTOE Core - Error Types
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

//! Error types for GUTOE operations

use thiserror::Error;

#[derive(Error, Debug)]
pub enum GutoeError {
    #[error("Invalid qubit index: {0}")]
    InvalidQubit(usize),

    #[error("Invalid state transition: {0}")]
    InvalidTransition(String),

    #[error("Coherence below threshold: {0}")]
    LowCoherence(f64),

    #[error("Invalid lattice coordinates")]
    InvalidCoordinates,

    #[error("Simulation overflow: {0}")]
    Overflow(String),

    #[error("Math error: {0}")]
    MathError(String),
}
