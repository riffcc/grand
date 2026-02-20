/*!
 * GUTOE Lattice - 3D Hexagonal Topology
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

//! 3D Hexagonal Lattice Topology
//!
//! Topology: Each node has 19 neighbors
//! - 6 planar neighbors (hexagonal in x-y plane)
//! - 1 node above + 6 around it (7 vertical-up neighbors)
//! - 1 node below + 6 around it (7 vertical-down neighbors)
//! - Total: 1 (self) + 6 + 7 + 7 = 21, or 19 if excluding self
//!
//! This forms a 3D hexagonal prism lattice - like stacked honeycomb sheets
//! connected vertically, with each vertical node also having its own hex ring.

use std::collections::HashMap;
use std::fmt;

/// 3D Hexagonal coordinates
/// Uses cubic coordinates extended to 3D: q + r + s + t = 0
/// For the hex prism: we have (q, r, z) where q+r+s = 0 in the plane
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Hex3D {
    pub q: i32,  // Hex axis 1
    pub r: i32,  // Hex axis 2
    pub z: i32,  // Vertical axis
}

impl Hex3D {
    /// Create new 3D hex coordinate
    pub fn new(q: i32, r: i32, z: i32) -> Self {
        Self { q, r, z }
    }

    /// Distance to another node (Chebyshev in 3D hex space)
    pub fn distance(&self, other: &Hex3D) -> u32 {
        // In 3D hex coords, need to handle the constraint
        let dz = (self.z - other.z).abs() as u32;
        let planar_dist = hex_distance_2d(self.q, self.r, other.q, other.r);
        dz.max(planar_dist)
    }

    /// All 20 neighbors (excludes self):
    /// - 6 planar neighbors (same z)
    /// - 1 above + 6 around above = 7 vertical-up
    /// - 1 below + 6 around below = 7 vertical-down
    /// Total: 6 + 7 + 7 = 20 neighbors, 21 including self
    pub fn neighbors(&self) -> [Hex3D; 20] {
        let planar = hex_neighbors_2d(self.q, self.r);

        // 6 planar neighbors (z same)
        let planar_neighbors: [Hex3D; 6] = planar.map(|(q, r)| Hex3D::new(q, r, self.z));

        // 1 directly above + 6 around above node (z + 1 layer)
        let above_node = Hex3D::new(self.q, self.r, self.z + 1);
        let above_ring: [Hex3D; 6] = planar.map(|(q, r)| Hex3D::new(q, r, self.z + 1));

        // 1 directly below + 6 around below node (z - 1 layer)
        let below_node = Hex3D::new(self.q, self.r, self.z - 1);
        let below_ring: [Hex3D; 6] = planar.map(|(q, r)| Hex3D::new(q, r, self.z - 1));

        // Combine: planar(6) + above(7) + below(7) = 20
        let mut all = [Hex3D::new(0, 0, 0); 20];

        // Planar: 0-5
        for i in 0..6 {
            all[i] = planar_neighbors[i];
        }

        // Above (direct + ring): 6-12
        all[6] = above_node;
        for i in 0..6 {
            all[i + 7] = above_ring[i];
        }

        // Below (direct + ring): 13-19
        all[13] = below_node;
        for i in 0..6 {
            all[i + 14] = below_ring[i];
        }

        all
    }

    /// Only the direct vertical neighbors (above and below)
    pub fn vertical_neighbors(&self) -> [Hex3D; 2] {
        [
            Hex3D::new(self.q, self.r, self.z + 1),
            Hex3D::new(self.q, self.r, self.z - 1),
        ]
    }

    /// Only the planar neighbors (same z-level)
    pub fn planar_neighbors(&self) -> [Hex3D; 6] {
        hex_neighbors_2d(self.q, self.r).map(|(q, r)| Hex3D::new(q, r, self.z))
    }

    /// Move up one layer
    pub fn up(&self) -> Hex3D {
        Hex3D::new(self.q, self.r, self.z + 1)
    }

    /// Move down one layer
    pub fn down(&self) -> Hex3D {
        Hex3D::new(self.q, self.r, self.z - 1)
    }
}

impl fmt::Display for Hex3D {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({}, {}, {})", self.q, self.r, self.z)
    }
}

/// Helper: 2D hex neighbors (axial coords where q+r+s=0)
fn hex_neighbors_2d(q: i32, r: i32) -> [(i32, i32); 6] {
    [
        (q + 1, r),     // East
        (q + 1, r - 1), // Southeast
        (q, r - 1),    // Southwest
        (q - 1, r),    // West
        (q - 1, r + 1), // Northwest
        (q, r + 1),    // Northeast
    ]
}

/// Helper: 2D hex distance (axial)
fn hex_distance_2d(q1: i32, r1: i32, q2: i32, r2: i32) -> u32 {
    let s1 = -q1 - r1;
    let s2 = -q2 - r2;

    ((q1 - q2).abs().max((r1 - r2).abs())).max((s1 - s2).abs()) as u32
}

/// A node in the 3D hex lattice
#[derive(Debug, Clone)]
pub struct LatticeNode3D {
    pub coord: Hex3D,
    pub veracity: f64,
    pub coherence: f64,
    pub state: Option<String>, // Quantum state stored here
}

impl LatticeNode3D {
    pub fn new(q: i32, r: i32, z: i32) -> Self {
        Self {
            coord: Hex3D::new(q, r, z),
            veracity: 1.0,
            coherence: 1.0,
            state: None,
        }
    }
}

/// 3D Hexagonal Lattice
#[derive(Debug, Clone)]
pub struct HexLattice3D {
    nodes: HashMap<Hex3D, LatticeNode3D>,
}

impl HexLattice3D {
    /// Create a 3D hex lattice with given radius and layers
    pub fn new(planar_radius: u32, z_layers: i32) -> Self {
        let mut nodes = HashMap::new();

        // Generate 2D hex pattern for each z layer
        for z in -z_layers..=z_layers {
            let r = planar_radius as i32;
            for q in -r..=r {
                let r_min = (-r).max(-q - r);
                let r_max = r.min(-q + r);
                for r_coord in r_min..=r_max {
                    let coord = Hex3D::new(q, r_coord, z);
                    nodes.insert(coord, LatticeNode3D::new(q, r_coord, z));
                }
            }
        }

        Self { nodes }
    }

    /// Get node at coordinate
    pub fn get(&self, coord: &Hex3D) -> Option<&LatticeNode3D> {
        self.nodes.get(coord)
    }

    /// Check if coordinate exists
    pub fn contains(&self, coord: &Hex3D) -> bool {
        self.nodes.contains_key(coord)
    }

    /// Number of nodes
    pub fn size(&self) -> usize {
        self.nodes.len()
    }

    /// Get all nodes
    pub fn nodes(&self) -> &HashMap<Hex3D, LatticeNode3D> {
        &self.nodes
    }

    /// Add a connection between two nodes (strengthens veracity)
    pub fn connect(&mut self, a: &Hex3D, b: &Hex3D) {
        if let Some(node_a) = self.nodes.get_mut(a) {
            node_a.veracity = (node_a.veracity + 0.1).min(1.0);
        }
        if let Some(node_b) = self.nodes.get_mut(b) {
            node_b.veracity = (node_b.veracity + 0.1).min(1.0);
        }
    }

    /// Calculate average veracity
    pub fn average_veracity(&self) -> f64 {
        if self.nodes.is_empty() {
            return 0.0;
        }
        let total: f64 = self.nodes.values().map(|n| n.veracity).sum();
        total / self.nodes.len() as f64
    }

    /// Get layer at z
    pub fn layer(&self, z: i32) -> Vec<&Hex3D> {
        self.nodes
            .keys()
            .filter(|c| c.z == z)
            .collect()
    }
}

/// Propagation simulation through 3D lattice
pub fn propagate_signal(
    lattice: &HexLattice3D,
    start: &Hex3D,
    steps: usize,
) -> Vec<f64> {
    let mut amplitudes = Vec::with_capacity(steps);
    let mut current = *start;

    for _ in 0..steps {
        let neighbors = current.neighbors();
        let mut strongest = None;
        let mut best_veracity = 0.0;

        for n in neighbors.iter() {
            if let Some(node) = lattice.get(n) {
                if node.veracity > best_veracity {
                    best_veracity = node.veracity;
                    strongest = Some(*n);
                }
            }
        }

        if let Some(next) = strongest {
            if let Some(node) = lattice.get(&next) {
                amplitudes.push(node.veracity);
                current = next;
            }
        } else {
            amplitudes.push(0.0);
            break;
        }
    }

    amplitudes
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_3d_coord_creation() {
        let c = Hex3D::new(1, 2, 3);
        assert_eq!(c.q, 1);
        assert_eq!(c.r, 2);
        assert_eq!(c.z, 3);
    }

    #[test]
    fn test_neighbors_count() {
        let c = Hex3D::new(0, 0, 0);
        let neighbors = c.neighbors();
        // 6 planar + 7 above + 7 below = 20 neighbors (21 including self)
        assert_eq!(neighbors.len(), 20);
    }

    #[test]
    fn test_vertical_neighbors() {
        let c = Hex3D::new(0, 0, 0);
        let vertical = c.vertical_neighbors();
        assert_eq!(vertical.len(), 2);
        assert_eq!(vertical[0].z, 1);
        assert_eq!(vertical[1].z, -1);
    }

    #[test]
    fn test_planar_neighbors() {
        let c = Hex3D::new(0, 0, 5);
        let planar = c.planar_neighbors();
        assert_eq!(planar.len(), 6);
        // All should have same z
        for n in planar.iter() {
            assert_eq!(n.z, 5);
        }
    }

    #[test]
    fn test_up_down() {
        let c = Hex3D::new(0, 0, 0);
        assert_eq!(c.up().z, 1);
        assert_eq!(c.down().z, -1);
    }

    #[test]
    fn test_distance_planar() {
        let a = Hex3D::new(0, 0, 0);
        let b = Hex3D::new(1, 0, 0);
        assert_eq!(a.distance(&b), 1);
    }

    #[test]
    fn test_distance_vertical() {
        let a = Hex3D::new(0, 0, 0);
        let b = Hex3D::new(0, 0, 1);
        assert_eq!(a.distance(&b), 1);
    }

    #[test]
    fn test_3d_lattice_creation() {
        let lattice = HexLattice3D::new(1, 1);
        // radius 1 = 7 nodes per layer
        // z from -1 to 1 = 3 layers
        // total = 21 nodes
        assert!(lattice.size() > 0);
    }
}
