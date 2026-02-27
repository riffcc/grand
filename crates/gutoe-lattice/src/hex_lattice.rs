/*!
 * GUTOE Lattice - Hexagonal Coordinate System
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

//! Hexagonal lattice using cubic coordinates
//!
//! From GUTOE.md:
//! "Uses cubic coordinates (q, r, s) where q+r+s=0"
//! "O(1) distance calculation: max(|q1-q2|, |r1-r2|, |s1-s2|)"
//! "6-fold symmetry with 60-degree rotations"

use std::collections::HashMap;
use std::fmt;

/// Cubic coordinates for hexagonal lattice
/// q + r + s = 0 invariant must hold
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct HexCoord {
    pub q: i32,
    pub r: i32,
    pub s: i32,
}

impl HexCoord {
    /// Create new coordinates, validating q+r+s=0
    pub fn new(q: i32, r: i32) -> Self {
        let s = -q - r;
        Self { q, r, s }
    }

    /// Distance to another hex (O(1) calculation)
    /// From GUTOE.md: max(|q1-q2|, |r1-r2|, |s1-s2|)
    pub fn distance(&self, other: &HexCoord) -> u32 {
        ((self.q - other.q).abs().max((self.r - other.r).abs())).max((self.s - other.s).abs())
            as u32
    }

    /// Get all 6 neighbors (6-fold symmetry)
    pub fn neighbors(&self) -> [HexCoord; 6] {
        [
            HexCoord::new(self.q + 1, self.r),
            HexCoord::new(self.q + 1, self.r - 1),
            HexCoord::new(self.q, self.r - 1),
            HexCoord::new(self.q - 1, self.r),
            HexCoord::new(self.q - 1, self.r + 1),
            HexCoord::new(self.q, self.r + 1),
        ]
    }

    /// Rotate by 60 degrees (π/3 radians)
    pub fn rotate60(&self) -> HexCoord {
        HexCoord::new(-self.r, -self.s)
    }

    /// Rotate by 120 degrees
    pub fn rotate120(&self) -> HexCoord {
        HexCoord::new(-self.s, -self.q)
    }

    /// Scale coordinates
    pub fn scale(&self, factor: i32) -> HexCoord {
        HexCoord::new(self.q * factor, self.r * factor)
    }
}

impl fmt::Display for HexCoord {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({}, {}, {})", self.q, self.r, self.s)
    }
}

/// A hexagonal lattice node that can hold quantum state
#[derive(Debug, Clone)]
pub struct LatticeNode {
    pub coord: HexCoord,
    pub veracity: f64,  // Connection strength (0.0 to 1.0)
    pub coherence: f64, // Quantum coherence
}

impl LatticeNode {
    pub fn new(q: i32, r: i32) -> Self {
        Self {
            coord: HexCoord::new(q, r),
            veracity: 1.0,
            coherence: 1.0,
        }
    }

    pub fn with_veracity(mut self, v: f64) -> Self {
        self.veracity = v.clamp(0.0, 1.0);
        self
    }

    pub fn with_coherence(mut self, c: f64) -> Self {
        self.coherence = c.clamp(0.0, 1.0);
        self
    }
}

/// The hexagonal lattice structure
#[derive(Debug, Clone)]
pub struct HexLattice {
    nodes: HashMap<HexCoord, LatticeNode>,
    radius: u32,
}

impl HexLattice {
    /// Create a hexagonal lattice of given radius
    /// Uses axial coordinates in a hexagonal pattern
    pub fn new(radius: u32) -> Self {
        let mut nodes = HashMap::new();
        let r = radius as i32;

        // Generate hexagonal ring pattern
        for q in -r..=r {
            let r_min = (-r).max(-q - r);
            let r_max = r.min(-q + r);
            for r_val in r_min..=r_max {
                let coord = HexCoord::new(q, r_val);
                nodes.insert(coord, LatticeNode::new(q, r_val));
            }
        }

        Self { nodes, radius }
    }

    /// Get a node at coordinates
    pub fn get(&self, coord: &HexCoord) -> Option<&LatticeNode> {
        self.nodes.get(coord)
    }

    /// Get a mutable node
    pub fn get_mut(&mut self, coord: &HexCoord) -> Option<&mut LatticeNode> {
        self.nodes.get_mut(coord)
    }

    /// Add a node at coordinates
    pub fn add_node(&mut self, coord: HexCoord) {
        self.nodes
            .entry(coord)
            .or_insert_with(|| LatticeNode::new(coord.q, coord.r));
    }

    /// Number of nodes in lattice
    pub fn size(&self) -> usize {
        self.nodes.len()
    }

    /// Get all nodes
    pub fn nodes(&self) -> &HashMap<HexCoord, LatticeNode> {
        &self.nodes
    }

    /// Check if coordinate exists
    pub fn contains(&self, coord: &HexCoord) -> bool {
        self.nodes.contains_key(coord)
    }

    /// All neighbors of a coordinate that exist in lattice
    pub fn existing_neighbors(&self, coord: &HexCoord) -> Vec<HexCoord> {
        coord
            .neighbors()
            .into_iter()
            .filter(|c| self.contains(c))
            .collect()
    }

    /// Calculate clustering coefficient
    /// From GUTOE.md experiments: "clustering coefficients of 0.72±0.05"
    pub fn clustering_coefficient(&self) -> f64 {
        let mut total = 0.0;
        let mut count = 0;

        for (coord, _) in &self.nodes {
            let neighbors = self.existing_neighbors(coord);
            if neighbors.len() < 2 {
                continue;
            }

            let mut edges = 0;
            for n1 in &neighbors {
                for n2 in &neighbors {
                    if n1 != n2 && self.contains(n1) && self.contains(n2) {
                        // Check if n1 and n2 are connected
                        if n1.neighbors().contains(n2) {
                            edges += 1;
                        }
                    }
                }
            }

            let possible = neighbors.len() * (neighbors.len() - 1);
            if possible > 0 {
                total += edges as f64 / possible as f64;
                count += 1;
            }
        }

        if count > 0 {
            total / count as f64
        } else {
            0.0
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hex_coord_creation() {
        let c = HexCoord::new(1, 2);
        assert_eq!(c.q, 1);
        assert_eq!(c.r, 2);
        assert_eq!(c.s, -3);
    }

    #[test]
    fn test_distance_calculation() {
        let a = HexCoord::new(0, 0);
        let b = HexCoord::new(1, 0);
        assert_eq!(a.distance(&b), 1);
    }

    #[test]
    fn test_rotation() {
        let c = HexCoord::new(1, 0);
        let rotated = c.rotate60();
        // After 60° rotation: q' = -r, r' = -s = q+r
        assert_eq!(rotated.q, 0);
        assert_eq!(rotated.r, 1);
    }

    #[test]
    fn test_lattice_creation() {
        let lattice = HexLattice::new(2);
        // Hex lattice of radius 2 has 1 + 6 + 12 = 19 nodes
        assert_eq!(lattice.size(), 19);
    }

    #[test]
    fn test_neighbors() {
        let c = HexCoord::new(0, 0);
        let neighbors = c.neighbors();
        assert_eq!(neighbors.len(), 6);
    }

    // ── Clustering coefficient ────────────────────────────────────────────────

    #[test]
    fn interior_node_has_6_planar_neighbors_with_6_ring_edges() {
        // For an interior hex node (q=0, r=0) in a radius≥2 lattice, all 6
        // hex neighbors exist AND are adjacent to their two ring-neighbours.
        // Edges among 6 ring nodes = 6 (the ring itself), not more.
        let lattice = HexLattice::new(3);
        let center = HexCoord::new(0, 0);
        let neighbors = lattice.existing_neighbors(&center);
        assert_eq!(
            neighbors.len(),
            6,
            "interior node should have exactly 6 neighbors"
        );

        // Count edges among neighbors (unordered pairs)
        let mut edge_count = 0;
        for i in 0..neighbors.len() {
            for j in (i + 1)..neighbors.len() {
                if neighbors[i].neighbors().contains(&neighbors[j]) {
                    edge_count += 1;
                }
            }
        }
        // In a hex lattice ring, consecutive neighbors are adjacent: exactly 6 edges.
        assert_eq!(
            edge_count, 6,
            "interior hex node: {edge_count} edges among neighbors, expected 6 (ring)"
        );
    }

    #[test]
    fn clustering_coefficient_is_near_0_4_not_0_72() {
        // GUTOE claims clustering coefficient of "0.72 ± 0.05" for hex lattice.
        //
        // Reality: interior nodes have k=6 neighbors with 6 ring edges.
        //   C = 6 / (6·5/2) = 6/15 = 0.4
        //
        // The code counts ordered pairs (each edge contributes 2) and divides by
        // k·(k-1) = 30, giving the same result: 12/30 = 0.4.
        //
        // The claimed 0.72 is not consistent with hex lattice geometry.
        let lattice = HexLattice::new(5); // Large enough for mostly-interior nodes
        let cc = lattice.clustering_coefficient();

        assert!(
            (cc - 0.4).abs() < 0.1,
            "clustering coefficient = {cc:.4}, expected ~0.4 \
             (GUTOE claims 0.72±0.05, which contradicts hex lattice geometry)"
        );
        assert!(
            cc < 0.67,
            "clustering coefficient = {cc:.4} is near the claimed 0.72±0.05, \
             but interior hex nodes give C = 6/15 = 0.4"
        );
    }
}
