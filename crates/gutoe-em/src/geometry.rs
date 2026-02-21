// GUTOE EM — Hex toroid geometry (intra-layer only, like Python mesh_neighbours)
// Copyright (C) 2026 Riff Labs, AGPL-3.0-or-later

use crate::config::LatticeConfig;

/// Decompose a flat site index into (row, col, layer).
pub fn site_coords(site: usize, cfg: &LatticeConfig) -> (usize, usize, usize) {
    let layer_sz = cfg.hex_rows * cfg.hex_cols;
    let z = site / layer_sz;
    let rem = site % layer_sz;
    (rem / cfg.hex_cols, rem % cfg.hex_cols, z)
}

fn flat_idx(r: usize, c: usize, z: usize, cfg: &LatticeConfig) -> usize {
    (z * cfg.hex_rows + r) * cfg.hex_cols + c
}

/// Six neighbours in the hex grid, wrapping on rows and cols.
/// Identical to Python `hex_neighbours` + `mesh_neighbours` — intra-layer only.
pub fn mesh_neighbours(r: usize, c: usize, z: usize, cfg: &LatticeConfig) -> Vec<usize> {
    let offsets: &[(i32, i32)] = if r % 2 == 0 {
        &[(-1, 0), (-1, 1), (0, -1), (0, 1), (1, 0), (1, 1)]
    } else {
        &[(-1, -1), (-1, 0), (0, -1), (0, 1), (1, -1), (1, 0)]
    };

    let rows = cfg.hex_rows as i32;
    let cols = cfg.hex_cols as i32;

    offsets
        .iter()
        .map(|&(dr, dc)| {
            let nr = ((r as i32 + dr).rem_euclid(rows)) as usize;
            let nc = ((c as i32 + dc).rem_euclid(cols)) as usize;
            flat_idx(nr, nc, z, cfg)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::LatticeConfig;

    #[test]
    fn site_coords_round_trip() {
        let cfg = LatticeConfig {
            hex_rows: 8,
            hex_cols: 8,
            layers: 1,
            ..Default::default()
        };
        let n = cfg.n_sites();
        for site in 0..n {
            let (r, c, z) = site_coords(site, &cfg);
            assert_eq!(flat_idx(r, c, z, &cfg), site, "round-trip failed for site {site}");
        }
    }

    #[test]
    fn each_site_has_six_neighbours() {
        let cfg = LatticeConfig {
            hex_rows: 8,
            hex_cols: 8,
            layers: 1,
            ..Default::default()
        };
        let n = cfg.n_sites();
        for site in 0..n {
            let (r, c, z) = site_coords(site, &cfg);
            let nbrs = mesh_neighbours(r, c, z, &cfg);
            assert_eq!(nbrs.len(), 6, "site {site} has {} neighbours, expected 6", nbrs.len());
        }
    }

    #[test]
    fn neighbours_are_symmetric() {
        let cfg = LatticeConfig {
            hex_rows: 8,
            hex_cols: 8,
            layers: 1,
            ..Default::default()
        };
        let n = cfg.n_sites();
        for site in 0..n {
            let (r, c, z) = site_coords(site, &cfg);
            let nbrs = mesh_neighbours(r, c, z, &cfg);
            for nb in nbrs {
                let (nr, nc, nz) = site_coords(nb, &cfg);
                let back = mesh_neighbours(nr, nc, nz, &cfg);
                assert!(
                    back.contains(&site),
                    "site {nb} does not list {site} as neighbour (asymmetry)"
                );
            }
        }
    }

    #[test]
    fn neighbours_are_intra_layer() {
        let cfg = LatticeConfig {
            hex_rows: 8,
            hex_cols: 8,
            layers: 3,
            ..Default::default()
        };
        let n = cfg.n_sites();
        for site in 0..n {
            let (r, c, z) = site_coords(site, &cfg);
            for nb in mesh_neighbours(r, c, z, &cfg) {
                let (_, _, zn) = site_coords(nb, &cfg);
                assert_eq!(z, zn, "neighbour of site {site} crosses layer boundary: z={z}, zn={zn}");
            }
        }
    }
}
