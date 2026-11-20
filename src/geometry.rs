//! Finding the atoms near one another: close pairs among positions, and contacts across the
//! edges of a periodic cell.

use std::collections::HashMap;

use cgmath::{InnerSpace, Matrix, Matrix3, Point3, SquareMatrix, Vector3};

use crate::Cell;

/// The pairs `(i, j, distance)` of `positions` with `i < j` and at most `cutoff` apart, ordered by
/// `i`, then `j`.
///
/// The positions are sorted into cubic cells with an edge of `cutoff`, and each is compared only
/// with those of its own and the 26 adjacent cells, so the work grows with the number of
/// positions and how many lie within `cutoff` of each, not with the number of all pairs.
pub fn pairs_within(positions: &[[f32; 3]], cutoff: f32) -> Vec<(usize, usize, f32)> {
    if !(cutoff.is_finite() && cutoff > 0.0) {
        return Vec::new();
    }

    let cell_of = |p: &[f32; 3]| p.map(|c| (c / cutoff).floor() as i64);

    let mut cells = HashMap::<[i64; 3], Vec<usize>>::new();
    for (index, position) in positions.iter().enumerate() {
        cells.entry(cell_of(position)).or_default().push(index);
    }

    let mut pairs = Vec::new();
    for (i, position) in positions.iter().enumerate() {
        let [x, y, z] = cell_of(position);
        for cx in x - 1..=x + 1 {
            for cy in y - 1..=y + 1 {
                for cz in z - 1..=z + 1 {
                    let Some(others) = cells.get(&[cx, cy, cz]) else {
                        continue;
                    };
                    for &j in others.iter().filter(|&&j| j > i) {
                        let other = positions[j];
                        let distance = ((position[0] - other[0]).powi(2)
                            + (position[1] - other[1]).powi(2)
                            + (position[2] - other[2]).powi(2))
                        .sqrt();
                        if distance <= cutoff {
                            pairs.push((i, j, distance));
                        }
                    }
                }
            }
        }
    }

    pairs.sort_unstable_by_key(|&(i, j, _)| (i, j));
    pairs
}

/// Every distance shorter than `cutoff` (Å) from each of the sites at `fractional` (coordinates
/// in `[0, 1)` of `cell`) to every atom of the periodic structure but itself: the site's index,
/// the other atom's site's index and the distance (Å). Each contact between two sites is listed
/// from both, so twice. `None` for a cell without volume.
///
/// The sites' copies within `cutoff` of the cell are searched as [`pairs_within`] searches, so
/// the work grows with the number of sites, not with its square.
pub fn periodic_contacts(
    cell: &Cell,
    fractional: &[Point3<f64>],
    cutoff: f32,
) -> Option<Vec<(usize, usize, f32)>> {
    let axis = |fractional: [f64; 3]| cell.cartesian(Vector3::from(fractional)).cast::<f64>();
    let lattice = Matrix3::from_cols(
        axis([1.0, 0.0, 0.0])?,
        axis([0.0, 1.0, 0.0])?,
        axis([0.0, 0.0, 1.0])?,
    );
    // How far (fractional) along each axis the cutoff reaches: the cutoff over the spacing of
    // the lattice planes across it.
    let inverse = lattice.invert()?;
    let reach = [0, 1, 2].map(|axis| f64::from(cutoff) * inverse.row(axis).magnitude());
    let shifts = |axis: usize| {
        let reach = reach[axis].ceil() as i32 + 1;
        -reach..=reach
    };

    // Each site's copies whose fractional coordinates lie within the reach of the cell, with
    // the site and whether the copy is the site itself.
    let mut positions = Vec::new();
    let mut copies = Vec::new();
    for (site, point) in fractional.iter().enumerate() {
        let point = Vector3::new(point.x, point.y, point.z);
        for a in shifts(0) {
            for b in shifts(1) {
                for c in shifts(2) {
                    let copy = point + Vector3::new(f64::from(a), f64::from(b), f64::from(c));
                    let near = (0..3)
                        .all(|axis| -reach[axis] <= copy[axis] && copy[axis] < 1.0 + reach[axis]);
                    if near {
                        let position = lattice * copy;
                        positions.push([position.x as f32, position.y as f32, position.z as f32]);
                        copies.push((site, (a, b, c) == (0, 0, 0)));
                    }
                }
            }
        }
    }

    let mut contacts = Vec::new();
    for (p, q, distance) in pairs_within(&positions, cutoff) {
        if distance >= cutoff {
            continue;
        }
        let ((site_p, home_p), (site_q, home_q)) = (copies[p], copies[q]);
        if home_p {
            contacts.push((site_p, site_q, distance));
        }
        if home_q {
            contacts.push((site_q, site_p, distance));
        }
    }
    Some(contacts)
}

/// `count` positions spread over a cube of edge `extent` around the origin, the same for every
/// call (a linear congruential sequence), plus a few on cell borders for a cutoff of 1.
#[cfg(test)]
pub(crate) fn scattered_positions(count: usize, extent: f32) -> Vec<[f32; 3]> {
    let mut state: u64 = 0x2545_f491_4f6c_dd1d;
    let mut next = || {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((state >> 40) as f32 / (1u64 << 24) as f32 - 0.5) * extent
    };

    let mut positions = (0..count)
        .map(|_| [next(), next(), next()])
        .collect::<Vec<_>>();
    positions.extend([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 1.0],
    ]);
    positions
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `count` values spread over `[-extent / 2, extent / 2)`, the same for every call (a linear
    /// congruential sequence).
    fn scattered(count: usize, extent: f64) -> Vec<f64> {
        let mut state: u64 = 0x2545_f491_4f6c_dd1d;
        (0..count)
            .map(|_| {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                ((state >> 40) as f64 / (1u64 << 24) as f64 - 0.5) * extent
            })
            .collect()
    }

    #[test]
    fn the_cell_search_finds_exactly_the_pairs_all_pairs_would() {
        let values = scattered(1200, 12.0);
        let mut positions: Vec<[f32; 3]> = values
            .chunks(3)
            .map(|c| [c[0] as f32, c[1] as f32, c[2] as f32])
            .collect();
        positions.extend([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 1.0],
        ]);

        for cutoff in [0.5, 1.0, 2.5] {
            let found: Vec<_> = pairs_within(&positions, cutoff)
                .into_iter()
                .map(|(i, j, _)| (i, j))
                .collect();

            let mut expected = Vec::new();
            for i in 0..positions.len() {
                for j in i + 1..positions.len() {
                    let [a, b] = [positions[i], positions[j]];
                    let distance =
                        ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2))
                            .sqrt();
                    if distance <= cutoff {
                        expected.push((i, j));
                    }
                }
            }
            assert_eq!(found, expected, "cutoff {cutoff}");
        }
    }

    #[test]
    fn no_pairs_without_a_positive_cutoff() {
        let positions = [[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]];

        assert!(pairs_within(&positions, 0.0).is_empty());
        assert!(pairs_within(&positions, f32::NAN).is_empty());
        assert!(pairs_within(&[], 1.0).is_empty());
    }

    #[test]
    fn periodic_contacts_are_those_every_copy_in_reach_gives() {
        // A triclinic cell, sites scattered over it and some on its faces and corners.
        let mut cell = Cell::default();
        (cell.a, cell.b, cell.c) = (5.1.into(), 6.3.into(), 4.4.into());
        (cell.alpha, cell.beta, cell.gamma) = (78.0.into(), 102.0.into(), 95.0.into());
        let values = scattered(90, 1.0);
        let mut fractional: Vec<Point3<f64>> = values
            .chunks(3)
            .map(|c| Point3::new(c[0] + 0.5, c[1] + 0.5, c[2] + 0.5))
            .collect();
        fractional.extend([Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.5, 0.5)]);

        let lattice = Matrix3::from_cols(
            cell.cartesian(Vector3::new(1.0, 0.0, 0.0))
                .cast::<f64>()
                .unwrap(),
            cell.cartesian(Vector3::new(0.0, 1.0, 0.0))
                .cast::<f64>()
                .unwrap(),
            cell.cartesian(Vector3::new(0.0, 0.0, 1.0))
                .cast::<f64>()
                .unwrap(),
        );
        let inverse = lattice.invert().unwrap();

        for cutoff in [1.0, 2.5, 4.0, 7.0] {
            let mut found = periodic_contacts(&cell, &fractional, cutoff).unwrap();

            // Every site against every copy of every site in the cells the cutoff reaches.
            let reach = [0, 1, 2]
                .map(|axis| (f64::from(cutoff) * inverse.row(axis).magnitude()).ceil() as i32);
            let mut expected = Vec::new();
            for (from, start) in fractional.iter().enumerate() {
                for (to, end) in fractional.iter().enumerate() {
                    let difference = (end - start).map(|d| d - d.round());
                    for a in -reach[0]..=reach[0] {
                        for b in -reach[1]..=reach[1] {
                            for c in -reach[2]..=reach[2] {
                                if from == to && (a, b, c) == (0, 0, 0) {
                                    continue;
                                }
                                let shift = Vector3::new(f64::from(a), f64::from(b), f64::from(c));
                                let distance = (lattice * (difference + shift)).magnitude() as f32;
                                if distance < cutoff {
                                    expected.push((from, to, distance));
                                }
                            }
                        }
                    }
                }
            }

            let order = |a: &(usize, usize, f32), b: &(usize, usize, f32)| {
                (a.0, a.1).cmp(&(b.0, b.1)).then(a.2.total_cmp(&b.2))
            };
            found.sort_by(order);
            expected.sort_by(order);
            assert_eq!(found.len(), expected.len(), "cutoff {cutoff}");
            for (found, expected) in found.iter().zip(&expected) {
                assert_eq!(
                    (found.0, found.1),
                    (expected.0, expected.1),
                    "cutoff {cutoff}"
                );
                assert!(
                    (found.2 - expected.2).abs() < 1e-3,
                    "{found:?} {expected:?}"
                );
            }
        }
    }
}
