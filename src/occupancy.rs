//! What partial occupancy tells about which atoms are present together. A site of occupancy
//! below 1 holds its atom in only part of the cells; two such sites whose occupancies add to at
//! most 1 may never be occupied in the same cell, and they are certain not to be where they lie
//! too close for two atoms (alternative positions of one atom: a split site, or two elements
//! sharing a site).

#[cfg(feature = "data")]
use crate::COVALENT_RADII;

/// How far a sum of CIF occupancies may be off from what they stand for, as CIFs round them
/// (1/3 written as 0.333).
pub const ROUNDING: f32 = 1e-3;

/// Two atoms closer than this times the sum of their covalent radii are never both present: in
/// 4,968 ICSD structures, 247 of 170,119 atoms filling their sites have a neighbour that close
/// (in 24 structures, mostly hydrogen placed by X-rays), while the shortest real bonds (W=O,
/// U=O) come to 0.68 of it. Partially occupied sites that close are alternative positions of one
/// atom (a split or shared site, or orientations of a disordered molecule).
#[cfg(feature = "data")]
pub const ALTERNATIVE_POSITION_FACTOR: f32 = 0.6;

/// Whether sites of the occupancies `a` and `b` must at times be occupied together: their
/// occupancies add to more than 1, beyond [`ROUNDING`].
pub fn together(a: f32, b: f32) -> bool {
    a + b > 1.0 + ROUNDING
}

/// Whether atoms of the elements `element_a` and `element_b` with the occupancies `occupancies`,
/// `distance` (Å) apart, are alternative positions of one atom: they needn't be present together
/// (see [`together`]) and lie closer than [`ALTERNATIVE_POSITION_FACTOR`] times the sum of their
/// covalent radii.
#[cfg(feature = "data")]
pub fn alternatives(
    (element_a, element_b): (&str, &str),
    occupancies: (f32, f32),
    distance: f32,
) -> bool {
    let radii = COVALENT_RADII.get(element_a) + COVALENT_RADII.get(element_b);
    !together(occupancies.0, occupancies.1) && distance < ALTERNATIVE_POSITION_FACTOR * radii
}

/// Whether `occupancy` is short of a full site, beyond [`ROUNDING`].
pub fn is_partial(occupancy: f32) -> bool {
    occupancy < 1.0 - ROUNDING
}

/// Positions closer than this (Å) are one place: sites atoms share, or a site listed twice.
pub const SAME_PLACE: f32 = 1.0e-3;

/// The atoms at `positions` (Cartesian, Å) with `occupancies` grouped by place (see
/// [`SAME_PLACE`]), as indices, in the order of each place's first atom; the atoms of a place
/// most occupied first.
pub fn places(positions: &[[f32; 3]], occupancies: &[f32]) -> Vec<Vec<usize>> {
    // The pairs come ordered by their first index, so each atom joins the place of the first one
    // it shares a place with.
    let mut place_of = (0..positions.len()).collect::<Vec<_>>();
    for (i, j, _) in crate::geometry::pairs_within(positions, SAME_PLACE) {
        if place_of[j] == j {
            place_of[j] = place_of[i];
        }
    }

    let mut places: Vec<Vec<usize>> = Vec::new();
    let mut index_of = vec![usize::MAX; positions.len()];
    for (atom, &first) in place_of.iter().enumerate() {
        if first == atom {
            index_of[atom] = places.len();
            places.push(vec![atom]);
        } else {
            places[index_of[first]].push(atom);
        }
    }
    for atoms in &mut places {
        atoms.sort_by(|&a, &b| occupancies[b].total_cmp(&occupancies[a]));
    }
    places
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sites_are_together_only_where_their_occupancies_add_to_more_than_one() {
        assert!(together(1.0, 0.25));
        assert!(together(0.6, 0.6));
        assert!(!together(0.5, 0.5));
        assert!(!together(0.667, 0.333), "a third rounded");
    }

    #[test]
    #[cfg(feature = "data")]
    fn close_partial_sites_are_alternatives_and_farther_ones_or_full_ones_are_not() {
        let split = (("C", "C"), (0.5, 0.5));
        assert!(alternatives(split.0, split.1, 0.4));
        assert!(
            !alternatives(split.0, split.1, 1.5),
            "bonded within one orientation"
        );
        assert!(!alternatives(("C", "C"), (1.0, 0.5), 0.4));
        assert!(alternatives(("Fe", "Mg"), (0.6, 0.4), 0.0), "a shared site");
    }

    #[test]
    fn atoms_at_one_place_are_grouped_most_occupied_first() {
        let positions = [[0.0; 3], [1.0, 0.0, 0.0], [0.0; 3], [0.0, 0.0, 0.0005]];
        let occupancies = [0.2, 1.0, 0.5, 0.3];
        assert_eq!(places(&positions, &occupancies), [vec![2, 3, 0], vec![1]]);
    }
}
