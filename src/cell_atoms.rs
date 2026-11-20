//! The atoms of a phase's unit cell: the copies of its sites that the space group's operations
//! make, with how often each one's place is occupied.

use cgmath::{EuclideanSpace, Point3};

use crate::{Atom, Phase, SpaceGroupSymmetryOperations, occupancy};

/// An atom of a phase's unit cell: a copy of one of its sites.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub struct CellAtom {
    /// The site it is a copy of, as its index in the phase's atoms.
    pub site: usize,
    /// Its fractional coordinates, in `[0, 1)`.
    pub fractional: Point3<f64>,
    /// How often its place is occupied, by it or by the atoms of other sites sharing its place:
    /// the sum of their occupancies (see [`occupancy::places`]).
    pub place_occupancy: f32,
}

/// The distinct positions (fractional, in `[0, 1)`) of `atom` in the unit cell of `phase`, as
/// many as the atom's stated multiplicity where its copies are close enough to merge (see
/// [`SpaceGroupSymmetryOperations::site_positions_in_cell`]).
pub fn site_positions(phase: &Phase, atom: &Atom) -> Option<Vec<Point3<f64>>> {
    let point = Point3::new(atom.x.value, atom.y.value, atom.z.value);
    let multiplicity = atom
        .multiplicity
        .filter(|multiplicity| *multiplicity >= 1.0)
        .map(|multiplicity| multiplicity.round() as usize);

    SpaceGroupSymmetryOperations::site_positions_in_cell(&phase.cell, point, multiplicity)
}

/// The atoms of `phase`'s unit cell, site by site. A site whose positions can't be made (its
/// space group unknown) has none, and a warning is logged.
pub fn cell_atoms(phase: &Phase) -> Vec<CellAtom> {
    let mut atoms = Vec::new();
    for (site, atom) in phase.atoms.iter().enumerate() {
        let Some(positions) = site_positions(phase, atom) else {
            log::warn!(
                "Failed to generate the symmetry positions of {} at ({}, {}, {})",
                atom.label,
                atom.x.value,
                atom.y.value,
                atom.z.value
            );
            continue;
        };
        atoms.extend(positions.into_iter().map(|fractional| CellAtom {
            site,
            fractional,
            place_occupancy: 0.0,
        }));
    }

    let positions: Vec<[f32; 3]> = atoms
        .iter()
        .map(|atom| phase.cell.cartesian(atom.fractional.to_vec()).into())
        .collect();
    let occupancies: Vec<f32> = atoms
        .iter()
        .map(|atom| phase.atoms[atom.site].occupancy.value as f32)
        .collect();
    for place in occupancy::places(&positions, &occupancies) {
        let occupancy = place.iter().map(|&atom| occupancies[atom]).sum();
        for atom in place {
            atoms[atom].place_occupancy = occupancy;
        }
    }
    atoms
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Atoms, Cell, Measured};

    #[test]
    fn a_place_is_occupied_by_all_the_atoms_sharing_it() {
        let site = |label: &str, [x, y, z]: [f64; 3], occupancy: f64| Atom {
            label: label.to_string(),
            type_: label.trim_end_matches(char::is_numeric).to_string(),
            x: x.into(),
            y: y.into(),
            z: z.into(),
            occupancy: Measured::new(occupancy),
            ..Atom::default()
        };
        let phase = Phase {
            cell: Cell {
                a: 5.0.into(),
                b: 5.0.into(),
                c: 5.0.into(),
                alpha: 90.0.into(),
                beta: 90.0.into(),
                gamma: 90.0.into(),
                space_group_number: 1,
                ..Cell::default()
            },
            atoms: Atoms(vec![
                site("Fe1", [0.0, 0.0, 0.0], 0.6),
                site("Mg1", [0.0, 0.0, 0.0], 0.4),
                site("O1", [0.5, 0.5, 0.5], 1.0),
                site("Na1", [0.25, 0.25, 0.25], 0.5),
            ]),
        };

        let atoms = cell_atoms(&phase);
        let place_occupancy = |site: usize| {
            atoms
                .iter()
                .find(|atom| atom.site == site)
                .unwrap()
                .place_occupancy
        };
        assert_eq!(atoms.len(), 4);
        assert!(
            (place_occupancy(0) - 1.0).abs() < 1e-6,
            "Fe and Mg share it"
        );
        assert!((place_occupancy(1) - 1.0).abs() < 1e-6);
        assert_eq!(place_occupancy(2), 1.0);
        assert_eq!(place_occupancy(3), 0.5);
    }
}
