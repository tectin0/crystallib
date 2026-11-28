//! The atoms of a phase's unit cell: the copies of its sites that the space group's operations
//! make, with how often each one's place is occupied.

use cgmath::{EuclideanSpace, InnerSpace, Point3};

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

/// The asymmetric unit as `phase` lists it: for each of its sites, the one of `atoms` (its cell
/// atoms, see [`cell_atoms`]) that is the copy of the site closest to the coordinates it states,
/// and the lattice translation (cells along a, b and c) from that atom to them, in the order of
/// the sites. Where a site's copies were merged (see [`site_positions`]), the closest copy lies
/// off the stated coordinates by up to the merged distance. A site without cell atoms has none.
pub fn asymmetric_unit(phase: &Phase, atoms: &[CellAtom]) -> Vec<(usize, [i32; 3])> {
    let mut closest: Vec<Option<(usize, [i32; 3], f32)>> = vec![None; phase.atoms.len()];
    for (index, atom) in atoms.iter().enumerate() {
        let site = &phase.atoms[atom.site];
        let stated = Point3::new(site.x.value, site.y.value, site.z.value);
        let translation = (stated - atom.fractional).map(f64::round);
        let distance = phase
            .cell
            .cartesian(stated - (atom.fractional + translation))
            .magnitude();

        if closest[atom.site].is_none_or(|(_, _, nearest)| distance < nearest) {
            let translation = translation.map(|c| c as i32).into();
            closest[atom.site] = Some((index, translation, distance));
        }
    }
    closest
        .into_iter()
        .flatten()
        .map(|(atom, translation, _)| (atom, translation))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Atoms, Cell, Measured};

    /// Writes for each CIF in the file `$ASYMMETRIC_UNIT_SAMPLE` (one path per line) a row of
    /// `$ASYMMETRIC_UNIT_OUT`: its sites, the atoms of its asymmetric unit, its sites without
    /// one, those with more than one, those whose atom lies off the stated coordinates (beyond
    /// 0.0001 Å) and the largest such distance (Å).
    #[test]
    #[cfg(feature = "cif")]
    #[ignore = "needs CIFs to read: ASYMMETRIC_UNIT_SAMPLE, ASYMMETRIC_UNIT_OUT"]
    fn asymmetric_unit_on_cifs() {
        use std::io::Write;

        let sample = std::env::var("ASYMMETRIC_UNIT_SAMPLE").expect("names the CIFs");
        let out = std::env::var("ASYMMETRIC_UNIT_OUT").expect("names the CSV to write");
        let paths: Vec<String> = std::fs::read_to_string(sample)
            .unwrap()
            .lines()
            .map(str::trim)
            .filter(|line| !line.is_empty())
            .map(str::to_string)
            .collect();
        let mut out = std::io::BufWriter::new(std::fs::File::create(out).unwrap());
        writeln!(
            out,
            "path,outcome,sites,unit,without,twice,off,largest_offset"
        )
        .unwrap();

        for path in &paths {
            let row = std::fs::read(path).ok().and_then(|bytes| {
                let (_, phase) = crate::cif_phase::try_phase_from_cif_bytes(&bytes)?;
                let atoms = cell_atoms(&phase);
                let unit = asymmetric_unit(&phase, &atoms);

                let mut per_site = vec![0; phase.atoms.len()];
                let (mut off, mut largest) = (0, 0.0_f32);
                for &(atom, translation) in &unit {
                    let site = &phase.atoms[atoms[atom].site];
                    per_site[atoms[atom].site] += 1;
                    let shown =
                        atoms[atom].fractional + cgmath::Vector3::from(translation.map(f64::from));
                    let stated = Point3::new(site.x.value, site.y.value, site.z.value);
                    let offset = phase.cell.cartesian(shown - stated).magnitude();
                    if offset > 1e-4 {
                        off += 1;
                        largest = largest.max(offset);
                    }
                }
                let without = per_site.iter().filter(|&&count| count == 0).count();
                let twice = per_site.iter().filter(|&&count| count > 1).count();
                Some(format!(
                    "ok,{},{},{without},{twice},{off},{largest:.4}",
                    phase.atoms.len(),
                    unit.len()
                ))
            });
            let row = row.unwrap_or_else(|| "unreadable,,,,,,".to_string());
            writeln!(out, "{path},{row}").unwrap();
        }
    }

    /// A phase of the carbons `sites` (label and stated coordinates) in a 5 Å cubic cell of the
    /// space group numbered `number`.
    fn phase_of(number: u8, sites: &[(&str, [f64; 3])]) -> Phase {
        let cell = Cell {
            a: 5.0.into(),
            b: 5.0.into(),
            c: 5.0.into(),
            alpha: 90.0.into(),
            beta: 90.0.into(),
            gamma: 90.0.into(),
            space_group_number: number,
            ..Cell::default()
        };
        let atoms = sites
            .iter()
            .map(|&(label, [x, y, z])| Atom {
                label: label.to_string(),
                type_: "C".to_string(),
                x: x.into(),
                y: y.into(),
                z: z.into(),
                ..Atom::default()
            })
            .collect();
        Phase {
            cell,
            atoms: Atoms(atoms),
        }
    }

    #[test]
    fn the_asymmetric_unit_is_each_site_where_it_is_stated() {
        // P−1: a general site has two copies in the cell, one at an inversion centre one; C3 is
        // stated outside the cell.
        let phase = phase_of(
            2,
            &[
                ("C1", [0.1, 0.2, 0.3]),
                ("C2", [0.5, 0.5, 0.5]),
                ("C3", [-0.1, 0.2, 1.3]),
            ],
        );
        let atoms = cell_atoms(&phase);
        assert_eq!(atoms.len(), 5);

        let unit = asymmetric_unit(&phase, &atoms);

        assert_eq!(unit.len(), 3);
        for (site, &(atom, translation)) in unit.iter().enumerate() {
            assert_eq!(atoms[atom].site, site);
            let shown = atoms[atom].fractional + cgmath::Vector3::from(translation.map(f64::from));
            let stated = &phase.atoms[site];
            let error =
                (shown - Point3::new(stated.x.value, stated.y.value, stated.z.value)).map(f64::abs);
            assert!(
                error.x.max(error.y).max(error.z) < 1e-9,
                "{site}: {shown:?}"
            );
        }
        assert_eq!(unit[2].1, [-1, 0, 1]);
    }

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
