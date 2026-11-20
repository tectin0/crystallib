//! Interactions between atoms the bond rule doesn't bond: hydrogen bonds and short contacts.
//!
//! Both are found geometrically, as Mercury finds them: atoms closer than the sum of their van
//! der Waals radii (crystallib's `VAN_DER_WAALS_RADII`) that aren't bonded and are more than
//! [`FEWEST_BONDS_APART`] bonds apart. A hydrogen bond is such a contact from a hydrogen bonded
//! to a donor to an acceptor, at a donor–H···acceptor angle of at least [`SMALLEST_ANGLE`].
//! Structures that state no hydrogen positions have no hydrogen bonds.
//!
//! A short contact (van der Waals or dipole–dipole) is any other such contact between two
//! non-metals of different fragments (molecules, layers, chains: atoms not joined by bonds).
//! Within a framework, the oxygens of neighbouring polyhedra come closer than their van der
//! Waals radii allow, and the radii of metals are large enough to put most cations of an ionic
//! structure within them; neither is a contact.
//!
//! Two atoms are in contact only where they must at times be present together, their places'
//! occupancies (by any atom) adding to more than 1 (see [`crate::occupancy::together`]): the
//! orientations of a disordered molecule (each copy on sites a fraction occupied) touch none of
//! the others. A disorder group (see [`crate::disorder`]) tells which atoms are present together
//! only within one piece of the structure, while neighbouring pieces take their configurations
//! independently, so it tells nothing of contacts between pieces.

use crate::{ELECTRONEGATIVITIES, VAN_DER_WAALS_RADII};
use cgmath::{InnerSpace, vec3};
use std::collections::{HashMap, HashSet};

use crate::bonding::BondAtom;
use crate::elements::{element_of, is_metal};
use crate::type_symbol::state_of;

/// The elements a hydrogen bonded to one of them gives a hydrogen bond from.
pub const DONORS: [&str; 4] = ["N", "O", "F", "S"];

/// The elements a hydrogen bond goes to, each with its valence electrons. An acceptor needs a
/// lone pair (IUPAC, Pure Appl. Chem. 83, 1637, 2011): its valence electrons less its oxidation
/// state are at least two. The centre of an oxyanion (S⁶⁺ in sulfate, N⁵⁺ in nitrate, Cl⁷⁺ in
/// perchlorate) has none and accepts no hydrogen bond, its oxygens do; N⁻ in hydroxylamine or
/// an oxime has one. An atom whose type symbol names no oxidation state counts as having none
/// where it is bonded to a more electronegative atom.
pub const ACCEPTORS: [(&str, f64); 7] = [
    ("N", 5.0),
    ("O", 6.0),
    ("F", 7.0),
    ("S", 6.0),
    ("Cl", 7.0),
    ("Br", 7.0),
    ("I", 7.0),
];

/// The smallest donor–H···acceptor angle (degrees) of a hydrogen bond.
pub const SMALLEST_ANGLE: f32 = 120.0;

/// The most bonds two atoms may be apart and still never be in contact: atoms of one octahedron
/// (two bonds apart through its centre) or of one small group stay out.
pub const FEWEST_BONDS_APART: usize = 3;

/// What kind of interaction a contact is.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ContactKind {
    /// From a hydrogen bonded to a donor to an acceptor.
    HydrogenBond,
    /// Any other contact closer than the van der Waals radii allow.
    Short,
}

/// Two atoms in contact.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Contact {
    /// The two atoms (indices into the atoms given): for a hydrogen bond, the hydrogen first.
    pub atoms: (usize, usize),
    /// Their distance (Å).
    pub distance: f32,
    /// What kind of contact it is.
    pub kind: ContactKind,
}

/// The contacts between the first `inner` of the atoms at `positions` (Cartesian, Å) whose
/// `bonds` (index pairs with their lengths, each listed once) the bond rule found in any of their
/// configurations (see [`crate::bonding::BondRule::bonds_of_any_configuration`]). The atoms past
/// `inner` only join others into fragments and count towards how many bonds apart they are, so
/// that atoms joined through atoms past them aren't taken for separate fragments; they should
/// surround the first `inner` by at least one cell.
pub fn contacts(
    positions: &[[f32; 3]],
    atoms: &[BondAtom],
    bonds: &[(usize, usize, f32)],
    inner: usize,
) -> Vec<Contact> {
    let elements: Vec<String> = atoms
        .iter()
        .map(|atom| element_of(atom.type_symbol))
        .collect();
    let radii: Vec<f32> = elements
        .iter()
        .map(|element| *VAN_DER_WAALS_RADII.get(element))
        .collect();

    let mut neighbours: Vec<Vec<usize>> = vec![Vec::new(); atoms.len()];
    for &(i, j, _) in bonds {
        neighbours[i].push(j);
        neighbours[j].push(i);
    }

    let position = |atom: usize| {
        let [x, y, z] = positions[atom];
        vec3(x, y, z)
    };
    let largest = radii.iter().copied().fold(0.0, f32::max);
    let inner = inner.min(positions.len());
    let touching = |i: usize, j: usize, distance: f32| {
        distance < radii[i] + radii[j]
            && crate::occupancy::together(atoms[i].place_occupancy, atoms[j].place_occupancy)
    };
    let close: Vec<(usize, usize, f32)> =
        crate::geometry::pairs_within(&positions[..inner], 2.0 * largest)
            .into_iter()
            .filter(|&(i, j, distance)| touching(i, j, distance))
            .collect();

    // Hydrogen bonds first, among all the atoms given: they join their atoms as bonds do, so that
    // a donor and its acceptor, or the hydrogens around them, aren't short contacts besides, also
    // where the hydrogen lies past the first `inner`.
    let largest_acceptor = ACCEPTORS
        .iter()
        .map(|(element, _)| *VAN_DER_WAALS_RADII.get(element))
        .fold(0.0, f32::max);
    let hydrogen_reach = *VAN_DER_WAALS_RADII.get("H") + largest_acceptor;
    let candidates = crate::geometry::pairs_within(positions, hydrogen_reach)
        .into_iter()
        .filter(|&(i, j, distance)| {
            (elements[i] == "H" || elements[j] == "H") && touching(i, j, distance)
        });

    let mut found = Vec::new();
    let electronegativity = |atom: usize| ELECTRONEGATIVITIES.get(&elements[atom]);
    let accepts = |atom: usize| {
        let Some(&(_, valence_electrons)) = ACCEPTORS
            .iter()
            .find(|(element, _)| *element == elements[atom])
        else {
            return false;
        };
        match state_of(atoms[atom].type_symbol) {
            Some(state) => valence_electrons - state >= 2.0,
            None => !neighbours[atom]
                .iter()
                .any(|&other| electronegativity(other) > electronegativity(atom)),
        }
    };

    let mut joined = neighbours.clone();
    for (i, j, distance) in candidates {
        let hydrogen_bond = [(i, j), (j, i)].into_iter().find(|&(hydrogen, acceptor)| {
            elements[hydrogen] == "H"
                && accepts(acceptor)
                && !within(&neighbours, hydrogen, FEWEST_BONDS_APART).contains(&acceptor)
                && neighbours[hydrogen].iter().any(|&donor| {
                    DONORS.contains(&elements[donor].as_str()) && {
                        let to_donor = position(donor) - position(hydrogen);
                        let to_acceptor = position(acceptor) - position(hydrogen);
                        to_donor.angle(to_acceptor).0.to_degrees() >= SMALLEST_ANGLE
                    }
                })
        });
        if let Some(atoms) = hydrogen_bond {
            joined[i].push(j);
            joined[j].push(i);
            if i < inner && j < inner {
                found.push(Contact {
                    atoms,
                    distance,
                    kind: ContactKind::HydrogenBond,
                });
            }
        }
    }

    let fragment = fragments(&neighbours);
    let mut near: HashMap<usize, HashSet<usize>> = HashMap::default();
    for &(i, j, distance) in &close {
        if fragment[i] == fragment[j] || is_metal(&elements[i]) || is_metal(&elements[j]) {
            continue;
        }
        let bonds_apart = near
            .entry(i)
            .or_insert_with(|| within(&joined, i, FEWEST_BONDS_APART));
        if !bonds_apart.contains(&j) {
            found.push(Contact {
                atoms: (i, j),
                distance,
                kind: ContactKind::Short,
            });
        }
    }
    found
}

/// Each atom's fragment: the same number for atoms joined by a path along `neighbours`.
fn fragments(neighbours: &[Vec<usize>]) -> Vec<usize> {
    let mut fragment = vec![usize::MAX; neighbours.len()];
    for start in 0..neighbours.len() {
        if fragment[start] != usize::MAX {
            continue;
        }
        fragment[start] = start;
        let mut stack = vec![start];
        while let Some(atom) = stack.pop() {
            for &next in &neighbours[atom] {
                if fragment[next] == usize::MAX {
                    fragment[next] = start;
                    stack.push(next);
                }
            }
        }
    }
    fragment
}

/// The atoms at most `most` bonds from `atom` along `neighbours`, `atom` included.
fn within(neighbours: &[Vec<usize>], atom: usize, most: usize) -> HashSet<usize> {
    let mut reached: HashSet<usize> = HashSet::default();
    reached.insert(atom);
    let mut frontier = vec![atom];
    for _ in 0..most {
        frontier = frontier
            .iter()
            .flat_map(|&atom| &neighbours[atom])
            .copied()
            .filter(|&next| reached.insert(next))
            .collect();
    }
    reached
}

#[cfg(test)]
mod tests {
    use super::*;

    fn atoms(symbols: &[&'static str]) -> Vec<BondAtom<'static>> {
        symbols
            .iter()
            .map(|&symbol| BondAtom::from(symbol))
            .collect()
    }

    #[test]
    fn a_hydrogen_between_two_waters_bonds_them_and_their_oxygens_touch_no_further() {
        // O–H···O, linear, 0.97 and 1.85 Å; the second water's hydrogens point away.
        let positions = [
            [0.0, 0.0, 0.0],
            [0.97, 0.0, 0.0],
            [2.82, 0.0, 0.0],
            [3.1, 0.9, 0.0],
            [3.1, -0.9, 0.0],
        ];
        let atoms = atoms(&["O2-", "H1+", "O2-", "H1+", "H1+"]);
        let bonds = [(0, 1, 0.97), (2, 3, 0.94), (2, 4, 0.94)];

        let found = contacts(&positions, &atoms, &bonds, positions.len());
        assert_eq!(found.len(), 1, "{found:?}");
        assert_eq!(found[0].kind, ContactKind::HydrogenBond);
        assert_eq!(found[0].atoms, (1, 2), "the hydrogen first");
        assert!(
            (found[0].distance - 1.85).abs() < 1e-4,
            "{}",
            found[0].distance
        );
    }

    #[test]
    fn a_bent_approach_is_a_short_contact_not_a_hydrogen_bond() {
        // The acceptor at 90° to O–H, 2.2 Å from the hydrogen.
        let positions = [[0.0, 0.0, 0.0], [0.97, 0.0, 0.0], [0.97, 2.2, 0.0]];
        let atoms = atoms(&["O2-", "H1+", "O2-"]);
        let bonds = [(0, 1, 0.97)];

        let found = contacts(&positions, &atoms, &bonds, positions.len());
        assert!(
            found
                .iter()
                .all(|contact| contact.kind == ContactKind::Short),
            "{found:?}"
        );
        assert!(!found.is_empty());
    }

    #[test]
    fn atoms_few_bonds_apart_are_never_in_contact() {
        // Two oxygens of one octahedron edge, 2.83 Å apart through their titanium.
        let positions = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]];
        let atoms = atoms(&["Ti4+", "O2-", "O2-"]);
        let bonds = [(0, 1, 2.0), (0, 2, 2.0)];
        assert!(contacts(&positions, &atoms, &bonds, positions.len()).is_empty());
    }

    #[test]
    fn atoms_that_may_never_be_present_together_are_never_in_contact() {
        // Two carbons of separate molecules, 3.2 Å apart.
        let positions = [[0.0, 0.0, 0.0], [3.2, 0.0, 0.0]];
        let found = |occupancies: [f32; 2]| {
            let atoms: Vec<BondAtom> = occupancies
                .iter()
                .map(|&place_occupancy| BondAtom {
                    place_occupancy,
                    disorder: None,
                    ..BondAtom::from("C")
                })
                .collect();
            contacts(&positions, &atoms, &[], positions.len()).len()
        };

        assert_eq!(found([1.0, 1.0]), 1);
        assert_eq!(found([1.0, 0.25]), 1, "the full one is always there");
        assert_eq!(
            found([0.5, 0.5]),
            0,
            "two orientations of a disordered molecule"
        );
    }

    /// The hydrogen bonds among the atoms of `symbols` at `positions` with `bonds`.
    fn hydrogen_bonds(
        positions: &[[f32; 3]],
        symbols: &[&'static str],
        bonds: &[(usize, usize, f32)],
    ) -> Vec<(usize, usize)> {
        contacts(positions, &atoms(symbols), bonds, positions.len())
            .into_iter()
            .filter(|contact| contact.kind == ContactKind::HydrogenBond)
            .map(|contact| contact.atoms)
            .collect()
    }

    #[test]
    fn a_hydrogen_bonds_to_a_sulfates_oxygen_not_to_its_sulfur() {
        // As in gypsum: O–H···O in a line, H···O 1.86 Å, the sulfur 2.91 Å from the hydrogen,
        // within 3.0 Å, the sum of their van der Waals radii.
        let positions = [
            [0.0, 0.0, 0.0],
            [0.96, 0.0, 0.0],
            [2.82, 0.0, 0.0],
            [3.585, 1.255, 0.0],
        ];
        let bonds = [(0, 1, 0.96), (2, 3, 1.47)];

        // S⁶⁺ has no lone pair.
        let stated = ["O2-", "H1+", "O2-", "S6+"];
        assert_eq!(hydrogen_bonds(&positions, &stated, &bonds), [(1, 2)]);

        // Without charges, it is bonded to the more electronegative oxygen.
        let unstated = ["O", "H", "O", "S"];
        assert_eq!(hydrogen_bonds(&positions, &unstated, &bonds), [(1, 2)]);

        let found = contacts(&positions, &atoms(&unstated), &bonds, positions.len());
        assert!(
            found.iter().all(|contact| contact.atoms != (1, 3)),
            "nor a short contact to the sulfur: {found:?}"
        );
    }

    #[test]
    fn hydroxylamines_nitrogen_bonded_to_oxygen_accepts_where_its_charge_says_it_has_a_lone_pair() {
        // O–H···N in a line, H···N 1.95 Å; the nitrogen bonded to an oxygen behind it.
        let positions = [
            [0.0, 0.0, 0.0],
            [0.96, 0.0, 0.0],
            [2.91, 0.0, 0.0],
            [3.6, 1.2, 0.0],
        ];
        let bonds = [(0, 1, 0.96), (2, 3, 1.45)];

        let stated = ["O2-", "H1+", "N1-", "O2-"];
        assert_eq!(hydrogen_bonds(&positions, &stated, &bonds), [(1, 2)]);

        let unstated = ["O", "H", "N", "O"];
        assert!(
            hydrogen_bonds(&positions, &unstated, &bonds).is_empty(),
            "without a charge, the electronegativity stand-in rules it out"
        );
    }

    #[test]
    fn layers_without_bonds_between_them_touch() {
        // Two sulfurs of neighbouring MoS2 layers, 3.5 Å apart, each bonded only in its layer.
        let positions = [
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.6],
            [0.0, 0.0, 5.1],
            [0.0, 0.0, 6.7],
        ];
        let atoms = atoms(&["Mo4+", "S2-", "S2-", "Mo4+"]);
        let bonds = [(0, 1, 1.6), (2, 3, 1.6)];

        let found = contacts(&positions, &atoms, &bonds, positions.len());
        assert_eq!(found.len(), 1, "{found:?}");
        assert_eq!(found[0].atoms, (1, 2));
        assert_eq!(found[0].kind, ContactKind::Short);
    }
}
