//! The fragments the bonds of a periodic structure join its sites into: molecules, which are
//! finite, and chains, layers and frameworks, which run on through the crystal.

use std::collections::VecDeque;

/// A bond of a periodic structure: from site `from` to the copy of site `to` shifted by `shift`
/// lattice translations (cells along a, b and c).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct PeriodicBond {
    /// The site it starts at.
    pub from: usize,
    /// The site whose copy it ends at.
    pub to: usize,
    /// The lattice translation from site `to` to the copy it ends at.
    pub shift: [i32; 3],
}

/// Sites joined by bonds.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Fragment {
    /// Its sites, each with the lattice translation from the site to its copy in the fragment,
    /// the first site's untranslated. A periodic fragment has a copy of its sites in every cell
    /// it runs through; these are the copies its bonds reach from the first site without
    /// returning to a site.
    pub sites: Vec<(usize, [i32; 3])>,
    /// In how many independent directions it runs on through the crystal: 0 for a molecule, 1
    /// for a chain, 2 for a layer, 3 for a framework. These are the lattice translations from a
    /// site to the other copies of it its bonds lead to.
    pub dimensions: u8,
}

impl Fragment {
    /// Whether it runs on through the crystal (a chain, layer or framework), not a molecule.
    pub fn is_periodic(&self) -> bool {
        self.dimensions > 0
    }
}

/// The independent directions among lattice translations: at most three, each added only if
/// it isn't a combination of those before.
#[derive(Debug, Default)]
struct Directions(Vec<[i64; 3]>);

impl Directions {
    fn add(&mut self, translation: [i32; 3]) {
        let v = translation.map(i64::from);
        let cross = |a: [i64; 3], b: [i64; 3]| {
            [
                a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0],
            ]
        };
        let independent = match self.0.as_slice() {
            [] => v != [0; 3],
            [a] => cross(*a, v) != [0; 3],
            [a, b] => {
                let normal = cross(*a, *b);
                normal[0] * v[0] + normal[1] * v[1] + normal[2] * v[2] != 0
            }
            _ => false,
        };
        if independent {
            self.0.push(v);
        }
    }
}

/// The fragments of a periodic structure's sites, each site in exactly one.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Fragments {
    fragments: Vec<Fragment>,
    of_site: Vec<(usize, [i32; 3])>,
}

impl Fragments {
    /// The fragments `bonds` join `site_count` sites into (a site without bonds is a fragment of
    /// its own). A bond may be listed once or from both of its ends.
    pub fn new(site_count: usize, bonds: &[PeriodicBond]) -> Self {
        let mut neighbours = vec![Vec::new(); site_count];
        for bond in bonds {
            debug_assert!(bond.from < site_count && bond.to < site_count, "{bond:?}");
            neighbours[bond.from].push((bond.to, bond.shift));
            neighbours[bond.to].push((bond.from, bond.shift.map(|c| -c)));
        }

        let mut placed: Vec<Option<(usize, [i32; 3])>> = vec![None; site_count];
        let mut fragments = Vec::new();
        for start in 0..site_count {
            if placed[start].is_some() {
                continue;
            }
            let index = fragments.len();
            let mut sites = Vec::new();
            let mut directions = Directions::default();

            placed[start] = Some((index, [0; 3]));
            let mut queue = VecDeque::from([start]);
            while let Some(site) = queue.pop_front() {
                let at = placed[site].map_or([0; 3], |(_, at)| at);
                sites.push((site, at));
                for &(to, shift) in &neighbours[site] {
                    let reached = [0, 1, 2].map(|axis| at[axis] + shift[axis]);
                    match placed[to] {
                        None => {
                            placed[to] = Some((index, reached));
                            queue.push_back(to);
                        }
                        Some((_, there)) => {
                            directions.add([0, 1, 2].map(|axis| reached[axis] - there[axis]));
                        }
                    }
                }
            }
            fragments.push(Fragment {
                sites,
                dimensions: directions.0.len() as u8,
            });
        }

        let of_site = placed
            .into_iter()
            .map(|place| place.unwrap_or((0, [0; 3])))
            .collect();
        Self { fragments, of_site }
    }

    /// Every fragment, in the order of their first sites.
    pub fn all(&self) -> &[Fragment] {
        &self.fragments
    }

    /// The fragment at `index` (see [`Self::of`]).
    pub fn get(&self, index: usize) -> &Fragment {
        &self.fragments[index]
    }

    /// The index of the fragment `site` belongs to, and the lattice translation from the site to
    /// its copy in it (see [`Fragment::sites`]).
    pub fn of(&self, site: usize) -> (usize, [i32; 3]) {
        self.of_site[site]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bond(from: usize, to: usize, shift: [i32; 3]) -> PeriodicBond {
        PeriodicBond { from, to, shift }
    }

    #[test]
    fn a_molecule_across_a_cell_face_is_finite_with_its_sites_placed_together() {
        // Site 1 is bonded to site 0 only through the copy of site 1 a cell along −a, and site 2
        // to that copy: the molecule is 0, 1 and 2 shifted by −a.
        let fragments = Fragments::new(3, &[bond(0, 1, [-1, 0, 0]), bond(1, 2, [0, 0, 0])]);

        assert_eq!(fragments.all().len(), 1);
        let molecule = fragments.get(0);
        assert_eq!(molecule.dimensions, 0);
        assert_eq!(
            molecule.sites,
            [(0, [0, 0, 0]), (1, [-1, 0, 0]), (2, [-1, 0, 0])]
        );
        assert_eq!(fragments.of(2), (0, [-1, 0, 0]));
    }

    #[test]
    fn a_bond_listed_from_both_ends_is_one_bond() {
        let fragments = Fragments::new(2, &[bond(0, 1, [0, 1, 0]), bond(1, 0, [0, -1, 0])]);

        assert_eq!(fragments.get(0).dimensions, 0);
        assert_eq!(fragments.of(1), (0, [0, 1, 0]));
    }

    #[test]
    fn a_site_bonded_to_its_own_copy_is_a_chain() {
        let fragments = Fragments::new(1, &[bond(0, 0, [0, 0, 1])]);

        assert_eq!(fragments.get(0).dimensions, 1);
    }

    #[test]
    fn bonds_along_two_and_three_directions_make_a_layer_and_a_framework() {
        let layer = Fragments::new(1, &[bond(0, 0, [1, 0, 0]), bond(0, 0, [1, 1, 0])]);
        let framework = Fragments::new(
            1,
            &[
                bond(0, 0, [1, 0, 0]),
                bond(0, 0, [0, 1, 0]),
                bond(0, 0, [0, 0, 1]),
            ],
        );

        assert_eq!(layer.get(0).dimensions, 2);
        assert_eq!(framework.get(0).dimensions, 3);
    }

    #[test]
    fn a_chain_with_bonds_in_both_directions_along_it_is_one_dimensional() {
        // Each site bonded to the next cell's copy and the previous cell's along a, and a zigzag
        // whose two bonds cross twice the same way.
        let straight = Fragments::new(1, &[bond(0, 0, [1, 0, 0]), bond(0, 0, [-1, 0, 0])]);
        let zigzag = Fragments::new(2, &[bond(0, 1, [0, 0, 0]), bond(1, 0, [2, 0, 0])]);

        assert_eq!(straight.get(0).dimensions, 1);
        assert_eq!(zigzag.get(0).dimensions, 1);
    }

    #[test]
    fn a_ring_of_bonds_that_ends_a_cell_further_on_runs_through_the_crystal() {
        // 0–1 in the cell, 1–0 into the next cell along b: …0–1–0–1… along b.
        let fragments = Fragments::new(2, &[bond(0, 1, [0, 0, 0]), bond(1, 0, [0, 1, 0])]);

        assert_eq!(fragments.all().len(), 1);
        assert_eq!(fragments.get(0).dimensions, 1);
    }

    #[test]
    fn a_ring_of_bonds_that_closes_in_one_copy_is_a_molecule() {
        // A three-membered ring whose bonds cross a face and come back.
        let fragments = Fragments::new(
            3,
            &[
                bond(0, 1, [1, 0, 0]),
                bond(1, 2, [0, 0, 0]),
                bond(2, 0, [-1, 0, 0]),
            ],
        );

        assert_eq!(fragments.get(0).dimensions, 0);
        assert_eq!(
            fragments.get(0).sites,
            [(0, [0, 0, 0]), (1, [1, 0, 0]), (2, [1, 0, 0])]
        );
    }

    #[test]
    fn sites_without_bonds_and_separate_molecules_are_fragments_of_their_own() {
        let fragments = Fragments::new(4, &[bond(1, 2, [0, 0, 0])]);

        assert_eq!(fragments.all().len(), 3);
        assert_eq!(fragments.of(0).0, 0);
        assert_eq!(fragments.of(1).0, fragments.of(2).0);
        assert_eq!(fragments.of(3).0, 2);
        assert!(
            fragments
                .all()
                .iter()
                .all(|fragment| !fragment.is_periodic())
        );
    }
}
