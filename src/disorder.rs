//! Which atoms of a disordered structure are present together: each partly occupied atom's
//! disorder assembly (a part of the structure in alternative configurations, one at a time) and
//! group (its configuration), so that atoms of different configurations aren't bonded.
//!
//! A CIF may state them (`_atom_site_disorder_assembly`, `_atom_site_disorder_group`, SHELX's
//! PART); a negative group's copies by symmetry are configurations of their own. Where it
//! doesn't, they are inferred as the Cambridge Structural Database does (Bruno et al., Acta
//! Cryst. B 67, 333, 2011): partly occupied atoms too close to be present together conflict,
//! those at bond length may be bonded, and both join them into an assembly. The bonds are taken
//! in order of how usual their length is, each joining two atoms' configurations unless that
//! would put conflicting atoms into one. An assembly whose configurations don't all hold the
//! same elements at one occupancy, or in which two elements share a partly occupied site (the
//! two ends of a molecule pointing either way), can't be told apart this way and stays
//! unresolved: its atoms are bonded only where their places must be occupied together.

use crate::{COVALENT_RADII, METALLIC_RADII, Phase};
use cgmath::{InnerSpace, Point3};
use std::collections::{HashMap, HashSet};

use crate::Atom;
use crate::bonding::{BondAtom, BondRule, LimitKind};
use crate::cell_atoms::CellAtom;
use crate::elements::element_of;
use crate::occupancy;
use crate::oxidation::SiteChemistry;

/// The configuration of a disordered part of a structure an atom belongs to.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Disorder {
    /// The disordered part (assembly), in alternative configurations one of which is present at
    /// a time.
    pub assembly: u32,
    /// Its configuration (group) within the assembly; `None` where they couldn't be told apart.
    pub group: Option<u32>,
}

/// Whether atoms `a` and `b` are never present together, so aren't bonded: they are of one
/// assembly and of different configurations; or one of them is of an assembly whose
/// configurations couldn't be told apart, which tells nothing of what it is present with, and
/// their places needn't be occupied together (see [`occupancy::together`]).
pub fn apart(a: &BondAtom, b: &BondAtom) -> bool {
    let unresolved = |disorder: Option<Disorder>| disorder.is_some_and(|d| d.group.is_none());
    match (a.disorder, b.disorder) {
        _ if unresolved(a.disorder) || unresolved(b.disorder) => {
            !occupancy::together(a.place_occupancy, b.place_occupancy)
        }
        (Some(a_disorder), Some(b_disorder)) if a_disorder.assembly == b_disorder.assembly => {
            a_disorder.group != b_disorder.group
        }
        _ => false,
    }
}

/// How far (fractional) a copy of a site may be from where a symmetry operation puts the site to
/// count as made by it.
const SAME_COPY: f64 = 1e-3;

/// How much the occupancies of one inferred configuration's atoms may differ.
const OCCUPANCY_SPREAD: f32 = 0.05;

/// Partly occupied atoms closer than this times the sum of their covalent radii are taken for
/// different configurations: closer than any bond between them usually is. In 4,968 ICSD
/// structures, the shortest bonds between atoms filling their sites (W=O, V=O, U=O) come to 0.68
/// to 0.8 of it; two carbons of a disordered molecule's configurations often lie 0.9 to 1.1 Å apart.
const CONFIGURATION_CONFLICT_FACTOR: f32 = 0.75;

/// What the disorder of an atom of a cell is decided by.
struct Site<'a> {
    /// The phase's atom it is a copy of.
    atom: &'a Atom,
    /// Its fractional coordinates.
    fractional: Point3<f64>,
    /// How often its place is occupied (see [`CellAtom::place_occupancy`]).
    place_occupancy: f32,
    /// Whether it has valence left to bond to other metals (see
    /// [`SiteChemistry::bonds_to_metals`]).
    bonds_to_metals: bool,
}

impl Site<'_> {
    /// The occupancy of its site.
    fn occupancy(&self) -> f32 {
        self.atom.occupancy.value as f32
    }

    /// The atom as the bond rule sees it, with no disorder.
    fn bond_atom(&self) -> BondAtom<'_> {
        BondAtom {
            type_symbol: &self.atom.type_,
            bonds_to_metals: self.bonds_to_metals,
            place_occupancy: self.place_occupancy,
            disorder: None,
        }
    }
}

/// Each of `atoms` of `phase`'s cell (see [`crate::cell_atoms::cell_atoms`]) its disorder, with
/// `chemistry` the phase's sites' (see [`crate::oxidation::site_chemistry`]): stated by the CIF
/// where it states groups that hold no atoms too close to be present together, else inferred for
/// the atoms on partly occupied places; `None` for the others.
pub fn of_atoms(
    phase: &Phase,
    atoms: &[CellAtom],
    chemistry: &[SiteChemistry],
) -> Vec<Option<Disorder>> {
    let sites: Vec<Site> = atoms
        .iter()
        .map(|atom| Site {
            atom: &phase.atoms[atom.site],
            fractional: atom.fractional,
            place_occupancy: atom.place_occupancy,
            bonds_to_metals: chemistry
                .get(atom.site)
                .is_none_or(|chemistry| chemistry.bonds_to_metals),
        })
        .collect();
    let sites = sites.as_slice();
    let (mut disorder, assemblies, groups) = stated(phase, sites);

    let involved: Vec<bool> = sites
        .iter()
        .zip(&disorder)
        .map(|(site, stated)| stated.is_some() || occupancy::is_partial(site.place_occupancy))
        .collect();
    let Some(pairs) = Pairs::of(phase, sites, &involved) else {
        return disorder;
    };

    // A stated group holding atoms too close to be present together isn't a configuration (the
    // copies by symmetry of a positive group around a symmetry element): its assembly is
    // inferred instead.
    let refuted: HashSet<u32> = pairs
        .conflicts
        .iter()
        .filter_map(|&(a, b)| match (disorder[a], disorder[b]) {
            (Some(a), Some(b)) if a == b => Some(a.assembly),
            _ => None,
        })
        .collect();
    for stated in &mut disorder {
        if stated.is_some_and(|stated| refuted.contains(&stated.assembly)) {
            *stated = None;
        }
    }

    let partial: Vec<bool> = sites
        .iter()
        .zip(&disorder)
        .map(|(site, stated)| stated.is_none() && occupancy::is_partial(site.place_occupancy))
        .collect();
    for (index, inferred) in infer(sites, &partial, &pairs) {
        disorder[index] = Some(Disorder {
            assembly: assemblies + inferred.assembly,
            group: inferred.group.map(|group| groups + group),
        });
    }
    disorder
}

/// The disorder the CIF of `phase` states for each of `sites`, and how many assemblies and
/// groups it numbers. Group `0` is none (SHELX's PART 0).
fn stated(phase: &Phase, sites: &[Site]) -> (Vec<Option<Disorder>>, u32, u32) {
    let mut disorder = vec![None; sites.len()];
    let mut assemblies: HashMap<String, u32> = HashMap::default();
    let mut groups: HashMap<(String, String, Option<usize>), u32> = HashMap::default();

    let operations = crate::SpaceGroupSymmetryOperations::of_cell(&phase.cell);
    for (index, site) in sites.iter().enumerate() {
        let atom = site.atom;
        let Some(group) = atom.disorder_group.as_deref().filter(|group| *group != "0") else {
            continue;
        };
        let assembly = atom.disorder_assembly.clone().unwrap_or_default();

        // A negative group's copy is a configuration of its own: told by the operation that made
        // it from the stated position.
        let copy = group.starts_with('-').then(|| {
            let stated = Point3::new(atom.x.value, atom.y.value, atom.z.value);
            operations.as_ref().and_then(|operations| {
                operations.iter().position(|operation| {
                    let made = Point3::from_homogeneous(operation * stated.to_homogeneous());
                    let offset = (made - site.fractional).map(|c| c - c.round());
                    offset.magnitude() < SAME_COPY
                })
            })
        });
        let copy = copy.flatten();

        let next = assemblies.len() as u32;
        let assembly_id = *assemblies.entry(assembly.clone()).or_insert(next);
        let next = groups.len() as u32;
        let group_id = *groups
            .entry((assembly, group.to_string(), copy))
            .or_insert(next);
        disorder[index] = Some(Disorder {
            assembly: assembly_id,
            group: Some(group_id),
        });
    }
    (disorder, assemblies.len() as u32, groups.len() as u32)
}

/// The pairs of a cell's sites that decide its disorder, each once.
struct Pairs {
    /// Atoms too close to be present together: on one place, or closer than
    /// [`CONFIGURATION_CONFLICT_FACTOR`] times the sum of their covalent radii.
    conflicts: HashSet<(usize, usize)>,
    /// Atoms at bond length, each pair with how unusual its length is (see [`unusualness`]).
    candidates: HashMap<(usize, usize), f32>,
    /// Atoms on a place they share with an atom of another element.
    shared_places: HashSet<usize>,
}

impl Pairs {
    /// The pairs among the `involved` of the cell's `sites` of `phase`; `None` for a cell without
    /// volume.
    fn of(phase: &Phase, sites: &[Site], involved: &[bool]) -> Option<Self> {
        let elements: Vec<String> = sites
            .iter()
            .map(|site| element_of(&site.atom.type_))
            .collect();
        let rule = BondRule::default();
        let involved_elements: HashSet<&str> = elements
            .iter()
            .zip(involved)
            .filter(|(_, involved)| **involved)
            .map(|(element, _)| element.as_str())
            .collect();
        let cutoff = involved_elements
            .iter()
            .flat_map(|a| involved_elements.iter().map(move |b| (a, b)))
            .map(|(a, b)| rule.limit(a, b).max)
            .fold(0.0, f32::max);

        let mut pairs = Self {
            conflicts: HashSet::default(),
            candidates: HashMap::default(),
            shared_places: HashSet::default(),
        };
        if cutoff <= 0.0 {
            return Some(pairs);
        }
        for (from, to, distance) in crate::geometry::periodic_contacts(
            &phase.cell,
            &sites.iter().map(|site| site.fractional).collect::<Vec<_>>(),
            cutoff,
        )? {
            if from >= to || !involved[from] || !involved[to] {
                continue;
            }
            let elements_of = (elements[from].as_str(), elements[to].as_str());
            if distance < occupancy::SAME_PLACE {
                if elements_of.0 != elements_of.1 {
                    pairs.shared_places.extend([from, to]);
                }
                pairs.conflicts.insert((from, to));
            } else if distance < CONFIGURATION_CONFLICT_FACTOR * covalent(elements_of) {
                pairs.conflicts.insert((from, to));
            } else {
                let atom = |index: usize| BondAtom {
                    place_occupancy: 1.0,
                    disorder: None,
                    ..sites[index].bond_atom()
                };
                if rule.are_bonded(atom(from), atom(to), distance) {
                    let score = unusualness(&rule, elements_of, distance);
                    pairs
                        .candidates
                        .entry((from, to))
                        .and_modify(|best| *best = best.min(score))
                        .or_insert(score);
                }
            }
        }
        Some(pairs)
    }
}

/// The disorder inferred for the `partial` ones of `sites`, from `pairs`, by index; numbered from
/// 0.
fn infer(sites: &[Site], partial: &[bool], pairs: &Pairs) -> Vec<(usize, Disorder)> {
    let elements: Vec<String> = sites
        .iter()
        .map(|site| element_of(&site.atom.type_))
        .collect();

    let among_partial = |&(a, b): &(usize, usize)| partial[a] && partial[b];
    let conflicts: HashSet<(usize, usize)> = pairs
        .conflicts
        .iter()
        .copied()
        .filter(among_partial)
        .collect();

    let candidates: HashMap<(usize, usize), f32> = pairs
        .candidates
        .iter()
        .filter(|(pair, _)| among_partial(pair))
        .map(|(&pair, &score)| (pair, score))
        .collect();

    let shared_places = &pairs.shared_places;

    // Assemblies: what conflicts and candidate bonds join.
    let mut assembly_of = UnionFind::new(sites.len());
    for &(a, b) in conflicts.iter().chain(candidates.keys()) {
        assembly_of.union(a, b);
    }

    // Configurations: candidate bonds between atoms of one occupancy (a configuration's atoms
    // share it), in order of how usual their length is, each joining two configurations unless
    // they hold conflicting atoms; hydrogen, whose small radius leaves its conflicts unseen, last,
    // each joining only the configuration of the atom it is bonded to most usually.
    let mut by_score: Vec<((usize, usize), f32)> = candidates
        .into_iter()
        .filter(|&((a, b), _)| {
            (sites[a].occupancy() - sites[b].occupancy()).abs() <= OCCUPANCY_SPREAD
        })
        .collect();

    by_score.sort_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)));
    let is_hydrogen = |index: usize| elements[index] == "H";
    let (hydrogen, heavy): (Vec<_>, Vec<_>) = by_score
        .into_iter()
        .partition(|&((a, b), _)| is_hydrogen(a) || is_hydrogen(b));

    let mut members: Vec<Vec<usize>> = (0..sites.len()).map(|index| vec![index]).collect();
    let mut group_of = UnionFind::new(sites.len());
    let conflicting = |a: usize, b: usize| conflicts.contains(&(a.min(b), a.max(b)));
    let mut join = |a: usize, b: usize| {
        let (root_a, root_b) = (group_of.find(a), group_of.find(b));
        if root_a == root_b
            || members[root_a]
                .iter()
                .any(|&x| members[root_b].iter().any(|&y| conflicting(x, y)))
        {
            return false;
        }
        let root = group_of.union(root_a, root_b);
        let other = if root == root_a { root_b } else { root_a };
        let moved = std::mem::take(&mut members[other]);
        members[root].extend(moved);
        true
    };
    for ((a, b), _) in heavy {
        join(a, b);
    }
    let mut placed: HashSet<usize> = HashSet::default();
    for ((a, b), _) in hydrogen {
        let (hydrogen, other) = if is_hydrogen(a) { (a, b) } else { (b, a) };
        if is_hydrogen(other) || placed.contains(&hydrogen) {
            continue;
        }
        if join(hydrogen, other) {
            placed.insert(hydrogen);
        }
    }

    // Each assembly of more than one atom, with its configurations; resolved where they hold the
    // same elements at one occupancy each and no partly occupied site is shared.
    let mut assembly_atoms: HashMap<usize, Vec<usize>> = HashMap::default();
    for index in (0..sites.len()).filter(|&index| partial[index]) {
        assembly_atoms
            .entry(assembly_of.find(index))
            .or_default()
            .push(index);
    }
    let mut assembly_roots: Vec<usize> = assembly_atoms
        .iter()
        .filter(|(_, atoms)| atoms.len() > 1)
        .map(|(&root, _)| root)
        .collect();
    assembly_roots.sort_unstable();

    let mut inferred = Vec::new();
    let mut next_group = 0;
    for (assembly, root) in assembly_roots.into_iter().enumerate() {
        let atoms = &assembly_atoms[&root];
        let mut configurations: HashMap<usize, Vec<usize>> = HashMap::default();
        for &index in atoms {
            configurations
                .entry(group_of.find(index))
                .or_default()
                .push(index);
        }

        let composition = |configuration: &[usize]| {
            let mut elements: Vec<&str> = configuration
                .iter()
                .map(|&index| elements[index].as_str())
                .collect();
            elements.sort_unstable();
            elements
        };

        let uniform = |configuration: &[usize]| {
            let occupancies = configuration.iter().map(|&index| sites[index].occupancy());
            let (low, high) = occupancies.fold((f32::MAX, f32::MIN), |(low, high), occupancy| {
                (low.min(occupancy), high.max(occupancy))
            });
            high - low <= OCCUPANCY_SPREAD
        };

        let mut configurations: Vec<Vec<usize>> = configurations.into_values().collect();
        configurations.sort_unstable();
        let first = composition(&configurations[0]);
        let resolved = !atoms.iter().any(|index| shared_places.contains(index))
            && configurations
                .iter()
                .all(|configuration| composition(configuration) == first && uniform(configuration));
        for configuration in &configurations {
            let group = resolved.then_some(next_group);
            next_group += 1;
            for &index in configuration {
                inferred.push((
                    index,
                    Disorder {
                        assembly: assembly as u32,
                        group,
                    },
                ));
            }
        }
    }
    inferred
}

/// How unusual a bond of `distance` (Å) between atoms of the elements `elements` is, for
/// ordering candidates: its relative departure from the sum of their covalent or metallic radii,
/// or for a cation–anion pair (whose length varies more) its share of the longest bond the rule
/// allows.
fn unusualness(rule: &BondRule, elements: (&str, &str), distance: f32) -> f32 {
    let limit = rule.limit(elements.0, elements.1);
    let usual = match limit.kind {
        LimitKind::Covalent => covalent(elements),
        LimitKind::Metallic => METALLIC_RADII
            .get(elements.0)
            .zip(METALLIC_RADII.get(elements.1))
            .map_or(0.0, |(a, b)| a + b),
        LimitKind::Ionic | LimitKind::Table => return distance / limit.max,
    };
    if usual > 0.0 {
        (distance / usual - 1.0).abs()
    } else {
        distance / limit.max
    }
}

/// The sum of the covalent radii (Å) of the elements `elements`.
fn covalent(elements: (&str, &str)) -> f32 {
    COVALENT_RADII.get(elements.0) + COVALENT_RADII.get(elements.1)
}

/// Disjoint sets of indices.
struct UnionFind(Vec<usize>);

impl UnionFind {
    fn new(count: usize) -> Self {
        Self((0..count).collect())
    }

    fn find(&mut self, index: usize) -> usize {
        let mut root = index;
        while self.0[root] != root {
            root = self.0[root];
        }
        let mut at = index;
        while self.0[at] != root {
            let next = self.0[at];
            self.0[at] = root;
            at = next;
        }
        root
    }

    /// Joins the sets of `a` and `b`; their root.
    fn union(&mut self, a: usize, b: usize) -> usize {
        let (a, b) = (self.find(a), self.find(b));
        self.0[b] = a;
        a
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A phase's cell with its atoms' disorder.
    struct Cell {
        phase: Phase,
        atoms: Vec<CellAtom>,
        disorder: Vec<Option<Disorder>>,
    }

    impl Cell {
        /// The label of the cell's atom `index`.
        fn label(&self, index: usize) -> &str {
            &self.phase.atoms[self.atoms[index].site].label
        }

        /// The cell's first atom labelled `label`.
        fn index(&self, label: &str) -> usize {
            (0..self.atoms.len())
                .find(|&index| self.label(index) == label)
                .unwrap()
        }

        /// The cell's atom `index` as the bond rule sees it.
        fn bond_atom(&self, index: usize) -> BondAtom<'_> {
            BondAtom {
                type_symbol: &self.phase.atoms[self.atoms[index].site].type_,
                bonds_to_metals: true,
                place_occupancy: self.atoms[index].place_occupancy,
                disorder: self.disorder[index],
            }
        }

        /// Whether the cell's atoms labelled `a` and `b` are bonded, by the default rule.
        fn bonded(&self, a: &str, b: &str) -> bool {
            let (a, b) = (self.index(a), self.index(b));
            let distance = (self.atoms[a].fractional - self.atoms[b].fractional).magnitude() * 10.0;
            BondRule::default().are_bonded(self.bond_atom(a), self.bond_atom(b), distance as f32)
        }
    }

    /// A 10 Å cubic cell of `space_group` with `sites` (`_atom_site_label`, type symbol,
    /// fractional x, y, z, occupancy, disorder group).
    fn cell(space_group: u8, operations: &[&str], sites: &str) -> Cell {
        let operations: String = operations.iter().map(|op| format!("'{op}'\n")).collect();
        let cif = format!(
            "data_test
_cell_length_a 10
_cell_length_b 10
_cell_length_c 10
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_cell_volume 1000
_space_group_IT_number {space_group}
loop_
_space_group_symop_operation_xyz
{operations}loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
_atom_site_disorder_group
{sites}"
        );
        let (_, phase) = crate::cif_phase::try_phase_from_cif_bytes(cif.as_bytes()).unwrap();
        let atoms = crate::cell_atoms::cell_atoms(&phase);
        let chemistry = crate::oxidation::site_chemistry(&phase);
        let disorder = of_atoms(&phase, &atoms, &chemistry);
        Cell {
            phase,
            atoms,
            disorder,
        }
    }

    #[test]
    fn two_configurations_are_inferred_from_atoms_too_close_to_be_present_together() {
        // C–O in two orientations, each half the time: the carbons 0.3 Å apart, each with its
        // oxygen 1.42 Å away. The second carbon is 1.45 Å from the first oxygen.
        let cell = cell(
            1,
            &["x, y, z"],
            "C1A C 0.5 0.5 0.5 0.5 .
O1A O 0.642 0.5 0.5 0.5 .
C1B C 0.5 0.53 0.5 0.5 .
O1B O 0.5 0.672 0.5 0.5 .
",
        );

        assert!(
            cell.disorder
                .iter()
                .all(|d| d.is_some_and(|d| d.group.is_some()))
        );
        assert!(cell.bonded("C1A", "O1A"));
        assert!(cell.bonded("C1B", "O1B"));
        assert!(!cell.bonded("C1B", "O1A"), "of the other configuration");
    }

    #[test]
    fn a_negative_groups_copy_by_symmetry_is_a_configuration_of_its_own() {
        // C1–C2 across an inversion centre: C1 and its inverted copy are 1.4 Å apart, as far as
        // C1 and C2.
        let molecule = |group: &str| {
            cell(
                2,
                &["x, y, z", "-x, -y, -z"],
                &format!(
                    "C1 C 0.07 0 0 0.5 {group}
C2 C 0.07 0.14 0 0.5 {group}
"
                ),
            )
        };
        let copy_of_c1 = |cell: &Cell| {
            (0..cell.atoms.len())
                .find(|&index| cell.label(index) == "C1" && cell.atoms[index].fractional.x > 0.5)
                .unwrap()
        };

        let negative = molecule("-1");
        let (c1, copy) = (negative.index("C1"), copy_of_c1(&negative));
        assert_ne!(negative.disorder[c1], negative.disorder[copy]);
        assert!(apart(&negative.bond_atom(c1), &negative.bond_atom(copy)));
        assert!(negative.bonded("C1", "C2"));

        let positive = molecule("1");
        let (c1, copy) = (positive.index("C1"), copy_of_c1(&positive));
        assert!(!apart(&positive.bond_atom(c1), &positive.bond_atom(copy)));
    }

    #[test]
    fn a_stated_group_of_atoms_too_close_to_be_present_together_is_not_taken() {
        // A C and an N stated to be of one group on one site, as the ends of a methylammonium
        // pointing either way are.
        let cell = cell(
            1,
            &["x, y, z"],
            "C1 C 0.5 0.5 0.5 0.25 1
N1 N 0.5 0.5 0.5 0.25 1
",
        );
        assert!(
            cell.disorder
                .iter()
                .all(|d| d.is_some_and(|d| d.group.is_none()))
        );
    }
}
