//! Whether two atoms are bonded ([`BondRule`]), and what decided it ([`BondRule::decide`]).
//!
//! A pair is decided in this order:
//! 1. Two ions whose [charges](Charge) have the same sign aren't bonded, unless they are of one
//!    element and [`BondRule::like_charges_of_one_element`] allows it. The charges are the ones
//!    the type symbols name: stated by the CIF, or guessed by
//!    [`crate::oxidation::with_guessed_oxidation_states`].
//! 2. Two [metals](crate::elements::is_metal) aren't bonded unless both have valence left for it
//!    ([`BondAtom::bonds_to_metals`], from [`crate::oxidation::site_chemistry`]).
//! 3. Otherwise their distance decides, against the pair's [limit](BondRule::limit): the bond
//!    table's if it lists the pair, else [`BondRule::default_limit`]'s, for a cation–anion pair
//!    by bond valence, for two metals by metallic radii, for any other pair by covalent radii.
//!    [`BondRule::tolerance`] widens the upper limit.
//!
//! What it knows about elements comes from [`crate::elements`]; what an atom's environment
//! adds (its charge, its valence left, its site's occupancy) comes with the atom, as a
//! [`BondAtom`].

use std::sync::Arc;

use crate::{BOND_VALENCE, BondValence, COVALENT_RADII, ELECTRONEGATIVITIES, METALLIC_RADII};
use std::collections::HashMap;

use cgmath::Point3;

use crate::Cell;
use crate::elements::{anion_state, bonds_covalently, element_of, is_anion_former, is_metal};
use crate::geometry::{copies_in_reach, pairs_within};
use crate::molecules::PeriodicBond;
use crate::type_symbol::Charge;

/// Bond limits per element pair: the shortest and longest distance (Å) at which two atoms of
/// those elements count as bonded.
#[derive(Default, Clone, PartialEq)]
pub struct BondTable {
    /// What the table is called, usually its file's name; empty where it has none.
    name: String,
    /// The CSV it was read from.
    csv: String,
    limits: HashMap<(String, String), (f32, f32)>,
}

impl std::fmt::Debug for BondTable {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("BondTable")
            .field("name", &self.name)
            .field("pairs", &self.limits.len())
            .finish()
    }
}

impl BondTable {
    /// The table called `name` (usually its file's name), read from the CSV `text` (see
    /// [`Self::from_csv`]).
    pub fn named(name: impl Into<String>, text: &str) -> Result<Self, String> {
        Ok(Self {
            name: name.into(),
            ..Self::from_csv(text)?
        })
    }

    /// What the table is called, usually its file's name; empty where it has none.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// The CSV the table was read from, to keep it and read it again.
    pub fn csv(&self) -> &str {
        &self.csv
    }

    /// Each element pair the table lists, with its limits (Å): `(element1, element2, min, max)`,
    /// the elements in alphabetical order.
    pub fn pairs(&self) -> impl Iterator<Item = (&str, &str, f32, f32)> {
        self.limits
            .iter()
            .map(|((a, b), (min, max))| (a.as_str(), b.as_str(), *min, *max))
    }

    /// Reads a bond table from CSV: a header naming the columns `element1`, `element2`, `min` and
    /// `max` (Å) in any order, then one row per element pair. Other columns are ignored.
    pub fn from_csv(text: &str) -> Result<Self, String> {
        let mut lines = text
            .lines()
            .enumerate()
            .map(|(index, line)| (index + 1, line.trim()))
            .filter(|(_, line)| !line.is_empty());

        let (_, header) = lines.next().ok_or("the bond table is empty")?;
        let header: Vec<&str> = header.split(',').map(str::trim).collect();
        let column = |name: &str| {
            header
                .iter()
                .position(|column| *column == name)
                .ok_or_else(|| format!("the bond table has no column {name}"))
        };
        let [element1, element2, min, max] = [
            column("element1")?,
            column("element2")?,
            column("min")?,
            column("max")?,
        ];

        let mut limits = HashMap::default();

        for (line_number, line) in lines {
            let fields: Vec<&str> = line.split(',').map(str::trim).collect();
            if fields.len() != header.len() {
                return Err(format!(
                    "line {line_number}: {} fields where the header names {}",
                    fields.len(),
                    header.len()
                ));
            }

            let number = |index: usize| {
                fields[index]
                    .parse::<f32>()
                    .ok()
                    .filter(|value| value.is_finite() && *value >= 0.0)
                    .ok_or_else(|| {
                        format!("line {line_number}: {} is not a distance", fields[index])
                    })
            };
            let (min, max) = (number(min)?, number(max)?);
            if min > max {
                return Err(format!("line {line_number}: min {min} is above max {max}"));
            }

            let pair = pair_key(fields[element1], fields[element2]);
            if limits.insert(pair, (min, max)).is_some() {
                return Err(format!(
                    "line {line_number}: {}-{} is listed twice",
                    fields[element1], fields[element2]
                ));
            }
        }

        if limits.is_empty() {
            return Err("the bond table has no element pairs".to_string());
        }

        Ok(Self {
            name: String::new(),
            csv: text.to_string(),
            limits,
        })
    }

    /// The limits of the pair `a`-`b`, in either order.
    pub fn limits(&self, a: &str, b: &str) -> Option<(f32, f32)> {
        self.limits.get(&pair_key(a, b)).copied()
    }

    /// The number of element pairs.
    pub fn len(&self) -> usize {
        self.limits.len()
    }

    /// Whether the table has no element pairs.
    pub fn is_empty(&self) -> bool {
        self.limits.is_empty()
    }
}

fn pair_key(a: &str, b: &str) -> (String, String) {
    if a <= b {
        (a.to_string(), b.to_string())
    } else {
        (b.to_string(), a.to_string())
    }
}

/// An atom as [`BondRule`] sees it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BondAtom<'a> {
    /// Its type symbol, with the charge it names if any (`Fe3+`).
    pub type_symbol: &'a str,
    /// Whether it has valence left to bond to other metals (see
    /// [`crate::oxidation::SiteChemistry::bonds_to_metals`]); only a metal's matters.
    pub bonds_to_metals: bool,
    /// How often its place is occupied, by it or by the atoms sharing its place (see
    /// [`crate::cell_atoms::CellAtom::place_occupancy`]).
    pub place_occupancy: f32,
    /// The configuration of a disordered part it belongs to, if any (see [`crate::disorder`]).
    pub disorder: Option<crate::disorder::Disorder>,
}

impl<'a> From<&'a str> for BondAtom<'a> {
    /// An atom of `type_symbol` on a fully occupied site that may bond to anything.
    fn from(type_symbol: &'a str) -> Self {
        Self {
            type_symbol,
            bonds_to_metals: true,
            place_occupancy: 1.0,
            disorder: None,
        }
    }
}

/// `atoms` of `phase`'s cell (see [`crate::cell_atoms::cell_atoms`]) as the bond rule sees
/// them, with `chemistry` the phase's sites' (see [`crate::oxidation::site_chemistry`]) and
/// `disorder` each atom's (see [`crate::disorder::of_atoms`]).
pub fn cell_bond_atoms<'a>(
    phase: &'a crate::Phase,
    atoms: &[crate::cell_atoms::CellAtom],
    chemistry: &[crate::oxidation::SiteChemistry],
    disorder: &[Option<crate::disorder::Disorder>],
) -> Vec<BondAtom<'a>> {
    atoms
        .iter()
        .zip(disorder)
        .map(|(atom, &disorder)| BondAtom {
            type_symbol: &phase.atoms[atom.site].type_,
            bonds_to_metals: chemistry
                .get(atom.site)
                .is_none_or(|chemistry| chemistry.bonds_to_metals),
            place_occupancy: atom.place_occupancy,
            disorder,
        })
        .collect()
}

/// The rule that gives a pair of elements its bond limits (see [`BondRule::limit`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LimitKind {
    /// The bond table's limits for the pair.
    Table,
    /// A cation–anion pair, bonded down to a small bond valence.
    Ionic,
    /// Two metals, bonded by their metallic radii.
    Metallic,
    /// Any other pair, bonded by their covalent radii.
    Covalent,
}

impl std::fmt::Display for LimitKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::Table => "table",
            Self::Ionic => "ionic",
            Self::Metallic => "metallic",
            Self::Covalent => "covalent",
        })
    }
}

/// The distances `[min, max)` (Å) at which a pair of elements is bonded, and the rule that gives
/// them.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BondLimit {
    /// The rule the limits come from.
    pub kind: LimitKind,
    /// The shortest distance (Å) bonded.
    pub min: f32,
    /// The distance (Å) from which on the pair isn't bonded, with the rule's tolerance.
    pub max: f32,
}

impl BondLimit {
    /// Whether `distance` (Å) is within the limits.
    pub fn contains(&self, distance: f32) -> bool {
        (self.min..self.max).contains(&distance)
    }
}

/// Whether two atoms are bonded, and what decided it (see [`BondRule::decide`]).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BondDecision {
    /// Not bonded: their charges have the same sign.
    LikeCharges,
    /// Not bonded: two metals of which one has no valence left to bond to metals.
    NoMetalValence,
    /// Not bonded: never present together, being of different configurations of a disordered
    /// part, or of one whose configurations couldn't be told apart (see
    /// [`crate::disorder::apart`]).
    NotTogether,
    /// Bonded if `distance` (Å) is within `limit`.
    ByDistance {
        /// The pair's limits.
        limit: BondLimit,
        /// The atoms' distance (Å).
        distance: f32,
        /// Whether `distance` is within `limit`.
        bonded: bool,
    },
}

impl BondDecision {
    /// Whether the atoms are bonded.
    pub fn is_bonded(&self) -> bool {
        matches!(self, Self::ByDistance { bonded: true, .. })
    }
}

impl std::fmt::Display for BondDecision {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::LikeCharges => f.write_str("not bonded: charges of the same sign"),
            Self::NoMetalValence => f.write_str("not bonded: a metal without valence left for it"),
            Self::NotTogether => f.write_str(
                "not bonded: different configurations of a disordered part, never present together",
            ),
            Self::ByDistance {
                limit,
                distance,
                bonded,
            } => write!(
                f,
                "{}: {distance:.3} Å, {} limit {:.3}–{:.3} Å",
                if *bonded { "bonded" } else { "not bonded" },
                limit.kind,
                limit.min,
                limit.max
            ),
        }
    }
}

/// Decides whether two atoms are bonded, by their type symbols (`Fe3+`, `O2-`, `Na`), whether
/// they have valence left for bonds to metals, their sites' occupancies and their distance.
///
/// Two ions whose type symbols name the same charge sign aren't bonded, unless they are of one
/// element and `like_charges_of_one_element` allows it (the O–O of a peroxide, the S–S of pyrite).
/// Two metals aren't bonded unless both have valence left for it (see
/// [`BondAtom::bonds_to_metals`]). Two atoms of different configurations of a disordered part
/// aren't bonded (see [`crate::disorder::apart`]): the orientations of a molecule disordered
/// over several lie at bond lengths from one another. Neither are alternative positions of one
/// atom (see [`crate::occupancy::alternatives`]). Otherwise, a pair of elements the table lists
/// is bonded within the table's limits, and any other pair by [`BondRule::default_limit`]. Either upper
/// limit is widened by `tolerance` (0.1 is 10 % longer).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct BondRule {
    /// The bond table, if one is loaded.
    pub table: Option<Arc<BondTable>>,
    /// How much longer than the upper limit a bond may be, as a fraction of it.
    pub tolerance: f32,
    /// Whether two ions of one element and the same charge sign may be bonded.
    pub like_charges_of_one_element: bool,
}

impl BondRule {
    /// The shortest distance (Å) at which a pair the table doesn't list is bonded. Two atoms
    /// closer than this are one atom placed twice or alternative positions of one atom (a split
    /// site), never two atoms bonded.
    pub const MIN_DISTANCE: f32 = 0.1;

    /// The smallest bond valence (v.u.) of a cation–anion bond the table doesn't list: such a
    /// pair is bonded up to the distance where its valence falls to this.
    pub const SMALLEST_BOND_VALENCE: f32 = 0.05;

    /// How much longer than the sum of their covalent radii two atoms of a pair that is neither a
    /// cation–anion nor a metal–metal pair (and the table doesn't list) may be apart to be bonded.
    pub const COVALENT_FACTOR: f32 = 1.15;

    /// How much longer than the sum of their metallic radii two metal atoms may be apart to be
    /// bonded. Every elemental metal gets its nearest shell (bcc 8, fcc and hcp 12) and no more:
    /// the shell lies 8–13 % inside the limit, bcc's second shell (1.15 times as far) at least
    /// 3.5 % outside. On 727 ordered ICSD intermetallics it bonds 95 % of each atom's first
    /// coordination shell (up to its neighbours' largest gap in distance), and leaves 0.5 % of
    /// atoms without a metal–metal bond.
    pub const METALLIC_FACTOR: f32 = 1.1;

    /// The distances at which atoms of the types `a` and `b` are bonded, by their elements alone,
    /// and the rule that gives them: the table's limits if it lists the pair, else
    /// [`Self::default_limit`]'s. The upper limit is widened by `tolerance`.
    pub fn limit(&self, a: &str, b: &str) -> BondLimit {
        let (a, b) = (element_of(a), element_of(b));

        let (kind, min, max) = match self.table.as_ref().and_then(|table| table.limits(&a, &b)) {
            Some((min, max)) => (LimitKind::Table, min, max),
            None => {
                let (kind, max) = Self::default_kind_and_limit(&a, &b);
                (kind, Self::MIN_DISTANCE, max)
            }
        };

        BondLimit {
            kind,
            min,
            max: max * (1.0 + self.tolerance),
        }
    }

    /// The longest distance (Å) at which atoms of the elements `a` and `b` are bonded where the
    /// table doesn't list them.
    ///
    /// A cation–anion pair (an anion former, see [`is_anion_former`], and a less electronegative
    /// element that doesn't [bond covalently](bonds_covalently)) is bonded while its bond valence
    /// `exp((R0 - R) / b)` is at least [`Self::SMALLEST_BOND_VALENCE`], up to
    /// `R0 - b ln(SMALLEST_BOND_VALENCE)` with the longest R0 any of the cation's oxidation states
    /// gives (crystallib's `BOND_VALENCE`). Two [metals](is_metal) are bonded up to
    /// [`Self::METALLIC_FACTOR`] times the sum of their metallic radii. Any other pair (C–C, C–H,
    /// O–H, Sn–Sn, a metal and a non-metal that aren't a cation–anion pair) is bonded up to
    /// [`Self::COVALENT_FACTOR`] times the sum of their covalent radii.
    pub fn default_limit(a: &str, b: &str) -> f32 {
        Self::default_kind_and_limit(a, b).1
    }

    /// [`Self::default_limit`] and the rule it comes from.
    fn default_kind_and_limit(a: &str, b: &str) -> (LimitKind, f32) {
        let ionic = Self::ionic_limit(a, b).or_else(|| Self::ionic_limit(b, a));
        let metallic =
            || Some((METALLIC_RADII.get(a)? + METALLIC_RADII.get(b)?) * Self::METALLIC_FACTOR);

        if let Some(limit) = ionic {
            (LimitKind::Ionic, limit)
        } else if let Some(limit) = metallic() {
            (LimitKind::Metallic, limit)
        } else {
            let radii = COVALENT_RADII.get(a) + COVALENT_RADII.get(b);
            (LimitKind::Covalent, radii * Self::COVALENT_FACTOR)
        }
    }

    /// The limit of `cation`–`anion` as a cation–anion pair, if it is one. This decides the
    /// distance at which two elements are bonded, by the elements alone; whether an atom of a
    /// structure is a cation or an anion is its role, which also looks at its neighbours (see
    /// [`crate::oxidation::SiteChemistry::role`]).
    fn ionic_limit(cation: &str, anion: &str) -> Option<f32> {
        if !is_cation_anion_pair(cation, anion) {
            return None;
        }

        let anion_state = anion_state(anion)?;
        BOND_VALENCE
            .cation_states(cation)
            .into_iter()
            .filter_map(|state| BOND_VALENCE.parameters(cation, state, anion, anion_state))
            .map(|parameters| BondValence::distance(parameters, Self::SMALLEST_BOND_VALENCE))
            .reduce(f32::max)
    }

    /// Whether atoms `a` and `b`, `distance` (Å) apart, are bonded.
    pub fn are_bonded<'a>(
        &self,
        a: impl Into<BondAtom<'a>>,
        b: impl Into<BondAtom<'a>>,
        distance: f32,
    ) -> bool {
        self.decide(a, b, distance).is_bonded()
    }

    /// Whether atoms `a` and `b`, `distance` (Å) apart, are bonded, and what decided it.
    pub fn decide<'a>(
        &self,
        a: impl Into<BondAtom<'a>>,
        b: impl Into<BondAtom<'a>>,
        distance: f32,
    ) -> BondDecision {
        let (a, b) = (a.into(), b.into());
        let (element_a, element_b) = (element_of(a.type_symbol), element_of(b.type_symbol));
        let charges = (Charge::of(a.type_symbol), Charge::of(b.type_symbol));
        if let Some(apart) = self.apart(&element_a, &element_b, charges, (a, b)) {
            return apart;
        }
        if alternatives((&element_a, &element_b), (a, b), distance) {
            return BondDecision::NotTogether;
        }

        let limit = self.limit(a.type_symbol, b.type_symbol);
        BondDecision::ByDistance {
            limit,
            distance,
            bonded: limit.contains(distance),
        }
    }

    /// The bonds [`Self::bonds`] finds with every atom present: those of every configuration of a
    /// disordered structure at once, which join each of its pieces as a whole (see
    /// [`crate::disorder`]).
    pub fn bonds_of_any_configuration(
        &self,
        positions: &[[f32; 3]],
        atoms: &[BondAtom],
    ) -> Vec<(usize, usize, f32)> {
        let present: Vec<BondAtom> = atoms
            .iter()
            .map(|atom| BondAtom {
                place_occupancy: 1.0,
                disorder: None,
                ..*atom
            })
            .collect();
        self.bonds(positions, &present)
    }

    /// The bonds `(i, j, distance)` between `atoms` at `positions`, `i < j`, ordered by `i`, then
    /// `j`: the pairs [`BondRule::are_bonded`] accepts, found without testing every pair (see
    /// [`pairs_within`]).
    pub fn bonds(&self, positions: &[[f32; 3]], atoms: &[BondAtom]) -> Vec<(usize, usize, f32)> {
        debug_assert_eq!(positions.len(), atoms.len());

        let ElementLimits {
            elements,
            indices,
            limits,
        } = self.element_limits(atoms);
        let charges = atoms
            .iter()
            .map(|atom| Charge::of(atom.type_symbol))
            .collect::<Vec<_>>();

        pairs_within(positions, longest(&limits))
            .into_iter()
            .filter(|&(i, j, distance)| {
                let (element_i, element_j) = (indices[i], indices[j]);
                let pair = (&elements[element_i], &elements[element_j]);
                self.apart(
                    pair.0,
                    pair.1,
                    (charges[i], charges[j]),
                    (atoms[i], atoms[j]),
                )
                .is_none()
                    && !alternatives((pair.0, pair.1), (atoms[i], atoms[j]), distance)
                    && limits[element_i * elements.len() + element_j].contains(distance)
            })
            .collect()
    }

    /// The bonds of the sites at `fractional` (coordinates in `[0, 1)` of `cell`), `atoms` as the
    /// rule sees them, in the periodic structure: [`Self::bonds`] among the copies of the sites
    /// a bond can reach (see [`copies_in_reach`]), each from a site to the copy of the other it
    /// is bonded to, and listed from both of its ends. `None` for a cell without volume.
    pub fn periodic_bonds(
        &self,
        cell: &Cell,
        fractional: &[Point3<f64>],
        atoms: &[BondAtom],
    ) -> Option<Vec<PeriodicBond>> {
        debug_assert_eq!(fractional.len(), atoms.len());

        let cutoff = longest(&self.element_limits(atoms).limits);
        let copies = copies_in_reach(cell, fractional, cutoff)?;
        let copy_atoms: Vec<BondAtom> =
            copies.copies.iter().map(|&(site, _)| atoms[site]).collect();

        let mut bonds = Vec::new();
        for (p, q, _) in self.bonds(&copies.positions, &copy_atoms) {
            let ((site_p, shift_p), (site_q, shift_q)) = (copies.copies[p], copies.copies[q]);
            if shift_p == [0; 3] {
                bonds.push(PeriodicBond {
                    from: site_p,
                    to: site_q,
                    shift: shift_q,
                });
            }
            if shift_q == [0; 3] {
                bonds.push(PeriodicBond {
                    from: site_q,
                    to: site_p,
                    shift: shift_p,
                });
            }
        }
        Some(bonds)
    }

    /// The elements of `atoms`, each atom's index among them, and the limits of every pair of
    /// them.
    fn element_limits(&self, atoms: &[BondAtom]) -> ElementLimits {
        let mut elements = Vec::<String>::new();
        let indices = atoms
            .iter()
            .map(|atom| {
                let element = element_of(atom.type_symbol);
                let index = elements.iter().position(|known| *known == element);
                index.unwrap_or_else(|| {
                    elements.push(element);
                    elements.len() - 1
                })
            })
            .collect::<Vec<_>>();

        let limits = elements
            .iter()
            .flat_map(|a| elements.iter().map(move |b| (a, b)))
            .map(|(a, b)| self.limit(a, b))
            .collect::<Vec<_>>();

        ElementLimits {
            elements,
            indices,
            limits,
        }
    }

    /// What keeps atoms `a` and `b`, of the elements `element_a` and `element_b` and with
    /// `charges`, from being bonded whatever their distance, if anything: charges of the same
    /// sign (unless they are of one element and `like_charges_of_one_element` allows it), two
    /// metals of which one has no valence left to bond to metals (see
    /// [`BondAtom::bonds_to_metals`]), or atoms never present together (see
    /// [`crate::disorder::apart`]).
    fn apart(
        &self,
        element_a: &str,
        element_b: &str,
        charges: (Option<Charge>, Option<Charge>),
        (a, b): (BondAtom, BondAtom),
    ) -> Option<BondDecision> {
        let one_element = element_a == element_b;
        if charges.0.is_some()
            && charges.0 == charges.1
            && !(self.like_charges_of_one_element && one_element)
        {
            return Some(BondDecision::LikeCharges);
        }
        if is_metal(element_a) && is_metal(element_b) && !(a.bonds_to_metals && b.bonds_to_metals) {
            return Some(BondDecision::NoMetalValence);
        }
        if crate::disorder::apart(&a, &b) {
            return Some(BondDecision::NotTogether);
        }
        None
    }
}

/// The limits of the element pairs among some atoms (see [`BondRule::element_limits`]).
struct ElementLimits {
    /// The elements, each once.
    elements: Vec<String>,
    /// Each atom's element, as its index in `elements`.
    indices: Vec<usize>,
    /// The limits of each pair of `elements`, by `first * elements.len() + second`.
    limits: Vec<BondLimit>,
}

/// The longest distance (Å) any of `limits` bonds at, 0 for none.
fn longest(limits: &[BondLimit]) -> f32 {
    limits.iter().map(|limit| limit.max).fold(0.0, f32::max)
}

/// Whether atoms `a` and `b`, of the elements `elements` and `distance` (Å) apart, are
/// alternative positions of one atom: too close to be present together, on places that needn't
/// be occupied together (see [`crate::occupancy::alternatives`]).
fn alternatives(elements: (&str, &str), (a, b): (BondAtom, BondAtom), distance: f32) -> bool {
    crate::occupancy::alternatives(elements, (a.place_occupancy, b.place_occupancy), distance)
}

/// Whether the elements `cation` and `anion` are a cation–anion pair: an anion former (see
/// [`is_anion_former`]) and another, less electronegative element that doesn't
/// [bond covalently](bonds_covalently).
fn is_cation_anion_pair(cation: &str, anion: &str) -> bool {
    let less_electronegative = matches!(
        (ELECTRONEGATIVITIES.get(cation), ELECTRONEGATIVITIES.get(anion)),
        (Some(cation), Some(anion)) if cation < anion
    );
    is_anion_former(anion) && !bonds_covalently(cation) && cation != anion && less_electronegative
}

/// Whether the bond between the elements `cation` and `anion` has a bond valence: an anion former
/// and a less electronegative element other than hydrogen. Carbon counts (C–O in carbonates)
/// although [`is_cation_anion_pair`] leaves it out; hydrogen doesn't, as an O–H bond would
/// outweigh every other bond of a hydrate.
fn has_bond_valence(cation: &str, anion: &str) -> bool {
    let less_electronegative = matches!(
        (ELECTRONEGATIVITIES.get(cation), ELECTRONEGATIVITIES.get(anion)),
        (Some(cation), Some(anion)) if cation < anion
    );
    is_anion_former(anion) && cation != "H" && cation != anion && less_electronegative
}

/// The bond valence (v.u., `exp((R0 - R) / b)`, crystallib's `BOND_VALENCE`) of a bond
/// `distance` (Å) long between atoms of the types `a` and `b` (type symbols), where it has one
/// (an anion former and a less electronegative element other than hydrogen): by the cation's
/// oxidation state where its type symbol names one the table has parameters for, else by the
/// state giving the largest valence. `None` for any other pair (covalent or metallic bonds,
/// bonds of hydrogen) and where the table has no
/// parameters.
pub fn bond_valence(a: &str, b: &str, distance: f32) -> Option<f32> {
    cation_bond_valence(a, b, distance).map(|(valence, _)| valence)
}

/// [`bond_valence`], with which of the two atoms is the cation: `true` for `a`.
pub fn cation_bond_valence(a: &str, b: &str, distance: f32) -> Option<(f32, bool)> {
    let (element_a, element_b) = (element_of(a), element_of(b));
    let (cation, cation_symbol, anion, a_is_cation) = if has_bond_valence(&element_a, &element_b) {
        (element_a, a, element_b, true)
    } else if has_bond_valence(&element_b, &element_a) {
        (element_b, b, element_a, false)
    } else {
        return None;
    };

    let anion_state = anion_state(&anion)?;
    let valence = |state: i8| {
        BOND_VALENCE
            .parameters(&cation, state, &anion, anion_state)
            .map(|parameters| BondValence::valence(parameters, distance))
    };

    let stated = crate::type_symbol::state_of(cation_symbol)
        .filter(|state| state.fract() == 0.0)
        .and_then(|state| valence(state as i8));

    stated
        .or_else(|| {
            BOND_VALENCE
                .cation_states(&cation)
                .into_iter()
                .filter_map(valence)
                .reduce(f32::max)
        })
        .map(|valence| (valence, a_is_cation))
}

#[cfg(test)]
mod tests {
    use super::*;

    const TABLE: &str = "element1,element2,min,max,count\nNa,O,1.9,2.9,491\nCl,Na,2.6,3.2,72\n";

    #[test]
    fn a_cation_anion_bond_has_a_bond_valence_and_others_none() {
        let titanium = bond_valence("Ti4+", "O2-", 2.0).unwrap();
        let barium = bond_valence("O2-", "Ba2+", 2.83).unwrap();
        assert!((0.5..0.75).contains(&titanium), "{titanium}");
        assert!((0.15..0.3).contains(&barium), "{barium}");
        assert!(
            bond_valence("Ba", "O", 2.83).is_some(),
            "without charges, by the cation's states"
        );

        let carbonate = bond_valence("C4+", "O2-", 1.28).unwrap();
        assert!((1.1..1.6).contains(&carbonate), "{carbonate}");
        assert_eq!(bond_valence("H1+", "O2-", 0.97), None);

        assert_eq!(bond_valence("C", "C", 1.54), None);
        assert_eq!(bond_valence("Fe", "Fe", 2.48), None);
        assert_eq!(bond_valence("O2-", "O2-", 2.8), None);
    }

    fn rule(tolerance: f32) -> BondRule {
        BondRule {
            table: Some(Arc::new(BondTable::from_csv(TABLE).expect("valid table"))),
            tolerance,
            like_charges_of_one_element: false,
        }
    }

    #[test]
    fn a_listed_pair_is_bonded_within_its_limits_in_either_order() {
        let rule = rule(0.0);

        assert!(rule.are_bonded("Na", "O", 2.4));
        assert!(rule.are_bonded("O", "Na", 2.4));
        assert!(!rule.are_bonded("Na", "O", 1.5));
        assert!(!rule.are_bonded("Na", "O", 3.0));
    }

    #[test]
    fn tolerance_widens_the_upper_limit() {
        assert!(!rule(0.0).are_bonded("Na", "Cl", 3.4));
        assert!(rule(0.1).are_bonded("Na", "Cl", 3.4));
        assert_eq!(rule(0.1).limit("Na", "Cl").min, 2.6);
    }

    #[test]
    fn a_pair_the_table_lacks_takes_the_default_limit() {
        assert_eq!(
            rule(0.0).limit("Hg", "S"),
            BondLimit {
                kind: LimitKind::Ionic,
                min: BondRule::MIN_DISTANCE,
                max: BondRule::default_limit("Hg", "S")
            }
        );
        assert_eq!(
            BondRule::default().limit("S2-", "Hg2+"),
            rule(0.0).limit("Hg", "S")
        );
        assert_eq!(rule(0.0).limit("Na", "Cl").kind, LimitKind::Table);
    }

    #[test]
    fn a_decision_names_what_decided_it() {
        let rule = BondRule::default();

        assert_eq!(rule.decide("Na1+", "Na1+", 3.0), BondDecision::LikeCharges);
        assert_eq!(
            rule.decide(
                BondAtom {
                    type_symbol: "Na",
                    bonds_to_metals: false,
                    place_occupancy: 1.0,
                    disorder: None,
                },
                "Na",
                3.0
            ),
            BondDecision::NoMetalValence
        );
        let kind = |a: &str, b: &str, distance: f32| match rule.decide(a, b, distance) {
            BondDecision::ByDistance { limit, bonded, .. } => Some((limit.kind, bonded)),
            _ => None,
        };
        assert_eq!(kind("Na", "Cl", 2.82), Some((LimitKind::Ionic, true)));
        assert_eq!(kind("Fe", "Fe", 2.48), Some((LimitKind::Metallic, true)));
        assert_eq!(kind("C", "C", 2.52), Some((LimitKind::Covalent, false)));
        assert_eq!(
            rule.decide("Na", "Cl", 2.82).to_string(),
            format!(
                "bonded: 2.820 Å, ionic limit 0.100–{:.3} Å",
                BondRule::default_limit("Na", "Cl")
            )
        );
    }

    #[test]
    fn atoms_of_different_configurations_or_alternative_positions_are_not_bonded() {
        use crate::disorder::Disorder;

        let rule = BondRule::default();
        let with = |type_symbol, place_occupancy, disorder| BondAtom {
            place_occupancy,
            disorder,
            ..BondAtom::from(type_symbol)
        };
        let configuration = |group| Some(Disorder { assembly: 0, group });

        // Carbonate in two configurations, each half the time.
        let carbon = with("C", 0.5, configuration(Some(0)));
        assert!(rule.are_bonded(carbon, with("O", 0.5, configuration(Some(0))), 1.28));
        assert_eq!(
            rule.decide(carbon, with("O", 0.5, configuration(Some(1))), 1.28),
            BondDecision::NotTogether
        );

        // Configurations that couldn't be told apart: bonded only where the places must be
        // occupied together.
        let unresolved =
            |type_symbol, place_occupancy| with(type_symbol, place_occupancy, configuration(None));
        assert_eq!(
            rule.decide(unresolved("C", 0.0834), unresolved("N", 0.0834), 1.42),
            BondDecision::NotTogether
        );
        assert!(rule.are_bonded(unresolved("C", 0.6), unresolved("N", 0.6), 1.48));

        // Alternative positions of one atom, whatever their configurations.
        assert_eq!(
            rule.decide(carbon, with("C", 0.5, configuration(Some(0))), 0.4),
            BondDecision::NotTogether
        );
        assert_eq!(
            rule.decide(with("C", 0.5, None), with("C", 0.5, None), 0.4),
            BondDecision::NotTogether
        );

        let positions = [[0.0, 0.0, 0.0], [1.28, 0.0, 0.0], [1.68, 0.0, 0.0]];
        let atoms = [
            carbon,
            with("O", 0.5, configuration(Some(0))),
            with("O", 0.5, configuration(Some(1))),
        ];
        assert_eq!(
            rule.bonds(&positions, &atoms)
                .iter()
                .map(|&(i, j, _)| (i, j))
                .collect::<Vec<_>>(),
            [(0, 1)]
        );
    }

    #[test]
    fn a_cation_anion_pair_is_bonded_down_to_a_small_bond_valence() {
        let rule = BondRule::default();

        // Rock salt: Na–Cl 2.82 Å, the next Na–Cl 4.88 Å.
        assert!(rule.are_bonded("Na", "Cl", 2.82));
        assert!(!rule.are_bonded("Na", "Cl", 4.88));
        // Cinnabar: Hg–S 2.37 Å.
        assert!(rule.are_bonded("Hg", "S", 2.37));
        // An iron oxide's first shell, Fe–O 2.0 Å, not its O at 3.4 Å.
        assert!(rule.are_bonded("Fe", "O", 2.0));
        assert!(!rule.are_bonded("Fe", "O", 3.4));
    }

    #[test]
    fn any_other_pair_is_bonded_by_covalent_radii() {
        let rule = BondRule::default();

        // Diamond: C–C 1.54 Å, the second neighbours 2.52 Å.
        assert!(rule.are_bonded("C", "C", 1.54));
        assert!(!rule.are_bonded("C", "C", 2.52));
        // A hydroxide's O–H 0.97 Å, not a hydrogen bond's H···O 1.8 Å.
        assert!(rule.are_bonded("H", "O", 0.97));
        assert!(!rule.are_bonded("D", "O", 1.8));
        // β-Sn, a metal without a metallic radius: its nearest Sn–Sn 3.02 Å.
        assert!(rule.are_bonded("Sn", "Sn", 3.02));
    }

    #[test]
    fn metals_are_bonded_to_their_nearest_shell_by_metallic_radii() {
        let rule = BondRule::default();

        // fcc Al: 12 at 2.86 Å, the next 6 at 4.05 Å.
        assert!(rule.are_bonded("Al", "Al", 2.86));
        assert!(!rule.are_bonded("Al", "Al", 4.05));
        // hcp Zn: 6 at 2.67 Å and 6 at 2.91 Å.
        assert!(rule.are_bonded("Zn", "Zn", 2.67));
        assert!(rule.are_bonded("Zn", "Zn", 2.91));
        // bcc Fe: 8 at 2.48 Å, not the 6 at 2.87 Å.
        assert!(rule.are_bonded("Fe", "Fe", 2.48));
        assert!(!rule.are_bonded("Fe", "Fe", 2.87));
        // Ni3Al: Ni–Ni and Ni–Al 2.52 Å.
        assert!(rule.are_bonded("Ni", "Ni", 2.52));
        assert!(rule.are_bonded("Ni", "Al", 2.52));
    }

    #[test]
    fn columns_are_found_by_their_header() {
        let table = BondTable::from_csv("max,element2,element1,min\n2.9,O,Na,1.9\n").unwrap();
        assert_eq!(table.limits("O", "Na"), Some((1.9, 2.9)));
    }

    #[test]
    fn a_malformed_table_names_its_line() {
        let error = |text: &str| BondTable::from_csv(text).unwrap_err();

        assert!(error("").contains("empty"));
        assert!(error("element1,element2,min\n").contains("max"));
        assert!(error("element1,element2,min,max\nNa,O,2.9,1.9\n").contains("line 2"));
        assert!(error("element1,element2,min,max\nNa,O,x,1.9\n").contains("line 2"));
        assert!(error("element1,element2,min,max\nNa,O,1,2\nO,Na,1,2\n").contains("twice"));
        assert!(error("element1,element2,min,max\n").contains("no element pairs"));
    }

    #[test]
    fn limits_are_looked_up_by_element_and_types_with_charges_find_them() {
        let rule = rule(0.0);

        assert!(rule.are_bonded("Na1+", "O2-", 2.4));
        assert_eq!(rule.limit("Na1+", "O2-"), rule.limit("Na", "O"));
    }

    #[test]
    fn ions_of_one_charge_sign_bond_only_within_an_element_where_allowed() {
        let table = BondTable::from_csv("element1,element2,min,max\nFe,Fe,2.6,3.8\nO,O,1.2,1.6\n");
        let mut rule = BondRule {
            table: Some(Arc::new(table.unwrap())),
            ..Default::default()
        };

        assert!(!rule.are_bonded("Fe2+", "Fe3+", 2.9));
        assert!(!rule.are_bonded("O1-", "O1-", 1.5));
        assert!(
            rule.are_bonded("Fe", "Fe3+", 2.9),
            "an uncharged type isn't filtered"
        );

        rule.like_charges_of_one_element = true;
        assert!(rule.are_bonded("O1-", "O1-", 1.5));
        assert!(rule.are_bonded("Fe2+", "Fe3+", 2.9));
        assert!(
            !rule.are_bonded("Na1+", "Fe3+", 2.9),
            "only within one element"
        );
    }

    #[test]
    fn bonds_are_the_pairs_are_bonded_accepts() {
        let positions = crate::geometry::scattered_positions(300, 10.0);
        let kinds = ["Na1+", "Cl1-", "O2-", "Fe3+", "Fe2+", "C", "Ne", "Al", "Ni"];
        // Every third atom has no valence left for bonds to metals.
        let atoms = (0..positions.len())
            .map(|index| BondAtom {
                type_symbol: kinds[index % kinds.len()],
                bonds_to_metals: index % 3 != 0,
                place_occupancy: 1.0,
                disorder: None,
            })
            .collect::<Vec<_>>();
        let table = "element1,element2,min,max\nCl,Na,2.2,3.4\nFe,O,1.6,2.4\nC,C,1.1,1.7\n\
                     Fe,Fe,2.4,3.0\nO,O,1.2,1.6\n";

        for like_charges_of_one_element in [false, true] {
            let rule = BondRule {
                table: Some(Arc::new(BondTable::from_csv(table).unwrap())),
                tolerance: 0.1,
                like_charges_of_one_element,
            };

            let mut expected = Vec::new();
            for i in 0..positions.len() {
                for j in i + 1..positions.len() {
                    let [a, b] = [positions[i], positions[j]];
                    let distance =
                        ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2))
                            .sqrt();
                    if rule.are_bonded(atoms[i], atoms[j], distance) {
                        expected.push((i, j));
                    }
                }
            }

            let found = rule
                .bonds(&positions, &atoms)
                .into_iter()
                .map(|(i, j, _)| (i, j))
                .collect::<Vec<_>>();

            assert!(!expected.is_empty());
            assert_eq!(
                found, expected,
                "like charges of one element: {like_charges_of_one_element}"
            );
        }
    }

    /// A cell with the edges `a`, `b` and `c` (Å) and right angles.
    fn orthorhombic(a: f64, b: f64, c: f64) -> Cell {
        let mut cell = Cell::default();
        (cell.a, cell.b, cell.c) = (a.into(), b.into(), c.into());
        (cell.alpha, cell.beta, cell.gamma) = (90.0.into(), 90.0.into(), 90.0.into());
        cell
    }

    #[test]
    fn periodic_bonds_reach_the_copies_in_the_neighbouring_cells() {
        // Two carbons 1.4 Å apart across the face at a = 0, 2.6 Å apart within the cell.
        let cell = orthorhombic(4.0, 4.0, 4.0);
        let fractional = [Point3::new(0.05, 0.5, 0.5), Point3::new(0.7, 0.5, 0.5)];
        let atoms = [BondAtom::from("C"), BondAtom::from("C")];

        let mut bonds = rule(0.0)
            .periodic_bonds(&cell, &fractional, &atoms)
            .unwrap();
        bonds.sort_by_key(|bond| (bond.from, bond.to));

        assert_eq!(
            bonds,
            [
                PeriodicBond {
                    from: 0,
                    to: 1,
                    shift: [-1, 0, 0]
                },
                PeriodicBond {
                    from: 1,
                    to: 0,
                    shift: [1, 0, 0]
                },
            ]
        );
    }

    #[test]
    fn periodic_bonds_reach_as_far_as_the_longest_bond_in_a_small_cell() {
        // Na+ and Cl− half a 1.2 Å cell apart along a: Cl− lies 0.6, 1.8 and 3.0 Å from Na+, and
        // only 3.0 Å is within the table's Na–Cl limits (2.6–3.2 Å), two and a half cells away.
        let cell = orthorhombic(1.2, 10.0, 10.0);
        let fractional = [Point3::new(0.0, 0.5, 0.5), Point3::new(0.5, 0.5, 0.5)];
        let atoms = [BondAtom::from("Na+"), BondAtom::from("Cl-")];

        let mut bonds = rule(0.0)
            .periodic_bonds(&cell, &fractional, &atoms)
            .unwrap();
        bonds.sort_by_key(|bond| (bond.from, bond.shift));

        let bond = |from, to, shift| PeriodicBond { from, to, shift };
        assert_eq!(
            bonds,
            [
                bond(0, 1, [-3, 0, 0]),
                bond(0, 1, [2, 0, 0]),
                bond(1, 0, [-2, 0, 0]),
                bond(1, 0, [3, 0, 0]),
            ]
        );
    }
}
