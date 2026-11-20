//! Oxidation states guessed from bond valence sums, for a phase whose type symbols name none,
//! and what the bond rule and the polyhedra take from each atom's environment.
//!
//! The bond valence of a cation–anion contact at distance `R` is `exp((R0 - R) / b)`, and the
//! valences around an ion sum to about its oxidation state (Brown's valence sum rule).
//! [`with_guessed_oxidation_states`] guesses in these steps:
//! 1. Roles: each atom is a cation, an anion or takes no part, by the neighbour of another
//!    element nearest to it (`role`, `nearest_partner`): hydrogen and carbon only next to their
//!    partners, a neighbour that is an alternative position of the atom (split or shared sites)
//!    or a channel anion doesn't count.
//! 2. Candidates: the states each site may take (`candidates`): an anion its state as an anion
//!    (−1 in an X₂ pair, two atoms of its kind as close as the bond rule bonds them), a cation
//!    the states it is known in and the midpoints between them.
//! 3. Contacts and sums: each site's contacts with the other role, down to a bond valence of
//!    `SMALLEST_CONTACT_VALENCE`, and its sum over them with crystallib's `BOND_VALENCE`
//!    parameters (`find_contacts`, `site_sum`).
//! 4. Assignment: of the combinations that leave the cell neutral (`CHARGE_TOLERANCE`), the one
//!    whose states fit the sums best, with rare states and midpoints costing extra
//!    (`assign_states`, `penalty`).
//! 5. Limits: the states are used only if the global instability index is at most
//!    `LARGEST_INSTABILITY`, so metallic and covalent structures get no charges (`guess`; the
//!    debug log says why a phase gets none).
//!
//! [`site_chemistry`] takes steps 1–4 also where the guess fails, for each atom's role (a
//! polyhedron center's charge where its type symbol names none) and whether it has valence left
//! to bond to other metals (the bond rule). The element facts the steps use come from
//! [`crate::elements`]; the results reach [`crate::bonding`] with each atom, as its type symbol's
//! charge and as [`crate::bonding::BondAtom::bonds_to_metals`].

use std::borrow::Cow;
use std::collections::BTreeMap;

use crate::{
    BOND_VALENCE, BondValence, COVALENT_RADII, ELECTRONEGATIVITIES, OXIDATION_STATES, Phase,
};
use cgmath::{EuclideanSpace, InnerSpace, Vector3};
use std::sync::Mutex;

use crate::bonding::BondRule;
use crate::elements::{
    AnionGroup, anion_group, anion_state, element_of, forms_covalent_hydrides, is_anion_former,
    is_carbon_partner, is_hydrogen_partner, is_metal, is_noble_gas, valence,
};
use crate::occupancy;
use crate::type_symbol::{Charge, with_state};

/// How far (Å) a site's nearest neighbour of another element may be for it to take a role.
const ROLE_REACH: f32 = 4.0;

/// How far the cell's charge may be from zero, per unit of weight (multiplicity times
/// occupancy) of its sites, for it to count as neutral.
const CHARGE_TOLERANCE: f32 = 0.01;

/// What a midpoint state (mixed valence, `Fe2.5+`) costs a site beyond its misfit, in v.u.²
/// per unit of weight (see [`penalty`]).
const MIDPOINT_PENALTY: f32 = 0.05;

/// What a state a site doesn't commonly take costs it beyond its misfit, in v.u.² per unit of
/// weight (see [`penalty`]).
const RARE_STATE_PENALTY: f32 = 0.5;

/// The smallest bond valence (v.u.) a contact adds to a site's sum with. It is smaller than
/// the one down to which a pair is bonded (`BondRule::SMALLEST_BOND_VALENCE`): weak contacts
/// still add to the sums, but aren't bonds.
const SMALLEST_CONTACT_VALENCE: f32 = 0.01;

/// The largest global instability index (v.u.) of a structure whose oxidation states are used.
/// On two ICSD samples of 5,000 structures, guesses up to it get the sign wrong on 0.1–0.3 % of
/// sites; up to 0.7 on up to 0.6 %, for 3 % more structures with charges per 0.05. A strained
/// structure reaches about 0.4 (cubic BaTiO3); beyond the limit, the sums miss bonds
/// cation–anion contacts don't show (the S–S of pyrite gives 2.1 with Brown's table, 0.64
/// without it),
/// or the structure is far from the conditions the parameters describe (cinnabar, HgS, gives
/// 0.08 at ambient pressure and 1.1 with its cell 6 % smaller, as at about 8 GPa).
const LARGEST_INSTABILITY: f32 = 0.55;

/// The oxidation state of each atom (`None` for one that takes no part), or `None` if no states
/// fit.
type Guess = Option<Vec<Option<f32>>>;

/// The phase last guessed, its result and the chemistry of its atoms: the same phase is often
/// asked for many times in a row.
static LAST_GUESS: Mutex<Option<(Phase, Guess, Vec<SiteChemistry>)>> = Mutex::new(None);

/// `phase`, with its sites' type symbols naming the oxidation states its bond valence sums give
/// (`Fe` → `Fe3+`). It is returned unchanged (`Cow::Borrowed`) if any type symbol names a
/// charge, if it has no cation or anion, or if no oxidation states fit its sums and neutrality
/// well enough; `Cow::Owned` exactly when its charges are guessed.
pub fn with_guessed_oxidation_states(phase: &Phase) -> Cow<'_, Phase> {
    if phase
        .atoms
        .iter()
        .any(|atom| Charge::of(&atom.type_).is_some())
    {
        return Cow::Borrowed(phase);
    }

    let (states, chemistry) = {
        let mut last = LAST_GUESS
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        match last.as_ref() {
            Some((known, states, chemistry)) if known == phase => {
                (states.clone(), chemistry.clone())
            }
            _ => {
                let analysis = analyse(phase);
                let states = guess(phase, analysis.fit);
                *last = Some((phase.clone(), states.clone(), analysis.chemistry.clone()));
                (states, analysis.chemistry)
            }
        }
    };

    let guessed = match states {
        Some(states) => {
            let mut phase = phase.clone();
            for (atom, state) in phase.atoms.0.iter_mut().zip(states) {
                if let Some(state) = state {
                    atom.type_ =
                        with_state(&crate::elements::letters(&atom.type_), f64::from(state));
                }
            }
            Cow::Owned(phase)
        }
        None => Cow::Borrowed(phase),
    };

    let mut last = LAST_SITE_CHEMISTRY
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
    if last.as_ref().is_none_or(|(known, _)| *known != *guessed) {
        *last = Some((guessed.clone().into_owned(), chemistry));
    }
    guessed
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum Role {
    Cation,
    Anion,
}

/// A site of the phase that takes part in the guess.
#[derive(Debug)]
struct Site {
    element: String,
    role: Role,
    /// The fraction of its positions an atom occupies.
    occupancy: f32,
    /// Multiplicity times occupancy: how much of the cell's charge one unit of its state is.
    weight: f32,
    /// The states it may take, in increasing magnitude; a cation's also with the midpoint
    /// between each two adjacent ones, for a site of mixed valence.
    candidates: Vec<f32>,
    /// The states it commonly takes: a cation's that Brown's table gives parameters for (all
    /// whole candidates if it lists none). Any other state costs [`RARE_STATE_PENALTY`].
    common: Vec<f32>,
    /// The state it starts from: an anion's only one, a cation's lowest common one (its lowest
    /// candidate where it has no common one).
    usual: f32,
    /// The sites of the other role it touches: their index and the distance (Å) of each contact.
    contacts: Vec<(usize, f32)>,
}

/// The oxidation state of each atom of `phase`, or `None` for an atom that takes no part, if its
/// `fit` is good enough (neutral, [`LARGEST_INSTABILITY`]); else `None`.
fn guess(phase: &Phase, fit: Result<Fit, NoFit>) -> Guess {
    let fit = match fit {
        Ok(fit) => fit,
        Err(no_fit) => {
            log::debug!("Oxidation states: {no_fit}");
            return None;
        }
    };

    if fit.charge.abs() > fit.charge_tolerance {
        log::debug!("Oxidation states: no neutral assignment fits the bond valence sums");
        return None;
    }
    if fit.instability > LARGEST_INSTABILITY {
        log::debug!(
            "Oxidation states: global instability index {:.2} v.u. is above {LARGEST_INSTABILITY}",
            fit.instability
        );
        return None;
    }

    log::debug!(
        "Oxidation states from bond valence sums, global instability index {:.2} v.u.: {}",
        fit.instability,
        phase
            .atoms
            .iter()
            .zip(&fit.states)
            .zip(&fit.sums)
            .filter_map(|((atom, state), sum)| Some(format!(
                "{} {:+} (sum {sum:.2})",
                atom.label,
                (*state)?
            )))
            .collect::<Vec<_>>()
            .join(", ")
    );
    Some(fit.states)
}

/// The oxidation states that fit a phase's bond valence sums best, and how well.
#[derive(Debug)]
struct Fit {
    /// Each atom's state, or `None` for one that takes no part.
    states: Vec<Option<f32>>,
    /// Each atom's bond valence sum (0 for one that takes no part).
    sums: Vec<f32>,
    /// The global instability index (v.u.): the root of the weighted mean square of each site's
    /// sum minus its state.
    instability: f32,
    /// The cell's charge under the states.
    charge: f32,
    /// How far from zero the charge may be for the cell to count as neutral.
    charge_tolerance: f32,
}

/// Why no states fit a phase.
#[derive(Debug, Clone, Copy, PartialEq)]
enum NoFit {
    /// Its positions in the cell can't be generated.
    NoCell,
    NoCation,
    NoAnion,
    /// No combination of candidate states leaves the cell neutral.
    NotNeutral,
}

impl std::fmt::Display for NoFit {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::NoCell => "the positions in the cell can't be generated",
            Self::NoCation => "no atom is a cation",
            Self::NoAnion => "no atom is an anion",
            Self::NotNeutral => "no combination of candidate states is neutral",
        })
    }
}

/// How well `states` fit the bond valence sums of `sites`.
fn fit_of(sites: &[Option<Site>], states: Vec<f32>) -> Fit {
    let total_weight: f32 = sites.iter().flatten().map(|site| site.weight).sum();

    Fit {
        sums: sums(sites, &states),
        instability: misfit(sites, &states).sqrt(),
        charge: charge(sites, &states),
        charge_tolerance: CHARGE_TOLERANCE * total_weight,
        states: sites
            .iter()
            .zip(states)
            .map(|(site, state)| site.as_ref().map(|_| state))
            .collect(),
    }
}

/// One pass over a phase's sites: the states that fit them best, for the guess, and what the
/// bond rule and the polyhedra take from each atom.
struct Analysis {
    fit: Result<Fit, NoFit>,
    chemistry: Vec<SiteChemistry>,
}

/// The [`Analysis`] of `phase`. A metal's valence left takes the sums with the states that fit
/// best, or with the usual ones where no combination is neutral.
fn analyse(phase: &Phase) -> Analysis {
    let Some(cell) = Cell::of(phase) else {
        return Analysis {
            fit: Err(NoFit::NoCell),
            chemistry: vec![SiteChemistry::UNKNOWN; phase.atoms.len()],
        };
    };
    let (elements, roles, partners) = roles(phase, &cell);

    let (fit, sums) = match sites_with_roles(phase, &cell, &elements, &roles) {
        Err(no_fit) => (Err(no_fit), None),
        Ok(sites) => match assign_states(&sites) {
            Some(states) => {
                let fit = fit_of(&sites, states);
                let sums = fit.sums.clone();
                (Ok(fit), Some(sums))
            }
            None => {
                let usual: Vec<f32> = sites
                    .iter()
                    .map(|site| site.as_ref().map_or(0.0, |site| site.usual))
                    .collect();
                (Err(NoFit::NotNeutral), Some(sums(&sites, &usual)))
            }
        },
    };

    let chemistry: Vec<SiteChemistry> = (0..elements.len())
        .map(|i| SiteChemistry {
            role: roles[i].map(|role| match role {
                Role::Cation => Charge::Positive,
                Role::Anion => Charge::Negative,
            }),
            bonds_to_metals: !is_metal(&elements[i]) || {
                let given = sums.as_ref().map_or(0.0, |sums| sums[i]);
                f32::from(valence(&elements[i])) - given > VALENCE_FOR_METAL_BONDS
            },
        })
        .collect();

    if log::log_enabled!(log::Level::Debug) {
        let label = |i: usize| phase.atoms.0[i].label.as_str();
        let sites = (0..elements.len()).map(|i| {
            let role = match roles[i] {
                Some(Role::Cation) => "cation",
                Some(Role::Anion) => "anion",
                None => "no role",
            };
            let partner = partners[i].map_or(String::new(), |(j, distance)| {
                format!(" by {} at {distance:.2} Å", label(j))
            });
            let metals = if chemistry[i].bonds_to_metals {
                ""
            } else {
                ", no valence left for metals"
            };
            format!("{} {role}{partner}{metals}", label(i))
        });
        log::debug!("Roles: {}", sites.collect::<Vec<_>>().join("; "));
    }

    Analysis { fit, chemistry }
}

/// What the bond rule and the polyhedra take from one atom's environment, whether or not
/// oxidation states fit the phase.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SiteChemistry {
    /// The sign of the atom's role (step 1 of the [module's](self) guess): positive for a
    /// cation, negative for an anion, `None` for one that takes no part.
    pub role: Option<Charge>,
    /// Whether the atom may bond to other metals. An atom's bonding electrons are shared out
    /// over its bonds (the valence sum rule): a metal's [`valence`] less the bond valence it
    /// gives to the anions around it is what it has left for bonds to metals, and it needs more
    /// than [`VALENCE_FOR_METAL_BONDS`] of it. So Na in rock salt has none left and no Na–Na
    /// bonds, while Mo in Mo6Cl12 or any atom of an intermetallic does. An atom that isn't a
    /// [metal](is_metal) is always `true`.
    pub bonds_to_metals: bool,
}

impl SiteChemistry {
    /// An atom whose environment is unknown: no role, free to bond to metals.
    const UNKNOWN: Self = Self {
        role: None,
        bonds_to_metals: true,
    };
}

/// The [`SiteChemistry`] of each atom of `phase`. A phase whose charges were just guessed (see
/// [`with_guessed_oxidation_states`]) takes the chemistry of that guess's pass; the charges
/// don't change it, as it is read from the elements alone.
pub fn site_chemistry(phase: &Phase) -> Vec<SiteChemistry> {
    let mut last = LAST_SITE_CHEMISTRY
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
    if let Some((known, chemistry)) = last.as_ref()
        && known == phase
    {
        return chemistry.clone();
    }

    let chemistry = analyse(phase).chemistry;
    *last = Some((phase.clone(), chemistry.clone()));
    chemistry
}

/// The valence (v.u.) a metal needs to have left after its bonds to anions to bond to other
/// metals (see [`SiteChemistry::bonds_to_metals`]). In ICSD structures with anions whose
/// oxidation states can't be guessed, it bonds 85 % of the metal–metal contacts (within 1.1
/// times the metallic radii) where the ICSD states a metal as 0, and 13 % of those where it
/// states both as cations.
pub const VALENCE_FOR_METAL_BONDS: f32 = 0.5;

/// The phase last given to [`site_chemistry`] and its result: the same phase is often asked for
/// many times in a row.
static LAST_SITE_CHEMISTRY: Mutex<Option<(Phase, Vec<SiteChemistry>)>> = Mutex::new(None);

/// The sites of `phase` that take part in the guess, with their contacts, given its cell and
/// each atom's element and role (see [`roles`]).
fn sites_with_roles(
    phase: &Phase,
    cell: &Cell,
    elements: &[String],
    roles: &[Option<Role>],
) -> Result<Vec<Option<Site>>, NoFit> {
    if !roles.contains(&Some(Role::Cation)) {
        return Err(NoFit::NoCation);
    }
    if !roles.contains(&Some(Role::Anion)) {
        return Err(NoFit::NoAnion);
    }

    let occupancies = occupancies(phase);

    let mut sites: Vec<Option<Site>> = (0..elements.len())
        .map(|i| {
            let role = roles[i]?;
            let candidates = site_candidates(i, elements, &occupancies, role, cell)?;
            let occupancy = occupancies[i];
            Some(Site {
                element: elements[i].clone(),
                role,
                occupancy,
                weight: cell.positions[i].len() as f32 * occupancy,
                candidates: candidates.states,
                common: candidates.common,
                usual: candidates.usual,
                contacts: Vec::new(),
            })
        })
        .collect();

    find_contacts(&mut sites, cell);
    Ok(sites)
}

/// The occupancy of each atom of `phase`.
fn occupancies(phase: &Phase) -> Vec<f32> {
    phase
        .atoms
        .iter()
        .map(|atom| atom.occupancy.value as f32)
        .collect()
}

/// The neighbour that decides an atom's role (see [`nearest_partner`]): its index and distance (Å).
type Partner = (usize, f32);

/// The element of each atom of `phase`, its role (see [`role`]) and the neighbour that decided
/// it.
fn roles(phase: &Phase, cell: &Cell) -> (Vec<String>, Vec<Option<Role>>, Vec<Option<Partner>>) {
    let elements: Vec<String> = phase
        .atoms
        .iter()
        .map(|atom| element_of(&atom.type_))
        .collect();
    let electronegativities: Vec<Option<f32>> = elements
        .iter()
        .map(|element| {
            (!is_noble_gas(element))
                .then(|| ELECTRONEGATIVITIES.get(element))
                .flatten()
        })
        .collect();

    let occupancies = occupancies(phase);

    let partners: Vec<Option<Partner>> = (0..elements.len())
        .map(|i| nearest_partner(i, &elements, &electronegativities, &occupancies, cell))
        .collect();
    let roles = (0..elements.len())
        .map(|i| role(i, &elements, &electronegativities, partners[i]))
        .collect();
    (elements, roles, partners)
}

/// The neighbour of atom `i` that decides its role (see [`role`]): its index and distance (Å).
/// It is the nearest one of another element with an electronegativity within [`ROLE_REACH`],
/// if atom `i` has an electronegativity itself. A neighbour whose occupancy adds to at most 1
/// with the atom's doesn't count if it is an alternative position of the atom (see
/// [`occupancy::alternatives`]), or if both are anion formers and they are farther apart
/// than the sum of their covalent radii: two partly occupied anions sharing a channel (Cl and F
/// in an apatite) aren't bonded, while the atoms of an oxyanion that is disordered as a whole
/// are closer.
fn nearest_partner(
    i: usize,
    elements: &[String],
    electronegativities: &[Option<f32>],
    occupancies: &[f32],
    cell: &Cell,
) -> Option<Partner> {
    electronegativities[i]?;

    let skipped = |j: usize, distance: f32| {
        let elements = (elements[i].as_str(), elements[j].as_str());
        let occupancies = (occupancies[i], occupancies[j]);
        let radii = COVALENT_RADII.get(elements.0) + COVALENT_RADII.get(elements.1);
        let both_anion_formers = is_anion_former(elements.0) && is_anion_former(elements.1);
        let sharing_a_channel = !occupancy::together(occupancies.0, occupancies.1)
            && both_anion_formers
            && distance > radii;
        occupancy::alternatives(elements, occupancies, distance) || sharing_a_channel
    };
    cell.neighbours(i, ROLE_REACH, |j| {
        electronegativities[j].is_some() && elements[j] != elements[i]
    })
    .into_iter()
    .filter(|&(j, distance)| !skipped(j, distance))
    .min_by(|a, b| a.1.total_cmp(&b.1))
}

/// The role of atom `i`, from its nearest neighbour of another element with an
/// electronegativity:
/// - hydrogen is a cation next to a [partner](is_hydrogen_partner), a hydride next to a less
///   electronegative element that doesn't [bond to it covalently](forms_covalent_hydrides), and
///   takes no part otherwise;
/// - carbon is a cation next to a [partner](is_carbon_partner) and takes no part otherwise;
/// - any other atom is a cation if less electronegative than that neighbour, an anion if more
///   and an [anion former](is_anion_former).
///
/// An atom without an electronegativity, without such a neighbour (see [`nearest_partner`]), or
/// as electronegative as it, has none.
fn role(
    i: usize,
    elements: &[String],
    electronegativities: &[Option<f32>],
    partner: Option<Partner>,
) -> Option<Role> {
    let own = electronegativities[i]?;
    let (nearest, _) = partner?;
    let partner = elements[nearest].as_str();
    let other = electronegativities[nearest]?;

    match elements[i].as_str() {
        "H" if is_hydrogen_partner(partner) => Some(Role::Cation),
        "H" => (other < own && !forms_covalent_hydrides(partner)).then_some(Role::Anion),
        "C" => is_carbon_partner(partner).then_some(Role::Cation),
        _ if own < other => Some(Role::Cation),
        element if own > other && is_anion_former(element) => Some(Role::Anion),
        _ => None,
    }
}

/// The candidate states of atom `i` in `role` (see [`candidates`]): an anion is in a dimer if
/// an atom of its own element is within the distance the bond rule bonds the two at
/// ([`BondRule::default_limit`]), on a site that can be occupied together with its own: their
/// occupancies add to more than 1. The positions of a split site (or of two sites sharing one
/// atom) lie that close and never are.
fn site_candidates(
    i: usize,
    elements: &[String],
    occupancies: &[f32],
    role: Role,
    cell: &Cell,
) -> Option<Candidates> {
    let element = &elements[i];
    let in_dimer = role == Role::Anion && {
        let reach = BondRule::default_limit(element, element);
        !cell
            .neighbours(i, reach, |j| {
                elements[j] == *element && occupancy::together(occupancies[i], occupancies[j])
            })
            .is_empty()
    };
    candidates(element, role, in_dimer)
}

/// The candidate states of a site of `element` in `role` and the one it starts from, or `None`
/// if it has none. An anion has one: its [state as an anion](anion_state), or `in_dimer` what
/// its [group](AnionGroup) takes in an X₂ pair (−1 for a chalcogen or halogen, none for a
/// pnictogen or a hydride's H). A cation takes every state it is known in or Brown's table
/// lists, and the midpoint between each two adjacent ones.
fn candidates(element: &str, role: Role, in_dimer: bool) -> Option<Candidates> {
    match role {
        Role::Anion => {
            let state = match (in_dimer, anion_group(element)) {
                (true, Some(AnionGroup::Chalcogen)) => -1,
                (true, Some(AnionGroup::Pnictogen) | None) => return None,
                _ => anion_state(element)?,
            };
            let state = f32::from(state);
            Some(Candidates {
                states: vec![state],
                common: vec![state],
                usual: state,
            })
        }
        Role::Cation => {
            // The states Brown's table gives parameters for, or all known ones where it lists
            // none of the element.
            let common: Vec<f32> = BOND_VALENCE
                .cation_states(element)
                .into_iter()
                .map(f32::from)
                .collect();

            let mut known: Vec<f32> = OXIDATION_STATES
                .get(element)
                .iter()
                .filter(|&&state| state > 0)
                .map(|&state| f32::from(state))
                .chain(common.iter().copied())
                .collect();
            known.sort_by(f32::total_cmp);
            known.dedup();

            let midpoints: Vec<f32> = known
                .windows(2)
                .filter(|pair| pair[1] - pair[0] == 1.0)
                .map(|pair| (pair[0] + pair[1]) / 2.0)
                .collect();
            let mut states = known;
            states.extend(midpoints);
            states.sort_by(f32::total_cmp);

            let usual = *common.first().or(states.first())?;
            Some(Candidates {
                states,
                common,
                usual,
            })
        }
    }
}

/// A site's candidate states (see [`Site`]).
struct Candidates {
    states: Vec<f32>,
    common: Vec<f32>,
    usual: f32,
}

/// A phase's atoms in its unit cell, with the geometry to find their periodic neighbours.
struct Cell {
    /// The distinct Cartesian positions (Å) in the unit cell of each atom of the phase.
    positions: Vec<Vec<Vector3<f32>>>,
    /// The Cartesian cell vectors a, b and c.
    axes: [Vector3<f32>; 3],
    /// The spacing (Å) of the lattice planes parallel to each pair of axes: bc, ca, ab.
    spacings: [f32; 3],
}

impl Cell {
    fn of(phase: &Phase) -> Option<Self> {
        let transformation = phase.cell.transformation_matrix_y_up();
        let cartesian = |f: Vector3<f64>| (transformation * crate::Cell::y_up_input(f)).truncate();

        let positions = phase
            .atoms
            .iter()
            .map(|atom| {
                Some(
                    crate::cell_atoms::site_positions(phase, atom)?
                        .into_iter()
                        .map(|point| cartesian(point.to_vec()))
                        .collect(),
                )
            })
            .collect::<Option<Vec<_>>>()?;

        let axes = [
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ]
        .map(cartesian);
        let volume = axes[0].dot(axes[1].cross(axes[2])).abs();
        let spacings = [
            axes[1].cross(axes[2]),
            axes[2].cross(axes[0]),
            axes[0].cross(axes[1]),
        ]
        .map(|face| volume / face.magnitude());

        Some(Self {
            positions,
            axes,
            spacings,
        })
    }

    /// The atoms `j` with `accept(j)` within `reach` (Å) of atom `i`'s first position, over as
    /// many neighbouring cells as `reach` spans: each with the distance of every such position.
    fn neighbours(
        &self,
        i: usize,
        reach: f32,
        accept: impl Fn(usize) -> bool,
    ) -> Vec<(usize, f32)> {
        let Some(&center) = self.positions[i].first() else {
            return Vec::new();
        };
        let images = self.spacings.map(|spacing| (reach / spacing).ceil() as i32);

        let mut found = Vec::new();
        for (j, positions) in self.positions.iter().enumerate() {
            if !accept(j) {
                continue;
            }
            for &position in positions {
                let shifts = (-images[0]..=images[0]).flat_map(|a| {
                    (-images[1]..=images[1])
                        .flat_map(move |b| (-images[2]..=images[2]).map(move |c| (a, b, c)))
                });
                for (a, b, c) in shifts {
                    let image = position
                        + self.axes[0] * a as f32
                        + self.axes[1] * b as f32
                        + self.axes[2] * c as f32;
                    let distance = (image - center).magnitude();
                    if distance > BondRule::MIN_DISTANCE && distance < reach {
                        found.push((j, distance));
                    }
                }
            }
        }
        found
    }
}

/// Fills in each site's contacts: the sites of the other role closer than the distance at which
/// their bond valence falls below [`SMALLEST_CONTACT_VALENCE`] for any of their states.
fn find_contacts(sites: &mut [Option<Site>], cell: &Cell) {
    let reach = |cation: &Site, anion: &Site| {
        cation
            .candidates
            .iter()
            .flat_map(|a| anion.candidates.iter().map(move |b| (a, b)))
            .filter_map(|(&state, &anion_state)| {
                let parameters = parameters(&cation.element, state, &anion.element, anion_state)?;
                Some(BondValence::distance(parameters, SMALLEST_CONTACT_VALENCE))
            })
            .fold(0.0, f32::max)
    };

    let contacts: Vec<Vec<(usize, f32)>> = (0..sites.len())
        .map(|i| {
            let Some(site) = &sites[i] else {
                return Vec::new();
            };
            let other_role = |j: usize| sites[j].as_ref().filter(|other| other.role != site.role);

            let largest_reach = (0..sites.len())
                .filter_map(&other_role)
                .map(|other| match site.role {
                    Role::Cation => reach(site, other),
                    Role::Anion => reach(other, site),
                })
                .fold(0.0, f32::max);

            cell.neighbours(i, largest_reach, |j| other_role(j).is_some())
                .into_iter()
                .filter(|&(j, distance)| {
                    let other = other_role(j).expect("only sites of the other role are found");
                    let limit = match site.role {
                        Role::Cation => reach(site, other),
                        Role::Anion => reach(other, site),
                    };
                    distance < limit
                })
                .collect()
        })
        .collect();

    for (site, contacts) in sites.iter_mut().zip(contacts) {
        if let Some(site) = site {
            site.contacts = contacts;
        }
    }
}

/// The bond valence parameters of a contact; a midpoint state takes those of the lower state.
fn parameters(cation: &str, state: f32, anion: &str, anion_state: f32) -> Option<(f32, f32)> {
    BOND_VALENCE.parameters(cation, state.floor() as i8, anion, anion_state as i8)
}

/// Each site's bond valence sum, with every site in the state `states` gives it.
fn sums(sites: &[Option<Site>], states: &[f32]) -> Vec<f32> {
    (0..sites.len())
        .map(|i| site_sum(sites, states, i, states[i]))
        .collect()
}

/// The bond valence sum of site `i` in `state`, with every other site in the state `states`
/// gives it; 0 for a site that takes no part. A contact counts as much as its neighbour's
/// occupancy: a position held half the time contributes half its valence.
fn site_sum(sites: &[Option<Site>], states: &[f32], i: usize, state: f32) -> f32 {
    let Some(site) = &sites[i] else {
        return 0.0;
    };

    site.contacts
        .iter()
        .filter_map(|&(j, distance)| {
            let other = sites[j].as_ref()?;
            let parameters = match site.role {
                Role::Cation => parameters(&site.element, state, &other.element, states[j])?,
                Role::Anion => parameters(&other.element, states[j], &site.element, state)?,
            };
            Some(other.occupancy * BondValence::valence(parameters, distance))
        })
        .sum()
}

/// How badly `states` fit: the weighted mean square of each site's sum minus its state.
fn misfit(sites: &[Option<Site>], states: &[f32]) -> f32 {
    let sums = sums(sites, states);
    let (squares, weights) = sites
        .iter()
        .enumerate()
        .filter_map(|(i, site)| Some((i, site.as_ref()?)))
        .fold((0.0, 0.0), |(squares, weights), (i, site)| {
            let deviation = sums[i] - states[i].abs();
            (
                squares + site.weight * deviation * deviation,
                weights + site.weight,
            )
        });
    squares / weights
}

/// What `state` costs `site` (v.u.², per unit of weight) beyond its misfit: a midpoint state
/// [`MIDPOINT_PENALTY`], a state the site doesn't commonly take (a midpoint next to one)
/// [`RARE_STATE_PENALTY`]. So the fit takes whole and common states unless the sums clearly
/// call for others.
fn penalty(site: &Site, state: f32) -> f32 {
    let is_common = |state: f32| site.common.contains(&state);
    let is_midpoint = state.fract() != 0.0;
    let is_rare = if is_midpoint {
        !(is_common(state.floor()) && is_common(state.ceil()))
    } else {
        !is_common(state)
    };

    let mut penalty = 0.0;
    if is_midpoint {
        penalty += MIDPOINT_PENALTY;
    }
    if is_rare {
        penalty += RARE_STATE_PENALTY;
    }
    penalty
}

/// The cell's charge under `states`.
fn charge(sites: &[Option<Site>], states: &[f32]) -> f32 {
    sites
        .iter()
        .zip(states)
        .filter_map(|(site, state)| Some(site.as_ref()?.weight * state))
        .sum()
}

/// The oxidation state of each site (0 for one that takes no part): of all combinations of
/// candidates that leave the cell neutral, the one whose states fit the sums best, or `None` if
/// no combination is neutral.
///
/// A site's sum depends on the states of its neighbours (through the parameters of each
/// contact), so the sums are taken with every neighbour in the state that fits it alone best,
/// settled over two rounds; the combination is then chosen on those sums. The search goes site
/// by site, keeping for each reachable charge of the cell (in steps of `CHARGE_STEP`) only the
/// cheapest way there, and dropping charges the remaining sites can't bring back to zero. A
/// state's cost is its misfit plus its penalty (see [`penalty`]).
fn assign_states(sites: &[Option<Site>]) -> Option<Vec<f32>> {
    let candidates = |site: &Site| -> Vec<f32> { site.candidates.clone() };

    let mut settled: Vec<f32> = sites
        .iter()
        .map(|site| site.as_ref().map_or(0.0, |site| site.usual))
        .collect();

    for _ in 0..2 {
        for i in 0..sites.len() {
            let Some(site) = &sites[i] else {
                continue;
            };
            let deviation = |state: f32| (site_sum(sites, &settled, i, state) - state.abs()).abs();
            settled[i] = candidates(site)
                .into_iter()
                .min_by(|&a, &b| deviation(a).total_cmp(&deviation(b)))?;
        }
    }

    const CHARGE_STEP: f32 = 0.01;
    let to_key = |charge: f32| (charge / CHARGE_STEP).round() as i64;

    let total_weight: f32 = sites.iter().flatten().map(|site| site.weight).sum();
    let tolerance = to_key(CHARGE_TOLERANCE * total_weight);

    let order: Vec<usize> = (0..sites.len()).filter(|&i| sites[i].is_some()).collect();

    // The least and greatest charge the sites from each point of `order` on can add.
    let mut remaining = vec![(0_i64, 0_i64); order.len() + 1];
    for (k, &i) in order.iter().enumerate().rev() {
        let site = sites[i].as_ref().expect("order holds sites only");
        let keys: Vec<i64> = candidates(site)
            .into_iter()
            .map(|state| to_key(site.weight * state))
            .collect();
        let (least, greatest) = (*keys.iter().min()?, *keys.iter().max()?);
        remaining[k] = (remaining[k + 1].0 + least, remaining[k + 1].1 + greatest);
    }

    // For each site in `order`, the cheapest way to each charge: its cost, the charge before
    // the site, and the site's state. Charges are visited in order, so of equally cheap ways the
    // one from the lowest charge and the earliest candidate wins, the same in every run.
    let mut layers: Vec<BTreeMap<i64, (f32, i64, f32)>> = Vec::with_capacity(order.len());
    let mut reached: BTreeMap<i64, f32> = BTreeMap::from([(0, 0.0)]);

    for (k, &i) in order.iter().enumerate() {
        let site = sites[i].as_ref().expect("order holds sites only");
        let options: Vec<(f32, i64, f32)> = candidates(site)
            .into_iter()
            .map(|state| {
                let deviation = site_sum(sites, &settled, i, state) - state.abs();
                (
                    state,
                    to_key(site.weight * state),
                    site.weight * (deviation * deviation + penalty(site, state)),
                )
            })
            .collect();

        let (least, greatest) = remaining[k + 1];
        let mut layer: BTreeMap<i64, (f32, i64, f32)> = BTreeMap::new();
        for (&key, &cost) in &reached {
            for &(state, step, state_cost) in &options {
                let next = key + step;
                if next + least > tolerance || next + greatest < -tolerance {
                    continue;
                }
                let cost = cost + state_cost;
                if layer.get(&next).is_none_or(|known| cost < known.0) {
                    layer.insert(next, (cost, key, state));
                }
            }
        }

        reached = layer
            .iter()
            .map(|(&key, &(cost, _, _))| (key, cost))
            .collect();
        layers.push(layer);
    }

    // The cheapest neutral charge; of equally cheap ones, the one nearest zero.
    let mut key = reached
        .into_iter()
        .filter(|(key, _)| key.abs() <= tolerance)
        .min_by(|a, b| a.1.total_cmp(&b.1).then(a.0.abs().cmp(&b.0.abs())))?
        .0;

    let mut states = vec![0.0; sites.len()];
    for (layer, &i) in layers.iter().zip(&order).rev() {
        let &(_, previous, state) = layer
            .get(&key)
            .expect("each reached charge has a way there");
        states[i] = state;
        key = previous;
    }

    Some(states)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::type_symbol::state_of;

    /// A cubic cell of edge `a` (Å) in space group `group`, with the sites `sites` (`label type
    /// x y z occupancy` per line).
    fn cubic(a: f32, group: u8, sites: &str) -> String {
        format!(
            "data_test
_cell_length_a {a}
_cell_length_b {a}
_cell_length_c {a}
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_cell_volume {}
_space_group_IT_number {group}
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
{sites}
",
            a * a * a
        )
    }

    fn cubic_batio3() -> String {
        cubic(
            4.0,
            221,
            "Ba1 Ba 0 0 0 1\nTi1 Ti 0.5 0.5 0.5 1\nO1 O 0.5 0.5 0 1",
        )
    }

    fn phase(text: &str) -> Phase {
        crate::cif_phase::try_phase_from_cif_bytes(text.as_bytes())
            .expect("the CIF describes a phase")
            .1
    }

    fn guessed_types(text: &str) -> Vec<String> {
        with_guessed_oxidation_states(&phase(text))
            .atoms
            .iter()
            .map(|atom| atom.type_.clone())
            .collect()
    }

    #[test]
    fn strained_cubic_batio3_gets_its_ions() {
        assert_eq!(guessed_types(&cubic_batio3()), ["Ba2+", "Ti4+", "O2-"]);
    }

    #[test]
    fn the_central_atom_of_an_oxyanion_is_a_cation() {
        // Sodium nitrate (calcite type, R-3c in hexagonal axes): N–O 1.23 Å, Na–O 2.40 Å.
        let text = "data_NaNO3
_cell_length_a 5.070
_cell_length_b 5.070
_cell_length_c 16.82
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 120
_cell_volume 374.4
_space_group_IT_number 167
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Na1 Na 0 0 0 1
N1 N 0 0 0.25 1
O1 O 0.2429 0 0.25 1
";

        assert_eq!(guessed_types(text), ["Na1+", "N5+", "O2-"]);
    }

    #[test]
    fn the_silicon_of_a_silicate_is_a_cation() {
        // Ideal β-cristobalite (Fd-3m, origin choice 1): Si on 8a, O on 16c, Si–O 1.60 Å. This
        // takes Brown's table (the feature bvparm2020): the element parameters alone put R0 of
        // Si–O at 1.70 Å instead of 1.62, the sums too high to fit.
        let text = cubic(7.40, 227, "Si1 Si 0 0 0 1\nO1 O 0.125 0.125 0.125 1");

        let expected = if cfg!(feature = "bvparm2020") {
            ["Si4+", "O2-"]
        } else {
            ["Si", "O"]
        };
        assert_eq!(guessed_types(&text), expected);
    }

    #[test]
    fn a_shared_site_counts_each_atom_by_its_occupancy() {
        // Rock salt with sodium and potassium sharing the cation site half and half.
        let text = cubic(
            6.0,
            225,
            "Na1 Na 0 0 0 0.5\nK1 K 0 0 0 0.5\nCl1 Cl 0.5 0.5 0.5 1",
        );

        assert_eq!(guessed_types(&text), ["Na1+", "K1+", "Cl1-"]);
    }

    #[test]
    fn each_oxygen_of_a_peroxide_holds_one_charge() {
        // ZnO2 in the pyrite type: O–O pairs 1.48 Å long.
        let text = cubic(4.871, 205, "Zn1 Zn 0 0 0 1\nO1 O 0.4125 0.4125 0.4125 1");

        assert_eq!(guessed_types(&text), ["Zn2+", "O1-"]);
    }

    fn bonds_to_metals(phase: &Phase) -> Vec<bool> {
        site_chemistry(phase)
            .iter()
            .map(|site| site.bonds_to_metals)
            .collect()
    }

    #[test]
    fn a_metal_bonds_to_metals_only_with_valence_left_after_its_anions() {
        // Rock salt NaCl: Na gives its one electron to Cl, none is left for Na–Na.
        let salt = phase(&cubic(5.64, 225, "Na1 Na 0 0 0 1\nCl1 Cl 0.5 0.5 0.5 1"));
        assert_eq!(bonds_to_metals(&salt), [false, true]);

        // Rock salt TiO, a metal: Ti gives O about 2 of its 4, Ti–Ti bonds take the rest.
        let monoxide = phase(&cubic(4.18, 225, "Ti1 Ti 0 0 0 1\nO1 O 0.5 0.5 0.5 1"));
        assert_eq!(bonds_to_metals(&monoxide), [true, true]);

        // Without anions all valence is left.
        let alloy = phase(&cubic(3.57, 221, "Ni1 Ni 0 0.5 0.5 1\nAl1 Al 0 0 0 1"));
        assert_eq!(bonds_to_metals(&alloy), [true, true]);
    }

    #[test]
    fn of_equally_good_assignments_the_one_reached_from_the_lowest_charge_wins() {
        let site = |element: &str, role: Role, candidates: Vec<f32>| Site {
            element: element.to_string(),
            role,
            occupancy: 1.0,
            weight: 1.0,
            common: candidates.clone(),
            usual: candidates[0],
            candidates,
            contacts: Vec::new(),
        };
        // Two alike cations of 1 or 2 beside an anion of −3: (1, 2) and (2, 1) cost the same.
        let sites = [
            Some(site("Tl", Role::Cation, vec![1.0, 2.0])),
            Some(site("Tl", Role::Cation, vec![1.0, 2.0])),
            Some(site("N", Role::Anion, vec![-3.0])),
        ];

        assert_eq!(assign_states(&sites), Some(vec![1.0, 2.0, -3.0]));
    }

    #[test]
    fn each_site_has_the_sign_of_its_role_also_where_no_states_fit() {
        let roles = |text: &str| {
            site_chemistry(&phase(text))
                .iter()
                .map(|site| site.role)
                .collect::<Vec<_>>()
        };

        assert_eq!(
            roles(&cubic(5.64, 225, "Na1 Na 0 0 0 1\nCl1 Cl 0.5 0.5 0.5 1")),
            [Some(Charge::Positive), Some(Charge::Negative)]
        );
        // Pyrite: no states fit, Fe is still the cation and S the anion.
        assert_eq!(
            roles(&cubic(
                5.418,
                205,
                "Fe1 Fe 0 0 0 1\nS1 S 0.385 0.385 0.385 1"
            )),
            [Some(Charge::Positive), Some(Charge::Negative)]
        );
        // A hydroxide's D, 1 Å from its O, is a cation as H would be.
        assert_eq!(
            roles(&cubic(
                4.0,
                221,
                "Na1 Na 0 0 0 1\nO1 O 0.5 0.5 0.5 1\nD1 D 0.5 0.5 0.75 0.166667"
            ))[2],
            Some(Charge::Positive)
        );
        // Ni3Al: Al, less electronegative than its Ni neighbours, is the cation; Ni no anion.
        assert_eq!(
            roles(&cubic(3.57, 221, "Ni1 Ni 0 0.5 0.5 1\nAl1 Al 0 0 0 1")),
            [None, Some(Charge::Positive)]
        );
    }

    #[test]
    fn positions_of_a_split_site_are_no_x2_pair() {
        // Rock salt MgO with each oxygen spread over six positions 0.25 Å apart, a sixth of
        // the time on each: they are never occupied together.
        let text = cubic(4.21, 225, "Mg1 Mg 0 0 0 1\nO1 O 0.53 0.5 0.5 0.166667");

        assert_eq!(guessed_types(&text), ["Mg2+", "O2-"]);
    }

    #[test]
    fn an_anion_sharing_a_channel_with_a_more_electronegative_one_stays_an_anion() {
        // Cl half the time at the cube centre, F a sixth of the time at each face centre, 2.1 Å
        // from it (more than the 1.59 Å of their covalent radii): never bonded to each other.
        let shared = phase(&cubic(
            4.2,
            221,
            "Na1 Na 0 0 0 1\nCl1 Cl 0.5 0.5 0.5 0.5\nF1 F 0.5 0.5 0 0.166667",
        ));
        let cell = Cell::of(&shared).unwrap();
        assert_eq!(roles(&shared, &cell).1[1], Some(Role::Anion));

        // Fully occupied, the same Cl and F would be present together: Cl is the cation there.
        let together = phase(&cubic(
            4.2,
            221,
            "Na1 Na 0 0 0 1\nCl1 Cl 0.5 0.5 0.5 1\nF1 F 0.5 0.5 0 1",
        ));
        let cell = Cell::of(&together).unwrap();
        assert_eq!(roles(&together, &cell).1[1], Some(Role::Cation));
    }

    #[test]
    fn an_alternative_position_of_another_element_decides_no_role() {
        // Rock salt with its cation site held by Li 60 % of the time and otherwise by Na,
        // 0.37 Å off: Na is less electronegative than Li, but never there together with it.
        let text = cubic(
            5.4,
            225,
            "Li1 Li 0 0 0 0.6\nNa1 Na 0.04 0.04 0.04 0.05\nCl1 Cl 0.5 0.5 0.5 1",
        );

        assert_eq!(guessed_types(&text), ["Li1+", "Na1+", "Cl1-"]);
    }

    #[test]
    fn an_alloy_gets_no_charges() {
        let text = cubic(3.57, 221, "Al1 Al 0 0 0 1\nNi1 Ni 0 0.5 0.5 1");

        assert_eq!(guessed_types(&text), ["Al", "Ni"]);
    }

    #[test]
    fn hydrogen_next_to_a_metal_is_a_hydride() {
        let text = cubic(4.88, 225, "Na1 Na 0 0 0 1\nH1 H 0.5 0.5 0.5 1");

        assert_eq!(guessed_types(&text), ["Na1+", "H1-"]);
    }

    #[test]
    fn carbon_and_hydrogen_take_no_part() {
        let text = cubic(
            4.0,
            221,
            "C1 C 0 0 0 1\nH1 D 0.5 0.5 0.5 1\nO1 O 0.5 0.5 0 1",
        );

        assert_eq!(guessed_types(&text), ["C", "D", "O"]);
    }

    #[test]
    fn pyrite_gets_no_charges_since_its_sums_miss_the_s_s_bond() {
        let text = cubic(5.418, 205, "Fe1 Fe 0 0 0 1\nS1 S 0.385 0.385 0.385 1");

        assert_eq!(guessed_types(&text), ["Fe", "S"]);
    }

    /// Guesses the oxidation states of ICSD structures with their stated ones removed and
    /// compares, site by site. Reads `$ICSD_DATA/data/icsdid_<id>.cif` for each id in the file
    /// `$CALIBRATION_SAMPLE` and writes one CSV row per structure to `$CALIBRATION_OUT`.
    #[test]
    #[ignore = "needs a local ICSD: ICSD_DATA, CALIBRATION_SAMPLE, CALIBRATION_OUT"]
    fn calibrate_on_icsd() {
        use std::io::Write;

        let folder = std::env::var("ICSD_DATA").expect("ICSD_DATA names the ICSD folder");
        let sample = std::env::var("CALIBRATION_SAMPLE").expect("CALIBRATION_SAMPLE names the ids");
        let out = std::env::var("CALIBRATION_OUT").expect("CALIBRATION_OUT names the CSV to write");

        let ids: Vec<String> = std::fs::read_to_string(sample)
            .unwrap()
            .lines()
            .map(str::trim)
            .filter(|line| !line.is_empty())
            .map(str::to_string)
            .collect();
        let mut out = std::io::BufWriter::new(std::fs::File::create(out).unwrap());
        writeln!(
            out,
            "id,outcome,cause,sites,stated_charged,mismatched,instability,charge,charge_tolerance,\
             correct,wrong_sign,wrong_magnitude,missed,milliseconds,stated_charge,stated_tolerance,misses"
        )
        .unwrap();

        for (n, id) in ids.iter().enumerate() {
            if n % 250 == 0 {
                eprintln!("{n} of {}", ids.len());
            }

            let path = format!("{folder}/data/icsdid_{id}.cif");
            let Ok(bytes) = std::fs::read(&path) else {
                writeln!(out, "{id},missing,,,,,,,,,,,,,,").unwrap();
                continue;
            };

            let started = std::time::Instant::now();
            let result = std::panic::catch_unwind(|| {
                let (_, phase) = crate::cif_phase::try_phase_from_cif_bytes(&bytes)?;
                // ICSD writes a mixed valence state's type symbol with a whole charge (`Cr3+`)
                // and its oxidation number with the fraction (3.75), so the number goes first.
                let numbers = crate::cif::read(&String::from_utf8_lossy(&bytes))
                    .ok()
                    .and_then(|document| {
                        document
                            .blocks
                            .iter()
                            .map(crate::cif_phase::oxidation_numbers)
                            .find(|numbers| !numbers.is_empty())
                    })
                    .unwrap_or_default();
                let stated: Vec<Option<f32>> = phase
                    .atoms
                    .iter()
                    .map(|atom| {
                        numbers
                            .iter()
                            .find(|(symbol, _)| *symbol == atom.type_)
                            .map(|&(_, number)| number as f32)
                            .or_else(|| state_of(&atom.type_).map(|state| state as f32))
                    })
                    .collect();
                let mut stripped = phase.clone();
                for atom in stripped.atoms.0.iter_mut() {
                    atom.type_ = crate::elements::letters(&atom.type_);
                }
                let labels: Vec<String> =
                    phase.atoms.iter().map(|atom| atom.label.clone()).collect();

                // Sites whose positions in the cell don't number the multiplicity the CIF states.
                let mismatched = Cell::of(&phase).map_or(usize::MAX, |cell| {
                    cell.positions
                        .iter()
                        .zip(phase.atoms.iter())
                        .filter(|(positions, atom)| {
                            atom.multiplicity
                                .is_some_and(|stated| (stated - positions.len() as f64).abs() > 0.5)
                        })
                        .count()
                });

                // The cell's charge under the stated states, and the tolerance a guess gets.
                let stated_charge =
                    Cell::of(&phase).map(|cell| {
                        let weights = cell.positions.iter().zip(phase.atoms.iter()).map(
                            |(positions, atom)| {
                                positions.len() as f32 * atom.occupancy.value as f32
                            },
                        );
                        let charge: f32 = weights
                            .clone()
                            .zip(&stated)
                            .map(|(weight, stated)| weight * stated.unwrap_or(0.0))
                            .sum();
                        (charge, CHARGE_TOLERANCE * weights.sum::<f32>())
                    });
                let stated_charge = stated_charge.map_or(String::new(), |(charge, tolerance)| {
                    format!("{charge:.4},{tolerance:.4}")
                });

                let fit = analyse(&stripped).fit;

                // For a structure without a fit: each site with a stated charge, the role it got
                // and whether its stated state was among its candidates.
                let notes = match (&fit, Cell::of(&stripped)) {
                    (Err(_), Some(cell)) => {
                        let (elements, roles, _) = roles(&stripped, &cell);
                        let occupancies = occupancies(&stripped);
                        stated
                            .iter()
                            .enumerate()
                            .filter_map(|(i, stated)| {
                                let stated = stated.filter(|state| *state != 0.0)?;
                                let note = match roles[i] {
                                    None => "norole".to_string(),
                                    Some(role) => {
                                        let candidates = site_candidates(
                                            i,
                                            &elements,
                                            &occupancies,
                                            role,
                                            &cell,
                                        );
                                        let role = match role {
                                            Role::Cation => "cation",
                                            Role::Anion => "anion",
                                        };
                                        match candidates {
                                            None => format!("{role}:nocandidates"),
                                            Some(candidates)
                                                if candidates
                                                    .states
                                                    .iter()
                                                    .any(|c| (c - stated).abs() < 0.01) =>
                                            {
                                                role.to_string()
                                            }
                                            Some(_) => format!("{role}:notcandidate"),
                                        }
                                    }
                                };
                                Some(format!("{}:{}:{stated:+}:{note}", labels[i], elements[i]))
                            })
                            .collect::<Vec<_>>()
                            .join(" ")
                    }
                    _ => String::new(),
                };

                Some((stated, labels, mismatched, stated_charge, fit, notes))
            });
            let milliseconds = started.elapsed().as_secs_f64() * 1e3;

            let (stated, labels, mismatched, stated_charge, fit, notes) = match result {
                Err(_) => {
                    writeln!(out, "{id},panic,,,,,,,,,,,,{milliseconds:.1},,,").unwrap();
                    continue;
                }
                Ok(None) => {
                    writeln!(out, "{id},unreadable,,,,,,,,,,,,{milliseconds:.1},,,").unwrap();
                    continue;
                }
                Ok(Some(parts)) => parts,
            };

            let stated_charged = stated.iter().flatten().any(|state| *state != 0.0);
            let fit = match fit {
                Ok(fit) => fit,
                Err(no_fit) => {
                    writeln!(
                        out,
                        "{id},no_fit,{no_fit:?},{},{stated_charged},{mismatched},,,,,,,,\
                         {milliseconds:.1},{stated_charge},{notes}",
                        stated.len()
                    )
                    .unwrap();
                    continue;
                }
            };

            let (mut correct, mut wrong_sign, mut wrong_magnitude, mut missed) = (0, 0, 0, 0);
            let mut misses = Vec::new();
            for ((stated, guessed), label) in stated.iter().zip(&fit.states).zip(&labels) {
                let stated = stated.unwrap_or(0.0);
                match guessed {
                    Some(guessed) if (guessed - stated).abs() < 0.01 => correct += 1,
                    Some(guessed) => {
                        if guessed.signum() != stated.signum() {
                            wrong_sign += 1;
                        } else {
                            wrong_magnitude += 1;
                        }
                        misses.push(format!("{label} {stated:+}>{guessed:+}"));
                    }
                    None if stated != 0.0 => {
                        missed += 1;
                        misses.push(format!("{label} {stated:+}>none"));
                    }
                    None => {}
                }
            }

            let accepted =
                fit.charge.abs() <= fit.charge_tolerance && fit.instability <= LARGEST_INSTABILITY;
            writeln!(
                out,
                "{id},{},,{},{stated_charged},{mismatched},{:.4},{:.4},{:.4},{correct},\
                 {wrong_sign},{wrong_magnitude},{missed},{milliseconds:.1},{stated_charge},{}",
                if accepted { "accepted" } else { "rejected" },
                stated.len(),
                fit.instability,
                fit.charge,
                fit.charge_tolerance,
                misses.iter().take(6).cloned().collect::<Vec<_>>().join(" ")
            )
            .unwrap();
        }
    }
}
