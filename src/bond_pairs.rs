//! How a bond rule decides each element pair of a structure: the rule that gives the pair its
//! limits, the bonds it makes, the contacts it keeps apart whatever their distance, and the
//! nearest contact it leaves unbonded.

use crate::Phase;
use std::collections::HashMap;

use crate::bonding::{BondDecision, BondLimit, BondRule};
use crate::elements::element_of;

/// How far (Å) beyond the longest limit of a structure's pairs contacts are looked for, so that
/// each pair's nearest unbonded contact is found where it lies just outside its limit.
pub const BEYOND_LIMITS: f32 = 1.0;

/// The longest distance (Å) contacts are looked for, whatever the limits.
pub const LONGEST_CONTACT: f32 = 8.0;

/// How close (Å) two bond lengths of a pair are to count as one length.
pub const SAME_LENGTH: f32 = 0.005;

/// How a bond rule decides the contacts between atoms of one pair of elements in one structure.
#[derive(Debug, Clone, PartialEq)]
pub struct PairBonds {
    /// The two elements, in alphabetical order.
    pub elements: (String, String),
    /// The pair's limits by the rule, with its tolerance: a bond table's where it lists the pair.
    pub limit: BondLimit,
    /// Where a bond table gives the limits, the ones the built-in rules would give instead.
    pub built_in: Option<BondLimit>,
    /// The bonds in the cell.
    pub bonds: usize,
    /// The bonds' lengths (Å), shortest first, each once: lengths closer than
    /// [`SAME_LENGTH`] to the previous one count as that one.
    pub bond_lengths: Vec<f32>,
    /// The contacts in the cell within the limits that aren't bonds all the same, and why.
    pub kept_apart: Vec<(KeptApart, usize)>,
    /// The shortest contact (Å) that isn't bonded, if there is one within the distance contacts
    /// are looked for.
    pub nearest_unbonded: Option<f32>,
}

/// Why two atoms within their pair's limits aren't bonded.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum KeptApart {
    /// Their charges have the same sign.
    LikeCharges,
    /// They are metals, and one has no valence left to bond to metals.
    NoMetalValence,
    /// They are never present together (see [`crate::disorder::apart`]).
    NotTogether,
}

impl std::fmt::Display for KeptApart {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::LikeCharges => "charges of the same sign",
            Self::NoMetalValence => "a metal without valence left for metal bonds",
            Self::NotTogether => "different configurations of a disordered part",
        })
    }
}

/// How `bond_rule` decides each pair of elements of `phase` (with its charges guessed where it
/// states none, see [`crate::oxidation::with_guessed_oxidation_states`]) that has a contact
/// within the distance looked at: the longest of the pairs' limits plus [`BEYOND_LIMITS`], at
/// most [`LONGEST_CONTACT`]. Sorted by elements.
pub fn pair_bonds(phase: &Phase, bond_rule: &BondRule) -> Vec<PairBonds> {
    let atoms = crate::cell_atoms::cell_atoms(phase);
    let chemistry = crate::oxidation::site_chemistry(phase);
    let disorder = crate::disorder::of_atoms(phase, &atoms, &chemistry);
    let bond_atoms = crate::bonding::cell_bond_atoms(phase, &atoms, &chemistry, &disorder);
    let elements: Vec<String> = bond_atoms
        .iter()
        .map(|atom| element_of(atom.type_symbol))
        .collect();

    let mut distinct = elements.clone();
    distinct.sort();
    distinct.dedup();
    let longest_limit = distinct
        .iter()
        .flat_map(|a| distinct.iter().map(move |b| (a, b)))
        .map(|(a, b)| bond_rule.limit(a, b).max)
        .fold(0.0, f32::max);
    let cutoff = (longest_limit + BEYOND_LIMITS).min(LONGEST_CONTACT);

    let Some(contacts) = crate::geometry::periodic_contacts(
        &phase.cell,
        &atoms.iter().map(|atom| atom.fractional).collect::<Vec<_>>(),
        cutoff,
    ) else {
        return Vec::new();
    };

    let built_in_rule = BondRule {
        table: None,
        ..bond_rule.clone()
    };

    // Each contact is listed from both of its atoms, so counts are halved at the end.
    let mut pairs: HashMap<(String, String), PairBonds> = HashMap::default();
    let mut kept_apart: HashMap<((String, String), KeptApart), usize> = HashMap::default();
    for (from, to, distance) in contacts {
        let key = if elements[from] <= elements[to] {
            (elements[from].clone(), elements[to].clone())
        } else {
            (elements[to].clone(), elements[from].clone())
        };
        let pair = pairs.entry(key.clone()).or_insert_with(|| {
            let limit = bond_rule.limit(&key.0, &key.1);
            let from_table = limit.kind == crate::bonding::LimitKind::Table;
            PairBonds {
                elements: key.clone(),
                limit,
                built_in: from_table.then(|| built_in_rule.limit(&key.0, &key.1)),
                bonds: 0,
                bond_lengths: Vec::new(),
                kept_apart: Vec::new(),
                nearest_unbonded: None,
            }
        });

        match bond_rule.decide(bond_atoms[from], bond_atoms[to], distance) {
            BondDecision::ByDistance { bonded: true, .. } => {
                pair.bonds += 1;
                pair.bond_lengths.push(distance);
            }
            decision => {
                let reason = match decision {
                    BondDecision::LikeCharges => Some(KeptApart::LikeCharges),
                    BondDecision::NoMetalValence => Some(KeptApart::NoMetalValence),
                    BondDecision::NotTogether => Some(KeptApart::NotTogether),
                    BondDecision::ByDistance { .. } => None,
                };
                if let Some(reason) = reason
                    && pair.limit.contains(distance)
                {
                    *kept_apart.entry((key, reason)).or_default() += 1;
                }
                pair.nearest_unbonded = Some(
                    pair.nearest_unbonded
                        .map_or(distance, |nearest| nearest.min(distance)),
                );
            }
        }
    }

    for ((key, reason), count) in kept_apart {
        if let Some(pair) = pairs.get_mut(&key) {
            pair.kept_apart.push((reason, count / 2));
        }
    }
    let mut pairs: Vec<PairBonds> = pairs
        .into_values()
        .map(|mut pair| {
            pair.bonds /= 2;
            pair.kept_apart.sort();

            pair.bond_lengths.sort_by(f32::total_cmp);
            let mut lengths: Vec<f32> = Vec::new();
            for length in pair.bond_lengths.drain(..) {
                if lengths
                    .last()
                    .is_none_or(|last| length - last >= SAME_LENGTH)
                {
                    lengths.push(length);
                }
            }
            pair.bond_lengths = lengths;

            pair
        })
        .collect();
    pairs.sort_by(|a, b| a.elements.cmp(&b.elements));
    pairs
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use super::*;
    use crate::bonding::{BondTable, LimitKind};

    /// Rock salt, Fm-3m, a 5.62 Å, with its charges.
    const ROCK_SALT: &str = "data_rock_salt
_cell_length_a 5.62
_cell_length_b 5.62
_cell_length_c 5.62
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_cell_volume 177.5
_space_group_IT_number 225
loop_
_atom_type_symbol
_atom_type_oxidation_number
Na1+ 1
Cl1- -1
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Na1 Na1+ 0 0 0
Cl1 Cl1- 0.5 0.5 0.5
";

    fn phase(cif: &[u8]) -> Phase {
        let (_, phase) = crate::cif_phase::try_phase_from_cif_bytes(cif).unwrap();
        crate::oxidation::with_guessed_oxidation_states(&phase).into_owned()
    }

    fn pair<'a>(pairs: &'a [PairBonds], a: &str, b: &str) -> &'a PairBonds {
        pairs
            .iter()
            .find(|pair| pair.elements == (a.to_string(), b.to_string()))
            .unwrap_or_else(|| panic!("no {a}–{b} in {pairs:?}"))
    }

    #[test]
    fn rock_salt_bonds_each_sodium_to_six_chlorines_and_keeps_like_ions_apart() {
        let pairs = pair_bonds(&phase(ROCK_SALT.as_bytes()), &BondRule::default());

        let sodium_chlorine = pair(&pairs, "Cl", "Na");
        assert_eq!(sodium_chlorine.limit.kind, LimitKind::Ionic);
        assert_eq!(sodium_chlorine.built_in, None);
        assert_eq!(
            sodium_chlorine.bonds,
            4 * 6,
            "four formula units, six bonds each"
        );
        assert_eq!(sodium_chlorine.bond_lengths.len(), 1, "{sodium_chlorine:?}");
        assert!((sodium_chlorine.bond_lengths[0] - 2.81).abs() < 0.01);
        assert!(
            sodium_chlorine.nearest_unbonded.unwrap() > 4.5,
            "the next Na–Cl shell, {sodium_chlorine:?}"
        );

        let sodium_sodium = pair(&pairs, "Na", "Na");
        assert_eq!(sodium_sodium.bonds, 0);
        assert!((sodium_sodium.nearest_unbonded.unwrap() - 3.99).abs() < 0.05);
    }

    #[test]
    fn a_tables_limits_replace_the_built_in_ones_for_its_pairs_only() {
        let rock_salt = phase(ROCK_SALT.as_bytes());
        let table = BondTable::named("short.csv", "element1,element2,min,max\nNa,Cl,2.0,2.5\n");
        let rule = BondRule {
            table: Some(Arc::new(table.unwrap())),
            ..BondRule::default()
        };
        let pairs = pair_bonds(&rock_salt, &rule);

        let sodium_chlorine = pair(&pairs, "Cl", "Na");
        assert_eq!(sodium_chlorine.limit.kind, LimitKind::Table);
        assert_eq!(sodium_chlorine.limit.max, 2.5);
        assert_eq!(sodium_chlorine.built_in.unwrap().kind, LimitKind::Ionic);
        assert_eq!(sodium_chlorine.bonds, 0);
        assert!(sodium_chlorine.bond_lengths.is_empty());
        assert!((sodium_chlorine.nearest_unbonded.unwrap() - 2.81).abs() < 0.01);

        assert_eq!(pair(&pairs, "Na", "Na").built_in, None);
    }
}
