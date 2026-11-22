//! What the bonding and oxidation-state rules know about elements: the element a type symbol
//! names, which elements count as metals, which form anions and how, and the few elements their
//! rules treat apart. Each fact is defined here once; the rules that use them live in
//! [`crate::bonding`] and [`crate::oxidation`].

use crate::{METALLIC_RADII, OXIDATION_STATES};

/// The letters of `text`, in order: a type symbol's element without its charge (`Fe3+` is `Fe`).
pub fn letters(text: &str) -> String {
    text.chars().filter(|c| c.is_alphabetic()).collect()
}

/// The element a type symbol names, written as elements are whatever case the type symbol uses
/// (`Fe3+`, `FE3+` and `fe` are Fe), with deuterium and tritium as hydrogen.
pub fn element_of(type_symbol: &str) -> String {
    let letters = letters(type_symbol);
    let mut chars = letters.chars();
    let element: String = match chars.next() {
        Some(first) => first
            .to_uppercase()
            .chain(chars.flat_map(char::to_lowercase))
            .collect(),
        None => String::new(),
    };

    match element.as_str() {
        "D" | "T" => "H".to_string(),
        _ => element,
    }
}

/// Whether `element` counts as a metal: it has a metallic radius (crystallib's
/// `METALLIC_RADII`). Sn, Pb, Sb and Bi have none there, so they don't.
pub fn is_metal(element: &str) -> bool {
    METALLIC_RADII.get(element).is_some()
}

/// The largest oxidation state `element` is known in (0 for one with no positive state): the
/// most bonding electrons it can give.
pub fn valence(element: &str) -> i8 {
    OXIDATION_STATES
        .get(element)
        .iter()
        .copied()
        .max()
        .unwrap_or(0)
        .max(0)
}

/// What the atoms of an anion former do when bonded to their own kind.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AnionGroup {
    /// N, P, As: bonded to their own kind (azide, polyphosphide) they get no guessed state.
    Pnictogen,
    /// O, S, Se, Te: in an X₂ pair each atom is one charge short of its usual state (O²⁻, but
    /// O⁻ in a peroxide).
    Chalcogen,
    /// F, Cl, Br, I: −1 in an X₂ pair too, their only negative state.
    Halogen,
}

/// The elements that may be anions, each with its group.
const ANION_FORMERS: [(&str, AnionGroup); 11] = [
    ("N", AnionGroup::Pnictogen),
    ("P", AnionGroup::Pnictogen),
    ("As", AnionGroup::Pnictogen),
    ("O", AnionGroup::Chalcogen),
    ("S", AnionGroup::Chalcogen),
    ("Se", AnionGroup::Chalcogen),
    ("Te", AnionGroup::Chalcogen),
    ("F", AnionGroup::Halogen),
    ("Cl", AnionGroup::Halogen),
    ("Br", AnionGroup::Halogen),
    ("I", AnionGroup::Halogen),
];

/// The group of `element` if it may be an anion.
pub fn anion_group(element: &str) -> Option<AnionGroup> {
    ANION_FORMERS
        .iter()
        .find(|(symbol, _)| *symbol == element)
        .map(|&(_, group)| group)
}

/// Whether `element` may be an anion: N, O, F, P, S, Cl, As, Se, Br, Te, I.
pub fn is_anion_former(element: &str) -> bool {
    anion_group(element).is_some()
}

/// The state `element` has as an anion: its most negative known state, or `None` if it has no
/// negative one.
pub fn anion_state(element: &str) -> Option<i8> {
    OXIDATION_STATES
        .get(element)
        .iter()
        .copied()
        .min()
        .filter(|state| *state < 0)
}

/// Whether `element` is a noble gas, which neither bonds ionically nor takes an oxidation state.
pub fn is_noble_gas(element: &str) -> bool {
    ["He", "Ne", "Ar", "Kr", "Xe", "Rn"].contains(&element)
}

/// Whether `element` is hydrogen or carbon, whose bonds to anion formers are covalent (O–H,
/// C–O): the bond rule bonds them by covalent radii rather than as cation–anion pairs, and the
/// oxidation guess makes them cations only next to their partners (see [`is_hydrogen_partner`]
/// and [`is_carbon_partner`]).
pub fn bonds_covalently(element: &str) -> bool {
    ["H", "C"].contains(&element)
}

/// Whether hydrogen (also as D and T) nearest to `element` is a cation: H in hydroxides,
/// hydrates and acids (N, O, F, Cl).
pub fn is_hydrogen_partner(element: &str) -> bool {
    ["N", "O", "F", "Cl"].contains(&element)
}

/// Whether carbon nearest to `element` is a cation: C in carbonates (O, F).
pub fn is_carbon_partner(element: &str) -> bool {
    ["O", "F"].contains(&element)
}

/// Whether `element` is less electronegative than hydrogen but bonds to it covalently (boranes,
/// silanes, phosphites: B, Si, Ge, P, As, Sb, Te), so hydrogen nearest to it is no hydride.
pub fn forms_covalent_hydrides(element: &str) -> bool {
    ["B", "Si", "Ge", "P", "As", "Sb", "Te"].contains(&element)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deuterium_and_tritium_are_hydrogen() {
        assert_eq!(element_of("D1+"), "H");
        assert_eq!(element_of("T"), "H");
        assert_eq!(element_of("Fe3+"), "Fe");
    }

    #[test]
    fn the_element_is_written_as_elements_are_whatever_case_the_symbol_uses() {
        assert_eq!(element_of("MG"), "Mg");
        assert_eq!(element_of("cl1-"), "Cl");
        assert_eq!(element_of("h"), "H");
        assert_eq!(element_of("d"), "H");
        assert_eq!(element_of(""), "");
    }

    #[test]
    fn metals_are_the_elements_with_a_metallic_radius() {
        assert!(is_metal("Fe"));
        assert!(is_metal("Al"));
        assert!(!is_metal("Sn"));
        assert!(!is_metal("O"));
    }

    #[test]
    fn anion_formers_have_a_group_and_a_negative_state() {
        assert_eq!(anion_group("O"), Some(AnionGroup::Chalcogen));
        assert_eq!(anion_group("P"), Some(AnionGroup::Pnictogen));
        assert_eq!(anion_group("I"), Some(AnionGroup::Halogen));
        assert_eq!(anion_group("Na"), None);
        for (element, _) in ANION_FORMERS {
            assert!(
                anion_state(element).is_some(),
                "{element} has no negative state"
            );
        }
        assert_eq!(anion_state("O"), Some(-2));
        assert_eq!(anion_state("He"), None);
    }

    #[test]
    fn a_metal_valence_is_its_largest_known_state() {
        assert_eq!(valence("Na"), 1);
        assert_eq!(valence("Ti"), 4);
        assert_eq!(valence("Xx"), 0);
    }
}
