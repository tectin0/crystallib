//! What kind of bond two elements make, from their electronegativities: where the pair lies in
//! the van Arkel–Ketelaar triangle of bond types (ionic, covalent, metallic).
//!
//! The classes follow G. Sproul, J. Chem. Educ. 78, 387 (2001), who sorted 312 binary compounds
//! of known bond type by their two electronegativities: a pair whose more electronegative
//! element lies below a first limit bonds metallically; otherwise one whose less electronegative
//! element lies below a second limit bonds ionically, and the rest covalently (2 % of the
//! covalent, 6 % of the ionic and none of the metallic compounds wrongly placed). His limits are
//! about 2.2 and 1.7 on Pauling's scale (2.195 between P and H), 13.5 and 9.5 eV of valence
//! electron energy: 2.28 and 1.61 on Allen's scale, which is that energy times 0.169. Covalent
//! bonds are nonpolar below [`NONPOLAR_DIFFERENCE`]. Sproul notes that for compounds of the
//! transition metals, whose electronegativity changes with their oxidation state, the classes
//! overlap: there they are a rough guide only.
//!
//! The ionic share is Pauling's, `1 - exp(-Δχ² / 4)` (The Nature of the Chemical Bond, 1960).

use crate::{ALLEN_ELECTRONEGATIVITIES, ELECTRONEGATIVITIES};

use crate::elements::element_of;

/// The largest electronegativity difference at which a covalent bond counts as nonpolar (the
/// textbook convention, on Pauling's scale).
pub const NONPOLAR_DIFFERENCE: f32 = 0.4;

/// A scale of electronegativity.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum ElectronegativityScale {
    /// Pauling's, from bond energies (as updated by Allred).
    #[default]
    Pauling,
    /// Allen's, the average energy of a free atom's valence electrons, scaled to Pauling's.
    Allen,
}

impl ElectronegativityScale {
    /// Every scale.
    pub const ALL: [Self; 2] = [Self::Pauling, Self::Allen];

    /// The electronegativity of `element` on this scale, if the scale gives one.
    pub fn of(self, element: &str) -> Option<f32> {
        match self {
            Self::Pauling => ELECTRONEGATIVITIES.get(element),
            Self::Allen => ALLEN_ELECTRONEGATIVITIES.get(element),
        }
    }

    /// The electronegativity below which a pair's more electronegative element makes the bond
    /// metallic.
    pub fn metallic_limit(self) -> f32 {
        match self {
            Self::Pauling => 2.195,
            Self::Allen => 2.28,
        }
    }

    /// The electronegativity below which a pair's less electronegative element makes a bond
    /// that isn't metallic ionic.
    pub fn ionic_limit(self) -> f32 {
        match self {
            Self::Pauling => 1.7,
            Self::Allen => 1.61,
        }
    }
}

impl std::fmt::Display for ElectronegativityScale {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::Pauling => "Pauling",
            Self::Allen => "Allen",
        })
    }
}

/// The kind of bond a pair of elements makes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BondClass {
    /// Electrons given from one atom to the other.
    Ionic,
    /// Electrons shared unequally.
    PolarCovalent,
    /// Electrons shared about equally.
    NonpolarCovalent,
    /// Electrons delocalised over the structure.
    Metallic,
}

impl std::fmt::Display for BondClass {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::Ionic => "ionic",
            Self::PolarCovalent => "polar covalent",
            Self::NonpolarCovalent => "nonpolar covalent",
            Self::Metallic => "metallic",
        })
    }
}

/// Where a pair of elements lies in the bond triangle.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BondCharacter {
    /// Its class, by Sproul's limits.
    pub class: BondClass,
    /// The difference of the two electronegativities.
    pub difference: f32,
    /// Their mean.
    pub mean: f32,
    /// Pauling's share of ionic character, 0 to 1.
    pub ionic_share: f32,
}

/// The character of a bond between atoms of the types `a` and `b` (type symbols or elements) on
/// `scale`; `None` where the scale gives one of the elements no electronegativity.
pub fn character(a: &str, b: &str, scale: ElectronegativityScale) -> Option<BondCharacter> {
    let (a, b) = (scale.of(&element_of(a))?, scale.of(&element_of(b))?);
    let (lower, higher) = (a.min(b), a.max(b));
    let difference = higher - lower;

    let class = if higher < scale.metallic_limit() {
        BondClass::Metallic
    } else if lower < scale.ionic_limit() {
        BondClass::Ionic
    } else if difference <= NONPOLAR_DIFFERENCE {
        BondClass::NonpolarCovalent
    } else {
        BondClass::PolarCovalent
    };

    Some(BondCharacter {
        class,
        difference,
        mean: (a + b) / 2.0,
        ionic_share: 1.0 - (-difference * difference / 4.0).exp(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn class(a: &str, b: &str, scale: ElectronegativityScale) -> BondClass {
        character(a, b, scale).unwrap().class
    }

    #[test]
    fn the_triangles_corners_and_sides_get_their_classes_on_both_scales() {
        for scale in ElectronegativityScale::ALL {
            assert_eq!(class("Cs", "F", scale), BondClass::Ionic, "{scale}");
            assert_eq!(class("Na+", "Cl-", scale), BondClass::Ionic, "{scale}");
            assert_eq!(class("Na", "H", scale), BondClass::Ionic, "{scale}");
            assert_eq!(class("Fe", "Fe", scale), BondClass::Metallic, "{scale}");
            assert_eq!(class("Cu", "Zn", scale), BondClass::Metallic, "{scale}");
            assert_eq!(
                class("C", "C", scale),
                BondClass::NonpolarCovalent,
                "{scale}"
            );
            assert_eq!(
                class("C", "H", scale),
                BondClass::NonpolarCovalent,
                "{scale}"
            );
            assert_eq!(
                class("Si4+", "O2-", scale),
                BondClass::PolarCovalent,
                "{scale}"
            );
            assert_eq!(class("B", "F", scale), BondClass::PolarCovalent, "{scale}");
        }
    }

    #[test]
    fn titanium_oxygen_is_ionic_but_less_so_than_barium_oxygen() {
        let pauling = ElectronegativityScale::Pauling;
        let titanium = character("Ti", "O", pauling).unwrap();
        let barium = character("Ba", "O", pauling).unwrap();

        assert_eq!(titanium.class, BondClass::Ionic);
        assert_eq!(barium.class, BondClass::Ionic);
        assert!((titanium.ionic_share - 0.59).abs() < 0.01, "{titanium:?}");
        assert!((barium.ionic_share - 0.80).abs() < 0.01, "{barium:?}");
    }

    #[test]
    fn an_element_without_a_value_on_the_scale_has_no_character() {
        assert!(character("La", "O", ElectronegativityScale::Allen).is_none());
        assert!(character("La", "O", ElectronegativityScale::Pauling).is_some());
    }
}
