use std::str::FromStr;

use crate::Measured;

#[cfg_attr(feature = "serde", derive(serde::Deserialize, serde::Serialize))]
#[cfg_attr(feature = "serde", serde(default))]
#[derive(Debug, Default, Clone, PartialEq)]
/// The atoms of a phase, one per site of its asymmetric unit.
pub struct Atoms(pub Vec<Atom>);

impl std::ops::Deref for Atoms {
    type Target = Vec<Atom>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl std::ops::DerefMut for Atoms {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

#[cfg_attr(feature = "serde", derive(serde::Deserialize, serde::Serialize))]
#[cfg_attr(feature = "get-size", derive(get_size2::GetSize))]
#[allow(clippy::upper_case_acronyms)]
#[derive(Debug, Default, Clone, PartialEq)]
/// How an atom's displacement is given, as CIF's `_atom_site_adp_type` names it.
pub enum AdpType {
    /// Isotropic, as U (Å²).
    #[default]
    Uiso,
    /// Anisotropic, as U (Å²).
    Uani,
    /// Overall, as U (Å²).
    Uovl,
    /// Multipole expansion, as U (Å²).
    Umpe,
    /// Anisotropic, as B (Å²).
    Bani,
    /// Isotropic, as B (Å²).
    Biso,
    /// Overall, as B (Å²).
    Bovl,
}

impl std::fmt::Display for AdpType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            AdpType::Uiso => write!(f, "Uiso"),
            AdpType::Uani => write!(f, "Uani"),
            AdpType::Uovl => write!(f, "Uovl"),
            AdpType::Umpe => write!(f, "Umpe"),
            AdpType::Bani => write!(f, "Bani"),
            AdpType::Biso => write!(f, "Biso"),
            AdpType::Bovl => write!(f, "Bovl"),
        }
    }
}

/// A text that names no [`AdpType`].
#[derive(Debug, PartialEq, Eq)]
pub struct ParseAdpTypeError(String);

impl ParseAdpTypeError {
    /// The error for the text `s`.
    pub fn new(s: impl AsRef<str>) -> Self {
        Self(format!("Invalid ADP type: {}", s.as_ref()))
    }
}

impl std::fmt::Display for ParseAdpTypeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for ParseAdpTypeError {}

/// Reads CIF's names (`Uiso`, `Uani`, ...) and GSAS-II's `I` (isotropic) and `A` (anisotropic).
impl FromStr for AdpType {
    type Err = ParseAdpTypeError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "I" => Ok(AdpType::Uiso),
            "A" => Ok(AdpType::Uani),
            "Uiso" => Ok(AdpType::Uiso),
            "Uani" => Ok(AdpType::Uani),
            "Uovl" => Ok(AdpType::Uovl),
            "Umpe" => Ok(AdpType::Umpe),
            "Bani" => Ok(AdpType::Bani),
            "Biso" => Ok(AdpType::Biso),
            "Bovl" => Ok(AdpType::Bovl),
            _ => Err(ParseAdpTypeError::new(s)),
        }
    }
}

#[cfg(test)]
mod test_parse_adp_type {
    #[test]
    fn test() {
        assert_eq!("Uiso".parse(), Ok(super::AdpType::Uiso));
        assert_eq!("Uani".parse(), Ok(super::AdpType::Uani));
        assert_eq!("Uovl".parse(), Ok(super::AdpType::Uovl));
        assert_eq!("Umpe".parse(), Ok(super::AdpType::Umpe));
        assert_eq!("Bani".parse(), Ok(super::AdpType::Bani));
        assert_eq!("Biso".parse(), Ok(super::AdpType::Biso));
        assert_eq!("Bovl".parse(), Ok(super::AdpType::Bovl));
        assert!("".parse::<super::AdpType>().is_err());
        assert!("Uiso1".parse::<super::AdpType>().is_err());
    }
}

#[cfg_attr(feature = "serde", derive(serde::Deserialize, serde::Serialize))]
#[cfg_attr(feature = "serde", serde(default))]
#[derive(Debug, Default, Clone, PartialEq)]
/// An atom on a site of a phase, as CIF's `_atom_site_*` and `_atom_site_aniso_*` give it. More
/// of the `_atom_site` items may join it without a breaking release, so it is built from
/// [`Atom::default`] outside the crate.
#[non_exhaustive]
pub struct Atom {
    /// The site's label (`Ba1`).
    pub label: String,
    /// The atom's type: its element, possibly with a charge (`Ba`, `O2-`).
    #[cfg_attr(feature = "serde", serde(rename = "type"))]
    pub type_: String,
    /// The fractional coordinate along a.
    pub x: Measured,
    /// The fractional coordinate along b.
    pub y: Measured,
    /// The fractional coordinate along c.
    pub z: Measured,
    /// The fraction of the site the atom occupies.
    pub occupancy: Measured,
    /// How many positions the site has in the unit cell, where known.
    pub multiplicity: Option<f64>,
    /// The disorder assembly the site belongs to (`_atom_site_disorder_assembly`, `A`): sites
    /// with alternative configurations, one of which is present at a time; `None` for a single
    /// assembly or none.
    pub disorder_assembly: Option<String>,
    /// The site's disorder group within its assembly (`_atom_site_disorder_group`, `1`, `-2`):
    /// the sites of one group are present together, never with those of another group of the
    /// assembly. A negative group's copies by symmetry are groups of their own (SHELX's negative
    /// PART). `None` where the site isn't disordered.
    pub disorder_group: Option<String>,
    /// How the displacement is given.
    pub adp_type: AdpType,
    /// The isotropic displacement, or the equivalent isotropic one of an anisotropic (Å²).
    pub u_iso_or_equiv: Measured,
    /// The anisotropic displacement U11 (Å²); likewise the five below.
    #[cfg_attr(feature = "serde", serde(rename = "U11"))]
    pub u11: Measured,
    /// U22 (Å²).
    #[cfg_attr(feature = "serde", serde(rename = "U22"))]
    pub u22: Measured,
    /// U33 (Å²).
    #[cfg_attr(feature = "serde", serde(rename = "U33"))]
    pub u33: Measured,
    /// U12 (Å²).
    #[cfg_attr(feature = "serde", serde(rename = "U12"))]
    pub u12: Measured,
    /// U13 (Å²).
    #[cfg_attr(feature = "serde", serde(rename = "U13"))]
    pub u13: Measured,
    /// U23 (Å²).
    #[cfg_attr(feature = "serde", serde(rename = "U23"))]
    pub u23: Measured,
}
