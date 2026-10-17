use std::str::FromStr;

use crate::Measured;

#[cfg_attr(feature = "serde", derive(serde::Deserialize, serde::Serialize))]
#[cfg_attr(feature = "serde", serde(default))]
#[derive(Debug, Default, Clone, PartialEq)]
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
pub enum AdpType {
    #[default]
    Uiso,
    Uani,
    Uovl,
    Umpe,
    Bani,
    Biso,
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

#[derive(Debug, PartialEq, Eq)]
pub struct ParseAdpTypeError(String);

impl ParseAdpTypeError {
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
pub struct Atom {
    pub label: String,
    #[cfg_attr(feature = "serde", serde(rename = "type"))]
    pub type_: String,
    pub x: Measured,
    pub y: Measured,
    pub z: Measured,
    pub occupancy: Measured,
    pub multiplicity: Option<f64>,
    pub adp_type: AdpType,
    pub u_iso_or_equiv: Measured,
    #[cfg_attr(feature = "serde", serde(rename = "U11"))]
    pub u11: Measured,
    #[cfg_attr(feature = "serde", serde(rename = "U22"))]
    pub u22: Measured,
    #[cfg_attr(feature = "serde", serde(rename = "U33"))]
    pub u33: Measured,
    #[cfg_attr(feature = "serde", serde(rename = "U12"))]
    pub u12: Measured,
    #[cfg_attr(feature = "serde", serde(rename = "U13"))]
    pub u13: Measured,
    #[cfg_attr(feature = "serde", serde(rename = "U23"))]
    pub u23: Measured,
}
