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
#[cfg_attr(feature = "serde", serde(default))]
#[derive(Debug, Default, Clone, PartialEq)]
pub struct Atom {
    pub label: String,
    #[cfg_attr(feature = "serde", serde(rename = "type"))]
    pub type_: String,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub occupancy: f64,
    pub multiplicity: f64,
    pub adp_type: String,
    pub u_iso_or_equiv: f64,
    #[cfg_attr(feature = "serde", serde(rename = "U11"))]
    pub u11: f64,
    #[cfg_attr(feature = "serde", serde(rename = "U22"))]
    pub u22: f64,
    #[cfg_attr(feature = "serde", serde(rename = "U33"))]
    pub u33: f64,
    #[cfg_attr(feature = "serde", serde(rename = "U12"))]
    pub u12: f64,
    #[cfg_attr(feature = "serde", serde(rename = "U13"))]
    pub u13: f64,
    #[cfg_attr(feature = "serde", serde(rename = "U23"))]
    pub u23: f64,
}
