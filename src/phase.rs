use crate::{atoms::Atoms, cell::Cell};

#[cfg_attr(feature = "serde", derive(serde::Deserialize, serde::Serialize))]
#[cfg_attr(feature = "serde", serde(default))]
#[derive(Default, Debug, Clone, PartialEq)]
pub struct Phase {
    pub cell: Cell,
    pub atoms: Atoms,
}
