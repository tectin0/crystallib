use crate::{atoms::Atoms, cell::Cell};

#[cfg_attr(feature = "serde", derive(serde::Deserialize, serde::Serialize))]
#[cfg_attr(feature = "serde", serde(default))]
#[derive(Default, Debug, Clone, PartialEq)]
/// A crystalline phase: its unit cell and the atoms on its sites.
pub struct Phase {
    /// The unit cell and space group.
    pub cell: Cell,
    /// The atoms, one per site of the asymmetric unit.
    pub atoms: Atoms,
}
