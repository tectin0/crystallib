//! Crystal structures (a phase's cell and atoms, values with their standard uncertainties),
//! space group symmetry (feature `symmetry`) and element data (feature `data`).
#![cfg_attr(all(feature = "symmetry", feature = "data"), doc = "")]
#![cfg_attr(all(feature = "symmetry", feature = "data"), doc = include_str!("../README.md"))]

/// The vector and matrix types of the API, so a user needs no `cgmath` of its own.
pub use cgmath;

mod atoms;
mod cell;
mod measured;
mod phase;

#[cfg(feature = "data")]
mod data;

#[cfg(feature = "data")]
pub use data::{
    BondValence, CovalentRadii, Electronegativities, IonicRadii, MetallicRadii, OxidationStates,
    VanDerWaalsRadii, ALLEN_ELECTRONEGATIVITIES, BOND_VALENCE, COVALENT_RADII, ELECTRONEGATIVITIES,
    IONIC_RADII, METALLIC_RADII, OXIDATION_STATES, VAN_DER_WAALS_RADII,
};

#[cfg(feature = "symmetry")]
mod symmetry;

#[cfg(feature = "symmetry")]
pub use symmetry::{
    independent_lattice_components_for_laue, independent_lattice_components_for_space_group_number,
    lattice_component_index, laue_and_unique_axis_from_sg_number, with_face_copies,
    IntoSpaceGroupNumber, IntoSpaceGroupSymbol, SpaceGroupGenerators, SpaceGroupSymmetryOperations,
    LARGEST_SITE_MERGE, SAME_POSITION, SPACEGROUP_NUMBERS, SPACEGROUP_SYMBOLS,
};

pub use atoms::Atom;
pub use atoms::Atoms;

pub use cell::Cell;
pub use measured::Measured;
pub use phase::Phase;

pub use atoms::{AdpType, ParseAdpTypeError};
