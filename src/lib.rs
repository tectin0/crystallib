//! Crystal structures (a phase's cell and atoms, values with their standard uncertainties),
//! space group symmetry (feature `symmetry`) and element data (feature `data`).
#![cfg_attr(all(feature = "symmetry", feature = "data"), doc = "")]
#![cfg_attr(all(feature = "symmetry", feature = "data"), doc = include_str!("../README.md"))]

/// The vector and matrix types of the API, so a user needs no `cgmath` of its own.
pub use cgmath;

mod atoms;
#[cfg(all(feature = "data", feature = "symmetry"))]
pub mod bond_character;
#[cfg(all(feature = "data", feature = "symmetry"))]
pub mod bond_pairs;
#[cfg(all(feature = "data", feature = "symmetry"))]
pub mod bonding;
mod cell;
#[cfg(feature = "symmetry")]
pub mod cell_atoms;
#[cfg(all(feature = "data", feature = "symmetry"))]
pub mod contacts;
#[cfg(all(feature = "data", feature = "symmetry"))]
pub mod disorder;
#[cfg(feature = "distance-distributions")]
pub mod distances;
#[cfg(all(feature = "data", feature = "symmetry"))]
pub mod elements;
pub mod geometry;
mod measured;
pub mod occupancy;
#[cfg(all(feature = "data", feature = "symmetry"))]
pub mod oxidation;
mod phase;
pub mod shapes;
pub mod type_symbol;

/// The CIF reader [`cif_phase`] reads with, so a user needs no `cif` of its own.
#[cfg(feature = "cif")]
pub use cif;
#[cfg(feature = "cif")]
pub mod cif_phase;

#[cfg(feature = "data")]
mod data;

#[cfg(feature = "data")]
pub use data::{
    ALLEN_ELECTRONEGATIVITIES, BOND_VALENCE, BondValence, COVALENT_RADII, CovalentRadii,
    ELECTRONEGATIVITIES, Electronegativities, IONIC_RADII, IonicRadii, METALLIC_RADII,
    MetallicRadii, OXIDATION_STATES, OxidationStates, VAN_DER_WAALS_RADII, VanDerWaalsRadii,
};

#[cfg(feature = "symmetry")]
mod symmetry;

#[cfg(feature = "symmetry")]
pub use symmetry::{
    IntoSpaceGroupNumber, IntoSpaceGroupSymbol, LARGEST_SITE_MERGE, SAME_POSITION,
    SPACEGROUP_NUMBERS, SPACEGROUP_SYMBOLS, SpaceGroupGenerators, SpaceGroupSymmetryOperations,
    independent_lattice_components_for_laue, independent_lattice_components_for_space_group_number,
    lattice_component_index, laue_and_unique_axis_from_sg_number, parse_symmetry_operation,
    with_face_copies,
};

pub use atoms::Atom;
pub use atoms::Atoms;

pub use cell::Cell;
pub use measured::Measured;
pub use phase::Phase;

pub use atoms::{AdpType, ParseAdpTypeError};
