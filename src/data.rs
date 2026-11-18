use std::{collections::HashMap, sync::LazyLock};

pub const VAN_DER_WAALS_RADII_RAW: &[u8] =
    include_bytes!("../assets/data/van-der-Waals-radii.toml");

pub struct VanDerWaalsRadii(LazyLock<HashMap<String, f32>>);

impl VanDerWaalsRadii {
    pub const fn load() -> Self {
        Self(LazyLock::new(|| {
            let data = std::str::from_utf8(VAN_DER_WAALS_RADII_RAW)
                .expect("Failed to convert van_der_Waals data to str");
            toml::from_str(data).expect("Failed to parse van_der_Waals data form toml")
        }))
    }

    pub fn get(&self, atom_type: &str) -> &f32 {
        self.0.get(atom_type).unwrap_or_else(|| {
            log::warn!(
                "Unknown atom type: {}. Using default radius of 1.0.",
                atom_type
            );
            &1.0
        })
    }
}

pub static VAN_DER_WAALS_RADII: VanDerWaalsRadii = VanDerWaalsRadii::load();

#[cfg(test)]
#[test]
fn test_van_der_waals_radii() {
    assert_eq!(VAN_DER_WAALS_RADII.get("H"), &1.2);
}

pub const COVALENT_RADII_RAW: &[u8] = include_bytes!("../assets/data/covalent-radii.toml");

pub struct CovalentRadii(LazyLock<HashMap<String, f32>>);

impl CovalentRadii {
    pub const fn load() -> Self {
        Self(LazyLock::new(|| {
            let data = std::str::from_utf8(COVALENT_RADII_RAW)
                .expect("Failed to convert covalent radii data to str");
            toml::from_str(data).expect("Failed to parse covalent radii data form toml")
        }))
    }

    pub fn get(&self, atom_type: &str) -> &f32 {
        self.0.get(atom_type).unwrap_or_else(|| {
            log::warn!(
                "Unknown atom type: {}. Using default radius of 1.0.",
                atom_type
            );
            &1.0
        })
    }
}

pub static COVALENT_RADII: CovalentRadii = CovalentRadii::load();

#[cfg(test)]
#[test]
fn test_covalent_radii() {
    assert_eq!(COVALENT_RADII.get("H"), &0.31);
}

pub const METALLIC_RADII_RAW: &[u8] = include_bytes!("../assets/data/metallic-radii.toml");

/// Metallic radii by element symbol (Å, mostly for 12-coordination).
pub struct MetallicRadii(LazyLock<HashMap<String, f32>>);

impl MetallicRadii {
    pub const fn load() -> Self {
        Self(LazyLock::new(|| {
            let data = std::str::from_utf8(METALLIC_RADII_RAW)
                .expect("Failed to convert metallic radii data to str");
            toml::from_str(data).expect("Failed to parse metallic radii data from toml")
        }))
    }

    /// The metallic radius of the element `symbol`, or `None` for a non-metal, a metal the
    /// source gives none (Sn, Pb, Sb, Bi) or an unknown symbol.
    pub fn get(&self, symbol: &str) -> Option<f32> {
        self.0.get(symbol).copied()
    }
}

pub static METALLIC_RADII: MetallicRadii = MetallicRadii::load();

#[cfg(test)]
#[test]
fn test_metallic_radii() {
    assert_eq!(METALLIC_RADII.get("Fe"), Some(1.26));
    assert_eq!(METALLIC_RADII.get("Al"), Some(1.43));
    assert_eq!(METALLIC_RADII.get("O"), None);
    assert_eq!(METALLIC_RADII.get("Xx"), None);
}

pub const ELECTRONEGATIVITIES_RAW: &[u8] = include_bytes!("../assets/data/electronegativity.toml");

/// Pauling electronegativities by element symbol.
pub struct Electronegativities(LazyLock<HashMap<String, f32>>);

impl Electronegativities {
    pub const fn load() -> Self {
        Self(LazyLock::new(|| {
            let data = std::str::from_utf8(ELECTRONEGATIVITIES_RAW)
                .expect("Failed to convert electronegativity data to str");
            toml::from_str(data).expect("Failed to parse electronegativity data from toml")
        }))
    }

    /// The Pauling electronegativity of the element `symbol`, or `None` for one the scale gives
    /// none (He, Ne, Ar) or an unknown symbol.
    pub fn get(&self, symbol: &str) -> Option<f32> {
        self.0.get(symbol).copied()
    }
}

pub static ELECTRONEGATIVITIES: Electronegativities = Electronegativities::load();

#[cfg(test)]
#[test]
fn test_electronegativities() {
    assert_eq!(ELECTRONEGATIVITIES.get("F"), Some(3.98));
    assert_eq!(ELECTRONEGATIVITIES.get("Na"), Some(0.93));
    assert_eq!(ELECTRONEGATIVITIES.get("He"), None);
    assert_eq!(ELECTRONEGATIVITIES.get("Xx"), None);
}

pub const OXIDATION_STATES_RAW: &[u8] = include_bytes!("../assets/data/oxidation-states.toml");

/// The known oxidation states of each element.
pub struct OxidationStates(LazyLock<HashMap<String, Vec<i8>>>);

impl OxidationStates {
    pub const fn load() -> Self {
        Self(LazyLock::new(|| {
            let data = std::str::from_utf8(OXIDATION_STATES_RAW)
                .expect("Failed to convert oxidation state data to str");
            toml::from_str(data).expect("Failed to parse oxidation state data from toml")
        }))
    }

    /// The known oxidation states of the element `symbol`, in increasing order; none for an
    /// unknown symbol.
    pub fn get(&self, symbol: &str) -> &[i8] {
        self.0.get(symbol).map_or(&[], Vec::as_slice)
    }
}

pub static OXIDATION_STATES: OxidationStates = OxidationStates::load();

pub const IONIC_RADII_RAW: &[u8] = include_bytes!("../assets/data/shannon-radii.csv");

/// Shannon's (1976) effective ionic radii (Å) by element, oxidation state and coordination
/// (`assets/data/shannon-radii.csv`, see `shannon-radii-NOTICE.md`).
pub struct IonicRadii(LazyLock<HashMap<(String, i8), Vec<IonicRadius>>>);

/// One of Shannon's radii of an element at an oxidation state.
struct IonicRadius {
    /// The number of neighbours.
    coordination: u8,
    /// Whether it is for one arrangement of them only (square planar, pyramidal).
    arranged: bool,
    /// Whether it is for the low spin state where the table gives both.
    low_spin: bool,
    /// The effective ionic radius (Å).
    radius: f32,
}

impl IonicRadii {
    pub const fn load() -> Self {
        Self(LazyLock::new(|| {
            let data = std::str::from_utf8(IONIC_RADII_RAW)
                .expect("Failed to convert ionic radii data to str");
            let mut radii: HashMap<(String, i8), Vec<IonicRadius>> = HashMap::new();
            for line in data.lines().skip(1) {
                let fields: Vec<&str> = line.split(',').collect();
                let [element, charge, coordination, spin, _crystal, ionic] = fields[..] else {
                    panic!("Failed to read the ionic radii row {line:?}");
                };
                let numeral = coordination
                    .strip_suffix("SQ")
                    .or_else(|| coordination.strip_suffix("PY"))
                    .unwrap_or(coordination);
                let radius = IonicRadius {
                    coordination: from_roman(numeral)
                        .unwrap_or_else(|| panic!("Unknown coordination {coordination:?}")),
                    arranged: numeral.len() != coordination.len(),
                    low_spin: spin == "Low Spin",
                    radius: ionic.parse().expect("Failed to read an ionic radius"),
                };
                let charge = charge.parse().expect("Failed to read an ionic charge");
                radii
                    .entry((element.to_string(), charge))
                    .or_default()
                    .push(radius);
            }
            radii
        }))
    }

    /// The effective ionic radius (Å) of the element `symbol` at the oxidation state `charge`
    /// with `coordination` neighbours: the table's for that many, else for the number nearest
    /// it (the larger of two as near). Where the table gives several for one number, the one
    /// for any arrangement comes before those for square planar or pyramidal, high spin before
    /// low. `None` for a state the table gives none for.
    pub fn get(&self, symbol: &str, charge: i8, coordination: usize) -> Option<f32> {
        let distance =
            |radius: &IonicRadius| usize::from(radius.coordination).abs_diff(coordination);
        self.0
            .get(&(symbol.to_string(), charge))?
            .iter()
            .min_by_key(|radius| {
                (
                    distance(radius),
                    std::cmp::Reverse(radius.coordination),
                    radius.arranged,
                    radius.low_spin,
                )
            })
            .map(|radius| radius.radius)
    }
}

/// The number a Roman numeral from I to XIV writes.
fn from_roman(numeral: &str) -> Option<u8> {
    const NUMERALS: [&str; 14] = [
        "I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV",
    ];
    let index = NUMERALS.iter().position(|known| *known == numeral)?;
    u8::try_from(index + 1).ok()
}

pub static IONIC_RADII: IonicRadii = IonicRadii::load();

#[cfg(test)]
#[test]
fn test_ionic_radii() {
    assert_eq!(IONIC_RADII.get("Na", 1, 6), Some(1.02));
    assert_eq!(IONIC_RADII.get("O", -2, 6), Some(1.4));
    // Fe3+ at 6 is high spin 0.645, low spin 0.55.
    assert_eq!(IONIC_RADII.get("Fe", 3, 6), Some(0.645));
    // Fe2+ at 4 is 0.63, square planar 0.64.
    assert_eq!(IONIC_RADII.get("Fe", 2, 4), Some(0.63));
    // Ag2+ has only square planar at 4.
    assert_eq!(IONIC_RADII.get("Ag", 2, 4), Some(0.79));
    // Al3+ has 4, 5 and 6; 7 takes 6, 3 takes 4.
    assert_eq!(IONIC_RADII.get("Al", 3, 7), IONIC_RADII.get("Al", 3, 6));
    assert_eq!(IONIC_RADII.get("Al", 3, 3), IONIC_RADII.get("Al", 3, 4));
    assert_eq!(IONIC_RADII.get("Fe", 7, 6), None);
    assert_eq!(IONIC_RADII.get("Xx", 1, 6), None);
}

#[cfg(test)]
#[test]
fn every_ionic_radius_row_is_read() {
    let rows = std::str::from_utf8(IONIC_RADII_RAW)
        .unwrap()
        .lines()
        .count()
        - 1;
    assert_eq!(IONIC_RADII.0.values().map(Vec::len).sum::<usize>(), rows);
}

pub const BOND_VALENCE_1991_RAW: &[u8] = include_bytes!("../assets/data/bond-valence-1991.toml");

#[cfg(feature = "bvparm2020")]
pub const BVPARM_2020_RAW: &[u8] = include_bytes!("../assets/data/bvparm2020.cif");

/// The parameters of the bond valence `s = exp((R0 - R) / b)` of a bond of length `R`.
///
/// With the `bvparm2020` feature, a pair Brown's table lists (for the two oxidation states, or
/// for any state where the table gives none) takes its first entry that isn't marked unchecked.
/// Any other pair takes R0 from the element parameters of O'Keeffe & Brese (1991), with
/// `b = 0.37`, the same for every oxidation state.
pub struct BondValence(LazyLock<BondValenceData>);

struct BondValenceData {
    elements: HashMap<String, (f32, f32)>,
    #[cfg(feature = "bvparm2020")]
    table: HashMap<(String, i8, String, i8), (f32, f32)>,
}

impl BondValence {
    /// The `b` of every pair without its own (Å).
    pub const B: f32 = 0.37;

    /// Brown's marker for an oxidation state the source of an entry didn't specify.
    #[cfg(feature = "bvparm2020")]
    const UNSPECIFIED: i8 = 9;

    pub const fn load() -> Self {
        Self(LazyLock::new(|| BondValenceData {
            elements: Self::read_elements(),
            #[cfg(feature = "bvparm2020")]
            table: Self::read_table(),
        }))
    }

    fn read_elements() -> HashMap<String, (f32, f32)> {
        let data = std::str::from_utf8(BOND_VALENCE_1991_RAW)
            .expect("Failed to convert bond valence data to str");
        let elements: HashMap<String, HashMap<String, f32>> =
            toml::from_str(data).expect("Failed to parse bond valence data from toml");

        elements
            .into_iter()
            .filter_map(|(element, values)| {
                let parameters = (*values.get("r")?, *values.get("c")?);
                Some((element, parameters))
            })
            .collect()
    }

    /// The table's rows: `cation state anion state R0 b reference [details]`, the most reliable
    /// entry of a pair first. Any other line (comments, the CIF header, the reference list) fails
    /// to read as one and is skipped.
    #[cfg(feature = "bvparm2020")]
    fn read_table() -> HashMap<(String, i8, String, i8), (f32, f32)> {
        let data =
            std::str::from_utf8(BVPARM_2020_RAW).expect("Failed to convert bvparm2020 to str");
        let mut table = HashMap::new();

        for line in data.split(['\r', '\n']) {
            if line.starts_with('#') {
                continue;
            }

            let fields: Vec<&str> = line.split_whitespace().collect();
            let [cation, state, anion, anion_state, r0, b, _reference, details @ ..] = &fields[..]
            else {
                continue;
            };
            if details.first() == Some(&"unchecked") {
                continue;
            }

            let (Ok(state), Ok(anion_state), Ok(r0), Ok(b)) =
                (state.parse(), anion_state.parse(), r0.parse(), b.parse())
            else {
                continue;
            };

            table
                .entry((cation.to_string(), state, anion.to_string(), anion_state))
                .or_insert((r0, b));
        }

        table
    }

    /// The bond valence (v.u.) of a bond `distance` (Å) long with the parameters `(R0, b)`:
    /// `exp((R0 - distance) / b)`.
    pub fn valence((r0, b): (f32, f32), distance: f32) -> f32 {
        ((r0 - distance) / b).exp()
    }

    /// The length (Å) at which a bond with the parameters `(R0, b)` has the bond valence
    /// `valence` (v.u.): `R0 - b ln(valence)`, the inverse of [`Self::valence`].
    pub fn distance((r0, b): (f32, f32), valence: f32) -> f32 {
        r0 - b * valence.ln()
    }

    /// `(R0, b)` (Å) of the bond between `cation` in oxidation state `state` and `anion` in
    /// `anion_state`, or `None` if neither source covers the pair.
    pub fn parameters(
        &self,
        cation: &str,
        state: i8,
        anion: &str,
        anion_state: i8,
    ) -> Option<(f32, f32)> {
        #[cfg(feature = "bvparm2020")]
        {
            let entry = |state: i8, anion_state: i8| {
                let key = (cation.to_string(), state, anion.to_string(), anion_state);
                self.0.table.get(&key).copied()
            };
            let found = entry(state, anion_state)
                .or_else(|| entry(Self::UNSPECIFIED, anion_state))
                .or_else(|| entry(state, Self::UNSPECIFIED));
            if found.is_some() {
                return found;
            }
        }

        #[cfg(not(feature = "bvparm2020"))]
        let _ = (state, anion_state);

        self.elements_only(cation, anion)
    }

    /// `(R0, b)` of the pair from the element parameters of O'Keeffe & Brese alone.
    fn elements_only(&self, cation: &str, anion: &str) -> Option<(f32, f32)> {
        let (ri, ci) = *self.0.elements.get(cation)?;
        let (rj, cj) = *self.0.elements.get(anion)?;
        let r0 = ri + rj - ri * rj * (ci.sqrt() - cj.sqrt()).powi(2) / (ci * ri + cj * rj);

        Some((r0, Self::B))
    }

    /// The positive oxidation states `cation` can take: the ones Brown's table gives
    /// parameters for, with the `bvparm2020` feature and where it lists the element, else its
    /// known ones (see [`OXIDATION_STATES`]).
    pub fn cation_states(&self, cation: &str) -> Vec<i8> {
        #[cfg(feature = "bvparm2020")]
        {
            let mut states: Vec<i8> = self
                .0
                .table
                .keys()
                .filter(|(element, state, _, _)| element == cation && *state != Self::UNSPECIFIED)
                .map(|(_, state, _, _)| *state)
                .filter(|state| *state > 0)
                .collect();
            states.sort_unstable();
            states.dedup();
            if !states.is_empty() {
                return states;
            }
        }

        OXIDATION_STATES
            .get(cation)
            .iter()
            .copied()
            .filter(|state| *state > 0)
            .collect()
    }
}

pub static BOND_VALENCE: BondValence = BondValence::load();

#[cfg(test)]
mod bond_valence_tests {
    use super::*;

    #[test]
    fn element_parameters_give_r0_by_o_keeffe_and_brese() {
        let (r0, b) = BondValence::load().elements_only("Ti", "O").unwrap();
        assert!((r0 - 1.814).abs() < 1e-3, "{r0}");
        assert_eq!(b, BondValence::B);
        assert_eq!(BondValence::load().elements_only("Xx", "O"), None);
    }

    #[test]
    fn oxidation_states_are_known_per_element() {
        assert_eq!(OXIDATION_STATES.get("Ba"), &[2]);
        assert!(OXIDATION_STATES.get("Fe").contains(&3));
        assert!(OXIDATION_STATES.get("Xx").is_empty());
    }

    #[test]
    fn a_bond_valence_and_its_distance_are_inverse() {
        let parameters = (1.759, 0.37);

        assert!((BondValence::valence(parameters, 1.759) - 1.0).abs() < 1e-6);
        let distance = BondValence::distance(parameters, 0.5);
        assert!((BondValence::valence(parameters, distance) - 0.5).abs() < 1e-6);
    }

    #[cfg(feature = "bvparm2020")]
    #[test]
    fn the_2020_table_gives_a_pair_its_first_reliable_entry() {
        assert_eq!(
            BOND_VALENCE.parameters("Fe", 3, "O", -2),
            Some((1.759, 0.37))
        );
        assert_eq!(
            BOND_VALENCE.parameters("Al", 3, "S", -2),
            Some((2.13, 0.37))
        );
        assert_eq!(
            BOND_VALENCE.parameters("Ag", 1, "Br", -1),
            Some((2.22, 0.37))
        );
        assert_eq!(BOND_VALENCE.cation_states("Fe"), vec![2, 3]);
    }

    #[cfg(feature = "bvparm2020")]
    #[test]
    fn the_2020_table_reads_only_its_parameter_rows() {
        let table = super::BondValence::read_table();

        assert_eq!(table.len(), 1168);
        assert!(
            table.keys().all(
                |(cation, _, anion, _)| [cation, anion].iter().all(|symbol| {
                    symbol.len() <= 2 && symbol.starts_with(|c: char| c.is_ascii_uppercase())
                })
            ),
            "a line that isn't a parameter row was read as one"
        );
    }

    #[cfg(not(feature = "bvparm2020"))]
    #[test]
    fn without_the_2020_table_every_state_takes_the_element_parameters() {
        assert_eq!(
            BOND_VALENCE.parameters("Fe", 2, "O", -2),
            BOND_VALENCE.parameters("Fe", 3, "O", -2)
        );
        assert_eq!(BOND_VALENCE.cation_states("Ba"), vec![2]);
    }
}
