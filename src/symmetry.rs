use std::{
    collections::{BTreeMap, HashMap},
    sync::LazyLock,
};

use cgmath::Transform;

use crate::Cell;

/// How close (in every fractional coordinate) two positions are to count as one. CIFs round
/// special positions (1/3 as 0.3333), so symmetry operations turn one position into copies that
/// differ in the fourth decimal; for a cell of 5–20 Å this is 0.005–0.02 Å.
pub const SAME_POSITION: f64 = 1e-3;

/// How far apart (Å) copies of one site may be to be merged down to the multiplicity the site
/// states (see [`SpaceGroupSymmetryOperations::site_positions_in_cell`]): closer than any two
/// atoms are bonded (the shortest bond, H–H, is 0.74 Å). CIFs give special positions a little
/// off (y = 0.249 for 1/4, 0.33 for 1/3, x = -0.0065 for 0), which puts their copies up to a
/// few tenths of an Ångström apart.
pub const LARGEST_SITE_MERGE: f64 = 0.5;

/// The distinct positions `operations` take `point` (fractional) to, in the unit cell
/// `[0, 1)`, with a copy on the opposite face of the cell for each coordinate on a face (so an
/// atom at a corner shows at all eight). Positions within [`SAME_POSITION`] of each other in
/// every coordinate (modulo 1) are one; a coordinate within it of a face is on it.
pub fn equivalent_points(
    operations: &[cgmath::Matrix4<f64>],
    point: cgmath::Point3<f64>,
) -> Vec<cgmath::Point3<f64>> {
    with_face_copies(distinct_points(operations, point, SAME_POSITION))
}

/// The positions `operations` take `point` (fractional) to, in the unit cell `[0, 1)`, those
/// within `tolerance` of each other in every coordinate (modulo 1) as one.
fn distinct_points(
    operations: &[cgmath::Matrix4<f64>],
    point: cgmath::Point3<f64>,
    tolerance: f64,
) -> Vec<cgmath::Point3<f64>> {
    let same = |a: &cgmath::Point3<f64>, b: &cgmath::Point3<f64>| {
        let delta = a - b;
        [delta.x, delta.y, delta.z]
            .iter()
            .all(|d| (d - d.round()).abs() < tolerance)
    };

    let mut distinct: Vec<cgmath::Point3<f64>> = Vec::new();
    for operation in operations {
        let new_point = operation.transform_point(point).map(|c| c.rem_euclid(1.0));
        if !distinct.iter().any(|known| same(known, &new_point)) {
            distinct.push(new_point);
        }
    }
    distinct
}

/// `points` (fractional, in `[0, 1)`) for drawing: a coordinate within [`SAME_POSITION`] of a
/// face of the cell is on it, and gets a copy on the opposite face, so an atom at a corner shows
/// at all eight.
pub fn with_face_copies(points: Vec<cgmath::Point3<f64>>) -> Vec<cgmath::Point3<f64>> {
    let on_face = |c: f64| {
        if c < SAME_POSITION || c > 1.0 - SAME_POSITION {
            0.0
        } else {
            c
        }
    };

    let mut distinct: Vec<cgmath::Point3<f64>> = Vec::with_capacity(points.len());
    for point in points.into_iter().map(|point| point.map(on_face)) {
        if !distinct.contains(&point) {
            distinct.push(point);
        }
    }

    let mut points = Vec::with_capacity(distinct.len());
    for point in distinct {
        let shifts = |c: f64| if c == 0.0 { vec![0.0, 1.0] } else { vec![0.0] };
        for dx in shifts(point.x) {
            for dy in shifts(point.y) {
                for dz in shifts(point.z) {
                    points.push(cgmath::Point3::new(
                        point.x + dx,
                        point.y + dy,
                        point.z + dz,
                    ));
                }
            }
        }
    }
    points
}

/// The symmetry operation `text` writes as CIFs do, each coordinate a sum of `x`, `y`, `z`
/// with signs and of a constant (`-x+y, -x, z+1/3`, `1/2+x, y, -z`, `x, y, z+0.5`), as the
/// matrix that transforms fractional coordinates. `None` if it doesn't read as three such sums.
pub fn parse_symmetry_operation(text: &str) -> Option<cgmath::Matrix4<f64>> {
    let text = text.trim().trim_matches(|c| c == '\'' || c == '"');
    let rows: Vec<&str> = text.split(',').collect();
    let [x, y, z] = rows[..] else {
        return None;
    };

    let mut matrix = cgmath::Matrix4::from_scale(0.0);
    matrix.w.w = 1.0;
    for (row, expression) in [x, y, z].into_iter().enumerate() {
        let (coefficients, constant) = parse_symmetry_component(expression)?;
        matrix.x[row] = coefficients[0];
        matrix.y[row] = coefficients[1];
        matrix.z[row] = coefficients[2];
        matrix.w[row] = constant;
    }
    Some(matrix)
}

/// The coefficients of `x`, `y` and `z` and the constant of one coordinate of a symmetry
/// operation (`-x+y`, `z+1/3`, `1/2-y`, `2x`).
fn parse_symmetry_component(expression: &str) -> Option<([f64; 3], f64)> {
    let expression: String = expression.chars().filter(|c| !c.is_whitespace()).collect();
    if expression.is_empty() {
        return None;
    }

    let mut coefficients = [0.0; 3];
    let mut constant = 0.0;

    let mut terms = Vec::new();
    let mut start = 0;
    for (index, c) in expression.char_indices() {
        if (c == '+' || c == '-') && index > start {
            terms.push(&expression[start..index]);
            start = index;
        }
    }
    terms.push(&expression[start..]);

    for term in terms {
        let (sign, body) = match term.as_bytes().first()? {
            b'-' => (-1.0, &term[1..]),
            b'+' => (1.0, &term[1..]),
            _ => (1.0, term),
        };
        let number = |text: &str| -> Option<f64> {
            if text.is_empty() {
                return Some(1.0);
            }
            match text.split_once('/') {
                Some((numerator, denominator)) => {
                    Some(numerator.parse::<f64>().ok()? / denominator.parse::<f64>().ok()?)
                }
                None => text.parse::<f64>().ok(),
            }
        };

        match body.chars().last()?.to_ascii_lowercase() {
            axis @ ('x' | 'y' | 'z') => {
                let factor = number(body[..body.len() - 1].trim_end_matches('*'))?;
                coefficients[(axis as u8 - b'x') as usize] += sign * factor;
            }
            _ => constant += sign * number(body)?,
        }
    }

    Some((coefficients, constant))
}

const SPACEGROUP_SETTINGS_RAW: &str = include_str!("../assets/spacegroup_settings.tsv");

/// One of the 530 settings of the space groups (one per Hall symbol): the axes and origin a
/// structure's coordinates are given in.
#[derive(Debug, Clone, PartialEq)]
pub struct SpaceGroupSetting {
    /// The space group's number (1–230).
    pub number: u8,
    /// Which of its group's settings this is, as spglib names it: the origin choice (`1`, `2`),
    /// hexagonal or rhombohedral axes (`H`, `R`), a monoclinic group's unique axis and cell
    /// choice (`b1`, `-c2`), an orthorhombic group's axes (`cab`, `2ba-c`); empty for a group's
    /// only setting.
    pub choice: &'static str,
    /// The short Hermann–Mauguin symbol, without spaces (`P2_1/c`, `Pbnm`).
    pub short_symbol: &'static str,
    /// The full Hermann–Mauguin symbol (`P 1 2_1/c 1`).
    pub full_symbol: &'static str,
    /// The Hall symbol (`-P 2ybc`).
    pub hall_symbol: &'static str,
    /// The symmetry operations, centring included, as matrices that transform fractional
    /// coordinates.
    pub operations: Vec<cgmath::Matrix4<f64>>,
}

/// The settings of the space groups, from spglib's database (BSD-3-Clause, see
/// `assets/spacegroup_settings-NOTICE.md`), in its order: by number, a group's origin choice 1
/// before 2 and hexagonal axes before rhombohedral (see [`SpaceGroupSetting::standard`]).
pub static SPACEGROUP_SETTINGS: LazyLock<Vec<SpaceGroupSetting>> = LazyLock::new(|| {
    SPACEGROUP_SETTINGS_RAW
        .lines()
        .map(|line| {
            let fields: Vec<&'static str> = line.split('\t').collect();
            let [
                number,
                choice,
                short_symbol,
                full_symbol,
                hall_symbol,
                operations,
            ] = fields[..]
            else {
                panic!("a space group setting has six fields: {line}");
            };

            SpaceGroupSetting {
                number: number.parse().expect("a space group setting's number"),
                choice,
                short_symbol,
                full_symbol,
                hall_symbol,
                operations: operations
                    .split(';')
                    .map(|operation| {
                        parse_symmetry_operation(operation)
                            .unwrap_or_else(|| panic!("`{operation}` of {hall_symbol} reads"))
                    })
                    .collect(),
            }
        })
        .collect()
});

/// The settings each name of a setting names. A setting's names are its short and full symbol
/// (see [`compact_symbol`]); for a monoclinic group (numbers 3–15), its full symbol without its
/// `1`s, the short symbol of each of its settings (spglib's short symbol is the standard
/// setting's for all of them: `P2_1/c` for `P 1 2_1/n 1`); and for a group with an e-glide
/// (numbers 39, 41, 64, 67 and 68), the older symbols that write each `e` as either of the two
/// glides it combines (`Cmca` and `Cmcb` for `Cmce`). No name names settings of two groups (a
/// test checks it).
static SETTINGS_BY_NAME: LazyLock<HashMap<String, Vec<&'static SpaceGroupSetting>>> =
    LazyLock::new(|| {
        let mut by_name: HashMap<String, Vec<&'static SpaceGroupSetting>> = HashMap::new();

        for setting in SPACEGROUP_SETTINGS.iter() {
            let mut names = vec![
                compact_symbol(setting.short_symbol),
                compact_symbol(setting.full_symbol),
            ];
            if (3..=15).contains(&setting.number) {
                let without_ones = setting
                    .full_symbol
                    .split_whitespace()
                    .filter(|&part| part != "1");
                names.push(compact_symbol(&without_ones.collect::<Vec<_>>().join(" ")));
            }
            names.extend(
                older_e_glide_symbols(setting.full_symbol)
                    .iter()
                    .map(|symbol| compact_symbol(symbol)),
            );
            names.sort();
            names.dedup();

            for name in names {
                by_name.entry(name).or_default().push(setting);
            }
        }

        by_name
    });

/// The full and short symbols an orthorhombic group's setting with an e-glide had before the
/// e-glide symbol was introduced: each `e` (in the part for the axis its plane is perpendicular
/// to) as either glide along the other two axes; the short symbol is each axis's part after its
/// `/`. None for a symbol without an `e`.
fn older_e_glide_symbols(full_symbol: &str) -> Vec<String> {
    let parts: Vec<&str> = full_symbol.split_whitespace().collect();
    let [lattice, axes @ ..] = &parts[..] else {
        return Vec::new();
    };
    if axes.len() != 3 || !axes.iter().any(|part| part.contains('e')) {
        return Vec::new();
    }

    let mut older: Vec<Vec<String>> = vec![Vec::new()];
    for (axis, part) in axes.iter().enumerate() {
        let ways: Vec<String> = match part.contains('e') {
            true => ["a", "b", "c"]
                .iter()
                .enumerate()
                .filter(|(other, _)| *other != axis)
                .map(|(_, glide)| part.replace('e', glide))
                .collect(),
            false => vec![part.to_string()],
        };
        older = older
            .iter()
            .flat_map(|written| {
                ways.iter()
                    .map(|way| [written.clone(), vec![way.clone()]].concat())
            })
            .collect();
    }

    older
        .iter()
        .flat_map(|axes| {
            let short: String = axes
                .iter()
                .map(|part| part.rsplit('/').next().unwrap_or(part))
                .collect();
            [
                format!("{lattice} {}", axes.join(" ")),
                format!("{lattice}{short}"),
            ]
        })
        .collect()
}

impl SpaceGroupSetting {
    /// The setting the Hermann–Mauguin symbol `symbol` names, or `None` for a symbol no setting
    /// has. Spaces are ignored, screw axes may be written as CIFs write them (`P 1 21/n 1`,
    /// `P63/mmc`), the glides of an e-glide group as they were written before (`C m c a`), and
    /// a cubic group's `-3` as it was written before (`F m 3 m`). A suffix chooses among the
    /// settings a symbol names (`:1`, `:2`, `:H`, `:R`, `:b1`; for an orthorhombic group with
    /// two origins and other axes, `:2` is origin choice 2 in those axes), and a symbol whose
    /// suffix names none of them names none. Without a suffix, a symbol names the first of them,
    /// except origin choice 2 of a group with two origins.
    pub fn named(symbol: &str) -> Option<&'static SpaceGroupSetting> {
        let (name, suffix) = match symbol.split_once(':') {
            Some((name, suffix)) => (name, suffix.trim()),
            None => (symbol, ""),
        };
        let settings = settings_named(name);

        if suffix.is_empty() {
            return preferred(settings);
        }

        settings.iter().copied().find(|setting| {
            setting.choice.eq_ignore_ascii_case(suffix)
                || matches!(suffix, "1" | "2") && setting.choice.starts_with(suffix)
        })
    }

    /// The standard setting of the space group numbered `number`, or `None` for a number out
    /// of 1–230: its first, except origin choice 2 of a group with two origins (as the Bilbao
    /// Crystallographic Server's standard), the setting its symbol names without a suffix (see
    /// [`Self::named`]).
    pub fn standard(number: u8) -> Option<&'static SpaceGroupSetting> {
        let settings: Vec<_> = SPACEGROUP_SETTINGS
            .iter()
            .filter(|setting| setting.number == number)
            .collect();
        preferred(&settings)
    }

    /// The setting `cell`'s coordinates are given in when it lists no symmetry operations: the
    /// one its symbol names (see [`Self::named`]), on rhombohedral axes where that names
    /// hexagonal ones but the cell's angles aren't those of hexagonal axes (α = β = 90°,
    /// γ = 120°). `None` where its symbol names no setting, or one of another group than its
    /// number (with a number of 0, any group).
    pub fn of_cell(cell: &Cell) -> Option<&'static SpaceGroupSetting> {
        let setting = Self::named(&cell.space_group)?;

        if cell.space_group_number != 0 && setting.number != cell.space_group_number {
            log::debug!(
                "Space group symbol {} names group {}, its number is {}",
                cell.space_group,
                setting.number,
                cell.space_group_number
            );
            return None;
        }

        const ANGLE_TOLERANCE: f64 = 0.1;
        let hexagonal_axes = (cell.alpha.value - 90.0).abs() < ANGLE_TOLERANCE
            && (cell.beta.value - 90.0).abs() < ANGLE_TOLERANCE
            && (cell.gamma.value - 120.0).abs() < ANGLE_TOLERANCE;
        if setting.choice == "H" && !hexagonal_axes {
            return SPACEGROUP_SETTINGS
                .iter()
                .find(|other| other.number == setting.number && other.choice == "R");
        }

        Some(setting)
    }
}

/// Of `settings` (one group's, in [`SPACEGROUP_SETTINGS`]' order), the one a symbol without a
/// suffix names: the first, except origin choice 2 where they hold both origin choices. CIFs that
/// write such a symbol without a suffix are in origin choice 2 (640 of the 668 COD CIFs that also
/// list their operations, 24 of the 28 without them that give a Hall symbol).
fn preferred(settings: &[&'static SpaceGroupSetting]) -> Option<&'static SpaceGroupSetting> {
    let has_origin_choice_1 = settings.iter().any(|s| s.choice.starts_with('1'));
    let origin_choice_2 = settings.iter().find(|s| s.choice.starts_with('2'));

    match origin_choice_2 {
        Some(setting) if has_origin_choice_1 => Some(setting),
        _ => settings.first().copied(),
    }
}

/// The settings `name` (a Hermann–Mauguin symbol without a suffix) names, in
/// [`SPACEGROUP_SETTINGS`]' order; none for a name no setting has.
fn settings_named(name: &str) -> &'static [&'static SpaceGroupSetting] {
    let compact = compact_symbol(name);
    if let Some(settings) = SETTINGS_BY_NAME.get(&compact) {
        return settings;
    }

    SETTINGS_BY_NAME
        .get(&with_older_cubic_threes(&compact))
        .filter(|settings| settings.iter().all(|setting| setting.number >= 195))
        .map_or(&[], Vec::as_slice)
}

/// `symbol` as names are kept: without spaces, its screw axes with `_` (see
/// [`screw_axes_as_in_the_table`]), its lattice letter upper case and the rest lower case.
fn compact_symbol(symbol: &str) -> String {
    let compact: String = symbol
        .split_whitespace()
        .map(screw_axes_as_in_the_table)
        .collect();
    let mut chars = compact.chars();

    match chars.next() {
        Some(first) => first
            .to_uppercase()
            .chain(chars.flat_map(char::to_lowercase))
            .collect(),
        None => String::new(),
    }
}

/// `compact` (see [`compact_symbol`]) with each `3` after a mirror or glide written `-3`: older
/// cubic symbols write `Fm3m` for `Fm-3m`, `Pa3` for `Pa-3`.
fn with_older_cubic_threes(compact: &str) -> String {
    let mut written = String::with_capacity(compact.len() + 2);
    let mut previous = None;

    for c in compact.chars() {
        if c == '3' && matches!(previous, Some('m' | 'n' | 'a' | 'd')) {
            written.push('-');
        }
        written.push(c);
        previous = Some(c);
    }

    written
}

/// The symmetry operations of each space group's standard setting (see
/// [`SpaceGroupSetting::standard`]), by number.
pub static SPACEGROUP_SYMMETRY_OPERATIONS: LazyLock<BTreeMap<u8, Vec<cgmath::Matrix4<f64>>>> =
    LazyLock::new(|| {
        (1..=230)
            .filter_map(|number| {
                Some((
                    number,
                    SpaceGroupSetting::standard(number)?.operations.clone(),
                ))
            })
            .collect()
    });

const SPACEGROUP_SYMBOLS_RAW: &str = include_str!("../assets/spacegroup_symbols.json");

/// The symbol of each space group by number (1–230), as pymatgen writes it: the full
/// Hermann–Mauguin symbol without spaces, a screw axis with `_` (`P12_1/c1`).
pub static SPACEGROUP_SYMBOLS: LazyLock<BTreeMap<u8, String>> = LazyLock::new(|| {
    let symbols: Vec<String> = serde_json::from_str(SPACEGROUP_SYMBOLS_RAW).unwrap();

    symbols
        .into_iter()
        .enumerate()
        .map(|(i, s)| ((i + 1) as u8, s))
        .collect()
});

/// The number of each space group by its symbol (see [`SPACEGROUP_SYMBOLS`]).
pub static SPACEGROUP_NUMBERS: LazyLock<BTreeMap<String, u8>> = LazyLock::new(|| {
    SPACEGROUP_SYMBOLS
        .iter()
        .map(|(i, s)| (s.clone(), *i))
        .collect()
});

const SPACEGROUP_GENERATORS_RAW: &str = include_str!("../assets/spacegroup_generators.json");

pub static SPACEGROUP_GENERATORS: LazyLock<BTreeMap<u8, Vec<cgmath::Matrix4<f64>>>> =
    LazyLock::new(|| serde_json::from_str(SPACEGROUP_GENERATORS_RAW).unwrap());

/// The generators of the space groups as pymatgen gives them, as matrices that transform
/// fractional coordinates: in their standard settings (see [`SpaceGroupSetting::standard`]),
/// except that a group with two origins is in origin choice 1.
pub struct SpaceGroupGenerators {}

impl SpaceGroupGenerators {
    /// The generators of `space_group_number`, or `None` for no space group.
    pub fn get<'a>(
        space_group_number: impl IntoSpaceGroupNumber,
    ) -> Option<&'a Vec<cgmath::Matrix4<f64>>> {
        let space_group_number = space_group_number.into_space_group_number()?;

        SPACEGROUP_GENERATORS.get(&space_group_number)
    }

    /// The generators of every space group, by number.
    pub fn get_all<'a>() -> &'a BTreeMap<u8, Vec<cgmath::Matrix4<f64>>> {
        &SPACEGROUP_GENERATORS
    }

    /// `point` (fractional) transformed by each generator of `space_group_number`, moved into the
    /// unit cell and translated by 0 or 1 along each of a, b and c, each position once; `None`
    /// for no space group.
    pub fn generate_symmetry_equivalent_points_from_point(
        space_group_number: impl IntoSpaceGroupNumber,
        point: cgmath::Point3<f64>,
    ) -> Option<Vec<cgmath::Point3<f64>>> {
        let generators = Self::get(space_group_number)?;

        let new_points: BTreeMap<[(u64, i16, i8); 3], cgmath::Point3<f64>> = generators
            .iter()
            .flat_map(|m| {
                let new_point = m.transform_point(point);
                let new_point = move_point_into_unit_cell(new_point);

                TRIVIAL_TRANSLATIONS
                    .iter()
                    .map(move |t| {
                        let new_point = t.transform_point(new_point);
                        let new_point = move_point_into_unit_cell(new_point);

                        let decoded = point3_decode(new_point);

                        (decoded, new_point)
                    })
                    .collect::<Vec<_>>()
            })
            .collect();

        Some(new_points.into_values().collect())
    }

    /// The positions of [`Self::generate_symmetry_equivalent_points_from_point`] for each of
    /// `points`, each position once.
    pub fn generate_symmetry_equivalent_points_from_points(
        space_group_number: impl IntoSpaceGroupNumber + Clone,
        points: &[cgmath::Point3<f64>],
    ) -> Option<Vec<cgmath::Point3<f64>>> {
        let space_group_number = space_group_number.into_space_group_number()?;

        let mut new_points = BTreeMap::new();

        for point in points {
            let new_points_ =
                Self::generate_symmetry_equivalent_points_from_point(space_group_number, *point)?;

            new_points.extend(
                new_points_
                    .into_iter()
                    .map(|point| (point3_decode(point), point)),
            );
        }

        Some(new_points.into_values().collect())
    }
}

/// The symmetry operations of the space groups in their standard settings (see
/// [`SpaceGroupSetting::standard`]), as matrices that transform fractional coordinates.
pub struct SpaceGroupSymmetryOperations {}

impl SpaceGroupSymmetryOperations {
    /// The symmetry operations of `space_group_number`, or `None` for no space group.
    pub fn get<'a>(
        space_group_number: impl IntoSpaceGroupNumber,
    ) -> Option<&'a Vec<cgmath::Matrix4<f64>>> {
        let space_group_number = space_group_number.into_space_group_number()?;

        SPACEGROUP_SYMMETRY_OPERATIONS.get(&space_group_number)
    }

    /// The symmetry operations of every space group, by number.
    pub fn get_all<'a>() -> &'a BTreeMap<u8, Vec<cgmath::Matrix4<f64>>> {
        &SPACEGROUP_SYMMETRY_OPERATIONS
    }

    /// The symmetry operations of `cell`: the ones it states (see [`Cell::symmetry_operations`])
    /// if it states any and all of them read, else those of the setting its symbol names (see
    /// [`SpaceGroupSetting::of_cell`]), else those of its space group number in the standard
    /// setting.
    pub fn of_cell(cell: &Cell) -> Option<Vec<cgmath::Matrix4<f64>>> {
        if !cell.symmetry_operations.is_empty() {
            let stated: Option<Vec<_>> = cell
                .symmetry_operations
                .iter()
                .map(|operation| parse_symmetry_operation(operation))
                .collect();
            match stated {
                Some(operations) => return Some(operations),
                None => log::warn!(
                    "A symmetry operation of space group {} doesn't read; using the setting its \
                     symbol or number names",
                    cell.space_group
                ),
            }
        }

        match SpaceGroupSetting::of_cell(cell) {
            Some(setting) => Some(setting.operations.clone()),
            None => Self::get(cell.space_group_number).cloned(),
        }
    }

    /// The positions in the unit cell `[0, 1]` equivalent to `point` (fractional) under the
    /// symmetry operations of `space_group_number` in its standard setting, with a copy on the
    /// opposite face for each coordinate on a face (see [`with_face_copies`]).
    pub fn generate_symmetry_equivalent_points_from_point(
        space_group_number: impl IntoSpaceGroupNumber,
        point: cgmath::Point3<f64>,
    ) -> Option<Vec<cgmath::Point3<f64>>> {
        Some(equivalent_points(Self::get(space_group_number)?, point))
    }

    /// The positions in the unit cell `[0, 1]` equivalent to `point` (fractional) under the
    /// symmetry operations of `cell` (see [`Self::of_cell`]), with a copy on the opposite face for
    /// each coordinate on a face (see [`with_face_copies`]).
    pub fn generate_symmetry_equivalent_points_in_cell(
        cell: &Cell,
        point: cgmath::Point3<f64>,
    ) -> Option<Vec<cgmath::Point3<f64>>> {
        Some(equivalent_points(&Self::of_cell(cell)?, point))
    }

    /// The distinct positions in the unit cell `[0, 1)` of a site at `point` (fractional) under
    /// the symmetry operations of `cell` (see [`Self::of_cell`]), without copies on opposite
    /// faces (see [`with_face_copies`]). Without a `multiplicity`, positions within
    /// [`SAME_POSITION`] of each other in every coordinate are one. With one, only coinciding
    /// positions are, and then the closest ones are merged until the site has that many, as
    /// long as they are within [`LARGEST_SITE_MERGE`] of each other: a split site's positions
    /// can be closer than [`SAME_POSITION`] and still distinct.
    pub fn site_positions_in_cell(
        cell: &Cell,
        point: cgmath::Point3<f64>,
        multiplicity: Option<usize>,
    ) -> Option<Vec<cgmath::Point3<f64>>> {
        let operations = Self::of_cell(cell)?;
        let Some(multiplicity) = multiplicity else {
            return Some(distinct_points(&operations, point, SAME_POSITION));
        };

        const COINCIDING: f64 = 1e-9;
        let mut positions = distinct_points(&operations, point, COINCIDING);
        let transformation = cell.transformation_matrix_y_up();
        let distance = |a: &cgmath::Point3<f64>, b: &cgmath::Point3<f64>| {
            let delta = (a - b).map(|d| d - d.round());
            let cartesian = transformation * Cell::y_up_input(delta);
            f64::from(cgmath::InnerSpace::magnitude(cartesian.truncate()))
        };

        while positions.len() > multiplicity {
            let closest = (0..positions.len())
                .flat_map(|i| (i + 1..positions.len()).map(move |j| (i, j)))
                .map(|(i, j)| (j, distance(&positions[i], &positions[j])))
                .min_by(|a, b| a.1.total_cmp(&b.1));
            match closest {
                Some((j, distance)) if distance < LARGEST_SITE_MERGE => {
                    positions.remove(j);
                }
                _ => break,
            }
        }
        Some(positions)
    }

    /// The positions of [`Self::generate_symmetry_equivalent_points_from_point`] for each of
    /// `points`, each position once.
    pub fn generate_symmetry_equivalent_points_from_points(
        space_group_number: impl IntoSpaceGroupNumber + Clone,
        points: &[cgmath::Point3<f64>],
    ) -> Option<Vec<cgmath::Point3<f64>>> {
        let space_group_number = space_group_number.into_space_group_number()?;

        let mut new_points = BTreeMap::new();

        for point in points {
            let new_points_ =
                Self::generate_symmetry_equivalent_points_from_point(space_group_number, *point)?;

            new_points.extend(
                new_points_
                    .into_iter()
                    .map(|point| (point3_decode(point), point)),
            );
        }

        Some(new_points.into_values().collect())
    }
}

/// A space group given by its number or its symbol.
pub trait IntoSpaceGroupNumber {
    /// The space group's number (1–230), or `None` for a number out of range or a symbol no
    /// space group has: a symbol of [`SPACEGROUP_NUMBERS`] or of any setting (see
    /// [`SpaceGroupSetting::named`]: `P 21 21 21`, `P 1 21/n 1`, `C m c a`), its suffix
    /// (`R -3 m :H`, `F d -3 m :2`) ignored.
    fn into_space_group_number(&self) -> Option<u8>;
}

macro_rules! impl_into_space_group_number {
    ($($t:ty),*) => {
        $(
            impl IntoSpaceGroupNumber for $t {
                fn into_space_group_number(&self) -> Option<u8> {
                    let id = u8::try_from(*self).ok()?;

                    if SPACEGROUP_SYMBOLS.contains_key(&id) {
                        Some(id)
                    } else {
                        None
                    }
                }
            }
        )*
    };
    () => {

    };
}

impl_into_space_group_number!(
    usize, u8, u16, u32, u64, u128, isize, i8, i16, i32, i64, i128
);

/// Maps a space-group number (1–230) to a Laue class label and optional unique-axis label.
///
/// Labels match the strings used by `SGLaue` / `SGUniq` in GSAS-II.
/// Returns `("-1", None)` for unrecognised or out-of-range numbers.
pub fn laue_and_unique_axis_from_sg_number(n: u32) -> (&'static str, Option<&'static str>) {
    match n {
        1 | 2 => ("-1", None),
        3..=15 => ("2/m", Some("b")),
        16..=74 => ("mmm", None),
        75..=88 => ("4/m", None),
        89..=142 => ("4/mmm", None),
        143..=148 => ("3", None),
        149 | 151 | 153 | 157 | 159 | 162 | 163 => ("31m", None),
        150 | 152 | 154 | 155 | 156 | 158 | 160 | 161 | 164 | 165 | 166 | 167 => ("3m1", None),
        168..=176 => ("6/m", None),
        177..=194 => ("6/mmm", None),
        195..=206 => ("m3", None),
        207..=230 => ("m3m", None),
        _ => ("-1", None),
    }
}

/// Returns symmetry-independent lattice component indices for a Laue class.
///
/// Component indices follow GSAS-II's reciprocal metric terms
/// `A = [G11, G22, G33, 2·G12, 2·G13, 2·G23]`: `0=a`, `1=b`, `2=c`, `3=gamma`, `4=beta`, `5=alpha`.
pub fn independent_lattice_components_for_laue(laue: &str, uniq: Option<&str>) -> Vec<u8> {
    match laue {
        "-1" => vec![0, 1, 2, 3, 4, 5],
        "2/m" => match uniq {
            Some("a") => vec![0, 1, 2, 5],
            Some("b") => vec![0, 1, 2, 4],
            _ => vec![0, 1, 2, 3],
        },
        "mmm" => vec![0, 1, 2],
        "4/m" | "4/mmm" => vec![0, 2],
        "6/m" | "6/mmm" | "3m1" | "31m" | "3" => vec![0, 2],
        "3R" | "3mR" => vec![0, 3],
        "m3m" | "m3" => vec![0],
        _ => vec![0, 1, 2, 3, 4, 5],
    }
}

/// Returns symmetry-independent lattice component indices for a space-group number.
pub fn independent_lattice_components_for_space_group_number(space_group_number: u32) -> Vec<u8> {
    let (laue, uniq) = laue_and_unique_axis_from_sg_number(space_group_number);
    independent_lattice_components_for_laue(laue, uniq)
}

/// Maps a unit-cell parameter name to its component index (0–5).
///
/// The index order matches GSAS-II convention and `independent_lattice_components_for_laue`:
/// `"a"=0`, `"b"=1`, `"c"=2`, `"gamma"=3`, `"beta"=4`, `"alpha"=5`.
///
/// Returns `None` for unrecognised names.
pub fn lattice_component_index(parameter_name: &str) -> Option<u8> {
    match parameter_name {
        "a" => Some(0),
        "b" => Some(1),
        "c" => Some(2),
        "gamma" => Some(3),
        "beta" => Some(4),
        "alpha" => Some(5),
        _ => None,
    }
}

impl IntoSpaceGroupNumber for String {
    fn into_space_group_number(&self) -> Option<u8> {
        let name = self.split(':').next().unwrap_or_default();

        SPACEGROUP_NUMBERS
            .get(&compact_symbol(name))
            .copied()
            .or_else(|| Some(settings_named(name).first()?.number))
    }
}

/// One space-separated part of a space group symbol with its screw axes as
/// [`SPACEGROUP_SYMBOLS`] writes them, `_` between their two digits (`2_1/c`): two digits that
/// are the whole part (`21`, `42`) or come before a `/` (`21/c`, `P63/mmc`). Elsewhere two digits
/// may be two axes (`P321`), so they are kept.
fn screw_axes_as_in_the_table(part: &str) -> std::borrow::Cow<'_, str> {
    let bytes = part.as_bytes();
    let screw_at = |i: usize| {
        matches!(bytes[i], b'2' | b'3' | b'4' | b'6')
            && bytes
                .get(i + 1)
                .is_some_and(|&next| (b'1'..bytes[i]).contains(&next))
            && (i == 0 || !bytes[i - 1].is_ascii_digit())
            && (bytes.len() == 2 || bytes.get(i + 2) == Some(&b'/'))
    };

    if !(0..bytes.len()).any(screw_at) {
        return part.into();
    }

    let mut written = String::with_capacity(part.len() + 2);
    for (i, c) in part.char_indices() {
        written.push(c);
        if screw_at(i) {
            written.push('_');
        }
    }
    written.into()
}

impl IntoSpaceGroupNumber for &str {
    fn into_space_group_number(&self) -> Option<u8> {
        self.to_string().into_space_group_number()
    }
}

/// A space group given by its number.
pub trait IntoSpaceGroupSymbol {
    /// The space group's symbol (see [`SPACEGROUP_SYMBOLS`]), or `None` for a number out of
    /// 1–230.
    fn into_space_group_symbol(&self) -> Option<&str>;
}

macro_rules! impl_into_space_group_symbol {
    ($($t:ty),*) => {
        $(
            impl IntoSpaceGroupSymbol for $t {
                fn into_space_group_symbol(&self) -> Option<&str> {
                    let id = u8::try_from(*self).ok()?;

                    SPACEGROUP_SYMBOLS.get(&id).map(|s| s.as_str())
                }
            }
        )*
    };
    () => {

    };
}

impl_into_space_group_symbol!(
    usize, u8, u16, u32, u64, u128, isize, i8, i16, i32, i64, i128
);

#[cfg(test)]
mod test_space_group_numbers {
    use super::{IntoSpaceGroupNumber, IntoSpaceGroupSymbol};

    #[test]
    fn a_number_out_of_range_names_no_space_group() {
        assert_eq!(225u16.into_space_group_symbol(), Some("Fm-3m"));
        assert_eq!(300u16.into_space_group_symbol(), None);
        assert_eq!(300u16.into_space_group_number(), None);
        assert_eq!((-31i32).into_space_group_number(), None);
        assert_eq!("F m -3 m".into_space_group_number(), Some(225));
    }

    #[test]
    fn a_symbol_may_write_screw_axes_and_settings_as_cifs_do() {
        for (symbol, number) in [
            ("P 21 21 21", 19),
            ("P2_12_12_1", 19),
            ("P 1 21/c 1", 14),
            ("P 2 2 2", 16),
            ("P 4 2 2", 89),
            ("P 42 2 2", 93),
            ("P 41", 76),
            ("P 61", 169),
            ("P 63/m m c", 194),
            ("R -3 m :H", 166),
            ("R -3 m:R", 166),
            ("P 1 21/n 1", 14),
            ("P 21/n", 14),
            ("P21/c", 14),
            ("P63/mmc", 194),
            ("P321", 150),
            ("P 21 c n", 33),
            ("C m c a", 64),
            ("A b a 2", 41),
            ("F m 3 m", 225),
            ("Fd3m", 227),
            ("p b n m", 62),
        ] {
            assert_eq!(symbol.into_space_group_number(), Some(number), "{symbol}");
        }
        for symbol in ["", "X4bm", "?P?", "P 4 3", "Pmn2"] {
            assert_eq!(symbol.into_space_group_number(), None, "{symbol}");
        }
    }

    #[test]
    fn every_tabulated_symbol_names_its_number() {
        for (number, symbol) in super::SPACEGROUP_SYMBOLS.iter() {
            assert_eq!(symbol.into_space_group_number(), Some(*number), "{symbol}");
        }
    }
}

#[cfg(test)]
mod test_space_group_settings {
    use std::collections::BTreeSet;

    use super::{
        SETTINGS_BY_NAME, SPACEGROUP_SETTINGS, SPACEGROUP_SYMMETRY_OPERATIONS, SpaceGroupSetting,
        SpaceGroupSymmetryOperations, parse_symmetry_operation,
    };
    use crate::{Cell, Measured};

    /// An operation's rotation and its translation in twelfths, modulo 1.
    fn key(operation: &cgmath::Matrix4<f64>) -> [i64; 12] {
        let mut key = [0; 12];
        for row in 0..3 {
            for column in 0..3 {
                key[3 * row + column] = operation[column][row].round() as i64;
            }
            key[9 + row] = (operation[3][row] * 12.0).round() as i64 % 12;
            key[9 + row] = (key[9 + row] + 12) % 12;
        }
        key
    }

    #[test]
    fn every_setting_is_a_group_and_each_number_has_its_settings() {
        assert_eq!(SPACEGROUP_SETTINGS.len(), 530);
        let numbers: BTreeSet<u8> = SPACEGROUP_SETTINGS.iter().map(|s| s.number).collect();
        assert_eq!(numbers, (1..=230).collect());

        for setting in SPACEGROUP_SETTINGS.iter() {
            let operations: BTreeSet<_> = setting.operations.iter().map(key).collect();
            assert_eq!(
                operations.len(),
                setting.operations.len(),
                "{}",
                setting.hall_symbol
            );
            for a in &setting.operations {
                for b in &setting.operations {
                    assert!(
                        operations.contains(&key(&(a * b))),
                        "{} is closed",
                        setting.hall_symbol
                    );
                }
            }
        }
    }

    #[test]
    fn a_number_and_its_tabulated_symbol_name_the_standard_setting() {
        assert_eq!(SPACEGROUP_SYMMETRY_OPERATIONS.len(), 230);
        for (number, operations) in SPACEGROUP_SYMMETRY_OPERATIONS.iter() {
            let standard = SpaceGroupSetting::standard(*number).unwrap();
            assert_eq!(operations, &standard.operations);

            let symbol = &super::SPACEGROUP_SYMBOLS[number];
            assert_eq!(SpaceGroupSetting::named(symbol), Some(standard), "{symbol}");
        }

        let choice = |number| SpaceGroupSetting::standard(number).map(|s| s.choice);
        assert_eq!(choice(227), Some("2"));
        assert_eq!(choice(14), Some("b1"));
        assert_eq!(choice(166), Some("H"));
        assert_eq!(choice(231), None);
    }

    #[test]
    fn no_name_names_settings_of_two_groups() {
        for (name, settings) in SETTINGS_BY_NAME.iter() {
            assert!(
                settings.iter().all(|s| s.number == settings[0].number),
                "{name} names settings of several groups"
            );
        }
    }

    #[test]
    fn a_symbol_names_its_setting() {
        let hall = |symbol: &str| SpaceGroupSetting::named(symbol).map(|s| s.hall_symbol);

        for (symbol, expected) in [
            ("P 1 21/n 1", "-P 2yn"),
            ("P 21/n", "-P 2yn"),
            ("P 1 21/c 1", "-P 2ybc"),
            ("P 2/m", "-P 2y"),
            ("P 1 1 2/m", "-P 2"),
            ("P b n m", "-P 2c 2ab"),
            ("P 21 c n", "P -2n 2a"),
            ("C m c a", "-C 2ac 2"),
            ("C m c b", "-C 2ac 2"),
            ("A b a 2", "A 2 -2ab"),
            ("F m 3 m", "-F 4 2 3"),
            ("F d -3 m", "-F 4vw 2vw 3"),
            ("F d -3 m :1", "F 4d 2 3 -1d"),
            ("F d -3 m:2", "-F 4vw 2vw 3"),
            ("R -3 m", "-R 3 2\""),
            ("R -3 m :R", "-P 3* 2"),
        ] {
            assert_eq!(hall(symbol), Some(expected), "{symbol}");
        }
        assert_eq!(hall("F d -3 m :3"), None);
        assert_eq!(hall("X4bm"), None);
    }

    fn cell(symbol: &str, number: u8, angles: [f64; 3]) -> Cell {
        Cell {
            alpha: Measured::new(angles[0]),
            beta: Measured::new(angles[1]),
            gamma: Measured::new(angles[2]),
            space_group: symbol.to_string(),
            space_group_number: number,
            ..Default::default()
        }
    }

    #[test]
    fn a_cell_without_operations_is_in_the_setting_its_symbol_names() {
        let operations = SpaceGroupSymmetryOperations::of_cell(&cell("P 21/n", 14, [90.0; 3]));
        let n_glide = parse_symmetry_operation("x+1/2,-y+1/2,z+1/2").unwrap();
        assert!(operations.unwrap().iter().any(|o| key(o) == key(&n_glide)));

        let choice = |symbol, number, angles| {
            SpaceGroupSetting::of_cell(&cell(symbol, number, angles)).map(|s| s.choice)
        };
        assert_eq!(choice("R -3 m", 166, [90.0, 90.0, 120.0]), Some("H"));
        assert_eq!(choice("R -3 m", 166, [57.3; 3]), Some("R"));
        assert_eq!(choice("R -3 m", 0, [57.3; 3]), Some("R"));
        assert_eq!(choice("P 21/n", 62, [90.0; 3]), None);
    }
}

#[cfg(test)]
mod test_lattice_components {
    use super::{independent_lattice_components_for_laue, lattice_component_index};

    #[test]
    fn the_unique_monoclinic_angle_is_independent() {
        for (uniq, angle) in [("a", "alpha"), ("b", "beta"), ("c", "gamma")] {
            let independent = independent_lattice_components_for_laue("2/m", Some(uniq));
            let component = lattice_component_index(angle).unwrap();
            assert!(
                independent.contains(&component),
                "{angle} ({component}) missing from {independent:?} for unique axis {uniq}"
            );
        }
    }
}

#[cfg(test)]
mod test_spacegroup_symmetry_operations {
    use std::collections::BTreeMap;

    use crate::symmetry::integer_decode;

    use super::SPACEGROUP_SYMMETRY_OPERATIONS;

    #[test]
    fn test_serde() {
        let symmetry_operations = SPACEGROUP_SYMMETRY_OPERATIONS.get(&1).unwrap();

        let spacegroup_1 = symmetry_operations[0];

        let expected = [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ]
        .into();

        assert_eq!(spacegroup_1, expected);

        let symbols = super::SPACEGROUP_SYMBOLS.get(&1).unwrap();

        assert_eq!(symbols, "P1");
    }

    #[test]
    #[should_panic]
    fn test_serde_fail() {
        super::SPACEGROUP_SYMBOLS.get(&0).unwrap();
    }

    #[test]
    fn test_dedup() {
        let x = [1.0, 2.0, 3.0, 3.0, 3.00002, 3.00002, 3.1];

        let y = x
            .iter()
            .map(|&x| {
                let (mantissa, exponent, sign) = integer_decode(x);

                ((mantissa, exponent, sign), x)
            })
            .collect::<BTreeMap<_, _>>()
            .into_values()
            .collect::<Vec<_>>();

        assert_eq!(y, [1.0, 2.0, 3.0, 3.00002, 3.1]);
    }

    #[test]
    #[ignore = "This test is for visualization purposes only"]
    fn plot_symmetry_equivalent_points() {
        let point = cgmath::Point3::new(0.0, 0.0, 0.0);

        let new_points =
            super::SpaceGroupSymmetryOperations::generate_symmetry_equivalent_points_from_point(
                227, point,
            )
            .unwrap();

        dbg!(&new_points);

        assert_eq!(new_points.len(), 38);

        let x = new_points.iter().map(|p| p.x).collect::<Vec<_>>();
        let y = new_points.iter().map(|p| p.y).collect::<Vec<_>>();
        let z = new_points.iter().map(|p| p.z).collect::<Vec<_>>();

        let mut figure = gnuplot::Figure::new();

        let axes = figure.axes3d();

        axes.points(x, y, z, &[]);

        figure.show().unwrap();
    }

    #[test]
    fn test_symmetry_equivalent_points() {
        // P1: an atom at the origin shows at all eight corners of the cell.
        let new_points =
            super::SpaceGroupSymmetryOperations::generate_symmetry_equivalent_points_from_point(
                1,
                cgmath::Point3::new(0.0, 0.0, 0.0),
            )
            .unwrap();
        assert_eq!(new_points.len(), 8);
        assert!(new_points.contains(&cgmath::Point3::new(1.0, 1.0, 1.0)));

        // P2 (unique axis b): (x, y, z) and (-x, y, -z), wrapped into the cell.
        let new_points =
            super::SpaceGroupSymmetryOperations::generate_symmetry_equivalent_points_from_point(
                3,
                cgmath::Point3::new(0.25, 0.5, 0.25),
            )
            .unwrap();
        assert_eq!(new_points.len(), 2);
        assert!(new_points.contains(&cgmath::Point3::new(0.75, 0.5, 0.75)));
    }

    #[test]
    fn rounded_special_positions_are_one_position() {
        // P6_3/mmc, Wyckoff 2c (1/3, 2/3, 1/4) as CIFs write it: two positions, not more.
        let new_points =
            super::SpaceGroupSymmetryOperations::generate_symmetry_equivalent_points_from_point(
                194,
                cgmath::Point3::new(0.3333, 0.6667, 0.25),
            )
            .unwrap();
        assert_eq!(new_points.len(), 2, "{new_points:?}");
    }

    #[test]
    fn a_site_a_little_off_its_special_position_has_its_stated_multiplicity() {
        // Pnma, Wyckoff 4c (x, 1/4, z) written with y = 0.249: its mirror copies are 0.015 Å
        // apart.
        let cell = crate::Cell {
            a: 5.4405.into(),
            b: 7.6889.into(),
            c: 5.4499.into(),
            alpha: 90.0.into(),
            beta: 90.0.into(),
            gamma: 90.0.into(),
            space_group_number: 62,
            ..Default::default()
        };
        let positions = |point, multiplicity| {
            super::SpaceGroupSymmetryOperations::site_positions_in_cell(&cell, point, multiplicity)
                .unwrap()
                .len()
        };

        let off = cgmath::Point3::new(0.486, 0.249, 0.064);
        assert_eq!(positions(off, None), 8);
        assert_eq!(positions(off, Some(4)), 4);

        // Further off, the copies are different atoms (0.77 Å apart): no merging.
        let far = cgmath::Point3::new(0.486, 0.2, 0.064);
        assert_eq!(positions(far, Some(4)), 8);
    }

    #[test]
    fn a_split_site_keeps_its_stated_multiplicity() {
        // P2: a site split about the 2-fold axis, its copies 0.011 Å apart and so within
        // SAME_POSITION of each other.
        let cell = crate::Cell {
            a: 10.0.into(),
            b: 10.0.into(),
            c: 10.0.into(),
            alpha: 90.0.into(),
            beta: 90.0.into(),
            gamma: 90.0.into(),
            space_group_number: 3,
            ..Default::default()
        };
        let point = cgmath::Point3::new(0.0004, 0.3, 0.0004);
        let positions = |multiplicity| {
            super::SpaceGroupSymmetryOperations::site_positions_in_cell(&cell, point, multiplicity)
                .unwrap()
                .len()
        };

        assert_eq!(positions(None), 1);
        assert_eq!(positions(Some(2)), 2);
    }

    #[test]
    fn symmetry_operations_read_as_cifs_write_them() {
        use cgmath::Transform;

        let operation = super::parse_symmetry_operation("'-x+y, -x, z+1/3'").unwrap();
        let point = operation.transform_point(cgmath::Point3::new(0.1, 0.4, 0.2));
        assert!((point.x - 0.3).abs() < 1e-12);
        assert!((point.y + 0.1).abs() < 1e-12);
        assert!((point.z - (0.2 + 1.0 / 3.0)).abs() < 1e-12);

        let operation = super::parse_symmetry_operation("1/2+X, y, -z+0.5").unwrap();
        let point = operation.transform_point(cgmath::Point3::new(0.1, 0.4, 0.2));
        assert_eq!(point, cgmath::Point3::new(0.6, 0.4, 0.3));

        assert!(super::parse_symmetry_operation("x, y").is_none());
        assert!(super::parse_symmetry_operation("x, y, q").is_none());
    }

    #[test]
    fn a_cell_with_its_own_operations_generates_in_its_setting() {
        // P4/nmm in origin choice 2: Al at (0, 0, 0) is on Wyckoff 4d there, but origin choice 1
        // puts that point on a 2-fold position.
        let operations = [
            "y+1/2, x+1/2, -z",
            "-y, x+1/2, -z",
            "y+1/2, -x, -z",
            "-y, -x, -z",
            "x+1/2, y+1/2, -z",
            "-x, y+1/2, -z",
            "x+1/2, -y, -z",
            "-x, -y, -z",
            "-y+1/2, -x+1/2, z",
            "y, -x+1/2, z",
            "-y+1/2, x, z",
            "y, x, z",
            "-x+1/2, -y+1/2, z",
            "x, -y+1/2, z",
            "-x+1/2, y, z",
            "x, y, z",
        ];
        let cell = crate::Cell {
            space_group: "P 4/n m m Z".to_string(),
            space_group_number: 129,
            symmetry_operations: operations.map(str::to_string).to_vec(),
            ..Default::default()
        };
        let distinct = |points: Vec<cgmath::Point3<f64>>| {
            points
                .into_iter()
                .filter(|p| p.x < 1.0 && p.y < 1.0 && p.z < 1.0)
                .count()
        };

        let in_cell =
            super::SpaceGroupSymmetryOperations::generate_symmetry_equivalent_points_in_cell(
                &cell,
                cgmath::Point3::new(0.0, 0.0, 0.0),
            )
            .unwrap();
        assert_eq!(distinct(in_cell), 4);

        let origin_choice_1 = super::SpaceGroupSetting::named("P 4/n m m :1").unwrap();
        let in_origin_choice_1 = super::equivalent_points(
            &origin_choice_1.operations,
            cgmath::Point3::new(0.0, 0.0, 0.0),
        );
        assert_ne!(distinct(in_origin_choice_1), 4);
    }

    #[test]
    fn a_cell_without_operations_uses_its_space_group_number() {
        let cell = crate::Cell {
            space_group_number: 194,
            ..Default::default()
        };

        assert_eq!(
            super::SpaceGroupSymmetryOperations::of_cell(&cell),
            super::SpaceGroupSymmetryOperations::get(194).cloned()
        );
    }
}

mod test_spacegroup_generators {

    #[test]
    #[ignore = "This test is for visualization purposes only"]
    fn plot_symmetry_equivalent_points() {
        let point = cgmath::Point3::new(0.0, 0.0, 0.0);

        let new_points =
            super::SpaceGroupGenerators::generate_symmetry_equivalent_points_from_point(227, point)
                .unwrap();

        dbg!(&new_points);

        assert_eq!(new_points.len(), 6);

        let x = new_points.iter().map(|p| p.x).collect::<Vec<_>>();

        let y = new_points.iter().map(|p| p.y).collect::<Vec<_>>();

        let z = new_points.iter().map(|p| p.z).collect::<Vec<_>>();

        let mut figure = gnuplot::Figure::new();

        let axes = figure.axes3d();

        axes.points(x, y, z, &[]);

        figure.show().unwrap();
    }
}

fn integer_decode(val: f64) -> (u64, i16, i8) {
    let bits = val.to_bits();
    let sign: i8 = if bits >> 63 == 0 { 1 } else { -1 };
    let mut exponent: i16 = ((bits >> 52) & 0x7ff) as i16;
    let mantissa = if exponent == 0 {
        (bits & 0xfffffffffffff) << 1
    } else {
        (bits & 0xfffffffffffff) | 0x10000000000000
    };

    exponent -= 1023 + 52;
    (mantissa, exponent, sign)
}

fn point3_decode(point: cgmath::Point3<f64>) -> [(u64, i16, i8); 3] {
    [
        integer_decode(point.x),
        integer_decode(point.y),
        integer_decode(point.z),
    ]
}

fn move_point_into_unit_cell(point: cgmath::Point3<f64>) -> cgmath::Point3<f64> {
    let x = match (point.x - 1.0).abs() < f64::EPSILON {
        true => 1.0,
        false => point.x.rem_euclid(1.0),
    };

    let y = match (point.y - 1.0).abs() < f64::EPSILON {
        true => 1.0,
        false => point.y.rem_euclid(1.0),
    };

    let z = match (point.z - 1.0).abs() < f64::EPSILON {
        true => 1.0,
        false => point.z.rem_euclid(1.0),
    };

    cgmath::Point3::new(x, y, z)
}

#[cfg(test)]
mod test_move_point_into_unit_cell {
    #[test]
    fn test_move_point_into_unit_cell() {
        let point = cgmath::Point3::new(0.0, 0.0, 0.0);
        let point = super::move_point_into_unit_cell(point);
        assert_eq!(point, cgmath::Point3::new(0.0, 0.0, 0.0));

        let point = cgmath::Point3::new(1.5, 1.5, 1.5);
        let point = super::move_point_into_unit_cell(point);
        assert_eq!(point, cgmath::Point3::new(0.5, 0.5, 0.5));

        let point = cgmath::Point3::new(-1.5, -1.5, -1.5);
        let point = super::move_point_into_unit_cell(point);
        assert_eq!(point, cgmath::Point3::new(0.5, 0.5, 0.5));

        let point = cgmath::Point3::new(2.5, 2.5, 2.5);
        let point = super::move_point_into_unit_cell(point);
        assert_eq!(point, cgmath::Point3::new(0.5, 0.5, 0.5));

        let point = cgmath::Point3::new(1.0, 1.0, 1.0);
        let point = super::move_point_into_unit_cell(point);
        assert_eq!(point, cgmath::Point3::new(1.0, 1.0, 1.0));

        let point = cgmath::Point3::new(-0.75, -0.75, -0.75);
        let point = super::move_point_into_unit_cell(point);
        assert_eq!(point, cgmath::Point3::new(0.25, 0.25, 0.25));

        let point = cgmath::Point3::new(-1.75, -0.75, -0.75);
        let point = super::move_point_into_unit_cell(point);
        assert_eq!(point, cgmath::Point3::new(0.25, 0.25, 0.25));
    }
}

static TRIVIAL_TRANSLATIONS: LazyLock<[cgmath::Matrix4<f64>; 8]> = LazyLock::new(|| {
    [
        cgmath::Matrix4::from_translation(cgmath::Vector3::new(0.0, 0.0, 0.0)),
        cgmath::Matrix4::from_translation(cgmath::Vector3::new(1.0, 0.0, 0.0)),
        cgmath::Matrix4::from_translation(cgmath::Vector3::new(0.0, 1.0, 0.0)),
        cgmath::Matrix4::from_translation(cgmath::Vector3::new(0.0, 0.0, 1.0)),
        cgmath::Matrix4::from_translation(cgmath::Vector3::new(1.0, 1.0, 0.0)),
        cgmath::Matrix4::from_translation(cgmath::Vector3::new(0.0, 1.0, 1.0)),
        cgmath::Matrix4::from_translation(cgmath::Vector3::new(1.0, 0.0, 1.0)),
        cgmath::Matrix4::from_translation(cgmath::Vector3::new(1.0, 1.0, 1.0)),
    ]
});
