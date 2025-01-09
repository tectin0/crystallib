use std::{collections::BTreeMap, sync::LazyLock};

use cgmath::Transform;

const SPACEGROUP_SYMMETRY_OPERATIONS_RAW: &str =
    include_str!("../assets/spacegroup_symmetry_operations.json");

pub static SPACEGROUP_SYMMETRY_OPERATIONS: LazyLock<BTreeMap<usize, Vec<cgmath::Matrix4<f64>>>> =
    LazyLock::new(|| serde_json::from_str(SPACEGROUP_SYMMETRY_OPERATIONS_RAW).unwrap());

const SPACEGROUP_SYMBOLS_RAW: &str = include_str!("../assets/spacegroup_symbols.json");

pub static SPACEGROUP_SYMBOLS: LazyLock<BTreeMap<usize, String>> = LazyLock::new(|| {
    let symbols: Vec<String> = serde_json::from_str(SPACEGROUP_SYMBOLS_RAW).unwrap();

    symbols
        .into_iter()
        .enumerate()
        .map(|(i, s)| (i + 1, s))
        .collect()
});

pub static SPACEGROUP_NUMBERS: LazyLock<BTreeMap<String, usize>> = LazyLock::new(|| {
    SPACEGROUP_SYMBOLS
        .iter()
        .map(|(i, s)| (s.clone(), *i))
        .collect()
});

pub struct SpaceGroupSymmetryOperations {}

impl SpaceGroupSymmetryOperations {
    pub fn get<'a>(
        space_group_number: impl IntoSpaceGroupNumber,
    ) -> Option<&'a Vec<cgmath::Matrix4<f64>>> {
        let space_group_number = space_group_number.into_space_group_number()?;

        SPACEGROUP_SYMMETRY_OPERATIONS.get(&space_group_number)
    }

    pub fn get_all<'a>() -> &'a BTreeMap<usize, Vec<cgmath::Matrix4<f64>>> {
        &SPACEGROUP_SYMMETRY_OPERATIONS
    }

    pub fn get_symbol<'a>(space_group_number: impl IntoSpaceGroupNumber) -> Option<&'a str> {
        let space_group_number = space_group_number.into_space_group_number()?;

        SPACEGROUP_SYMBOLS
            .get(&space_group_number)
            .map(|s| s.as_str())
    }

    pub fn get_number<'a>(space_group_symbol: impl IntoSpaceGroupSymbol) -> Option<usize> {
        let space_group_symbol = space_group_symbol.into_space_group_symbol()?;

        SPACEGROUP_NUMBERS.get(space_group_symbol).copied()
    }

    pub fn generate_symmetry_equivalent_points_from_point(
        space_group_number: impl IntoSpaceGroupNumber,
        point: cgmath::Point3<f64>,
    ) -> Option<Vec<cgmath::Point3<f64>>> {
        let symmetry_operations = Self::get(space_group_number)?;

        let new_points: BTreeMap<[(u64, i16, i8); 3], cgmath::Point3<f64>> = symmetry_operations
            .iter()
            .map(|m| {
                let new_point = m.transform_point(point);

                (
                    [
                        integer_decode(new_point.x),
                        integer_decode(new_point.y),
                        integer_decode(new_point.z),
                    ],
                    new_point,
                )
            })
            .collect();

        Some(new_points.into_values().collect())
    }

    pub fn generate_symmetry_equivalent_points_from_points(
        space_group_number: impl IntoSpaceGroupNumber,
        points: &[cgmath::Point3<f64>],
    ) -> Option<Vec<Vec<cgmath::Point3<f64>>>> {
        let symmetry_operations = Self::get(space_group_number)?;

        let new_points: BTreeMap<[(u64, i16, i8); 3], Vec<cgmath::Point3<f64>>> =
            symmetry_operations
                .iter()
                .map(|m| {
                    let new_points: Vec<cgmath::Point3<f64>> = points
                        .iter()
                        .map(|point| {
                            

                            m.transform_point(*point)
                        })
                        .collect();

                    (
                        [
                            integer_decode(new_points[0].x),
                            integer_decode(new_points[0].y),
                            integer_decode(new_points[0].z),
                        ],
                        new_points,
                    )
                })
                .collect();

        Some(new_points.into_values().collect())
    }
}

pub trait IntoSpaceGroupNumber {
    fn into_space_group_number(&self) -> Option<usize>;
}

macro_rules! impl_into_space_group_number {
    ($($t:ty),*) => {
        $(
            impl IntoSpaceGroupNumber for $t {
                fn into_space_group_number(&self) -> Option<usize> {
                    let id = *self as usize;

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

impl_into_space_group_number!(usize, u8, u16, u32, u64, u128, isize, i8, i16, i32, i64, i128);

impl IntoSpaceGroupNumber for String {
    fn into_space_group_number(&self) -> Option<usize> {
        let split = self.split_whitespace();
        let number_of_splits = split.clone().count();

        let space_group = match number_of_splits {
            0 => None,
            _ => Some(split.collect::<Vec<_>>().join("")),
        };

        let space_group = space_group.as_ref().unwrap_or(self);

        SPACEGROUP_NUMBERS.get(space_group).copied()
    }
}

impl IntoSpaceGroupNumber for &str {
    fn into_space_group_number(&self) -> Option<usize> {
        self.to_string().into_space_group_number()
    }
}

pub trait IntoSpaceGroupSymbol {
    fn into_space_group_symbol(&self) -> Option<&str>;
}

macro_rules! impl_into_space_group_symbol {
    ($($t:ty),*) => {
        $(
            impl IntoSpaceGroupSymbol for $t {
                fn into_space_group_symbol(&self) -> Option<&str> {
                    let id = *self as usize;

                    SPACEGROUP_SYMBOLS.get(&id).map(|s| s.as_str())
                }
            }
        )*
    };
    () => {

    };
}

impl_into_space_group_symbol!(usize, u8, u16, u32, u64, u128, isize, i8, i16, i32, i64, i128);

#[cfg(test)]
mod test_spacegroup_symmetry_operations {
    use std::collections::BTreeMap;

    use cgmath::Matrix;
    use gnuplot::AxesCommon;

    use crate::{symmetry::integer_decode, SpaceGroupSymmetryOperations};

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
            .into_iter()
            .map(|(_, v)| v)
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
        let point = cgmath::Point3::new(0.0, 0.0, 0.0);

        let new_points =
            super::SpaceGroupSymmetryOperations::generate_symmetry_equivalent_points_from_point(
                1, point,
            )
            .unwrap();

        let expected = vec![point];

        assert_eq!(new_points, expected);

        let point = cgmath::Point3::new(0.5, 0.5, 0.5);

        let new_points =
            super::SpaceGroupSymmetryOperations::generate_symmetry_equivalent_points_from_point(
                3, point,
            )
            .unwrap();

        let expected = vec![
            cgmath::Point3::new(-0.5, 0.5, -0.5),
            cgmath::Point3::new(0.5, 0.5, 0.5),
        ];

        assert_eq!(new_points, expected);
    }
}

fn integer_decode(val: f64) -> (u64, i16, i8) {
    let bits: u64 = unsafe { std::mem::transmute(val) };
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
