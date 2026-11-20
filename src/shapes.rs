//! The shape of the atoms bonded to a centre, among the known coordination polyhedra
//! (tetrahedron, octahedron, …), by the angles between their directions from it.

use std::f32::consts::FRAC_1_SQRT_2;
use std::fmt::{Display, Formatter};
use std::sync::LazyLock;

use cgmath::{InnerSpace, Vector3, vec3};

/// The largest root-mean-square difference between a polyhedron's signature (the cosines of the
/// angles between each two of its bonds, see [`classify`]) and a reference shape's at which it
/// counts as that shape; a polyhedron farther from every shape has none.
pub const MAX_SHAPE_ERROR: f32 = 0.1;

/// The shape the atoms at `offsets` (from a centre) form, where it is one of the known ones within
/// [`MAX_SHAPE_ERROR`].
pub fn shape_of(offsets: &[Vector3<f32>]) -> Option<Shape> {
    classify(offsets, MAX_SHAPE_ERROR).map(|classification| classification.shape)
}

#[derive(Clone, Copy, Hash, Eq, PartialEq, Debug)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Shape {
    /// Linear
    Linear,
    /// TrigonalPlanar
    TrigonalPlanar,
    /// Tetrahedron
    Tetrahedron,
    /// SquarePlanar
    SquarePlanar,
    /// TrigonalBipyramidal
    TrigonalBipyramidal,
    /// SquarePyramidal
    SquarePyramidal,
    /// Octahedron
    Octahedron,
    /// TrigonalPrismatic
    TrigonalPrismatic,
    /// PentagonalBipyramidal
    PentagonalBipyramidal,
    /// CappedOctahedral
    CappedOctahedral,
    /// CappedTrigonalPrismatic
    CappedTrigonalPrismatic,
    /// SquareAntiprismatic
    SquareAntiprismatic,
    /// Bisdisphenoid
    Bisdisphenoid,
    /// BicappedTrigonalPrismatic
    BicappedTrigonalPrismatic,
    /// Cubic
    Cubic,
    /// HexagonalBipyramidal
    HexagonalBipyramidal,
    /// TransBicappedOctahedral
    TransBicappedOctahedral,
    /// TriangularFaceBicappedTrigonalPrismatic
    TriangularFaceBicappedTrigonalPrismatic,
    /// TricappedTrigonalPrismatic
    TricappedTrigonalPrismatic,
    /// CappedSquareAntiprismatic
    CappedSquareAntiprismatic,
    /// BicappedSquareAntiprismatic
    BicappedSquareAntiprismatic,
    /// BicappedCubic
    BicappedCubic,
    /// Icosahedron
    Icosahedron,
    /// Cuboctahedron
    Cuboctahedron,
    /// Anticuboctahedron
    Anticuboctahedron,
    /// BicappedHexagonalAntiprismatic
    BicappedHexagonalAntiprismatic,
}

impl Display for Shape {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let name = format!("{self:?}");
        for (index, character) in name.char_indices() {
            if index > 0 && character.is_uppercase() {
                write!(f, " {}", character.to_ascii_lowercase())?;
            } else {
                write!(f, "{character}")?;
            }
        }
        Ok(())
    }
}

impl Shape {
    /// A list of all shapes
    pub const ALL_SHAPES: &[Shape] = &[
        Shape::Linear,
        Shape::TrigonalPlanar,
        Shape::Tetrahedron,
        Shape::SquarePlanar,
        Shape::TrigonalBipyramidal,
        Shape::SquarePyramidal,
        Shape::Octahedron,
        Shape::TrigonalPrismatic,
        Shape::PentagonalBipyramidal,
        Shape::CappedOctahedral,
        Shape::CappedTrigonalPrismatic,
        Shape::SquareAntiprismatic,
        Shape::Bisdisphenoid,
        Shape::BicappedTrigonalPrismatic,
        Shape::Cubic,
        Shape::HexagonalBipyramidal,
        Shape::TransBicappedOctahedral,
        Shape::TriangularFaceBicappedTrigonalPrismatic,
        Shape::TricappedTrigonalPrismatic,
        Shape::CappedSquareAntiprismatic,
        Shape::BicappedSquareAntiprismatic,
        Shape::BicappedCubic,
        Shape::Icosahedron,
        Shape::Cuboctahedron,
        Shape::Anticuboctahedron,
        Shape::BicappedHexagonalAntiprismatic,
    ];

    // Canonical representation of all shapes (used to construct the shape reference)
    fn canonical_points(&self) -> Option<Vec<Vector3<f32>>> {
        // The golden ratio, which places an icosahedron's vertices.
        let phi = std::f32::consts::GOLDEN_RATIO;

        match self {
            Shape::TrigonalPlanar => Some(vec![
                vec3(1.0, 0.0, 0.0),
                vec3(-0.50, 0.866025, 0.0),
                vec3(-0.50, -0.866025, 0.0),
            ]),
            Shape::Tetrahedron => Some(vec![
                vec3(1.0, 1.0, 1.0),
                vec3(1.0, -1.0, -1.0),
                vec3(-1.0, 1.0, -1.0),
                vec3(-1.0, -1.0, 1.0),
            ]),
            Shape::SquarePlanar => Some(vec![
                vec3(FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0),
                vec3(-FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0),
                vec3(-FRAC_1_SQRT_2, -FRAC_1_SQRT_2, 0.0),
                vec3(FRAC_1_SQRT_2, -FRAC_1_SQRT_2, 0.0),
            ]),
            Shape::TrigonalBipyramidal => Some(vec![
                vec3(1.0, 0.0, 0.0),
                vec3(-0.50, 0.866025, 0.0),
                vec3(-0.50, -0.866025, 0.0),
                vec3(0.0, 0.0, 1.0),
                vec3(0.0, 0.0, -1.0),
            ]),
            Shape::SquarePyramidal => Some(vec![
                vec3(FRAC_1_SQRT_2, FRAC_1_SQRT_2, -0.350),
                vec3(-FRAC_1_SQRT_2, FRAC_1_SQRT_2, -0.350),
                vec3(-FRAC_1_SQRT_2, -FRAC_1_SQRT_2, -0.350),
                vec3(FRAC_1_SQRT_2, -FRAC_1_SQRT_2, -0.350),
                vec3(0.0, 0.0, 0.950),
            ]),
            Shape::Octahedron => Some(vec![
                vec3(1.0, 0.0, 0.0),
                vec3(-1.0, 0.0, 0.0),
                vec3(0.0, 1.0, 0.0),
                vec3(0.0, -1.0, 0.0),
                vec3(0.0, 0.0, 1.0),
                vec3(0.0, 0.0, -1.0),
            ]),
            Shape::TrigonalPrismatic => Some(vec![
                vec3(1.0, 0.0, 0.50),
                vec3(-0.50, 0.866025, 0.50),
                vec3(-0.50, -0.866025, 0.50),
                vec3(1.0, 0.0, -0.50),
                vec3(-0.50, 0.866025, -0.50),
                vec3(-0.50, -0.866025, -0.50),
            ]),
            Shape::PentagonalBipyramidal => Some(vec![
                vec3(1.0, 0.0, 0.0),
                vec3(0.309017, 0.951057, 0.0),
                vec3(-0.809017, 0.587785, 0.0),
                vec3(-0.809017, -0.587785, 0.0),
                vec3(0.309017, -0.951057, 0.0),
                vec3(0.0, 0.0, 1.0),
                vec3(0.0, 0.0, -1.0),
            ]),
            Shape::CappedOctahedral => Some(vec![
                vec3(1.0, 0.0, 0.0),
                vec3(-1.0, 0.0, 0.0),
                vec3(0.0, 1.0, 0.0),
                vec3(0.0, -1.0, 0.0),
                vec3(0.0, 0.0, 1.0),
                vec3(0.0, 0.0, -1.0),
                vec3(0.850, 0.850, 0.850),
            ]),
            Shape::CappedTrigonalPrismatic => Some(vec![
                vec3(1.0, 0.0, 0.50),
                vec3(-0.50, 0.866025, 0.50),
                vec3(-0.50, -0.866025, 0.50),
                vec3(1.0, 0.0, -0.50),
                vec3(-0.50, 0.866025, -0.50),
                vec3(-0.50, -0.866025, -0.50),
                vec3(0.750, 1.299038, 0.0),
            ]),
            Shape::SquareAntiprismatic => Some(vec![
                vec3(1.0, 0.0, 0.50),
                vec3(0.0, 1.0, 0.50),
                vec3(-1.0, 0.0, 0.50),
                vec3(0.0, -1.0, 0.50),
                vec3(FRAC_1_SQRT_2, FRAC_1_SQRT_2, -0.50),
                vec3(-FRAC_1_SQRT_2, FRAC_1_SQRT_2, -0.50),
                vec3(-FRAC_1_SQRT_2, -FRAC_1_SQRT_2, -0.50),
                vec3(FRAC_1_SQRT_2, -FRAC_1_SQRT_2, -0.50),
            ]),
            Shape::Bisdisphenoid => Some(vec![
                vec3(0.720, 0.0, 1.0),
                vec3(-0.720, 0.0, 1.0),
                vec3(0.0, 0.720, -1.0),
                vec3(0.0, -0.720, -1.0),
                vec3(1.0, 0.0, -0.720),
                vec3(-1.0, 0.0, -0.720),
                vec3(0.0, 1.0, 0.720),
                vec3(0.0, -1.0, 0.720),
            ]),
            Shape::BicappedTrigonalPrismatic => Some(vec![
                vec3(1.0, 0.0, 0.50),
                vec3(-0.50, 0.866025, 0.50),
                vec3(-0.50, -0.866025, 0.50),
                vec3(1.0, 0.0, -0.50),
                vec3(-0.50, 0.866025, -0.50),
                vec3(-0.50, -0.866025, -0.50),
                vec3(0.750, 1.299038, 0.0),
                vec3(-0.750, -1.299038, 0.0),
            ]),
            Shape::Cubic => Some(vec![
                vec3(-1.0, -1.0, -1.0),
                vec3(-1.0, -1.0, 1.0),
                vec3(-1.0, 1.0, -1.0),
                vec3(-1.0, 1.0, 1.0),
                vec3(1.0, -1.0, -1.0),
                vec3(1.0, -1.0, 1.0),
                vec3(1.0, 1.0, -1.0),
                vec3(1.0, 1.0, 1.0),
            ]),
            Shape::HexagonalBipyramidal => Some(vec![
                vec3(1.0, 0.0, 0.0),
                vec3(0.50, 0.866025, 0.0),
                vec3(-0.50, 0.866025, 0.0),
                vec3(-1.0, 0.0, 0.0),
                vec3(-0.50, -0.866025, 0.0),
                vec3(0.50, -0.866025, 0.0),
                vec3(0.0, 0.0, 1.0),
                vec3(0.0, 0.0, -1.0),
            ]),
            Shape::TransBicappedOctahedral => Some(vec![
                vec3(1.0, 0.0, 0.0),
                vec3(-1.0, 0.0, 0.0),
                vec3(0.0, 1.0, 0.0),
                vec3(0.0, -1.0, 0.0),
                vec3(0.0, 0.0, 1.0),
                vec3(0.0, 0.0, -1.0),
                vec3(0.923760, 0.923760, 0.923760),
                vec3(-0.923760, -0.923760, -0.923760),
            ]),
            Shape::TriangularFaceBicappedTrigonalPrismatic => Some(vec![
                vec3(1.0, 0.0, 0.50),
                vec3(-0.50, 0.866025, 0.50),
                vec3(-0.50, -0.866025, 0.50),
                vec3(1.0, 0.0, -0.50),
                vec3(-0.50, 0.866025, -0.50),
                vec3(-0.50, -0.866025, -0.50),
                vec3(0.0, 0.0, 1.30),
                vec3(0.0, 0.0, -1.30),
            ]),
            Shape::TricappedTrigonalPrismatic => Some(vec![
                vec3(1.0, 0.0, 0.50),
                vec3(-0.50, 0.866025, 0.50),
                vec3(-0.50, -0.866025, 0.50),
                vec3(1.0, 0.0, -0.50),
                vec3(-0.50, 0.866025, -0.50),
                vec3(-0.50, -0.866025, -0.50),
                vec3(0.750, 1.299038, 0.0),
                vec3(-1.50, 0.0, 0.0),
                vec3(0.750, -1.299038, 0.0),
            ]),
            Shape::CappedSquareAntiprismatic => Some(vec![
                vec3(1.0, 0.0, 0.50),
                vec3(0.0, 1.0, 0.50),
                vec3(-1.0, 0.0, 0.50),
                vec3(0.0, -1.0, 0.50),
                vec3(FRAC_1_SQRT_2, FRAC_1_SQRT_2, -0.50),
                vec3(-FRAC_1_SQRT_2, FRAC_1_SQRT_2, -0.50),
                vec3(-FRAC_1_SQRT_2, -FRAC_1_SQRT_2, -0.50),
                vec3(FRAC_1_SQRT_2, -FRAC_1_SQRT_2, -0.50),
                vec3(0.0, 0.0, 1.40),
            ]),
            Shape::BicappedSquareAntiprismatic => Some(vec![
                vec3(1.0, 0.0, 0.50),
                vec3(0.0, 1.0, 0.50),
                vec3(-1.0, 0.0, 0.50),
                vec3(0.0, -1.0, 0.50),
                vec3(FRAC_1_SQRT_2, FRAC_1_SQRT_2, -0.50),
                vec3(-FRAC_1_SQRT_2, FRAC_1_SQRT_2, -0.50),
                vec3(-FRAC_1_SQRT_2, -FRAC_1_SQRT_2, -0.50),
                vec3(FRAC_1_SQRT_2, -FRAC_1_SQRT_2, -0.50),
                vec3(0.0, 0.0, 1.40),
                vec3(0.0, 0.0, -1.40),
            ]),
            Shape::BicappedCubic => Some(vec![
                vec3(-1.0, -1.0, -1.0),
                vec3(-1.0, -1.0, 1.0),
                vec3(-1.0, 1.0, -1.0),
                vec3(-1.0, 1.0, 1.0),
                vec3(1.0, -1.0, -1.0),
                vec3(1.0, -1.0, 1.0),
                vec3(1.0, 1.0, -1.0),
                vec3(1.0, 1.0, 1.0),
                vec3(0.0, 0.0, 1.70),
                vec3(0.0, 0.0, -1.70),
            ]),
            Shape::Icosahedron => Some(vec![
                vec3(0.0, 1.0, phi),
                vec3(0.0, -1.0, phi),
                vec3(0.0, 1.0, -phi),
                vec3(0.0, -1.0, -phi),
                vec3(1.0, phi, 0.0),
                vec3(-1.0, phi, 0.0),
                vec3(1.0, -phi, 0.0),
                vec3(-1.0, -phi, 0.0),
                vec3(phi, 0.0, 1.0),
                vec3(-phi, 0.0, 1.0),
                vec3(phi, 0.0, -1.0),
                vec3(-phi, 0.0, -1.0),
            ]),
            Shape::Cuboctahedron => Some(vec![
                vec3(0.0, -1.0, -1.0),
                vec3(0.0, -1.0, 1.0),
                vec3(0.0, 1.0, -1.0),
                vec3(0.0, 1.0, 1.0),
                vec3(-1.0, 0.0, -1.0),
                vec3(1.0, 0.0, -1.0),
                vec3(-1.0, 0.0, 1.0),
                vec3(1.0, 0.0, 1.0),
                vec3(-1.0, -1.0, 0.0),
                vec3(-1.0, 1.0, 0.0),
                vec3(1.0, -1.0, 0.0),
                vec3(1.0, 1.0, 0.0),
            ]),
            Shape::Anticuboctahedron => Some(vec![
                vec3(1.0, 0.0, 1.0),
                vec3(-0.50, 0.866025, 1.0),
                vec3(-0.50, -0.866025, 1.0),
                vec3(0.995929, 0.5750, 0.0),
                vec3(0.0, 1.150, 0.0),
                vec3(-0.995929, 0.5750, 0.0),
                vec3(-0.995929, -0.5750, 0.0),
                vec3(0.0, -1.150, 0.0),
                vec3(0.995929, -0.5750, 0.0),
                vec3(-1.0, 0.0, -1.0),
                vec3(0.50, -0.866025, -1.0),
                vec3(0.50, 0.866025, -1.0),
            ]),
            Shape::BicappedHexagonalAntiprismatic => Some(vec![
                vec3(1.0, 0.0, 0.450),
                vec3(0.309017, 0.951057, 0.450),
                vec3(-0.809017, 0.587785, 0.450),
                vec3(-0.809017, -0.587785, 0.450),
                vec3(0.309017, -0.951057, 0.450),
                vec3(0.809017, 0.587785, -0.450),
                vec3(-0.309017, 0.951057, -0.450),
                vec3(-1.0, 0.0, -0.450),
                vec3(-0.309017, -0.951057, -0.450),
                vec3(0.809017, -0.587785, -0.450),
                vec3(0.0, 0.0, 1.40),
                vec3(0.0, 0.0, -1.40),
            ]),
            Shape::Linear => Some(vec![vec3(-1.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0)]),
        }
    }
}
/// The shape a set of directions from a centre forms, and how far it is from it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Classification {
    /// The shape.
    pub shape: Shape,
    /// The root-mean-square difference between the directions' signature and the shape's (see
    /// [`classify`]).
    pub error: f32,
}
/// The signatures of all shapes. Only calculated once.
static SHAPE_REFERENCES: LazyLock<Vec<ShapeReference>> = LazyLock::new(|| {
    Shape::ALL_SHAPES
        .iter()
        .filter_map(|&shape| {
            let signature = shape_signature(&shape.canonical_points()?);
            Some(ShapeReference {
                geometry: shape,
                signature,
            })
        })
        .collect()
});

/// Representation of a canonical shape
struct ShapeReference {
    geometry: Shape,
    signature: Vec<f32>,
}

/// Calculates the signature of a shape
/// The signature is the ordered list of all dot products of normalized vectors from the center to
/// the vertices.  
fn shape_signature(points: &[Vector3<f32>]) -> Vec<f32> {
    let directions: Vec<_> = points.iter().map(|p| p.normalize()).collect();

    let mut signature = Vec::with_capacity((directions.len() * (directions.len() - 1)) / 2);

    for i in 0..directions.len() {
        for j in (i + 1)..directions.len() {
            signature.push(directions[i].dot(directions[j]));
        }
    }
    signature.sort_unstable_by(f32::total_cmp);
    signature
}

/// The shape `vertices` (directions from a centre) form: the one whose signature (the cosines of
/// the angles between each two of them, sorted) is closest to theirs, by root-mean-square
/// difference, where that is at most `max_error`; none for fewer than two vertices, which no shape
/// has.
pub fn classify(vertices: &[Vector3<f32>], max_error: f32) -> Option<Classification> {
    if vertices.len() < 2 {
        return None;
    }

    let signature = shape_signature(vertices);
    SHAPE_REFERENCES
        .iter()
        .filter(|reference| signature.len() == reference.signature.len())
        .map(|reference| {
            let error = signature
                .iter()
                .zip(&reference.signature)
                .map(|(a, b)| (a - b).powi(2))
                .sum::<f32>()
                / signature.len() as f32;

            Classification {
                shape: reference.geometry,
                error: error.sqrt(),
            }
        })
        .filter(|class_result| class_result.error <= max_error)
        .min_by(|a, b| a.error.total_cmp(&b.error))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn an_atom_with_fewer_than_two_neighbours_forms_no_shape() {
        assert_eq!(shape_of(&[]), None);
        assert_eq!(shape_of(&[vec3(1.0, 0.0, 0.0)]), None);
        assert_eq!(
            shape_of(&[vec3(1.0, 0.0, 0.0), vec3(-1.0, 0.0, 0.0)]),
            Some(Shape::Linear)
        );
    }

    #[test]
    fn every_shape_is_classified_as_itself() {
        for &shape in Shape::ALL_SHAPES {
            let points = shape.canonical_points().unwrap();
            assert_eq!(shape_of(&points), Some(shape), "{shape}");
        }
    }
}
