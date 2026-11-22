use cgmath::{Matrix3, Matrix4, Vector3, Vector4, vec3, vec4};

use crate::Measured;

#[cfg_attr(feature = "serde", derive(serde::Deserialize, serde::Serialize))]
#[cfg_attr(feature = "serde", serde(default))]
#[derive(Default, Debug, Clone, PartialEq)]
/// A phase's unit cell and space group. More of what a CIF tells may join it without a breaking
/// release, so it is built from [`Cell::default`] outside the crate.
#[non_exhaustive]
pub struct Cell {
    /// The length of a (Å).
    pub a: Measured,
    /// The length of b (Å).
    pub b: Measured,
    /// The length of c (Å).
    pub c: Measured,
    /// The angle between b and c (°).
    pub alpha: Measured,
    /// The angle between a and c (°).
    pub beta: Measured,
    /// The angle between a and b (°).
    pub gamma: Measured,
    /// The cell's volume (Å³).
    pub volume: Measured,
    /// The space group's Hermann–Mauguin symbol, as the source writes it.
    pub space_group: String,
    /// The space group's number (1–230); 0 where unknown.
    pub space_group_number: u8,
    /// The symmetry operations as the CIF writes them (`-x+y, -x, z+1/3`). They fix the setting
    /// the atoms' coordinates are given in (origin choice, axes); empty for the setting
    /// `space_group` names (see `SpaceGroupSetting::of_cell`), or the standard setting of
    /// `space_group_number` where it names none.
    pub symmetry_operations: Vec<String>,
}

impl Cell {
    /// The Cartesian vector (Å, Y-up) of the fractional vector `fractional` (`x`, `y`, `z` along
    /// a, b, c).
    pub fn cartesian(&self, fractional: Vector3<f64>) -> Vector3<f32> {
        (self.transformation_matrix_y_up() * Self::y_up_input(fractional)).truncate()
    }

    /// The fractional vector `fractional` (`x`, `y`, `z` along a, b, c) in the order
    /// [`Self::transformation_matrix_y_up`] takes it: `(x, z, y, 0)`, c up.
    pub fn y_up_input(fractional: Vector3<f64>) -> Vector4<f32> {
        vec4(
            fractional.x as f32,
            fractional.z as f32,
            fractional.y as f32,
            0.0,
        )
    }

    /// Returns the transformation matrix for the unit cell with Y-up coordinates. It takes a
    /// fractional vector in the order a, c, b (see [`Self::y_up_input`]). a lies along x, c in the
    /// xy plane (up for an orthogonal cell), and b towards −z, so a, b and c form a right-handed
    /// frame, as crystallographic axes do: a structure isn't drawn as its mirror image.
    pub fn transformation_matrix_y_up(&self) -> Matrix4<f32> {
        let x = self.a.value as f32;
        let y = self.c.value as f32;
        let z = self.b.value as f32;
        let angle_yz = (self.alpha.value as f32).to_radians();
        let angle_xz = (self.gamma.value as f32).to_radians();
        let angle_yx = (self.beta.value as f32).to_radians();

        let x0 = x;
        let x1 = 0.0;
        let x2 = 0.0;

        let y0 = y * (angle_yx.cos());
        let y1 = y * (angle_yx.sin());
        let y2 = 0.0;

        let z0 = z * (angle_xz.cos());
        let z1 = z * (angle_yz.cos() - angle_xz.cos() * angle_yx.cos()) / angle_yx.sin();
        let z2 = -z
            * (angle_xz.sin().powi(2)
                - (angle_yz.cos() - angle_xz.cos() * angle_yx.cos()).powi(2)
                    / angle_yx.sin().powi(2))
            .sqrt();

        let c0 = vec3(x0, x1, x2);
        let c1 = vec3(y0, y1, y2);
        let c2 = vec3(z0, z1, z2);

        let m = Matrix3::from_cols(c0, c1, c2);

        Matrix4::from(m)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_b_and_c_form_a_right_handed_frame_of_the_cells_volume() {
        // Triclinic, so every angle enters.
        let cell = Cell {
            a: 5.0.into(),
            b: 6.0.into(),
            c: 7.0.into(),
            alpha: 80.0.into(),
            beta: 95.0.into(),
            gamma: 110.0.into(),
            ..Default::default()
        };
        let [a, b, c] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
            .map(|axis| cell.cartesian(Vector3::from(axis)));

        let triple_product = cgmath::InnerSpace::dot(a, b.cross(c));

        let [alpha, beta, gamma] = [80f32, 95.0, 110.0].map(f32::to_radians);
        let volume = 5.0
            * 6.0
            * 7.0
            * (1.0 - alpha.cos().powi(2) - beta.cos().powi(2) - gamma.cos().powi(2)
                + 2.0 * alpha.cos() * beta.cos() * gamma.cos())
            .sqrt();
        assert!(
            (triple_product - volume).abs() < 1e-3,
            "a · (b × c) = {triple_product}, the cell's volume is {volume}"
        );
    }

    #[test]
    fn c_points_up_and_b_along_negative_z() {
        let cell = Cell {
            a: 2.0.into(),
            b: 3.0.into(),
            c: 4.0.into(),
            alpha: 90.0.into(),
            beta: 90.0.into(),
            gamma: 90.0.into(),
            ..Default::default()
        };
        let close = |v: Vector3<f32>, expected: [f32; 3]| {
            cgmath::InnerSpace::magnitude(v - Vector3::from(expected)) < 1e-5
        };

        assert!(close(
            cell.cartesian(Vector3::new(1.0, 0.0, 0.0)),
            [2.0, 0.0, 0.0]
        ));
        assert!(close(
            cell.cartesian(Vector3::new(0.0, 1.0, 0.0)),
            [0.0, 0.0, -3.0]
        ));
        assert!(close(
            cell.cartesian(Vector3::new(0.0, 0.0, 1.0)),
            [0.0, 4.0, 0.0]
        ));
    }
}
